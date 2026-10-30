#!/usr/bin/env python3
"""
bispec_biphase.py
-----------------
Models whose bispectrum has structured biphases on either side of the QPO,
compared with the sideband predictions in biphase_models.md.

A slow flux component s(t) appears in the light curve and modulates the
complex envelope of a QPO at f0 through a response H(nu). For a coupling
frequency nu = f1 (the flux leg of the bispectrum) the biphase is

    upper side  (f1, f2 ~ f0,      f3 = f2 + f1 ~ f0 + f1):  -arg H(+nu)
    lower side  (f1, f2 ~ f0 - f1, f3 = f2 + f1 ~ f0):       +arg H(-nu)

Implemented models
    fm_response   QPO frequency f(t) = f0 + kappa r(t), where r is s passed
                  through a first-order low-pass of response time tau_r
                  (tau_r dr/dt = s - r). H(nu) = (kappa/nu) / (1 + 2 pi i nu tau_r),
                  so upper = atan(2 pi nu tau_r) and lower = pi + atan(2 pi nu tau_r):
                  the pair stays a pi apart and rotates with nu at a rate set
                  by tau_r.

The side biphases are measured for FWHM < f1 < f0/2 - FWHM/2, where the two
sides do not overlap and each triple has a single flux leg: above f0/2 the
lower-side triple (f1, f0 - f1, f0) is dominated by the term whose flux leg is
f0 - f1 instead.

Light curves are Gaussian rates (no counting noise unless --poisson).
Bispectrum maps use bispec_sims.compute_bispec (stingray, Hann-tapered
segments). Side biphases are measured with the same segments and taper, with
bootstrap errors over segments; they are checked against the stingray maps.

Example
-------
    python bispec_biphase.py fm_response --outdir figures/biphase
"""
from __future__ import annotations

import argparse
import os
import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.signal import lfilter

from stingray import Lightcurve

import qpo_sims as qs
import bispec_sims as bsm

warnings.filterwarnings("ignore", message="The stingray bispectrum interface")

F0 = 1.5            # QPO frequency [Hz]
Q0 = 10.0           # QPO quality factor
FWHM = F0 / Q0
_PI_TICKS = bsm._PI_TICKS


# =========================================================================
# Signal generators
# =========================================================================

def ou_unit(n, dt, tau, rng):
    """Unit-rms Ornstein-Uhlenbeck process (Lorentzian at zero frequency)."""
    a = np.exp(-dt / tau)
    x = lfilter([np.sqrt(1 - a * a)], [1, -a], rng.standard_normal(n))
    return x / x.std()


def lowpass(x, dt, tau):
    """First-order low-pass, tau dr/dt = x - r (identity for tau = 0)."""
    if tau <= 0:
        return x.copy()
    a = np.exp(-dt / tau)
    return lfilter([1 - a], [1, -a], x)


def qpo(phase_mod, dt, rng, fwhm=FWHM):
    """Gaussian QPO with an imposed phase modulation.

    z = w exp(i (2 pi f0 t + phase_mod)), with w a zero-frequency complex OU
    process of FWHM `fwhm`. This has the same distribution as the stochastic
    oscillator z[n+1] = exp((-g + 2 pi i f[n]) dt) z[n] + xi[n] of
    bispec_nonstat.py for f = f0 + d(phase_mod)/dt / (2 pi), because
    exp(-i phi) xi is again circular Gaussian noise, but it is vectorised.
    Returns the complex z scaled to unit rms of its real part.
    """
    n = phase_mod.size
    g = np.pi * fwhm
    a = np.exp(-g * dt)
    xi = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    w = lfilter([1.0], [1, -a], xi)
    t = np.arange(n) * dt
    z = w * np.exp(1j * (2 * np.pi * F0 * t + phase_mod))
    return z / np.sqrt(np.mean(z.real ** 2))


def to_lc(x, dt, mean=1e3, rng=None, poisson=False):
    t = np.arange(x.size) * dt
    rate = mean * (1 + x)
    if poisson:
        rate = rng.poisson(np.clip(rate * dt, 0, None)) / dt
    return Lightcurve(t, rate, dt=dt, skip_checks=True, err_dist="gauss")


# =========================================================================
# Model 1: frequency coupled to the flux with a response time
# =========================================================================

def model_fm_response(args, rng):
    """Returns a list of (label, sweep value, x1, x2, prediction function)."""
    n = int(round(args.nseg * args.segment / args.dt))
    dt = args.dt
    s = ou_unit(n, dt, args.tau_s, rng)               # flux component in the light curve
    cases = []
    for tau_r in args.tau_r:
        r = lowpass(s, dt, tau_r)
        r /= r.std()                                   # frequency rms = kappa for every tau_r
        phase_mod = 2 * np.pi * args.kappa * np.cumsum(r) * dt
        z = qpo(phase_mod, dt, rng)
        x1 = args.c_s1 * s + args.a_q1 * z.real
        x2 = args.c_s2 * s + args.a_q2 * z.real
        cases.append(dict(label=f"{tau_r:g} s", sweep=tau_r, x1=x1, x2=x2,
                          pred=lambda nu, tr=tau_r: _pred_fm(nu, tr), fit_tau=True))
    return cases


def _pred_fm(nu, tau_r):
    th = np.arctan(2 * np.pi * nu * tau_r)
    return th + np.pi, th          # lower, upper


def _fit_tau(f1, sides):
    """Weighted least-squares tau_r from the mean side biphase, plus mean Delta."""
    lo, up = sides["lower"]["phase"], sides["upper"]["phase"]
    err = 0.5 * np.hypot(sides["lower"]["err"], sides["upper"]["err"])
    taus = np.concatenate([[0.0], np.geomspace(1e-3, 30, 2000)])
    chi2 = []
    for tr in taus:
        lp, upp = _pred_fm(f1, tr)
        sig = 0.5 * (_wrap_near(lo, lp) + _wrap_near(up, upp))
        chi2.append(np.sum(((sig - 0.5 * (lp + upp)) / err) ** 2))
    chi2 = np.array(chi2)
    k = int(np.argmin(chi2))
    ok = taus[chi2 <= chi2[k] + 1.0]
    lp, upp = _pred_fm(f1, taus[k])
    half = 0.5 * (_wrap_near(lo, lp) - _wrap_near(up, upp))
    w = 1 / err ** 2
    print(f"    fit: tau_r = {taus[k]:.3f} s (1 sigma {ok.min():.3f}-{ok.max():.3f}), "
          f"chi2 = {chi2[k]:.1f}/{f1.size - 1}; mean Delta = "
          f"{np.degrees(np.sum(w * half) / np.sum(w)):.1f} +/- "
          f"{np.degrees(1 / np.sqrt(np.sum(w))):.1f} deg (pred 90)")


MODELS = {"fm_response": (model_fm_response, r"$\tau_r$ (s)")}


# =========================================================================
# Side biphases with bootstrap errors
# =========================================================================

def segment_ffts(lc, segment, window):
    segs = bsm._segments(lc, segment, window)
    y = np.array([sg.counts - sg.counts.mean() for sg in segs])
    return np.fft.rfft(y, axis=1), np.fft.rfftfreq(y.shape[1], lc.dt)


def side_biphase(X, freq, f1_edges, half_width, n_boot=200, rng=None):
    """Biphase on the lower and upper side of the QPO per f1 bin.

    Lower: f3 = f1 + f2 within f0 +/- half_width.  Upper: f2 within f0 +/- half_width.
    For each f1 bin, the triple products X(f1) X(f2) X*(f1 + f2) of each segment
    are summed over the region; the biphase is the phase of the mean over
    segments, and its error the circular spread of segment-bootstrap means.
    """
    df = freq[1] - freq[0]
    nf = freq.size
    k0 = np.arange(nf)
    out = {side: dict(phase=[], err=[], b2=[]) for side in ("lower", "upper")}
    centres = []
    for lo, hi in zip(f1_edges[:-1], f1_edges[1:]):
        i1 = k0[(freq >= lo) & (freq < hi)]
        centres.append(freq[i1].mean())
        for side in ("lower", "upper"):
            seg_sum = np.zeros(X.shape[0], complex)
            norm1 = np.zeros(X.shape[0])
            norm2 = np.zeros(X.shape[0])
            for i in i1:
                if side == "upper":
                    j = k0[np.abs(freq - F0) <= half_width]
                else:
                    j = k0[np.abs(freq + freq[i] - F0) <= half_width]
                j = j[(j > 0) & (i + j < nf)]
                t = X[:, [i]] * X[:, j] * np.conj(X[:, i + j])
                seg_sum += t.sum(1)
                norm1 += (np.abs(X[:, [i]] * X[:, j]) ** 2).sum(1)
                norm2 += (np.abs(X[:, i + j]) ** 2).sum(1)
            m = seg_sum.mean()
            b2 = np.abs(seg_sum.sum()) ** 2 / (norm1.sum() * norm2.sum())
            idx = rng.integers(0, X.shape[0], (n_boot, X.shape[0]))
            boot = seg_sum[idx].mean(1)
            dphi = np.angle(boot * np.conj(m))
            out[side]["phase"].append(np.angle(m))
            out[side]["err"].append(np.sqrt(np.mean(dphi ** 2)))
            out[side]["b2"].append(b2)
    return np.array(centres), {k: {q: np.array(v) for q, v in d.items()} for k, d in out.items()}


def _wrap_near(phi, ref):
    """Wrap phi into (ref - pi, ref + pi]."""
    return ref + np.angle(np.exp(1j * (phi - ref)))


# =========================================================================
# Plots
# =========================================================================

def plot_side_biphase(results, sweep_label, savepath, f_lim):
    """Top: lower (open) and upper (filled) side biphase vs f1 with predictions.
    Middle: mean of the two (symmetric part).  Bottom: half difference."""
    fig, axes = qs._setup_xspec_axes(figsize=(7, 8.5), nrows=3)
    ax1, ax2, ax3 = axes
    colours = plt.cm.viridis(np.linspace(0.1, 0.85, len(results)))
    nu = np.logspace(np.log10(f_lim[0]), np.log10(f_lim[1]), 300)
    for res, col in zip(results, colours):
        f1, sides, pred = res["f1"], res["sides"], res["pred"]
        lo_p, up_p = pred(nu)
        lo_pf, up_pf = pred(f1)
        lo = _wrap_near(sides["lower"]["phase"], lo_pf)
        up = _wrap_near(sides["upper"]["phase"], up_pf)
        ax1.plot(nu, up_p, color=col, lw=1.2)
        ax1.plot(nu, lo_p, color=col, lw=1.2, ls="--")
        ax1.errorbar(f1, up, sides["upper"]["err"], fmt="o", color=col, ms=5, capsize=0,
                     label=res["label"])
        ax1.errorbar(f1, lo, sides["lower"]["err"], fmt="o", mfc="white", color=col, ms=5,
                     capsize=0)
        mean_p = 0.5 * (lo_p + up_p)
        half_p = 0.5 * (lo_p - up_p)
        mean_m = 0.5 * (lo + up)
        half_m = 0.5 * (lo - up)
        err = 0.5 * np.hypot(sides["lower"]["err"], sides["upper"]["err"])
        ax2.plot(nu, mean_p, color=col, lw=1.2)
        ax2.errorbar(f1, mean_m, err, fmt="o", color=col, ms=5, capsize=0)
        ax3.plot(nu, half_p, color=col, lw=1.2)
        ax3.errorbar(f1, half_m, err, fmt="o", color=col, ms=5, capsize=0)

    ax1.set_ylabel("Side biphase (rad)")
    ax2.set_ylabel("Mean (rad)")
    ax3.set_ylabel("Half difference (rad)")
    ax1.set_ylim(-np.pi / 2, 3 * np.pi / 2)
    ax1.set_yticks([0, np.pi / 2, np.pi])
    ax1.set_yticklabels(["0", r"$\pi/2$", r"$\pi$"])
    ax2.set_ylim(0, np.pi)
    ax2.set_yticks([np.pi / 4, np.pi / 2, 3 * np.pi / 4])
    ax2.set_yticklabels([r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$"])
    ax3.set_ylim(0, np.pi)
    ax3.set_yticks([np.pi / 4, np.pi / 2, 3 * np.pi / 4])
    ax3.set_yticklabels([r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$"])
    for ax in axes:
        ax.axvline(FWHM, color="0.6", ls=":", lw=0.8)
    ax3.legend(handles=[Line2D([0], [0], color=col, marker="o", label=res["label"])
                        for res, col in zip(results, colours)],
               title=sweep_label, fontsize=8, loc="lower right", ncol=4, title_fontsize=8)
    ax1.legend(handles=[
        Line2D([0], [0], color="gray", marker="o", ls="-", label="upper side"),
        Line2D([0], [0], color="gray", marker="o", mfc="white", ls="--", label="lower side"),
    ], fontsize=8, loc="lower right")
    ax3.set_xlabel(r"$f_1$ (Hz)")
    from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter
    ticks = [t for t in (0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0) if f_lim[0] <= t <= f_lim[1]]
    ax3.xaxis.set_major_locator(FixedLocator(ticks))
    ax3.xaxis.set_major_formatter(ScalarFormatter())
    ax3.xaxis.set_minor_formatter(NullFormatter())
    axes[0].set_xlim(*f_lim)
    fig.align_ylabels(axes)
    fig.savefig(savepath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {savepath}")


# =========================================================================
# Driver
# =========================================================================

def run(name, args):
    fn, sweep_label = MODELS[name]
    rng = np.random.default_rng(args.seed)
    print(f"[{name}] {args.nseg} segments of {args.segment:g} s, dt = {args.dt:g} s")
    cases = fn(args, rng)

    f1_edges = np.geomspace(args.f1_min, args.f1_max, args.n_f1 + 1)
    results, sets = [], []
    for c in cases:
        lc1, lc2 = to_lc(c["x1"], args.dt, rng=rng, poisson=args.poisson), \
            to_lc(c["x2"], args.dt, rng=rng, poisson=args.poisson)
        X, freq = segment_ffts(lc1, args.segment, args.window)
        f1, sides = side_biphase(X, freq, f1_edges, FWHM / 2, rng=rng)
        results.append(dict(label=c["label"], f1=f1, sides=sides, pred=c["pred"]))
        lo_p, up_p = c["pred"](f1)
        print(f"  {c['label']}")
        for k in range(f1.size):
            lo = _wrap_near(sides["lower"]["phase"][k], lo_p[k])
            up = _wrap_near(sides["upper"]["phase"][k], up_p[k])
            print(f"    f1 = {f1[k]:.2f} Hz (b2 {sides['lower']['b2'][k]:.3f}/{sides['upper']['b2'][k]:.3f}): lower {np.degrees(lo):7.1f} "
                  f"(pred {np.degrees(lo_p[k]):6.1f}) +/- {np.degrees(sides['lower']['err'][k]):4.1f}"
                  f"   upper {np.degrees(up):7.1f} (pred {np.degrees(up_p[k]):6.1f}) "
                  f"+/- {np.degrees(sides['upper']['err'][k]):4.1f}")
        if c.get("fit_tau"):
            _fit_tau(f1, sides)
        if args.maps:
            # maps (and the stingray cross-check) on the first map_nseg segments
            nm = int(round(args.map_nseg * args.segment / args.dt))
            m1 = to_lc(c["x1"][:nm], args.dt, rng=rng, poisson=args.poisson)
            m2 = to_lc(c["x2"][:nm], args.dt, rng=rng, poisson=args.poisson)
            bs = bsm.compute_bispec(m1, m2, args.dt, args.segment, args.window,
                                    label=f"{c['sweep']:g}")
            Xm, fm = segment_ffts(m1, args.segment, args.window)
            _, sm = side_biphase(Xm, fm, f1_edges, FWHM / 2, n_boot=2, rng=rng)
            bs.check_phase = sm["upper"]["phase"][0]
            sets.append(bs)

    f_lim = (args.f1_min, args.f1_max)
    plot_side_biphase(results, sweep_label, os.path.join(args.outdir, f"{name}_side_biphase.{args.ext}"),
                      f_lim)

    if args.maps:
        grid = bsm.LogGrid(sets[0].freq, args.rebin)
        null = bsm.NullModel(args.dt, args.segment, sets[0].m, grid, window=args.window,
                             n_real=4)
        keep = bsm._subsample(len(sets), args.max_cols)
        labels = [f"{sweep_label} = {sets[i].label}" for i in keep]
        map_lim = (0.1, 4.0)
        for phase, tag in ((False, "bicoherence"), (True, "biphase")):
            bsm.plot_maps([sets[i] for i in keep], grid, null,
                          os.path.join(args.outdir, f"{name}_{tag}_maps.{args.ext}"),
                          map_lim, (F0,), col_labels=labels, phase=phase)
        # f1 band slice across the QPO, from the stingray estimator
        f_ref = np.sqrt(args.f1_min * args.f1_max)
        bsm.plot_slices([sets[i] for i in keep], grid, null,
                        os.path.join(args.outdir, f"{name}_bicoherence_slices.{args.ext}"),
                        (0.1, 4.0), (F0,), f_ref, 2.0, legend_title=sweep_label)
        # consistency check: stingray vs segment estimator, upper side, first f1 bin
        for bs, res in zip(sets, results):
            fr = bs.freq
            i1 = (fr >= f1_edges[0]) & (fr < f1_edges[1])
            j = np.abs(fr - F0) <= FWHM / 2
            ph = np.angle(bs.bispec["11"][np.ix_(i1, j)].sum())
            d = np.degrees(np.angle(np.exp(1j * (ph - bs.check_phase))))
            print(f"  stingray vs segment estimator ({res['label']}): {d:+.2f} deg")


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("models", nargs="+", choices=list(MODELS))
    p.add_argument("--outdir", default="figures/biphase")
    p.add_argument("--ext", default="png", choices=["png", "pdf"])
    p.add_argument("--dt", type=bsm._parse_dt, default=1 / 64)
    p.add_argument("--segment", type=float, default=16.0)
    p.add_argument("--nseg", type=int, default=512)
    p.add_argument("--window", default="hann", choices=["hann", "none"])
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--poisson", action="store_true")
    # flux component and amplitudes (fractional rms)
    p.add_argument("--tau_s", type=float, default=0.25, help="OU timescale of the flux component (s)")
    p.add_argument("--c_s1", type=float, default=0.20)
    p.add_argument("--c_s2", type=float, default=0.10)
    p.add_argument("--a_q1", type=float, default=0.10)
    p.add_argument("--a_q2", type=float, default=0.10)
    # fm_response
    p.add_argument("--kappa", type=float, default=0.12, help="QPO frequency rms (Hz)")
    p.add_argument("--tau_r", type=float, nargs="+", default=[0.0, 0.1, 0.3, 1.0])
    # side-biphase measurement
    p.add_argument("--f1_min", type=float, default=0.2)
    p.add_argument("--f1_max", type=float, default=0.65,
                   help="Keep below f0/2 - FWHM/2: above it the lower-side triple "
                        "(f1, f0 - f1, f0) is dominated by the flux leg at f0 - f1.")
    p.add_argument("--n_f1", type=int, default=4)
    # maps
    p.add_argument("--no-maps", dest="maps", action="store_false")
    p.add_argument("--map-nseg", type=int, default=256, dest="map_nseg",
                   help="Segments used for the 2D maps (the side biphases use all).")
    p.add_argument("--rebin", type=float, default=0.05)
    p.add_argument("--max-cols", type=int, default=6, dest="max_cols")
    return p


def main():
    args = build_parser().parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    for m in args.models:
        run(m, args)


if __name__ == "__main__":
    main()

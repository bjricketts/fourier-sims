#!/usr/bin/env python3
"""
bispec_sims.py
--------------
Bispectral diagnostics for the simulated scenarios in qpo_sims.py.

Every model is run exactly as in qpo_sims.py (same defaults, seeds and
light curves). The two output light curves are rebinned to a coarser time
resolution and passed to the stingray bispectrum implementation (direct
Fourier estimator, Maccarone 2013; stingray main branch), after splitting
them into Hann-tapered segments (see _segments):

    band 1 auto    B11(f1, f2) = < X1(f1) X1(f2) X1*(f1 + f2) >
    band 2 auto    B22(f1, f2) = < X2(f1) X2(f2) X2*(f1 + f2) >
    cross          B112(f1, f2) = < X1(f1) X1(f2) X2*(f1 + f2) >

The bicoherence is the squared Kim & Powers (1979) normalisation, without the
~1/M bias subtraction. For plotting, b^2 is averaged over log-spaced cells of
the (f1, f2) plane (see LogGrid), and each cell is compared with a 99% null
threshold calibrated on Gaussian white-noise light curves segmented
and tapered identically (see NullModel): cells above it show coupling that a
linear Gaussian process would not produce. No Poisson-noise correction is
applied (the Wirnitzer 1985 correction assumes untapered counts); Poisson white
noise at the count rates of these simulations (down to well below 40 counts
per 1/512 s bin) passes this null test at the nominal 1% rate.

Each run produces two figures:

    {name}_bicoherence_maps.{ext}     2D log10 b^2 (and biphase above the null)
    {name}_biphase_maps.{ext}         sweeps only: 2D biphase above the null
    {name}_bicoherence_slices.{ext}   b^2 along f1 = f2 and in a band of f2 around f_ref,
                                      plus the biphase in that band

The cross bispectrum lives on the signed frequency plane; only the
sum-frequency quadrant (f1, f2 > 0) is plotted.

Example
-------
    python bispec_sims.py all --outdir figures/bispec
    python bispec_sims.py damping fm_phase --outdir figures/bispec --dt 1/128
"""
from __future__ import annotations

import argparse
import gc
import os
import warnings
from dataclasses import dataclass, field
from fractions import Fraction

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from scipy.signal.windows import hann

from stingray import Lightcurve
from stingray.bispectrum import AveragedBispectrum, AveragedCrossBispectrum

import qpo_sims as qs

warnings.filterwarnings("ignore", message="The stingray bispectrum interface")


# =========================================================================
# Paper runs: name -> (qpo_sims model, CLI overrides, reference frequency)
# =========================================================================

@dataclass
class Run:
    model: str
    cli: list = field(default_factory=list)
    f_ref: float = 1.5          # centre of the f2 slice band [Hz]
    q_ref: float = 10.0         # slice band is f_ref * (1 +/- 1/(2 q_ref))
    guides: tuple = None        # guide frequencies if the model sets none
    fig: str = ""               # paper figure, for the log only


PAPER_RUNS = {
    "damping":              Run("damping", f_ref=5.0, fig="2"),
    "resonance":            Run("resonance", f_ref=5.0, fig="3"),
    "coherent_sum":         Run("coherent_sum", f_ref=3.0, q_ref=5.0, fig="4"),
    "hidden":               Run("hidden", ["--seed", "0"], fig="5"),
    "multi_coherent":       Run("multi_coherent", fig="6"),
    "multi_incoherent":     Run("multi_incoherent", fig="7"),
    # slice on the swept mode (mode 2, 1.6 Hz)
    "multi_multiplicative": Run("multi_multiplicative", f_ref=1.6, fig="8"),
    # FM: peak width set by the wander (sigma_f = 0.15 Hz) -> Q ~ 5
    "fm":                   Run("fm", q_ref=5.0, fig="9"),
    "fm_phase":             Run("fm", ["--fm_sigma_phase", f"{np.pi / 2}"],
                                q_ref=5.0, fig="10"),
    "fm_offset":            Run("fm", ["--fm_df", "0.15"], q_ref=5.0, fig="11"),
    "count_rate_fm":        Run("count_rate_fm", q_ref=5.0, guides=(1.5,), fig="12"),
}


# =========================================================================
# Bispectrum computation
# =========================================================================

@dataclass
class BispecSet:
    """Bicoherence / biphase for one pair of light curves, on the positive
    (f1, f2) grid shared by the auto and cross estimators."""
    freq: np.ndarray
    b2: dict            # key in {"11", "22", "112"} -> (nf, nf) float32
    bispec: dict        # same keys -> (nf, nf) complex64 bispectrum
    m: int
    label: str = ""


def _parse_dt(s):
    return float(Fraction(s))


def _rebin(lc, dt_new):
    if np.isclose(lc.dt, dt_new):
        return lc
    return lc.rebin(dt_new)


def _segments(lc, segment, window):
    """Split a light curve into contiguous segments and taper each one.

    Tapering matters here: above a high-Q QPO the PSD falls ~f^-4, steeper than
    the sidelobes of a rectangular window, so without a taper the high-frequency
    Fourier amplitudes of each segment are dominated by phase-locked leakage
    from the strong low-frequency variability, which shows up as large
    spurious bicoherence. The Hann taper is applied to the mean-subtracted
    counts and rescaled to preserve the variance (sqrt(8/3)); the mean is
    added back so the segments remain valid count light curves.
    """
    n = int(round(segment / lc.dt))
    if window == "hann":
        w = hann(n, sym=False) * np.sqrt(8.0 / 3.0)
    elif window == "none":
        w = np.ones(n)
    else:
        raise ValueError(f"Unknown window {window!r}")
    out = []
    for i in range(lc.n // n):
        sl = slice(i * n, (i + 1) * n)
        y = np.asarray(lc.counts[sl], dtype=float)
        mu = y.mean()
        out.append(Lightcurve(lc.time[sl], mu + (y - mu) * w, dt=lc.dt,
                              skip_checks=True, err_dist="gauss"))
    return out


def _positive_quadrant(cb):
    """Restrict a signed-grid cross bispectrum to f1, f2 > 0."""
    pos = cb.freq > 0
    ix = np.ix_(pos, pos)
    return cb.freq[pos], cb.bicoherence[ix], cb.bispec[ix]


def compute_bispec(lc1, lc2, dt, segment, window="hann", label=""):
    s1 = _segments(_rebin(lc1, dt), segment, window)
    s2 = _segments(_rebin(lc2, dt), segment, window)
    kw = dict(segment_size=segment, silent=True)
    # No Poisson subtraction: the Wirnitzer correction assumes untapered
    # counts. The Poisson contribution is instead part of the null test
    # (see NullModel / the README note).
    b11 = AveragedBispectrum(s1, **kw)
    b22 = AveragedBispectrum(s2, **kw)
    b112 = AveragedCrossBispectrum(s1, s1, s2, channels_overlap=False, **kw)
    if not (b11.m == b22.m == b112.m == len(s1)):
        raise RuntimeError("Unexpected number of averaged segments.")

    freq = b11.freq
    f_c, bc, pc = _positive_quadrant(b112)
    if not np.allclose(f_c, freq):
        raise RuntimeError("Auto and cross bispectrum frequency grids differ.")

    out = BispecSet(
        freq=freq,
        b2={"11": b11.bicoherence.astype(np.float32),
            "22": b22.bicoherence.astype(np.float32),
            "112": bc.astype(np.float32)},
        bispec={"11": b11.bispec.astype(np.complex64),
                "22": b22.bispec.astype(np.complex64),
                "112": pc.astype(np.complex64)},
        m=b11.m, label=label)
    del b11, b22, b112
    gc.collect()
    return out


class LogGrid:
    """Logarithmic 2D cells on the (f1, f2) plane.

    The native bispectral grid is linear in frequency, so at high frequency a
    log-scaled map is dominated by thousands of independent, noise-limited
    pixels. Here each native frequency is assigned to a log-spaced bin of
    fractional width `rebin`, and maps are reduced over the resulting cells:

      * b^2: incoherent mean of the native b^2 over the valid pixels of a
        cell (the bicoherence of each pixel is kept; only its variance drops);
      * biphase: phase of the summed complex bispectrum over the cell.

    With rebin = 0 every native frequency is its own bin (no rebinning).
    """

    def __init__(self, freq, rebin):
        self.freq_native = freq
        df = freq[1] - freq[0]
        if rebin <= 0:
            edges = np.concatenate([freq - df / 2, [freq[-1] + df / 2]])
        else:
            edges = [freq[0] - df / 2]
            while edges[-1] < freq[-1] + df / 2:
                edges.append(max(edges[-1] * (1 + rebin), edges[-1] + df))
            edges = np.asarray(edges)
        self.idx = np.digitize(freq, edges) - 1
        used = np.unique(self.idx)
        # Re-index to the occupied bins only, so no empty rows/columns
        remap = -np.ones(len(edges), int)
        remap[used] = np.arange(used.size)
        self.idx = remap[self.idx]
        self.nb = used.size
        self.edges = np.concatenate([edges[used], [edges[used[-1] + 1]]])
        self.centre = np.array([freq[self.idx == i].mean() for i in range(self.nb)])

    def _reduce(self, arr):
        """Sum of `arr` and number of finite pixels per 2D cell."""
        ok = np.isfinite(arr)
        I, J = np.meshgrid(self.idx, self.idx, indexing="ij")
        tot = np.zeros((self.nb, self.nb), dtype=arr.dtype)
        cnt = np.zeros((self.nb, self.nb))
        np.add.at(tot, (I[ok], J[ok]), arr[ok])
        np.add.at(cnt, (I[ok], J[ok]), 1)
        return tot, cnt

    def b2(self, b2):
        tot, cnt = self._reduce(b2)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(cnt > 0, tot / cnt, np.nan), cnt

    def phase(self, bispec):
        tot, cnt = self._reduce(bispec)
        return np.where(cnt > 0, np.angle(tot), np.nan)

    def band_mask(self, f_lo, f_hi):
        """Native frequencies within [f_lo, f_hi] (at least the nearest one)."""
        m = (self.freq_native >= f_lo) & (self.freq_native <= f_hi)
        if not m.any():
            m[np.argmin(np.abs(self.freq_native - 0.5 * (f_lo + f_hi)))] = True
        return m

    def _reduce_slice(self, arr, mask):
        """Sum over the native f2 in `mask`, then over each log f1 cell."""
        sub = arr[:, mask]
        ok = np.isfinite(sub)
        tot = np.zeros(self.nb, dtype=arr.dtype)
        cnt = np.zeros(self.nb)
        np.add.at(tot, self.idx, np.where(ok, sub, 0).sum(axis=1))
        np.add.at(cnt, self.idx, ok.sum(axis=1))
        return tot, cnt

    def slice_b2(self, b2, mask):
        tot, cnt = self._reduce_slice(b2, mask)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(cnt > 0, tot / cnt, np.nan), cnt

    def slice_phase(self, bispec, mask):
        tot, cnt = self._reduce_slice(bispec, mask)
        return np.where(cnt > 0, np.angle(tot), np.nan)


class NullModel:
    """99% threshold on the cell-averaged b^2 of a linear Gaussian process.

    A native-resolution b^2 under the null is close to exponential with mean
    mu0 ~ 1/M, so a mean over k pixels is close to Gamma(k, mu0/k). mu0 and a
    single multiplicative correction kappa (for the correlation between
    neighbouring bispectral pixels, which share Fourier amplitudes) are
    calibrated on `n_real` Gaussian white-noise light curves segmented exactly
    like the data, so that the per-cell threshold gives a 1% false-positive
    rate on the null realisations.
    """

    def __init__(self, dt, segment, n_seg, grid: LogGrid, window="hann", q=0.99,
                 n_real=10, seed=1234):
        from scipy.stats import gamma
        self._gamma, self.q, self.grid = gamma, q, grid
        self._natives = []
        n = int(round(segment / dt))
        rng = np.random.default_rng(seed)
        t = np.arange(n * n_seg) * dt
        maps, natives = [], []
        for _ in range(n_real):
            lc = Lightcurve(t, rng.normal(1e3, 30.0, t.size), dt=dt,
                            skip_checks=True, err_dist="gauss")
            b = AveragedBispectrum(_segments(lc, segment, window),
                                   segment_size=segment, silent=True)
            natives.append(b.bicoherence[np.isfinite(b.bicoherence)])
            maps.append(grid.b2(b.bicoherence))
            self._natives.append(b.bicoherence.astype(np.float32))
        self.mu0 = float(np.mean(np.concatenate(natives)))
        cnt = maps[0][1]
        base = self._base(cnt)
        ratios = np.concatenate([(m / base)[cnt > 0] for m, _ in maps])
        self.kappa = float(np.percentile(ratios, 100 * q))
        self.threshold = self.kappa * base          # (nb, nb) per cell
        self.native_level = float(np.percentile(np.concatenate(natives), 100 * q))

    def _base(self, cnt):
        k = np.maximum(cnt, 1)
        with np.errstate(invalid="ignore"):
            return np.where(cnt > 0, self._gamma.ppf(self.q, a=k, scale=self.mu0 / k),
                            np.nan)

    def slice_threshold(self, mask):
        """Per-f1-cell 99% threshold for `LogGrid.slice_b2` with this mask,
        with its own kappa calibrated on the stored null realisations."""
        vals = [self.grid.slice_b2(m, mask) for m in self._natives]
        cnt = vals[0][1]
        base = self._base(cnt)
        ratios = np.concatenate([(v / base)[cnt > 0] for v, _ in vals])
        kappa = float(np.percentile(ratios, 100 * self.q))
        return kappa * base, kappa


# =========================================================================
# Plotting
# =========================================================================

_C_GUIDE = "#b2abd2"
_CHANNELS = [("11", "Band 1"), ("22", "Band 2"), ("112", r"Cross ($1\,1\,2^{*}$)")]
_PI_TICKS = ([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi],
             [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])


def _guides(ax, f_guides, f_lim):
    for f in f_guides:
        ax.axvline(f, color=_C_GUIDE, ls="--", lw=1.2)
        ax.axhline(f, color=_C_GUIDE, ls="--", lw=1.2)
    ax.plot(f_lim, f_lim, color="0.6", ls=":", lw=0.8)


def _style_map_ax(ax, f_lim, left, bottom):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*f_lim)
    ax.set_ylim(*f_lim)
    ax.set_aspect("equal")
    if bottom:
        ax.set_xlabel(r"$f_1$ (Hz)")
    else:
        ax.tick_params(labelbottom=False)
    if left:
        ax.set_ylabel(r"$f_2$ (Hz)")
    else:
        ax.tick_params(labelleft=False)


def plot_maps(sets, grid, null, savepath, f_lim, f_guides, col_labels=None,
              vmin=-2.5, vmax=0.0, phase=False):
    """Grid of cell-averaged log10 b^2 or biphase maps.

    One BispecSet -> 2 rows (b^2, biphase where b^2 > 99% null) x 3 channels.
    Several (a sweep) -> 3 rows (channels) x N sweep points, showing b^2, or
    the biphase where b^2 > 99% null if `phase` is True.
    In b^2 maps, cells above the 99% null threshold are outlined in white; in
    biphase maps, cells below it are grey.
    """
    single = len(sets) == 1
    nrows, ncols = (2, 3) if single else (3, len(sets))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols + 1.2, 3.2 * nrows),
                             squeeze=False, gridspec_kw=dict(wspace=0.06, hspace=0.06))
    e = grid.edges
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("white")
    pcmap = plt.get_cmap("twilight").copy()
    pcmap.set_bad("0.55")   # distinct from the light +/-pi ends of twilight
    im_b = im_p = None
    for c in range(ncols):
        for r in range(nrows):
            ax = axes[r, c]
            if single:
                bs, (key, name) = sets[0], _CHANNELS[c]
                show_phase = r == 1
            else:
                bs, (key, name) = sets[c], _CHANNELS[r]
                show_phase = phase

            b2, _ = grid.b2(bs.b2[key])
            sig = b2 > null.threshold
            if show_phase:
                ph = np.where(sig, grid.phase(bs.bispec[key]), np.nan)
                im_p = ax.pcolormesh(e, e, ph.T, cmap=pcmap, vmin=-np.pi,
                                     vmax=np.pi, rasterized=True)
            else:
                with np.errstate(divide="ignore", invalid="ignore"):
                    lb = np.log10(b2)
                im_b = ax.pcolormesh(e, e, lb.T, cmap=cmap, vmin=vmin, vmax=vmax,
                                     rasterized=True)
                if sig.any():
                    ax.contour(grid.centre, grid.centre, sig.T.astype(float),
                               levels=[0.5], colors="white", linewidths=0.7)
            _style_map_ax(ax, f_lim, left=(c == 0), bottom=(r == nrows - 1))
            _guides(ax, f_guides, f_lim)

            txt = ""
            if single and r == 0:
                txt = name
            elif not single:
                txt = name if r else f"{name}\n{col_labels[c]}"
            if txt:
                ax.text(0.04, 0.96, txt, transform=ax.transAxes, va="top",
                        ha="left", fontsize=7, color="white",
                        bbox=dict(fc="black", alpha=0.45, lw=0, pad=2))

    if im_b is not None:
        b_axes = axes[0, :] if single else axes.ravel()
        fig.colorbar(im_b, ax=b_axes.tolist(), label=r"$\log_{10} b^2$",
                     shrink=0.9, pad=0.02)
    if im_p is not None:
        p_axes = axes[1, :] if single else axes.ravel()
        cb = fig.colorbar(im_p, ax=p_axes.tolist(), shrink=0.9, pad=0.02,
                          label=r"Biphase (rad), $b^2 >$ 99% null")
        cb.set_ticks(_PI_TICKS[0])
        cb.set_ticklabels(_PI_TICKS[1])

    fig.savefig(savepath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {savepath}")


def plot_slices(sets, grid, null, savepath, f_lim, f_guides, f_ref, q_ref,
                legend_title=None):
    """Stacked panels on a shared log-frequency axis (qpo_sims layout):
    b^2(f, f) / b^2(f, f2 in band) / biphase(f, f2 in band).

    The f2 band is f_ref * (1 +/- 1/(2 q_ref)), i.e. the FWHM of the reference
    QPO; b^2 is averaged over it (and over each log f1 cell), and the biphase
    is that of the summed bispectrum.

    Single realisation: band 1, band 2 and cross in dark / light / black-dashed.
    Sweep: band 1 only, colour = sweep point (all channels are in the maps).
    The grey dash-dot curves are the 99% null thresholds.
    """
    fig, axes = qs._setup_xspec_axes(figsize=(7, 8.5), nrows=3)
    ax1, ax2, ax3 = axes
    single = len(sets) == 1
    f = grid.centre
    diag = np.arange(grid.nb)
    hw = 1.0 / (2.0 * q_ref)
    f_lo, f_hi = f_ref * (1 - hw), f_ref * (1 + hw)
    mask = grid.band_mask(f_lo, f_hi)
    thr_slice, _ = null.slice_threshold(mask)

    if single:
        channels = _CHANNELS
        styles = {"11": dict(c=qs._C_DARK, ls="-"), "22": dict(c=qs._C_LIGHT, ls="-"),
                  "112": dict(c="black", ls="--")}
        colours = [None]
    else:
        channels = _CHANNELS[:1]
        styles = {"11": dict(ls="-")}
        colours = plt.cm.viridis(np.linspace(0.15, 0.9, len(sets)))
    markers = {"11": "o", "22": "s", "112": "^"}

    for bs, col in zip(sets, colours):
        for key, name in channels:
            st = dict(styles[key])
            if col is not None:
                st["c"] = col
            lab = name if single else bs.label
            b2, _ = grid.b2(bs.b2[key])
            ax1.plot(f, b2[diag, diag], lw=1.3, alpha=0.9, label=lab,
                     drawstyle="steps-mid", **st)
            sl, _ = grid.slice_b2(bs.b2[key], mask)
            ph = grid.slice_phase(bs.bispec[key], mask)
            ax2.plot(f, sl, lw=1.3, alpha=0.9, drawstyle="steps-mid", **st)
            sig = sl > thr_slice
            ax3.plot(f[sig], ph[sig], marker=markers[key], ms=4, lw=0, alpha=0.9,
                     c=st["c"])

    thr_kw = dict(color="0.35", ls="-.", lw=1.0, drawstyle="steps-mid")
    ax1.plot(f, null.threshold[diag, diag], label="99% null", **thr_kw)
    ax2.plot(f, thr_slice, **thr_kw)
    for ax in (ax1, ax2):
        ax.set_yscale("log")
        ax.set_ylim(3e-3, 0.7)
    ax1.set_ylabel(r"$b^2(f,\,f)$")
    band_txt = rf"[{f_lo:.2f},\,{f_hi:.2f}]"
    ax2.set_ylabel(rf"$b^2(f,\,f_2 \in {band_txt}\,\mathrm{{Hz}})$")
    ax3.set_ylabel("Biphase (rad)")
    ax3.set_ylim(-np.pi * 1.08, np.pi * 1.08)
    ax3.set_yticks(_PI_TICKS[0])
    ax3.set_yticklabels(_PI_TICKS[1])
    ax3.axhline(0, color="gray", lw=0.5)
    ax2.axvspan(f_lo, f_hi, color=_C_GUIDE, alpha=0.25, lw=0)

    if single:
        ax1.legend(loc="upper right", fontsize=8)
    else:
        ax1.legend(title=f"Band 1, {legend_title}", fontsize=7, ncol=2,
                   loc="upper right", title_fontsize=8)

    for ax in axes:
        qs._add_vlines(ax, f_guides)
    axes[0].set_xlim(*f_lim)
    fig.align_ylabels(axes)
    fig.savefig(savepath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {savepath}")


# =========================================================================
# Driver
# =========================================================================

def _sweep_label(entry):
    if entry.get("label"):
        return entry["label"]
    return f'{entry["frac_qpo"] * 100:.2f}%'


def _subsample(n, k):
    """Indices of at most k roughly evenly spaced sweep points (keeps ends)."""
    if n <= k:
        return list(range(n))
    return sorted(set(np.round(np.linspace(0, n - 1, k)).astype(int).tolist()))


_NULL_CACHE = {}


def run(name, args):
    spec = PAPER_RUNS[name]
    print(f"[{name}]  (paper Fig. {spec.fig}; qpo_sims model '{spec.model}')")
    model_args = qs.build_parser().parse_args([spec.model] + spec.cli)
    res = qs.MODELS[spec.model](model_args)

    f_lim = tuple(v if v is not None else d for v, d in zip(res.f_lim, (0.1, 20.0)))
    f_lim = (max(f_lim[0], 0.5 / args.segment), min(f_lim[1], 0.5 / args.dt))
    f_guides = tuple(v for v in res.vlines if v is not None) or (spec.guides or ())

    if isinstance(res, qs.SweepResult):
        sets = []
        legend_title = res.legend_title
        while res.entries:
            e = res.entries.pop(0)
            sets.append(compute_bispec(e["lc1"], e["lc2"], args.dt, args.segment,
                                       args.window, label=_sweep_label(e)))
            del e
            gc.collect()
    else:
        sets = [compute_bispec(res.lc1, res.lc2, args.dt, args.segment,
                               args.window)]
        legend_title = None
    del res
    gc.collect()

    grid = LogGrid(sets[0].freq, args.rebin)
    key = (args.dt, args.segment, sets[0].m, args.rebin, args.window)
    if key not in _NULL_CACHE:
        _NULL_CACHE[key] = NullModel(args.dt, args.segment, sets[0].m, grid,
                                     window=args.window)
    null = _NULL_CACHE[key]
    print(f"  M = {sets[0].m}, null mean b^2 = {null.mu0:.4f}, "
          f"native 99% = {null.native_level:.4f}, kappa = {null.kappa:.2f}")

    in_lim = (grid.centre >= f_lim[0]) & (grid.centre <= f_lim[1])
    win = np.ix_(in_lim, in_lim)
    for bs in sets:
        for k, _ in _CHANNELS:
            b2, _ = grid.b2(bs.b2[k])
            exc = (b2 / null.threshold)[win]
            frac = np.nanmean(exc > 1)
            print(f"    {bs.label:>10s} {k:>4s}: max b^2/thr = {np.nanmax(exc):5.2f}, "
                  f"frac cells > 99% null = {frac:.3f}")

    base = os.path.join(args.outdir, name)
    if len(sets) > 1:
        keep = _subsample(len(sets), args.max_cols)
        plot_maps([sets[i] for i in keep], grid, null,
                  f"{base}_bicoherence_maps.{args.ext}", f_lim, f_guides,
                  col_labels=[f"{legend_title} = {sets[i].label}" for i in keep])
        plot_maps([sets[i] for i in keep], grid, null,
                  f"{base}_biphase_maps.{args.ext}", f_lim, f_guides,
                  col_labels=[f"{legend_title} = {sets[i].label}" for i in keep],
                  phase=True)
    else:
        plot_maps(sets, grid, null, f"{base}_bicoherence_maps.{args.ext}",
                  f_lim, f_guides)
    plot_slices(sets, grid, null, f"{base}_bicoherence_slices.{args.ext}", f_lim,
                f_guides, spec.f_ref, spec.q_ref, legend_title=legend_title)

    if args.save_npz:
        np.savez_compressed(
            f"{base}_bispec.npz", freq=sets[0].freq, m=sets[0].m,
            native_null_99=null.native_level,
            labels=np.array([bs.label for bs in sets]),
            **{f"b2_{k}": np.stack([bs.b2[k] for bs in sets]) for k, _ in _CHANNELS},
            **{f"bispec_{k}": np.stack([bs.bispec[k] for bs in sets])
               for k, _ in _CHANNELS})
    del sets
    gc.collect()


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("runs", nargs="+", choices=list(PAPER_RUNS) + ["all"])
    p.add_argument("--outdir", default="figures/bispec")
    p.add_argument("--ext", default="png", choices=["png", "pdf"])
    p.add_argument("--dt", type=_parse_dt, default=1 / 64,
                   help="Time resolution the light curves are rebinned to before "
                        "the bispectrum (default 1/64 s, Nyquist 32 Hz). Accepts "
                        "fractions, e.g. 1/128.")
    p.add_argument("--segment", type=float, default=16.0,
                   help="Segment length in s (default 16, as in the paper).")
    p.add_argument("--window", default="hann", choices=["hann", "none"],
                   help="Per-segment taper (default hann; see _segments).")
    p.add_argument("--rebin", type=float, default=0.05,
                   help="Fractional width of the log frequency cells used for "
                        "plotting (default 0.05; 0 = native resolution).")
    p.add_argument("--max-cols", type=int, default=6, dest="max_cols",
                   help="Max sweep points shown in the 2D map figure (all points "
                        "are always shown in the slice figure).")
    p.add_argument("--save-npz", action="store_true", dest="save_npz",
                   help="Also save native-resolution b^2 and complex bispectra "
                        "to {name}_bispec.npz.")
    return p


def main():
    args = build_parser().parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    names = list(PAPER_RUNS) if "all" in args.runs else args.runs
    for name in names:
        run(name, args)


if __name__ == "__main__":
    main()

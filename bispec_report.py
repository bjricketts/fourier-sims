"""
bispec_report.py
----------------
Second-order Fourier products and per-model reports, shared by bispec_sims.py
and bispec_biphase.py.

For each case (one realisation, or one point of a sweep) `summarise_case`
keeps everything needed later, so the full light curves can be freed:

    * the first `t_snip` seconds of both light curves;
    * power spectra of both bands, the cross spectrum, phase lag and coherence
      (stingray AveragedPowerspectrum / AveragedCrossspectrum, 'frac' norm,
      untapered segments of the bispectrum segment length, as in qpo_sims);
    * dynamical power spectra of both bands (stingray DynamicalPowerspectrum),
      averaged to at most `n_time` time bins and rebinned in frequency.

The cross spectrum is AveragedCrossspectrum(band 1, band 2), the order used by
qpo_sims.make_spectra.

Outputs
    {name}_fourier.{ext}   PSDs / cross spectrum (Re, Im) / phase lag / coherence,
                           all cases overlaid (stacked qpo_sims layout).
    {name}_report.pdf      Page 1: the Fourier overview of all cases, and the
                           model-specific summary figure if one is passed.
                           Then one page per case: light curves, dynamical
                           power spectra, Fourier products, and the band-1 and
                           cross (1 1 2*) bicoherence and biphase maps.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

from stingray import AveragedCrossspectrum, AveragedPowerspectrum, DynamicalPowerspectrum

import qpo_sims as qs
import bispec_sims as bsm


@dataclass
class CaseSummary:
    label: str
    t: np.ndarray
    c1: np.ndarray
    c2: np.ndarray
    fourier: dict
    dyn: dict
    intrinsic: bool


# =========================================================================
# Computation
# =========================================================================

def _rebin_spec(spec, scale, segment):
    if scale == "log":
        return spec.rebin_log(f=0.03)
    return spec.rebin(df=2.0 / segment)


def summarise_case(lc1, lc2, segment, label, intrinsic, scale="log", t_snip=10.0,
                   n_time=64, f_max=None):
    """Light-curve snippet, second-order Fourier products and dynamical PSDs.

    intrinsic: use the intrinsic coherence (Poisson noise subtracted) rather
    than the raw coherence; only meaningful for counting-noise light curves.
    """
    kw = dict(segment_size=segment, norm="frac", silent=True)
    ps1 = _rebin_spec(AveragedPowerspectrum(lc1, **kw), scale, segment)
    ps2 = _rebin_spec(AveragedPowerspectrum(lc2, **kw), scale, segment)
    cs = _rebin_spec(AveragedCrossspectrum(lc1, lc2, **kw), scale, segment)
    lag, lag_e = cs.phase_lag()
    coh = cs.intrinsic_coherence()[0] if intrinsic else cs.coherence()[0]
    fourier = dict(f_ps=ps1.freq, ps1=ps1.power, ps1_e=ps1.power_err, ps2=ps2.power,
                   ps2_e=ps2.power_err, f_cs=cs.freq, cs=cs.power, lag=lag,
                   lag_e=lag_e, coh=np.real(coh))

    dyn = {}
    for b, lc in ((1, lc1), (2, lc2)):
        d = DynamicalPowerspectrum(lc, segment_size=segment, norm="frac")
        P, tt = d.dyn_ps, d.time            # (n_freq, n_time)
        k = int(np.ceil(P.shape[1] / n_time))
        if k > 1:
            m = (P.shape[1] // k) * k
            P = P[:, :m].reshape(P.shape[0], -1, k).mean(2)
            tt = tt[:m].reshape(-1, k).mean(1)
        freq = d.freq
        if f_max is not None:
            keep = freq <= f_max
            freq, P = freq[keep], P[keep]
        grid = bsm.LogGrid(freq, 0.05 if scale == "log" else 0.0)
        tot = np.zeros((grid.nb, P.shape[1]))
        cnt = np.zeros(grid.nb)
        np.add.at(tot, grid.idx, np.real(P))
        np.add.at(cnt, grid.idx, 1)
        Pb = tot / cnt[:, None]
        dt_bin = tt[1] - tt[0] if tt.size > 1 else segment
        dyn[b] = dict(t_edges=np.concatenate([tt - dt_bin / 2, [tt[-1] + dt_bin / 2]]),
                      f_edges=grid.edges, f=grid.centre, fP=grid.centre[:, None] * Pb)

    n = int(round(t_snip / lc1.dt))
    return CaseSummary(label, lc1.time[:n] - lc1.time[0], np.asarray(lc1.counts[:n]),
                       np.asarray(lc2.counts[:n]), fourier, dyn, intrinsic)


# =========================================================================
# Plotting
# =========================================================================

def draw_fourier(axes, sums, f_lim, f_guides, legend_title=None, scale="log"):
    """PSDs / cross spectrum / phase lag / coherence into four stacked axes.

    One case: band 1 dark, band 2 light; Re orange, Im blue; lag and coherence
    black. Several: colour = case; band 2 and Im dotted/dashed.
    """
    ax1, ax2, ax3, ax4 = axes
    single = len(sums) == 1
    colours = [None] if single else plt.cm.viridis(np.linspace(0.15, 0.9, len(sums)))
    re_all, im_all = [], []
    for s, col in zip(sums, colours):
        F = s.fourier
        if single:
            ax1.errorbar(F["f_ps"], F["f_ps"] * F["ps1"], F["f_ps"] * F["ps1_e"],
                         c=qs._C_DARK, lw=1.5, label="band 1")
            ax1.errorbar(F["f_ps"], F["f_ps"] * F["ps2"], F["f_ps"] * F["ps2_e"],
                         c=qs._C_LIGHT, lw=1.5, label="band 2")
            ax2.plot(F["f_cs"], F["cs"].real, c=qs._C_REAL, lw=1.5, label="Real")
            ax2.plot(F["f_cs"], F["cs"].imag, c=qs._C_IMAG, lw=1.5, label="Imaginary")
            ax3.errorbar(F["f_cs"], F["lag"], F["lag_e"], c="black", lw=1.5)
            ax4.plot(F["f_cs"], F["coh"], c="black", lw=1.5, drawstyle="steps-mid")
        else:
            ax1.plot(F["f_ps"], F["f_ps"] * F["ps1"], c=col, lw=1.5, label=s.label)
            ax1.plot(F["f_ps"], F["f_ps"] * F["ps2"], c=col, lw=1.2, ls="--", alpha=0.7)
            ax2.plot(F["f_cs"], F["cs"].real, c=col, lw=1.5)
            ax2.plot(F["f_cs"], F["cs"].imag, c=col, lw=1.2, ls=":")
            ax3.errorbar(F["f_cs"], F["lag"], F["lag_e"], c=col, lw=1.2, alpha=0.85)
            ax4.plot(F["f_cs"], F["coh"], c=col, lw=1.2, drawstyle="steps-mid")
        re_all.append(F["cs"].real)
        im_all.append(F["cs"].imag)

    ax1.set_yscale("log")
    ax1.set_ylabel(r"$\nu P_\nu$ [(rms/mean)$^2$]")
    ax2.axhline(0, color="gray", lw=0.5)
    ax2.set_yscale("symlog", linthresh=qs._symlog_threshold(re_all + im_all))
    ax2.set_ylabel("Cross spectrum")
    ax3.axhline(0, color="gray", lw=0.5)
    ax3.set_ylim(-np.pi, np.pi)
    ax3.set_yticks([-np.pi / 2, 0, np.pi / 2])
    ax3.set_yticklabels([r"$-\pi/2$", "0", r"$\pi/2$"])
    ax3.set_ylabel("Phase lag (rad)")
    ax4.set_ylim(0, 1.05)
    ax4.set_ylabel("Intrinsic coherence" if sums[0].intrinsic else "Coherence")
    if single:
        ax1.legend(loc="lower left", fontsize=8)
        ax2.legend(loc="lower left", fontsize=8)
    else:
        ax1.legend(title=legend_title, fontsize=7, ncol=2, loc="lower left",
                   title_fontsize=8)
        ax2.legend(handles=[
            Line2D([0], [0], color="gray", label="Real / band 1"),
            Line2D([0], [0], color="gray", ls=":", label="Imaginary"),
            Line2D([0], [0], color="gray", ls="--", label="band 2 (PSD)"),
        ], fontsize=7, loc="lower left")
    for ax in axes:
        ax.set_xscale(scale)
        qs._add_vlines(ax, f_guides)
    axes[0].set_xlim(0.0 if scale == "linear" else f_lim[0], f_lim[1])
    axes[-1].set_xlabel("Frequency (Hz)")


def plot_fourier(sums, savepath, f_lim, f_guides, legend_title=None, scale="log"):
    fig, axes = qs._setup_xspec_axes(figsize=(7, 10), nrows=4)
    draw_fourier(axes, sums, f_lim, f_guides, legend_title, scale)
    fig.align_ylabels(axes)
    fig.savefig(savepath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {savepath}")


def _draw_dyn(ax, dyn, norm, f_lim, f_guides, scale, text):
    im = ax.pcolormesh(dyn["t_edges"], dyn["f_edges"], dyn["fP"], norm=norm,
                       cmap="magma", rasterized=True)
    ax.set_yscale(scale)
    ax.set_ylim(0.0 if scale == "linear" else f_lim[0], f_lim[1])
    for f in f_guides:
        ax.axhline(f, color="#b2abd2", ls="--", lw=1.0)
    ax.set_xlabel("Time (s)")
    ax.text(0.02, 0.95, text, transform=ax.transAxes, va="top", ha="left", fontsize=8,
            color="white", bbox=dict(fc="black", alpha=0.45, lw=0, pad=2))
    return im


def _case_page(s, bs, grid, null, f_lim, map_lim, f_guides, scale, header):
    fig = plt.figure(figsize=(15, 19))
    outer = fig.add_gridspec(3, 2, height_ratios=[0.7, 1.0, 2.6], hspace=0.28,
                             wspace=0.22)

    # Light curves
    ax = fig.add_subplot(outer[0, :])
    ax.plot(s.t, s.c1, c=qs._C_DARK, lw=1.0, label="band 1")
    ax.plot(s.t, s.c2, c=qs._C_LIGHT, lw=1.0, label="band 2")
    ax.set_xlim(s.t[0], s.t[-1])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Rate")
    ax.legend(loc="upper right", fontsize=8)
    ax.text(0.01, 0.95, header, transform=ax.transAxes, va="top", ha="left",
            fontsize=10, bbox=dict(fc="white", alpha=0.8, lw=0, pad=3))

    # Dynamical power spectra, shared colour scale
    vals = np.concatenate([s.dyn[b]["fP"].ravel() for b in (1, 2)])
    vals = vals[np.isfinite(vals) & (vals > 0)]
    norm = LogNorm(*np.percentile(vals, [2, 99.5]))
    for col, b in enumerate((1, 2)):
        ax = fig.add_subplot(outer[1, col])
        im = _draw_dyn(ax, s.dyn[b], norm, f_lim, f_guides, scale, f"band {b}")
        if col == 0:
            ax.set_ylabel("Frequency (Hz)")
    fig.colorbar(im, ax=fig.axes[-2:], label=r"$\nu P_\nu$ [(rms/mean)$^2$]",
                 shrink=0.9, pad=0.02)

    # Fourier products
    inner = outer[2, 0].subgridspec(4, 1, hspace=0)
    fa = [fig.add_subplot(inner[i]) for i in range(4)]
    for a in fa[1:]:
        a.sharex(fa[0])
    for a in fa[:-1]:
        a.tick_params(labelbottom=False)
    draw_fourier(fa, [s], f_lim, f_guides, scale=scale)
    fig.align_ylabels(fa)

    # Bispectrum maps: band 1 and cross, b^2 and biphase
    if bs is not None:
        mlim = (0.0, map_lim[1]) if scale == "linear" else map_lim
        inner = outer[2, 1].subgridspec(2, 2, hspace=0.08, wspace=0.08)
        maps = []
        for r, (key, name) in enumerate((("11", "Band 1"), ("112", r"Cross ($1\,1\,2^{*}$)"))):
            for c, phase in enumerate((False, True)):
                ax = fig.add_subplot(inner[r, c])
                im = bsm.draw_map_panel(ax, bs, key, grid, null, mlim, f_guides,
                                        phase=phase, scale=scale, left=(c == 0),
                                        bottom=(r == 1), text=name)
                maps.append((ax, im, phase))
        fig.colorbar([m[1] for m in maps if not m[2]][0],
                     ax=[m[0] for m in maps if not m[2]], label=r"$\log_{10} b^2$",
                     orientation="horizontal", pad=0.06, shrink=0.8)
        bsm.phase_colorbar(fig, [m[1] for m in maps if m[2]][0],
                           [m[0] for m in maps if m[2]], orientation="horizontal",
                           pad=0.06, shrink=0.8)
    return fig


def write_report(path, name, sums, sets, grid, null, f_lim, map_lim, f_guides,
                 legend_title=None, scale="log", summary_png=None):
    """Multi-page PDF: overview, then one page per case.

    sets: BispecSets aligned with sums (or None to omit the maps).
    summary_png: optional path of a model-specific summary figure for page 1.
    """
    with PdfPages(path) as pdf:
        fig = plt.figure(figsize=(15, 11))
        gs = fig.add_gridspec(1, 2 if summary_png else 1, wspace=0.15)
        inner = gs[0, 0].subgridspec(4, 1, hspace=0)
        fa = [fig.add_subplot(inner[i]) for i in range(4)]
        for a in fa[1:]:
            a.sharex(fa[0])
        for a in fa[:-1]:
            a.tick_params(labelbottom=False)
        draw_fourier(fa, sums, f_lim, f_guides, legend_title, scale)
        fa[0].text(0.02, 1.02, f"{name}: all cases", transform=fa[0].transAxes,
                   fontsize=10)
        fig.align_ylabels(fa)
        if summary_png:
            ax = fig.add_subplot(gs[0, 1])
            ax.imshow(plt.imread(summary_png))
            ax.axis("off")
        pdf.savefig(fig, dpi=120)
        plt.close(fig)

        for i, s in enumerate(sums):
            header = f"{name}" + (f": {legend_title} = {s.label}" if legend_title
                                  else (f": {s.label}" if s.label else ""))
            bs = sets[i] if sets is not None else None
            fig = _case_page(s, bs, grid, null, f_lim, map_lim, f_guides, scale, header)
            pdf.savefig(fig, dpi=120)
            plt.close(fig)
    print(f"  -> {path}")

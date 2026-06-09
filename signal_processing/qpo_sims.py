#!/usr/bin/env python3
"""
qpo_sims.py
-----------
Run the simulated scenarios from the QPO Paper A (Simulations) introduction
and save lightcurve + cross-spectral diagnostic figures.

Each model below produces two output files:

    {model}_lightcurve.{ext}
    {model}_crossspectral.{ext}

Models (taken from the paper introduction):

  Dampened sine waves
    damping              - Two sine waves, identical freq, different damping factors
    coherent_sum         - Two coherent sine waves added together (dual-signal case)
    resonance            - Different resonant frequencies   [PLACEHOLDER]

  Oscillator interaction
    hidden               - "Hidden" oscillator in broadband (Konig+24-like) [PLACEHOLDER]
    multi_incoherent     - Multiple incoherent oscillators                  [PLACEHOLDER]
    multi_multiplicative - Multiplicative oscillators                       [PLACEHOLDER]

  Non-stationarity
    fm                   - Frequency-modulated QPO (time-varying omega)
    count_rate_fm        - GRS1915-style count-rate / frequency coupling
    envelope_qpo         - Non-stationary envelope multiplying a QPO         [PLACEHOLDER]

  Non-minimum phase transfer functions
    allpass              - All-pass filter
    echo                 - Direct + delayed echo
    coupled              - Coupled oscillators

Example
-------
    python qpo_sims.py damping --outdir figures/
    python qpo_sims.py damping --outdir figures/ --omega1 5 --omega2 5 \\
                               --dampen1 0.125 --dampen2 0.25 --segment 16
    python qpo_sims.py allpass --outdir figures/ --omega_z 1.0 --gamma_z 0.1
"""
from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve, lfilter

from stingray import Lightcurve, AveragedCrossspectrum, AveragedPowerspectrum
from stingray import simulator


# =========================================================================
# Core signal-generation utilities (extracted from fourier_sims.ipynb)
# =========================================================================

def make_decaying_sine(constant, factor, nu, Dnu, phi, t, phase=0):
    """Decaying sine wave: C + A*C * exp(-pi*Dnu*|t|) * sin(2 pi nu t - 2 pi phi + phase)."""
    alpha = np.pi * Dnu
    return (constant
            + factor * constant
              * np.exp(-alpha * np.abs(t))
              * np.sin(2 * np.pi * nu * t - 2 * np.pi * phi + phase))


def lorentz(freq, peak_f, q, rms):
    """Lorentzian model from ndspec (imported lazily so the rest of the script
    still runs if ndspec is not installed)."""
    import ndspec
    par_array = np.array([peak_f, q, rms])
    return ndspec.Models.lorentz(freq, par_array)


def lorentzian(f, f0, q, norm):
    """Closed-form zero-centred Lorentzian (no ndspec dependency)."""
    df = f0 / q
    return norm * (df / np.pi) / ((f - f0)**2 + df**2)


def dho_filter(x, dt, f0, zeta, gain=1.0):
    """Pass a real signal through a damped harmonic oscillator H(f).

    H(f) = w0^2 / (w0^2 - w^2 + 2j zeta w0 w),  w = 2*pi*f, w0 = 2*pi*f0.
    Q-factor is 1/(2*zeta).
    """
    N = len(x)
    freq = np.fft.rfftfreq(N, d=dt)
    w_ = 2 * np.pi * freq
    w0 = 2 * np.pi * f0
    H = w0**2 / (w0**2 - w_**2 + 2j * zeta * w0 * w_)
    X = np.fft.rfft(x - x.mean())
    return np.fft.irfft(gain * H * X, n=N)


def timmer_koenig_two_lorentz(t, mean, dt, omega1, omega2, q1, q2, rms1, rms2):
    """Generate two Timmer-Koenig realisations from a two-Lorentzian PSD."""
    rms = rms1 + rms2
    sim = simulator.Simulator(N=len(t), mean=mean, dt=dt, rms=rms, poisson=True)
    w = np.fft.rfftfreq(sim.N, d=sim.dt)[1:]
    l1 = lorentz(w, omega1, q1, rms1)
    l2 = lorentz(w, omega2, q2, rms2)
    psd = l1 + l2
    lc1 = sim.simulate(l1)
    lc2 = sim.simulate(l2)
    return lc1, lc2, w, l1, l2, psd


def timmer_koenig_single(t, mean, dt, omega, q, rms):
    """Single-Lorentzian Timmer-Koenig realisation."""
    sim = simulator.Simulator(N=len(t), mean=mean, dt=dt, rms=rms, poisson=True)
    w = np.fft.rfftfreq(sim.N, d=sim.dt)[1:]
    lor = lorentz(w, omega, q, rms)
    lc = sim.simulate(lor)
    return lc, w, lor


# ---------------- Non-minimum-phase transfer functions ----------------

def delay_tf(freqs, tau):
    omega = 2 * np.pi * freqs
    return np.exp(-1j * omega * tau)


def allpass_tf(freqs, omega_z, gamma_z):
    """All-pass: unit magnitude, frequency-dependent phase (UHP zero / LHP pole)."""
    omega = 2 * np.pi * freqs
    z = omega_z + 1j * gamma_z / 2
    z_conj = omega_z - 1j * gamma_z / 2
    return (omega - z_conj) / (omega - z)


def echo_tf(freqs, alpha, tau):
    omega = 2 * np.pi * freqs
    return 1 + alpha * np.exp(-1j * omega * tau)


def coupled_oscillators_tf(freqs, omega1, gamma1, omega2, gamma2, kappa,
                           zero_type='minimum'):
    omega = 2 * np.pi * freqs
    denom1 = omega1**2 - omega**2 + 1j * gamma1 * omega
    denom2 = omega2**2 - omega**2 + 1j * gamma2 * omega
    denom = denom1 * denom2 - kappa**2

    omega_z = (omega1 + omega2) / 2
    gamma_z = (gamma1 + gamma2) / 2
    sign = +1 if zero_type == 'minimum' else -1
    numer = omega_z**2 - omega**2 + sign * 1j * gamma_z * omega
    return numer / denom


def frequency_to_impulse_response(H_func, t_max, dt, **tf_kwargs):
    """Convert a frequency-domain transfer function to a real impulse response."""
    N = int(t_max / dt)
    if N % 2 == 1:
        N += 1
    f_fft = np.fft.fftfreq(N, d=dt)
    H_fft = np.zeros(N, dtype=complex)

    pos = f_fft > 0
    neg = f_fft < 0
    H_fft[pos] = H_func(f_fft[pos], **tf_kwargs)
    H_fft[neg] = np.conj(H_func(-f_fft[neg], **tf_kwargs))
    H_fft[0] = np.real(H_func(np.array([1e-10]), **tf_kwargs)[0])

    h = np.real(np.fft.ifft(H_fft)) / dt
    t = np.arange(N) * dt
    return t, h


def apply_transfer_function(input_signal, h, dt):
    return fftconvolve(input_signal, h, mode='full')[:len(input_signal)] * dt


# ---------------- Time-varying-omega damped oscillator (for FM) ----------------

def damped_oscillator_convolve(noise, omega, gamma, dt=1.0):
    """Drive a damped oscillator with time-varying omega (state-space step)."""
    noise = np.asarray(noise, dtype=float)
    N = noise.size
    omega = np.broadcast_to(np.asarray(omega, dtype=float), (N,))
    gamma = np.broadcast_to(np.asarray(gamma, dtype=float), (N,))

    x = np.zeros(N)
    v = np.zeros(N)
    for n in range(1, N):
        w, g = omega[n - 1], gamma[n - 1]
        wd_sq = w * w - g * g
        wd = np.sqrt(wd_sq) if wd_sq > 0 else 1e-12
        e = np.exp(-g * dt)
        c = np.cos(wd * dt)
        s = np.sin(wd * dt)
        a11 = e * (c + (g / wd) * s)
        a12 = e * s / wd
        a21 = -e * (w * w / wd) * s
        a22 = e * (c - (g / wd) * s)
        x[n] = a11 * x[n - 1] + a12 * v[n - 1]
        v[n] = a21 * x[n - 1] + a22 * v[n - 1] + noise[n - 1] * dt
    return x


def frequency_from_noise(noise, f_base, f_scale, smooth_tau=None, dt=1.0):
    """sqrt(counts) -> instantaneous frequency map (with optional 1-pole smoothing)."""
    counts = np.clip(np.asarray(noise, dtype=float), 0.0, None)
    if smooth_tau is not None and smooth_tau > 0:
        alpha = np.exp(-dt / smooth_tau)
        smoothed = np.empty_like(counts)
        smoothed[0] = counts[0]
        for n in range(1, counts.size):
            smoothed[n] = alpha * smoothed[n - 1] + (1 - alpha) * counts[n]
        counts = smoothed
    root = np.sqrt(counts)
    norm = np.sqrt(np.mean(counts)) + 1e-12
    return f_base + f_scale * (root / norm)


def ou_process(N, dt, mean, tau, sigma, rng=None):
    """Discrete-time Ornstein-Uhlenbeck process.

    Update:  x[n+1] = mean + alpha*(x[n] - mean) + eta * xi[n]
    with alpha = exp(-dt/tau), eta = sigma * sqrt(1 - alpha**2).
    Stationary distribution: N(mean, sigma**2). Power spectrum is a
    Lorentzian centred at 0 with HWHM 1/(2*pi*tau). Implemented via an
    AR(1) lfilter for speed at large N.
    """
    if rng is None:
        rng = np.random.default_rng()
    alpha = np.exp(-dt / tau)
    eta = sigma * np.sqrt(1.0 - alpha**2)
    y = lfilter([eta], [1.0, -alpha], rng.standard_normal(N))
    return mean + y


# =========================================================================
# Generic plotting helpers
# =========================================================================

@dataclass
class SimResult:
    """Container for whatever a model produces, passed to the plotters."""
    lc1: Lightcurve
    lc2: Lightcurve
    ps1: AveragedPowerspectrum
    ps2: AveragedPowerspectrum
    avgcs12: AveragedCrossspectrum
    vlines: tuple = ()            # vertical guide frequencies on the spectral plot
    t_lim: tuple = (0, 5)         # lightcurve x-range
    f_lim: tuple = (None, None)   # spectral x-range
    label1: str = r"$\phi_1$"
    label2: str = r"$\phi_2$"


@dataclass
class SweepResult:
    """Container for a rms-sweep model (cf. paper Figs. 4/5/6).

    `entries` is a list of dicts with keys frac_qpo, ps1, ps2, cs, coh, lag, lag_e.
    `lc1`/`lc2` are the highest-rms-iteration lightcurves, used by `plot_lightcurves`.
    """
    entries: list
    lc1: Lightcurve
    lc2: Lightcurve
    vlines: tuple = ()
    t_lim: tuple = (0, 20)
    f_lim: tuple = (0.1, 20)
    lag_ylim: tuple = (-50, 100)        # in ms
    coh_ylim: tuple = (0.5, 1.02)
    title: str = ""
    legend_title: str = "QPO frac rms"
    label1: str = "band 1"
    label2: str = "band 2"


def _add_vlines(ax, vlines):
    for v in vlines:
        if v is not None:
            ax.axvline(v, color="#b2abd2", linestyle="--", lw=0.8)


def plot_lightcurves(res: SimResult, savepath: str):
    """Plot the two lightcurves over `t_lim`."""
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(res.lc1.time, res.lc1.counts, c="#e66101", label=res.label1, lw=0.8)
    ax.plot(res.lc2.time, res.lc2.counts, c="#5e3c99", label=res.label2, lw=0.8,
            alpha=0.85)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.set_xlim(*res.t_lim)
    ax.legend()
    fig.tight_layout()
    fig.savefig(savepath, dpi=150)
    plt.close(fig)
    print(f"  -> {savepath}")


def plot_crossspectral(res: SimResult, savepath: str):
    """Power spectra, cross-spectrum (re/im), phase lag, coherence."""
    fig, axes = plt.subplots(1, 4, figsize=(20, 4.5))
    ax1, ax2, ax3, ax4 = axes

    # PSD
    ax1.errorbar(res.ps1.freq, res.ps1.freq * res.ps1.power,
                 yerr=res.ps1.freq * res.ps1.power_err,
                 c="#e66101", label=f"{res.label1} PSD", lw=0.8)
    ax1.errorbar(res.ps2.freq, res.ps2.freq * res.ps2.power,
                 yerr=res.ps2.freq * res.ps2.power_err,
                 c="#5e3c99", label=f"{res.label2} PSD", lw=0.8)
    _add_vlines(ax1, res.vlines)
    ax1.set(xscale="log", yscale="log",
            xlabel="Frequency (Hz)", ylabel=r"Frequency * Power (rms$^2$)")
    ax1.set_xlim(*res.f_lim)
    ax1.legend(loc="lower left", fontsize=9)

    # Real / Imag cross-spectrum
    ax2.plot(res.avgcs12.freq, np.real(res.avgcs12.power),
             c="#e66101", label="Real", lw=0.8)
    ax2.plot(res.avgcs12.freq, np.imag(res.avgcs12.power),
             c="#5e3c99", label="Imaginary", lw=0.8)
    _add_vlines(ax2, res.vlines)
    ax2.set(xscale="log", xlabel="Frequency (Hz)", ylabel="Cross-spectrum")
    ax2.set_yscale("symlog", linthresh=1e-4)
    ax2.set_xlim(*res.f_lim)
    ax2.legend(fontsize=9)

    # Phase lag
    phase, phase_err = res.avgcs12.phase_lag()
    ax3.errorbar(res.avgcs12.freq, phase, yerr=res.avgcs12.lag_err,
                 c="#fdb863", lw=0.8)
    _add_vlines(ax3, res.vlines)
    ax3.set(xscale="log", xlabel="Frequency (Hz)", ylabel="Phase (radians)")
    ax3.set_xlim(*res.f_lim)
    ax3.set_ylim(-np.pi, np.pi)

    # Coherence
    coh = res.avgcs12.coherence()[0]
    ax4.plot(res.avgcs12.freq, coh, c="#fdb863", lw=0.8)
    _add_vlines(ax4, res.vlines)
    ax4.set(xscale="log", xlabel="Frequency (Hz)", ylabel="Coherence")
    ax4.set_xlim(*res.f_lim)
    ax4.set_ylim(0, 1.01)

    fig.tight_layout()
    fig.savefig(savepath, dpi=150)
    plt.close(fig)
    print(f"  -> {savepath}")


def plot_crossspectral_sweep(res: SweepResult, savepath: str):
    """3-panel rms-sweep figure (PSD / time lag / coherence), as in paper Figs 4-6."""
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    cmap = plt.cm.viridis(np.linspace(0.15, 0.9, len(res.entries)))
    for ax in axes:
        _add_vlines(ax, res.vlines)

    for r, c in zip(res.entries, cmap):
        lbl = f'{r["frac_qpo"] * 100:.2f}%'
        axes[0].loglog(r['ps1'].freq, r['ps1'].freq * r['ps1'].power,
                       color=c, alpha=0.85, label=lbl)
        axes[0].loglog(r['ps2'].freq, r['ps2'].freq * r['ps2'].power,
                       color=c, alpha=0.5, ls='--')
        axes[1].errorbar(r['cs'].freq, r['lag'] * 1e3, yerr=r['lag_e'] * 1e3,
                         color=c, alpha=0.8, ms=3, lw=1)
        axes[2].semilogx(r['cs'].freq, r['coh'], color=c, alpha=0.85,
                         drawstyle='steps-mid')

    axes[0].set(ylabel=r'PSD $\times \nu$ [(rms/mean)$^2$]', title=res.title)
    axes[0].legend(title=res.legend_title, fontsize=8, ncol=2)
    axes[1].axhline(0.0, ls=':', color='gray', alpha=0.5)
    axes[1].set(ylabel='Time lag [ms]', ylim=res.lag_ylim)
    axes[2].set(xscale='log', xlabel='Frequency [Hz]', ylabel='Coherence',
                ylim=res.coh_ylim, xlim=res.f_lim)

    fig.tight_layout()
    fig.savefig(savepath, dpi=150)
    plt.close(fig)
    print(f"  -> {savepath}")


def make_spectra(lc1, lc2, segment_size):
    ps1 = AveragedPowerspectrum.from_lightcurve(lc1, norm="frac",
                                                segment_size=segment_size)
    ps2 = AveragedPowerspectrum.from_lightcurve(lc2, norm="frac",
                                                segment_size=segment_size)
    cs  = AveragedCrossspectrum.from_lightcurve(lc1, lc2, norm="frac",
                                                segment_size=segment_size)
    return ps1, ps2, cs


# =========================================================================
# Models
# =========================================================================

def _opt(value, default):
    """Return `value` if not None, else `default`. Lets each model define
    its own paper-faithful defaults while CLI args still act as overrides."""
    return value if value is not None else default


def _dho_band(driver, t, dt, dho_specs, mean, frac, rng):
    """Build a single light curve: mean + frac * sum(DHO responses) on a
    shared driver, Poisson-sampled.

    `dho_specs` is a list of (f0, zeta) pairs; their responses are summed
    before injection so multiple peaks per band share the driver coherently.
    The amplitude factor `frac` is applied relative to `driver.std()`, so
    different-ζ DHOs on the same driver keep their relative amplitudes.
    """
    resp = sum(dho_filter(driver, dt, f0=f, zeta=z, gain=1.0)
               for f, z in dho_specs)
    amp = frac * mean / driver.std()
    rate = np.clip(mean + amp * resp, 0, None)
    return Lightcurve(t, rng.poisson(rate), dt=dt, skip_checks=True)


# ---- 1. Differing damping ------------------------------------------------
def model_damping(args) -> SimResult:
    """Paper Fig. 1: one broadband driver into two DHOs at the *same*
    resonant frequency but with different damping ratios."""
    dt = 1 / 512
    T_bins = int(1024 / dt)
    driver, t = _build_broadband(dt, T_bins, seed=1)

    f1 = _opt(args.omega1, 5.0)
    f2 = _opt(args.omega2, 5.0)
    z1 = _opt(args.dampen1, 0.05)     # Q = 10
    z2 = _opt(args.dampen2, 0.15)     # Q ~ 3.3
    mean = _opt(args.constant, 1e4)

    rng = np.random.default_rng(11)
    lc1 = _dho_band(driver, t, dt, [(f1, z1)], mean, frac=0.05, rng=rng)
    lc2 = _dho_band(driver, t, dt, [(f2, z2)], mean, frac=0.05, rng=rng)

    ps1, ps2, cs = make_spectra(lc1, lc2, args.segment)
    return SimResult(lc1, lc2, ps1, ps2, cs,
                     vlines=(f1, f2), t_lim=(0, 5), f_lim=(0.5, 20))


# ---- 2. Two coherent peaks per band (dual-signal case) -------------------
def model_coherent_sum(args) -> SimResult:
    """Paper Fig. 2: each band is the sum of two DHO responses on a shared
    driver. The low-frequency response is identical between bands; the
    high-frequency one sits at slightly different center frequencies."""
    dt = 1 / 512
    T_bins = int(1024 / dt)
    driver, t = _build_broadband(dt, T_bins, seed=1)

    f_shared = _opt(args.omega1, 1.0)
    f_sec_a  = _opt(args.omega2, 3.0)    # band 1's high-freq peak
    f_sec_b  = _opt(args.omega3, 2.7)    # band 2's high-freq peak
    z_shared = _opt(args.dampen1, 0.1)
    z_sec_a  = _opt(args.dampen2, 0.1)
    z_sec_b  = _opt(args.dampen3, 0.15)
    mean = _opt(args.constant, 1e4)

    rng = np.random.default_rng(11)
    lc1 = _dho_band(driver, t, dt,
                    [(f_shared, z_shared), (f_sec_a, z_sec_a)],
                    mean, frac=0.05, rng=rng)
    lc2 = _dho_band(driver, t, dt,
                    [(f_shared, z_shared), (f_sec_b, z_sec_b)],
                    mean, frac=0.05, rng=rng)

    ps1, ps2, cs = make_spectra(lc1, lc2, args.segment)
    return SimResult(lc1, lc2, ps1, ps2, cs,
                     vlines=(f_shared, f_sec_a, f_sec_b),
                     t_lim=(0, 5), f_lim=(0.1, 12))


# ---- 3. Different resonant frequencies -----------------------------------
def model_resonance(args) -> SimResult:
    """Paper Section 2.1.2 (no figure in the draft): one broadband driver
    into two DHOs at different resonant frequencies, identical damping."""
    dt = 1 / 512
    T_bins = int(1024 / dt)
    driver, t = _build_broadband(dt, T_bins, seed=1)

    f1 = _opt(args.omega1, 5.0)
    f2 = _opt(args.omega2, 6.0)
    z = _opt(args.dampen1, 0.05)
    mean = _opt(args.constant, 1e4)

    rng = np.random.default_rng(11)
    lc1 = _dho_band(driver, t, dt, [(f1, z)], mean, frac=0.05, rng=rng)
    lc2 = _dho_band(driver, t, dt, [(f2, z)], mean, frac=0.05, rng=rng)

    ps1, ps2, cs = make_spectra(lc1, lc2, args.segment)
    return SimResult(lc1, lc2, ps1, ps2, cs,
                     vlines=(f1, f2),
                     t_lim=(0, 5), f_lim=(0.5, 20))


# ---- 4. "Hidden" QPO in broadband (paper Fig. 4) -------------------------

def _build_broadband(dt, T_bins, rms=0.5, mean=100, seed=None):
    """Two-Lorentzian broadband Timmer-Koenig realisation. Returns (broadband, t).

    `broadband` is mean-subtracted.
    """
    kwargs = {} if seed is None else {"random_state": seed}
    sim = simulator.Simulator(N=int(T_bins), mean=mean, dt=dt, rms=rms,
                              poisson=True, **kwargs)
    w = np.fft.rfftfreq(sim.N, d=sim.dt)[1:]
    psd_model = lorentzian(w, 0.4, 0.3, 1.0) + lorentzian(w, 2.0, 0.5, 1.2)
    lc_bb = sim.simulate(psd_model)
    bb = lc_bb.counts - lc_bb.counts.mean()
    t = np.arange(lc_bb.n) * dt
    return bb, t


def _sweep_qpo_rms(broadband_pri, broadband_sec, qpo_pri_unit, qpo_sec_unit,
                   t, dt, segment, frac_qpo_vals, frac_qpo_sec_vals,
                   mean1=100.0, mean2=40.0, frac_bb_1=0.25, frac_bb_2=0.12,
                   rebin=None, seed=11):
    """Generic two-band, two-QPO rms sweep -> list of result dicts.

    For the single-QPO case, pass `qpo_sec_unit=None` and frac_qpo_sec_vals=0.
    Independent vs shared broadband driver is controlled by which broadband
    arrays the caller passes for `broadband_pri` / `broadband_sec`.
    """
    bb1 = broadband_pri * (frac_bb_1 * mean1 / broadband_pri.std())
    bb2 = broadband_sec * (frac_bb_2 * mean2 / broadband_sec.std())
    rng = np.random.default_rng(seed)

    out = []
    for frac_qpo, frac_qpo_sec in zip(frac_qpo_vals, frac_qpo_sec_vals):
        qpo_pri_1 = qpo_pri_unit * (frac_qpo * mean1)
        qpo_pri_2 = qpo_pri_unit * (frac_qpo * mean2)
        if qpo_sec_unit is not None:
            qpo_sec_1 = qpo_sec_unit * (frac_qpo_sec * mean1)
            qpo_sec_2 = qpo_sec_unit * (frac_qpo_sec * mean2)
        else:
            qpo_sec_1 = qpo_sec_2 = 0.0

        rate1 = np.clip(mean1 + bb1 + qpo_pri_1 + qpo_sec_1, 0, None)
        rate2 = np.clip(mean2 + bb2 + qpo_pri_2 + qpo_sec_2, 0, None)
        lc1 = Lightcurve(t, rng.poisson(rate1), dt=dt, skip_checks=True)
        lc2 = Lightcurve(t, rng.poisson(rate2), dt=dt, skip_checks=True)

        ps1 = AveragedPowerspectrum(lc1, segment_size=segment, norm='frac')
        ps2 = AveragedPowerspectrum(lc2, segment_size=segment, norm='frac')
        cs = AveragedCrossspectrum(lc2, lc1, segment_size=segment, norm='frac')
        if rebin is not None:
            ps1.rebin_log(f=rebin)
            ps2.rebin_log(f=rebin)
            cs.rebin_log(f=rebin)

        try:
            coh, _ = cs.intrinsic_coherence()
        except Exception:
            coh, _ = cs.coherence()
        lag, lag_e = cs.time_lag()

        out.append(dict(frac_qpo=frac_qpo, ps1=ps1, ps2=ps2, cs=cs,
                        coh=coh, lag=lag, lag_e=lag_e, lc1=lc1, lc2=lc2))
    return out


def model_hidden(args) -> SweepResult:
    """Single 'hidden' QPO sweep -- paper Fig. 4."""
    dt = 1 / 512
    T_bins = int(1024 / dt)
    broadband, t = _build_broadband(dt, T_bins, seed=args.seed)

    f0_qpo, zeta = 1.5, 0.05      # Q = 10
    test_resp = dho_filter(broadband, dt, f0=f0_qpo, zeta=zeta, gain=1.0)
    qpo_unit = test_resp / test_resp.std()

    frac_vals = np.geomspace(1e-3, 5e-2, 6)
    entries = _sweep_qpo_rms(
        broadband, broadband, qpo_unit, None,
        t, dt, segment=16.0,
        frac_qpo_vals=frac_vals, frac_qpo_sec_vals=np.zeros_like(frac_vals),
        seed=11)

    last = entries[-1]
    return SweepResult(
        entries=entries, lc1=last['lc1'], lc2=last['lc2'],
        vlines=(f0_qpo,), t_lim=(0, 20), f_lim=(0.1, 20),
        lag_ylim=(-50, 100), coh_ylim=(0.5, 1.02),
        title=f"Hidden-QPO rms sweep (Q={1/(2*zeta):.1f}, "
              "solid: band 1, dashed: band 2)")


# ---- 5. Multiple incoherent oscillators (paper Fig. 6a) ------------------

def _setup_two_qpo_drivers(dt, T_bins, shared_driver: bool, seed_a=1, seed_b=2):
    """Return (broadband_obs, broadband_pri, broadband_sec, t).

    `broadband_obs` is the broadband that goes into the observed lightcurves
    of both bands (always realisation A). `broadband_pri` / `broadband_sec`
    are the drivers of the primary / secondary QPO -- the same array for the
    'coherent' case, independent realisations for the 'incoherent' case.
    """
    bb_a, t = _build_broadband(dt, T_bins, seed=seed_a)
    if shared_driver:
        return bb_a, bb_a, bb_a, t
    bb_b, _ = _build_broadband(dt, T_bins, seed=seed_b)
    return bb_a, bb_a, bb_b, t


def _model_multi_qpo(args, shared_driver: bool) -> SweepResult:
    dt = 1 / 512
    T_bins = int(1024 / dt)
    bb_obs, bb_pri, bb_sec, t = _setup_two_qpo_drivers(
        dt, T_bins, shared_driver=shared_driver)

    f0_pri, zeta_pri = 1.5, 0.05    # Q = 10
    f0_sec, zeta_sec = 1.6, 0.05    # Q = 10

    resp_pri = dho_filter(bb_pri, dt, f0=f0_pri, zeta=zeta_pri, gain=1.0)
    resp_sec = dho_filter(bb_sec, dt, f0=f0_sec, zeta=zeta_sec, gain=1.0)
    qpo_pri_unit = resp_pri / resp_pri.std()
    qpo_sec_unit = resp_sec / resp_sec.std()

    frac_vals = np.geomspace(1e-3, 1e-1, 11)
    entries = _sweep_qpo_rms(
        bb_obs, bb_obs, qpo_pri_unit, qpo_sec_unit,
        t, dt, segment=16.0,
        frac_qpo_vals=frac_vals, frac_qpo_sec_vals=frac_vals * 0.5,
        mean1=100.0, mean2=100.0, frac_bb_1=0.12, frac_bb_2=0.06,
        rebin=0.05, seed=11)

    last = entries[-1]
    kind = "Coherent" if shared_driver else "Incoherent"
    return SweepResult(
        entries=entries, lc1=last['lc1'], lc2=last['lc2'],
        vlines=(f0_pri, f0_sec), t_lim=(0, 20), f_lim=(0.1, 20),
        lag_ylim=(-120, 120), coh_ylim=(0.0, 1.02),
        title=f"{kind} secondary QPO  "
              rf"($Q_1$={1/(2*zeta_pri):.1f}, $Q_2$={1/(2*zeta_sec):.1f})",
        legend_title="Primary QPO frac rms")


def model_multi_incoherent(args) -> SweepResult:
    return _model_multi_qpo(args, shared_driver=False)


def model_multi_coherent(args) -> SweepResult:
    """Multiple oscillators sharing a common driver -- paper Fig. 6b."""
    return _model_multi_qpo(args, shared_driver=True)


# ---- 6. Multiplicative oscillators (paper §3.3) --------------------------
#
# Two-stage signal:
#   1) broadband B drives a first oscillator -> O1 (additive QPO, as in §3.2)
#   2) a second oscillator O2 *multiplies* either O1 alone or the whole (B+O1)
#
# mult_mode = 'osc1_only':  observed = B + O1 * (1 + a * O2_unit)
# mult_mode = 'combined':   observed = (B + O1) * (1 + a * O2_unit)
#
# Both O1 and O2 are DHO responses driven by the same broadband realisation
# (matches the §3.2 "coherent secondary" setup with multiplicativity added).
# To make O2 driven by O1, or by an independent broadband, swap the `o2_drv`
# line below.

def model_multi_multiplicative(args) -> SweepResult:
    dt = 1 / 512
    T_bins = int(1024 / dt)
    bb, t = _build_broadband(dt, T_bins, seed=args.seed)

    f1, zeta1 = 1.5, 0.05      # primary  (Q = 10)
    f2, zeta2 = 1.6, 0.05      # secondary (Q = 10)

    # Both oscillators driven by the shared broadband. Swap `bb` for the
    # output of the first dho_filter call (or for an independent broadband
    # via _build_broadband(..., seed=other)) to change the driver.
    o1 = dho_filter(bb, dt, f0=f1, zeta=zeta1, gain=1.0)
    o2_drv = bb
    o2 = dho_filter(o2_drv, dt, f0=f2, zeta=zeta2, gain=1.0)
    o1_unit = o1 / o1.std()
    o2_unit = o2 / o2.std()

    # Fixed broadband and primary-QPO contributions; sweep the secondary
    # modulation depth `a`.
    mean1, mean2 = 100.0, 40.0
    frac_bb_1, frac_bb_2 = 0.25, 0.12
    frac_o1 = 0.05             # primary at a comfortably-detectable 5%
    bb1 = bb * (frac_bb_1 * mean1 / bb.std())
    bb2 = bb * (frac_bb_2 * mean2 / bb.std())
    o1_1 = o1_unit * (frac_o1 * mean1)
    o1_2 = o1_unit * (frac_o1 * mean2)

    a_vals = np.geomspace(1e-2, 5e-1, 6)
    rng = np.random.default_rng(11)
    mode = getattr(args, "mult_mode", "osc1_only")

    entries = []
    for a in a_vals:
        mod = 1.0 + a * o2_unit
        if mode == "osc1_only":
            rate1 = mean1 + bb1 + o1_1 * mod
            rate2 = mean2 + bb2 + o1_2 * mod
        elif mode == "combined":
            rate1 = mean1 + (bb1 + o1_1) * mod
            rate2 = mean2 + (bb2 + o1_2) * mod
        else:
            raise ValueError(f"unknown --mult_mode {mode!r}")

        lc1 = Lightcurve(t, rng.poisson(np.clip(rate1, 0, None)),
                         dt=dt, skip_checks=True)
        lc2 = Lightcurve(t, rng.poisson(np.clip(rate2, 0, None)),
                         dt=dt, skip_checks=True)

        ps1 = AveragedPowerspectrum(lc1, segment_size=16.0, norm='frac')
        ps2 = AveragedPowerspectrum(lc2, segment_size=16.0, norm='frac')
        cs = AveragedCrossspectrum(lc2, lc1, segment_size=16.0, norm='frac')
        ps1.rebin_log(f=0.05)
        ps2.rebin_log(f=0.05)
        cs.rebin_log(f=0.05)
        try:
            coh, _ = cs.intrinsic_coherence()
        except Exception:
            coh, _ = cs.coherence()
        lag, lag_e = cs.time_lag()

        # `frac_qpo` here is the modulation depth a -- the sweep variable,
        # which the plotter renders as a percentage in the legend.
        entries.append(dict(frac_qpo=a, ps1=ps1, ps2=ps2, cs=cs,
                            coh=coh, lag=lag, lag_e=lag_e,
                            lc1=lc1, lc2=lc2))

    last = entries[-1]
    return SweepResult(
        entries=entries, lc1=last['lc1'], lc2=last['lc2'],
        vlines=(f1, f2), t_lim=(0, 20), f_lim=(0.1, 20),
        lag_ylim=(-120, 120), coh_ylim=(0.0, 1.02),
        title=(rf"Multiplicative secondary QPO, mode={mode}  "
               rf"($Q_1$={1/(2*zeta1):.1f}, $Q_2$={1/(2*zeta2):.1f})"),
        legend_title="O2 mod. depth $a$")


# ---- 7. Frequency-modulated QPO ------------------------------------------
def model_fm(args) -> SimResult:
    """QPO from frequency wandering rather than damping.

    The instantaneous frequency f(t) is an Ornstein-Uhlenbeck process with
    mean f0, std `fm_sigma` and correlation time `fm_tau`. The QPO itself
    is sin(integral of 2*pi*f(t)dt), so it is undamped and unit-amplitude;
    the width of the resulting PSD peak comes purely from f(t) wandering.

    Two-band setup: each band sees the same shared wander but can also
    differ in three independent, optional ways (all default to identity,
    so the two bands are identical out of the box):
      * fm_sigma_band/fm_tau_band - independent OU perturbations to f(t).
        Even small values strongly decohere because phase accumulates the
        frequency-difference integral over the segment length.
      * fm_sigma_phase/fm_tau_phase - independent OU perturbations applied
        directly to the phase before sin(). Linear, intuitive knob in rad:
        ~0.1 mild dip, ~pi total decoherence; independent of segment size.
      * fm_df, fm_wander_scale - band 2 sees (f0 + fm_df) as its central
        frequency and a wander amplitude scaled by fm_wander_scale.

    Band 2 also lags band 1 by a constant `fm_delay`.
    """
    dt = 1 / 512
    T = 1024
    N = int(T / dt)
    t = np.arange(N) * dt
    rng = np.random.default_rng(11)

    f0 = _opt(args.omega1, 1.5)               # central QPO freq [Hz]
    sigma_f = _opt(args.fm_sigma, 0.15)       # std of common f-wander [Hz]
    tau_f = _opt(args.fm_tau, 2.0)            # correlation time of common wander [s]
    sigma_band = _opt(args.fm_sigma_band, 0.0)   # std of per-band f-perturbation [Hz]
    tau_band = _opt(args.fm_tau_band, tau_f)     # correlation time of perturbation [s]
    sigma_phase = _opt(args.fm_sigma_phase, 0.0) # std of per-band phase noise [rad]
    tau_phase = _opt(args.fm_tau_phase, tau_f)   # correlation time of phase noise [s]
    df_offset = _opt(args.fm_df, 0.0)         # offset added to band-2 central f [Hz]
    wander_scale = _opt(args.fm_wander_scale, 1.0)  # band-2 wander amplitude scaling
    delay_s = _opt(args.fm_delay, 5e-3)       # constant band-2 delay [s]

    # Shared wander (deviation from f0) + optional per-band frequency
    # perturbations (sigma_band) and phase noise (sigma_phase). Band 2 can
    # also have an offset central frequency and a rescaled wander amplitude.
    # All extra knobs default to the identity, so bands stay identical.
    df_common = ou_process(N, dt, mean=0.0, tau=tau_f, sigma=sigma_f, rng=rng)
    if sigma_band > 0:
        rng_b1 = np.random.default_rng(21)
        rng_b2 = np.random.default_rng(22)
        df1 = ou_process(N, dt, mean=0.0, tau=tau_band, sigma=sigma_band, rng=rng_b1)
        df2 = ou_process(N, dt, mean=0.0, tau=tau_band, sigma=sigma_band, rng=rng_b2)
    else:
        df1 = df2 = 0.0
    if sigma_phase > 0:
        rng_p1 = np.random.default_rng(31)
        rng_p2 = np.random.default_rng(32)
        dphi1 = ou_process(N, dt, mean=0.0, tau=tau_phase, sigma=sigma_phase, rng=rng_p1)
        dphi2 = ou_process(N, dt, mean=0.0, tau=tau_phase, sigma=sigma_phase, rng=rng_p2)
    else:
        dphi1 = dphi2 = 0.0

    f1_t = f0 + df_common + df1
    f2_t = f0 + df_offset + wander_scale * df_common + df2
    phase1 = 2 * np.pi * np.cumsum(f1_t) * dt + dphi1
    phase2 = 2 * np.pi * np.cumsum(f2_t) * dt + dphi2
    qpo1 = np.sin(phase1)
    qpo2 = np.sin(phase2)

    # Broadband floor (shared between bands)
    bb, _ = _build_broadband(dt, N, seed=1)

    # Inject into two bands; band 2 lags band 1 by `delay_s`
    mean1, mean2 = 100.0, 40.0
    frac_bb_1, frac_bb_2 = 0.25, 0.12
    frac_qpo = 0.10
    bb1 = bb * (frac_bb_1 * mean1 / bb.std())
    bb2 = bb * (frac_bb_2 * mean2 / bb.std())

    delay_n = int(round(delay_s / dt))
    qpo2_delayed = np.roll(qpo2, delay_n)

    rate1 = np.clip(mean1 + bb1 + frac_qpo * mean1 * qpo1,         0, None)
    rate2 = np.clip(mean2 + bb2 + frac_qpo * mean2 * qpo2_delayed, 0, None)
    lc1 = Lightcurve(t, rng.poisson(rate1), dt=dt, skip_checks=True)
    lc2 = Lightcurve(t, rng.poisson(rate2), dt=dt, skip_checks=True)

    ps1, ps2, cs = make_spectra(lc1, lc2, args.segment)
    return SimResult(lc1, lc2, ps1, ps2, cs,
                     vlines=(f0,),
                     t_lim=(0, 20), f_lim=(0.1, 20))


# ---- 8. Count-rate / frequency coupling (GRS1915-style) ------------------
def model_count_rate_fm(args) -> SimResult:
    dt = 0.01
    N = 400000
    t = np.arange(N) * dt
    sim = simulator.Simulator(N=N, mean=100000, dt=dt, rms=0.25, poisson=False)
    drive = sim.simulate(2)
    f_t = frequency_from_noise(drive.counts, f_base=1.0, f_scale=3.0,
                               smooth_tau=0.5, dt=dt)
    omega_t = 2 * np.pi * f_t
    x1 = damped_oscillator_convolve(drive.counts, omega_t, args.gamma1, dt=dt)
    x2 = damped_oscillator_convolve(drive.counts, omega_t, args.gamma2, dt=dt)
    C = 100
    lc1 = Lightcurve(t, x1 * C * 0.1 + C * 10)
    lc2 = Lightcurve(t, x2 * C * 0.1 + C * 10)
    ps1, ps2, cs = make_spectra(lc1, lc2, args.segment)
    return SimResult(lc1, lc2, ps1, ps2, cs,
                     vlines=(1.0, 4.0),
                     t_lim=(0, 20), f_lim=(1 / args.segment, 20))


# ---- 9. Non-stationary envelope * QPO (PLACEHOLDER) ----------------------
def model_envelope_qpo(args) -> SimResult:
    print("[PLACEHOLDER] model_envelope_qpo: GRB-like envelope multiplying a "
          "QPO is not implemented yet. Stub returns a simple Gaussian envelope.")
    dt = args.dt
    T = args.segment * 32
    t = np.arange(0, T, dt)
    env = np.exp(-((t - T / 2) ** 2) / (2 * (T / 8) ** 2))
    base, _, _ = timmer_koenig_single(t, mean=500, dt=dt,
                                      omega=_opt(args.omega1, 1.5), q=10, rms=0.5)
    sig = base.counts * env
    lc1 = Lightcurve(t, sig + 100)
    lc2 = Lightcurve(t, sig * 0.9 + 100)
    ps1, ps2, cs = make_spectra(lc1, lc2, args.segment)
    return SimResult(lc1, lc2, ps1, ps2, cs,
                     vlines=(_opt(args.omega1, 1.5),),
                     t_lim=(0, T), f_lim=(1 / args.segment, 20))

# =========================================================================
# CLI
# =========================================================================

MODELS: dict[str, Callable[[argparse.Namespace], SimResult]] = {
    "damping":              model_damping,
    "coherent_sum":         model_coherent_sum,
    "resonance":            model_resonance,
    "hidden":               model_hidden,
    "multi_incoherent":     model_multi_incoherent,
    "multi_coherent":       model_multi_coherent,
    "multi_multiplicative": model_multi_multiplicative,
    "fm":                   model_fm,
    "count_rate_fm":        model_count_rate_fm,
    "envelope_qpo":         model_envelope_qpo,
}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("model", choices=list(MODELS.keys()) + ["all"],
                   help="Which scenario to run ('all' = run every model).")
    p.add_argument("--outdir", default="figures",
                   help="Directory to write the two output figures into.")
    p.add_argument("--ext", default="png", choices=["png", "pdf"],
                   help="Output file extension.")
    p.add_argument("--seed", type=int, default=None,
                   help="RNG seed for reproducibility.")

    # Generic time-domain options
    p.add_argument("--dt", type=float, default=1 / 512,
                   help="Sample spacing in seconds (default 1/512).")
    p.add_argument("--segment", type=int, default=16,
                   help="Segment size for averaging (s, default 16).")
    # Frequencies / damping (used by several models; each model has its own
    # paper-faithful defaults that are used when these are not passed).
    p.add_argument("--omega1", type=float, default=None)
    p.add_argument("--omega2", type=float, default=None)
    p.add_argument("--omega3", type=float, default=None)
    p.add_argument("--dampen1", type=float, default=None,
                   help="DHO damping ratio (zeta); Q = 1/(2*zeta).")
    p.add_argument("--dampen2", type=float, default=None)
    p.add_argument("--dampen3", type=float, default=None)
    p.add_argument("--constant", type=float, default=None,
                   help="Mean count rate for sine-wave / DHO models.")
    p.add_argument("--gamma1", type=float, default=0.3,
                   help="Damping rate (used by count_rate_fm / coupled-oscillator models).")
    p.add_argument("--gamma2", type=float, default=0.5)

    # FM-model knobs
    p.add_argument("--fm_sigma", type=float, default=None,
                   help="fm: std of frequency wandering [Hz] (default 0.15).")
    p.add_argument("--fm_tau", type=float, default=None,
                   help="fm: correlation time of frequency wandering [s] "
                        "(default 2.0).")
    p.add_argument("--fm_sigma_band", type=float, default=None,
                   help="fm: std of per-band independent frequency "
                        "perturbation [Hz] (default 0 = identical bands).")
    p.add_argument("--fm_tau_band", type=float, default=None,
                   help="fm: correlation time of per-band perturbation [s] "
                        "(default same as fm_tau).")
    p.add_argument("--fm_sigma_phase", type=float, default=None,
                   help="fm: std of per-band phase noise [rad] applied directly "
                        "before sin() (default 0). Linear decoherence knob: "
                        "~0.1 mild, ~pi total.")
    p.add_argument("--fm_tau_phase", type=float, default=None,
                   help="fm: correlation time of per-band phase noise [s] "
                        "(default same as fm_tau).")
    p.add_argument("--fm_df", type=float, default=None,
                   help="fm: offset added to band-2 central frequency [Hz] "
                        "(default 0 = same central freq).")
    p.add_argument("--fm_wander_scale", type=float, default=None,
                   help="fm: scaling of band-2 wander amplitude relative to "
                        "band 1 (default 1.0 = same wander).")
    p.add_argument("--fm_delay", type=float, default=None,
                   help="fm: constant time delay between bands [s] (default 5e-3).")

    # Non-min phase TF options
    p.add_argument("--omega_z", type=float, default=1.0)
    p.add_argument("--gamma_z", type=float, default=0.1)
    p.add_argument("--alpha", type=float, default=1.0,
                   help="Echo amplitude.")
    p.add_argument("--tau", type=float, default=0.1,
                   help="Echo delay (s).")
    p.add_argument("--kappa", type=float, default=0.1,
                   help="Coupling strength for coupled oscillators.")
    p.add_argument("--zero_type", default="minimum",
                   choices=["minimum", "nonminimum"])
    p.add_argument("--mult_mode", default="osc1_only",
                   choices=["osc1_only", "combined"],
                   help="multi_multiplicative: O2 modulates O1 alone "
                        "('osc1_only') or the full B+O1 signal ('combined').")

    return p


def run_one(name: str, args) -> None:
    print(f"[{name}]")
    res = MODELS[name](args)
    os.makedirs(args.outdir, exist_ok=True)
    base = os.path.join(args.outdir, name)
    plot_lightcurves(res, f"{base}_lightcurve.{args.ext}")
    if isinstance(res, SweepResult):
        plot_crossspectral_sweep(res, f"{base}_crossspectral.{args.ext}")
    else:
        plot_crossspectral(res, f"{base}_crossspectral.{args.ext}")


def main():
    args = build_parser().parse_args()
    if args.seed is not None:
        np.random.seed(args.seed)

    targets = list(MODELS) if args.model == "all" else [args.model]
    for name in targets:
        run_one(name, args)


if __name__ == "__main__":
    main()
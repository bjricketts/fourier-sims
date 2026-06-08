"""
PROPFLUC — Propagating Fluctuations model.

Python port of the Fortran 90 code by Adam Ingram, which implements the
analytic propagating-fluctuations power-/lag-spectrum model described in:

    Ingram & Done (2011, MNRAS, 415, 2323)
    Ingram & Done (2012, MNRAS, 419, 2369)
    Ingram & van der Klis (2013, MNRAS, 434, 1476)

The model treats the inner hot flow of an accreting black hole as a series
of concentric rings, each generating multiplicative mass-accretion-rate
fluctuations at its local viscous frequency. Fluctuations propagate inward,
so each ring's variability is the convolution (in linear space) of the
fluctuations from all rings outside it. Hard and soft X-ray emissivity
profiles let us synthesise band-limited power spectra and the
cross-spectrum between bands; a Lorentzian QPO can be added or multiplied
on top.

The original Fortran exposes a subroutine `propfluc` that XSPEC calls to
fit data; the program/`geobins` wrapper at the top is just for standalone
testing. The translation keeps function boundaries and naming close to the
original where the original names are descriptive, and renames a handful of
single-letter Fortran variables when a longer name is materially clearer.

The non-trivial pieces are the FFT-convolution routine `fftconv` and the
ring-loop in `dopropfluc`, which together implement equations 20–22 of
Ingram & van der Klis (2013).
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Parameter index reference (matches the 1-based Fortran indexing in comments)
# ---------------------------------------------------------------------------
# param[ 0] / par(1)  Sigma0   — surface-density normalisation
# param[ 1] / par(2)  rbw      — bending-wave radius (Rg)
# param[ 2] / par(3)  kappa    — surface-density profile shape
# param[ 3] / par(4)  lambda   — surface-density profile shape
# param[ 4] / par(5)  zeta     — surface-density profile shape
# param[ 5] / par(6)  Fvar     — fractional variability per radial decade
# param[ 6] / par(7)  ro       — outer (truncation) radius (Rg)
# param[ 7] / par(8)  ri       — inner radius of the hot flow (Rg)
# param[ 8] / par(9)  Q        — quality factor of the QPO harmonics
# param[ 9] / par(10) Qsub     — quality factor of the sub-harmonic
# param[10] / par(11) n_qpo    — normalisation of fundamental
# param[11] / par(12) n_2qpo   — normalisation of 2nd harmonic
# param[12] / par(13) n_3qpo   — normalisation of 3rd harmonic
# param[13] / par(14) n_05qpo  — normalisation of sub-harmonic
# param[14] / par(15) gamma_h  — hard-band emissivity index
# param[15] / par(16) gamma_s  — soft-band emissivity index
# param[16] / par(17) M        — black-hole mass (Msun)
# param[17] / par(18) a        — dimensionless spin
# param[18] / par(19) Ndec     — number of rings per radial decade
# param[19] / par(20) mode     — 1: hard PSD, 2: soft PSD, 3: lag spectrum
# param[20] / par(21) conv     — 0: add QPO, 1: multiplicatively convolve QPO

# Internal hard cap on the working FFT grid size, mirroring the Fortran nmax.
NMAX_DEFAULT = 2**20


# ---------------------------------------------------------------------------
# Geometric frequency bins
# ---------------------------------------------------------------------------
def geobins(dt: float, n: int, c: float) -> np.ndarray:
    """
    Build a set of geometrically-spaced frequency bin edges.

    Mirrors the Fortran ``geobins`` subroutine. The first bin edge sits at
    ``df/2`` (so that the DC element is centred on f=0 in the natural FFT
    convention) and successive bin widths grow by factor ``c``. The walk
    stops as soon as adding the next bin would exceed the Nyquist index.

    Parameters
    ----------
    dt : float
        Time step of the underlying FFT grid (seconds).
    n : int
        Number of points in the underlying FFT grid (power of 2).
    c : float
        Geometric growth factor between successive bins (e.g. 1.05).

    Returns
    -------
    far : np.ndarray, shape (nf+1,)
        Bin edges in Hz. ``far[0]`` is the lower edge of bin 1; the centre
        and width of bin i are ``(far[i]+far[i-1])/2`` and
        ``far[i]-far[i-1]`` respectively.
    """
    df = 1.0 / (n * dt)

    # Start the walk: far[0] is the low edge of the first output bin.
    far_list = [df / 2.0]
    n_per_bin_float = 1.0   # the running "np" in the Fortran
    n_accumulated = 0       # cumulative number of fine FFT bins consumed

    while True:
        # Geometric growth of the (real-valued) bin width in fine-bin units.
        n_per_bin_float = c * n_per_bin_float
        n_accumulated += int(n_per_bin_float)

        if n_accumulated > n // 2:
            # Final bin — clip so the last edge sits at Nyquist.
            n_accumulated = n // 2
            far_list.append((n_accumulated + 0.5) * df)
            break

        far_list.append((n_accumulated + 0.5) * df)

    return np.asarray(far_list, dtype=np.float64)


# ---------------------------------------------------------------------------
# Main entry point — XSPEC-style wrapper
# ---------------------------------------------------------------------------
def propfluc(
    far: np.ndarray,
    param: np.ndarray | list,
    nmax: int = NMAX_DEFAULT,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Top-level wrapper that XSPEC would call.

    Picks an internal FFT grid that adequately samples the requested output
    frequency range, runs the model on that working grid, and rebins the
    result onto the user-supplied bin edges ``far``.

    Parameters
    ----------
    far : np.ndarray, shape (nf+1,)
        Output frequency bin edges in Hz (e.g. from :func:`geobins`).
    param : array-like, length 21
        Model parameters; see the header of this file for the index map.
    nmax : int, optional
        Cap on the internal FFT length. Defaults to 2**20 to match Fortran.

    Returns
    -------
    photar : np.ndarray, shape (nf,)
        XSPEC-style "photon rate" per bin, i.e. ``P_bin * df_bin`` for the
        PSD modes and the (unscaled) lag for mode 3.
    photer : np.ndarray, shape (nf,)
        Statistical errors. Always zero — this is a deterministic model.
    """
    param = np.asarray(param, dtype=np.float64)
    far = np.asarray(far, dtype=np.float64)
    nf = far.size - 1

    # Choose the working FFT grid: dt resolves the highest output frequency,
    # and n is large enough to cover at least the lowest output frequency.
    fmax = far[-1]
    fmin = far[0]
    dt = 1.0 / (2.0 * fmax)
    n = 2 ** int(np.ceil(-np.log(fmin * dt) / np.log(2.0)))
    if n > 2**15:
        n = 2**15

    # Run the core model on the working grid.
    P = dopropfluc(param, n, nmax, dt)

    # Bin / interpolate onto the requested output grid.
    photar = np.zeros(nf, dtype=np.float64)
    photer = np.zeros(nf, dtype=np.float64)
    for i in range(nf):
        f_centre = 0.5 * (far[i + 1] + far[i])
        df_bin = far[i + 1] - far[i]
        Pin = modbin(dt, n, P, f_centre, df_bin)
        photar[i] = Pin * df_bin

    return photar, photer


# ---------------------------------------------------------------------------
# The core calculation
# ---------------------------------------------------------------------------
def dopropfluc(
    param: np.ndarray,
    n: int,
    nmax: int,
    dt: float,
) -> np.ndarray:
    """
    Compute the model PSD or lag spectrum on the working FFT grid.

    Implements the analytic ring-by-ring sum in the propagating-fluctuations
    formalism. For each ring k from outermost to innermost, we:
      1. Compute the local viscous frequency ``fv`` and the hard/soft
         emissivity weights ``h_k`` and ``s_k``.
      2. Build the input Lorentzian PSD ``p_k`` of mass-accretion-rate
         fluctuations injected at this ring.
      3. Convolve ``p_k`` (in frequency space) onto the running propagated
         spectrum ``cp`` — mathematically this is the multiplicative
         combination of fluctuations in the time domain.
      4. Add the diagonal contribution ``h_k**2 * mf[k]`` (and similarly
         for soft and cross terms) to the band-PSDs and cross-spectrum.
      5. Add the off-diagonal cross-pair terms between ring k and every
         outer ring l < k, which carry the propagation lag through their
         phase factor exp(i * 2π f τ_{l→k}).

    See equations 20–22 of Ingram & van der Klis (2013).

    Returns
    -------
    pout : np.ndarray, shape (n//2,)
        Either the band-limited PSD (modes 1, 2) or the time-lag spectrum
        (mode 3), defined on FFT frequency bins j = 1..n//2.
    """
    # Radial grid setup
    m, a_ratio, dronr = rgrid(param)

    df = 1.0 / (n * dt)

    # Mean and standard deviation of the per-ring fluctuations.
    mu = 1.0
    sigma = param[5] / np.sqrt(param[18])      # Fvar / sqrt(Ndec)

    # Per-ring storage. All length n//2+1 so index 0 holds the DC element.
    half = n // 2 + 1

    # Hard- and soft-band PSDs and the real/imag parts of the cross-spectrum.
    ph = np.zeros(half)
    ps = np.zeros(half)
    rc = np.zeros(half)
    ic = np.zeros(half)

    # Per-ring weights and propagated PSDs.
    h = np.zeros(m + 1)              # hard emissivity, h[0] unused
    s = np.zeros(m + 1)              # soft emissivity, s[0] unused
    il = np.zeros(m + 1, dtype=int)  # propagation lag in fine FFT bins
    muk = np.zeros(m + 1)            # sqrt(p_k(0))

    # mf[k] is the running propagated PSD up to and including ring k.
    # pk[k] is the locally-injected PSD at ring k alone.
    mf = np.zeros((m + 1, half))
    pk = np.zeros((m + 1, half))

    # Running totals used for the QPO frequency.
    Td = 0.0  # total LT torque on the disc
    Jd = 0.0  # total angular momentum of the precessing region

    # Running propagated PSD: cp at iteration k holds mf[k].
    cp = np.zeros(half)

    # Loop from outermost (k=1) to innermost (k=m) ring.
    for k in range(1, m + 1):
        # Physical properties of this ring.
        fv, hard, soft, dTd, dJd = ring(k, param, a_ratio, dronr)
        Td += dTd
        Jd += dJd
        h[k] = hard
        s[k] = soft

        # Propagation lag from ring k+1 to ring k, expressed in *fine FFT
        # bins* of the working grid. The Fortran stores this with index l+1
        # in the cross-term phase below.
        il[k] = int(dronr / (fv * dt))

        # Lorentzian PSD of fluctuations injected at this ring (zero-centred).
        p = mkpwr(n, dt, mu, sigma, fv, 0.0)

        # Multiplicative propagation: each new ring's PSD is convolved with
        # the running propagated PSD. Done in the time domain via FFTs.
        if k == 1:
            cp = p.copy()
        else:
            prev = cp.copy()
            cp = fftconv(dt, n, p, prev)

        muk[k] = np.sqrt(p[0])
        pk[k] = p
        mf[k] = cp

        # Diagonal contributions (l == k terms) to band PSDs and cross.
        ph += h[k] ** 2 * mf[k]
        ps += s[k] ** 2 * mf[k]
        rc += h[k] * s[k] * mf[k]

        # Off-diagonal cross-pair contributions (Ingram & van der Klis 2013,
        # eqs. 20-22): for every outer ring l < k, the propagated PSD up to
        # the *outer* ring is mf[l], and the time delay between l and k
        # accumulates the il[l+1] sample shifts.
        phase = 0.0
        for i in range(1, k):
            l = k - i
            # Accumulate the propagation phase increment for ring l+1.
            phase += 2.0 * np.pi * il[l + 1] / n

            # Vectorised over frequency bin j = 0..n//2.
            j = np.arange(half)
            x = np.mod(phase * j, 2.0 * np.pi)
            cosx = np.cos(x)
            sinx = np.sin(x)

            ph += 2.0 * h[l] * h[k] * cosx * mf[l]
            ps += 2.0 * s[l] * s[k] * cosx * mf[l]
            rc += (h[k] * s[l] + h[l] * s[k]) * cosx * mf[l]
            ic += (h[k] * s[l] - h[l] * s[k]) * sinx * mf[l]

    # QPO centroid follows from the global precession frequency Td/Jd.
    fqpo = Td / Jd
    mode = int(param[19])
    conv = int(param[20])

    if mode == 1:
        ptot = incQPO(dt, param, conv, fqpo, n, ph)
        pout = ptot[1 : n // 2 + 1]
    elif mode == 2:
        ptot = incQPO(dt, param, conv, fqpo, n, ps)
        pout = ptot[1 : n // 2 + 1]
    elif mode == 3:
        tlag = lincQPO(dt, param, conv, fqpo, n, rc, ic)
        pout = tlag[1 : n // 2 + 1]
    else:
        raise ValueError(f"mode must be 1, 2 or 3 (got {mode})")

    return pout


# ---------------------------------------------------------------------------
# Physics: per-ring properties
# ---------------------------------------------------------------------------
def rgrid(par: np.ndarray) -> tuple[int, float, float]:
    """
    Set up the logarithmic radial grid.

    Returns
    -------
    m : int
        Number of rings.
    a_ratio : float
        Logarithmic spacing factor between successive rings (r_{k+1}/r_k).
    dronr : float
        Constant value of dr/r used by every ring.
    """
    ro = par[6]
    ri = par[7]
    Ndec = par[18]
    m = int(Ndec * (np.log10(ro) - np.log10(ri)))
    a_ratio = np.exp(np.log(ri / ro) / m)
    dronr = 2.0 * (1.0 - a_ratio) / (1.0 + a_ratio)
    return m, a_ratio, dronr


def ring(
    k: int,
    par: np.ndarray,
    a_ratio: float,
    dronr: float,
) -> tuple[float, float, float, float, float]:
    """
    Compute physical properties of ring ``k``.

    Returns
    -------
    fv : float
        Viscous frequency at the ring (Hz).
    h, s : float
        Hard- and soft-band emissivity weights for the ring.
    dTd, dJd : float
        Contributions of this ring to the total Lense-Thirring torque and
        angular momentum, used to set the QPO centroid frequency.
    """
    pi = np.pi
    c_speed = 3e8

    Sigma0 = par[0]
    rbw = par[1]
    kappa = par[2]
    lam = par[3]
    zeta = par[4]
    ro = par[6]
    gammaH = par[14]
    gammaS = par[15]
    Rg = par[16] * 1474.8           # gravitational radius in metres
    spin = par[17]

    # Ring radius (geometric mean of edges in log space).
    r = ro * (a_ratio ** (k - 1) + a_ratio ** k) / 2.0

    # Surface density profile (Ingram & Done 2012, eq. 6).
    x = r / rbw
    Sigma = x ** lam / ((1.0 + x ** kappa) ** ((zeta + lam) / kappa))
    Sigma *= Sigma0

    # Viscous frequency from the surface-density profile.
    fv = (1.0 + x ** kappa) ** ((lam + zeta) / kappa) / (x ** (lam + 2.0))
    fv = fv / (2.0 * pi * Sigma0 * rbw ** 2.0)
    fv = fv * c_speed / Rg

    # Emissivity-weighted contributions. Power-law emissivity in radius,
    # weighted by surface density and the dimensionless ring width.
    h = dronr * r ** (2.0 - gammaH) * Sigma
    s = dronr * r ** (2.0 - gammaS) * Sigma

    # Keplerian and Lense-Thirring frequencies at this ring.
    fk = r ** (-1.5) / (2.0 * pi) * c_speed / Rg
    flt = 1.0 - np.sqrt(1.0 - 4.0 * spin * r ** (-1.5) + 3.0 * spin ** 2 * r ** (-2.0))
    flt *= fk

    # Local LT torque and angular-momentum contributions.
    dTd = flt * fk * Sigma * r ** 4.0 * dronr
    dJd = fk * Sigma * r ** 4.0 * dronr

    return fv, h, s, dTd, dJd


# ---------------------------------------------------------------------------
# QPO inclusion
# ---------------------------------------------------------------------------
def incQPO(
    dt: float,
    param: np.ndarray,
    conv: int,
    fqpo: float,
    n: int,
    psum: np.ndarray,
) -> np.ndarray:
    """
    Add or convolve the QPO into the broadband PSD ``psum``.

    Normalises ``psum`` by its DC value before combining, then either
    adds a Lorentzian QPO sum (``conv==0``) or multiplicatively convolves
    it (``conv==1``).
    """
    psum = psum.copy()
    half = n // 2 + 1

    # Normalise so DC == 1; the absolute scale is set by the QPO normalisations.
    psum[1:half] = psum[1:half] / psum[0]
    psum[0] = 1.0

    if conv == 0:
        mu = 0.0
        pqpo = QPOpower(dt, param, mu, fqpo, n)
        ptot = psum + pqpo
    else:
        mu = 1.0
        pqpo = QPOpower(dt, param, mu, fqpo, n)
        ptot = fftconv(dt, n, psum, pqpo)

    return ptot


def lincQPO(
    dt: float,
    par: np.ndarray,
    conv: int,
    fqpo: float,
    n: int,
    rc: np.ndarray,
    ic: np.ndarray,
) -> np.ndarray:
    """
    Combine the QPO with the cross-spectrum and convert it to a phase lag
    spectrum (mode 3).

    Normalises the cross-spectrum by ``|C(0)|`` first, then either adds the
    QPO into the real part (``conv==0``) or convolves both real and
    imaginary parts with the QPO PSD (``conv==1``). Finally the phase lag
    is converted to a time lag via ``tlag = arg(C) / (2π f)``, with the
    correct branch choice for negative real parts.
    """
    rc = rc.copy()
    ic = ic.copy()

    half = n // 2 + 1
    df = 1.0 / (n * dt)

    cov0 = rc[0] ** 2 + ic[0] ** 2
    rc /= np.sqrt(cov0)
    ic /= np.sqrt(cov0)

    tlag = np.zeros(half)

    if conv == 0:
        mu = 0.0
        pqpo = QPOpower(dt, par, mu, fqpo, n)
        rc[1:half] = rc[1:half] + pqpo[1:half]
    else:
        mu = 1.0
        pqpo = QPOpower(dt, par, mu, fqpo, n)
        rcqpo = fftconv(dt, n, rc, pqpo)
        icqpo = fftconv(dt, n, ic, pqpo)
        rc[1:half] = rcqpo[1:half]
        ic[1:half] = icqpo[1:half]

    # Branch-corrected phase, then convert to time lag in seconds.
    j = np.arange(1, half)
    f = j * df
    phase = np.arctan(ic[1:half] / rc[1:half])
    # arctan returns the principal value in (-pi/2, pi/2); fix the branches.
    mask_pos = (ic[1:half] > 0.0) & (rc[1:half] < 0.0)
    mask_neg = (ic[1:half] < 0.0) & (rc[1:half] < 0.0)
    phase = np.where(mask_pos, phase + np.pi, phase)
    phase = np.where(mask_neg, phase - np.pi, phase)
    tlag[1:half] = phase / (2.0 * np.pi * f)

    return tlag


def QPOpower(
    dt: float,
    param: np.ndarray,
    mu: float,
    fqpo: float,
    n: int,
) -> np.ndarray:
    """
    Build the four-Lorentzian QPO PSD: fundamental, 2nd, 3rd harmonics and
    the half-frequency sub-harmonic (which uses its own quality factor).

    Each Lorentzian has FWHM = ``f_centre / (2 * Q)``.
    """
    half = n // 2 + 1
    Pqpo = np.zeros(half)
    Pqpo[0] = mu ** 2

    Q = param[8]            # Q for harmonics 1, 2, 3
    Qsub = param[9]         # Q for the sub-harmonic

    # Loop over k = 1, 2, 3, sub. The Fortran uses a single loop with
    # k=4 redirected to fac=0.5 and Q->Qsub; we keep the same shape.
    for k in range(1, 5):
        if k != 4:
            fac = float(k)
            Q_use = Q
        else:
            fac = 0.5
            Q_use = Qsub

        fcen = fac * fqpo
        width = fcen / (2.0 * Q_use)
        Norm = param[10 + (k - 1)]    # param(11..14) in 1-based indexing
        P = mkpwr(n, dt, mu, Norm, width, fcen)
        Pqpo[1:half] += P[1:half]

    return Pqpo


# ---------------------------------------------------------------------------
# Power-spectrum primitive
# ---------------------------------------------------------------------------
def mkpwr(
    np_grid: int,
    dt: float,
    mu: float,
    r: float,
    Delta: float,
    f0: float,
) -> np.ndarray:
    """
    Return a single Lorentzian PSD on the FFT grid.

    The PSD shape is ``Delta / (Delta**2 + (f - f0)**2)``, normalised so
    that integrating over [0, +inf) gives ``r**2`` (i.e. the variance is
    ``r**2``). The DC component is set to ``mu**2``.
    """
    df = 1.0 / (np_grid * dt)
    half = np_grid // 2 + 1

    ap = np.zeros(half)
    ap[0] = mu ** 2

    # One-sided normalisation: integral from 0 to infty of Delta/(Delta^2+(f-f0)^2) df
    # equals (pi/2 + atan(f0/Delta)).
    norm = r ** 2 / (np.pi / 2.0 + np.arctan(f0 / Delta))

    j = np.arange(1, half)
    f = j * df
    ap[1:half] = norm * Delta / (Delta ** 2 + (f - f0) ** 2)

    return ap


# ---------------------------------------------------------------------------
# FFT-based frequency-space convolution
# ---------------------------------------------------------------------------
def fftconv(
    dt: float,
    n: int,
    ap: np.ndarray,
    bp: np.ndarray,
) -> np.ndarray:
    """
    Convolve two real one-sided PSDs ``ap`` and ``bp`` in frequency space.

    The Fortran routine builds Hermitian-symmetric complex spectra of length
    n, IFFTs each to the (real) autocovariance, multiplies them in the time
    domain, then forward-FFTs back. We do the same with numpy's ``irfft`` /
    ``rfft`` pair, taking some care with the normalisation so that the
    output matches the original Fortran ``cp`` in both the DC and AC bins.

    The Fortran convention used `four1` (Numerical Recipes style), with
    sign=-1 the inverse and sign=+1 the forward. The DC element of the
    *output* is divided by an extra factor of ``n*dt`` relative to the AC
    bins because the Fortran stuffs the input ``ap(0)`` and ``bp(0)`` into
    the DC slot of the data arrays as ``ap(0)*n*dt`` (and likewise for bp),
    so multiplying gives an extra ``(n*dt)**2`` that we have to scale away
    at the end.
    """
    half = n // 2 + 1

    # Build full Hermitian-symmetric complex spectra of length n.
    # The Fortran code packs:
    #   data[1] = ap(0) * n * dt  (DC)
    #   data[2j+1] = ap(j) for j=1..n/2 with conjugate mirror in data[2n+2-(2j+1)].
    # In numpy's rfft layout, a length-n real signal has rfft length n//2+1
    # with the DC at index 0 and the Nyquist at index n//2.
    A = np.zeros(half, dtype=np.complex128)
    B = np.zeros(half, dtype=np.complex128)
    A[0] = ap[0] * n * dt
    B[0] = bp[0] * n * dt
    A[1:half] = ap[1:half]
    B[1:half] = bp[1:half]

    # IFFT to the time domain. numpy.fft.irfft(X, n) divides by n internally
    # to match the convention X = sum_k x_k exp(-2π i k j / n).
    # The Fortran four1(-1) does NOT divide by n, so we multiply back.
    a_t = np.fft.irfft(A, n=n) * n
    b_t = np.fft.irfft(B, n=n) * n

    # Multiply in the time domain (this is the convolution in frequency).
    c_t = a_t * b_t

    # Forward FFT back to frequency. numpy.fft.rfft has no 1/n; Fortran's
    # four1(+1) also has no 1/n, so the conventions match here.
    C = np.fft.rfft(c_t)

    # The Fortran normalisation: cp(j) = cdata(2j+1) / (n^2 * dt) for j>=1,
    # and cp(0) gets an extra division by (n*dt).
    cp = np.zeros(half)
    cp[1:half] = C[1:half].real / (n ** 2 * dt)
    cp[0] = C[0].real / (n ** 2 * dt) / (n * dt)

    return cp


# ---------------------------------------------------------------------------
# Binning / interpolation onto the output grid
# ---------------------------------------------------------------------------
def modbin(
    dt: float,
    n: int,
    P: np.ndarray,
    fbin: float,
    dfbin: float,
) -> float:
    """
    Average ``P`` over an output bin centred on ``fbin`` of width ``dfbin``,
    or interpolate if the bin is narrower than the FFT spacing.

    ``P`` is the model on FFT bins j = 1..n//2; here we follow the Fortran
    convention of indexing it 1-based via ``P[j-1]`` in Python.
    """
    df = 1.0 / (n * dt)
    T = n * dt

    jmax = int(np.floor((fbin + dfbin / 2.0) / df))
    jmin = int(np.ceil((fbin - dfbin / 2.0) / df))

    if jmax > jmin:
        # The FFT bins fully covered by this output bin: average them.
        # P here is 0-indexed but represents Fortran indices 1..n/2 mapped to
        # P[0..n/2-1]. Convert (jmin..jmax) Fortran-style to Python slicing.
        # Note: P came in as the slice pout[1:n//2+1], so P[j-1] is FFT bin j.
        Pbin = float(np.sum(P[jmin - 1 : jmax])) / (jmax - jmin + 1)
    else:
        # Output bin is narrower than the FFT bin width — interpolate.
        Pbin = interpolate(fbin, T, df, P, n)

    return Pbin


def interpolate(
    fin: float,
    T: float,
    df: float,
    P: np.ndarray,
    np_grid: int,
) -> float:
    """
    Log-linear-style interpolation of ``P`` at frequency ``fin``.

    Reproduces the slightly idiosyncratic Fortran formula: it builds the
    line through the two adjacent FFT samples in (f, P) coordinates but
    uses the integer-frequency ratio ``j/js`` to set the intercept rather
    than the usual ``f_j - f_js`` form. The end result is identical for
    well-separated samples, and we keep it verbatim to preserve the
    original output exactly at very narrow output bins.
    """
    j = int(np.ceil(fin * T))
    js = int(np.floor(fin * T))

    if js < 1:
        j = 2
        js = 1
    elif j > np_grid // 2:
        j = np_grid // 2
        js = np_grid // 2 - 1

    if j == js:
        # Single sample falls exactly on the bin centre.
        Pin = float(P[j - 1])
    else:
        # P[j-1] in Python is FFT bin j in Fortran.
        Pj = float(P[j - 1])
        Pjs = float(P[js - 1])
        m = (Pj - Pjs) / df
        fr = j / js
        c = (fr * Pjs - Pj) / (fr - 1.0)
        Pin = m * fin + c

    return Pin


# ---------------------------------------------------------------------------
# Standalone test driver, replacing the Fortran `program test`
# ---------------------------------------------------------------------------
def _default_params() -> np.ndarray:
    """Return the same parameter set the Fortran test driver uses."""
    param = np.zeros(21)
    param[0] = 13.91     # Sigma0
    param[1] = 7.36      # rbw
    param[2] = 3.0       # kappa
    param[3] = 0.9       # lambda
    param[4] = 0.0       # zeta
    param[5] = 0.165     # Fvar
    param[6] = 25.63     # ro
    param[7] = 3.3       # ri
    param[8] = 10.0      # Q
    param[9] = 10.0      # Qsub
    param[10] = 0.1638   # n_qpo
    param[11] = 0.049    # n_2qpo
    param[12] = 0.0249   # n_3qpo
    param[13] = 0.0      # n_05qpo
    param[14] = 3.738    # gamma_h
    param[15] = 3.0      # gamma_s
    param[16] = 10.0     # M
    param[17] = 0.5      # a
    param[18] = 35.0     # Ndec
    param[19] = 1.0      # mode (1 = hard PSD)
    param[20] = 1.0      # conv (1 = multiply QPO)
    return param


def main(out_path: str = "fort99.txt") -> None:
    """Reproduce the Fortran test driver: dump (f, P/df_bin) to a text file."""
    n = 2 ** 15
    dt = 1.0 / 256.0
    c = 1.05

    far = geobins(dt, n, c)
    param = _default_params()
    photar, _photer = propfluc(far, param)

    with open(out_path, "w") as f:
        for i in range(far.size - 1):
            f_centre = 0.5 * (far[i + 1] + far[i])
            df_bin = far[i + 1] - far[i]
            f.write(f"{f_centre:.7e} {photar[i] / df_bin:.7e}\n")
        f.write("log\n")
        f.write("li s on\n")


if __name__ == "__main__":
    main()

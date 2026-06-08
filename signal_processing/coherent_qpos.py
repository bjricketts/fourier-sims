import numpy as np
import matplotlib.pyplot as plt
from stingray import (simulator, Lightcurve, AveragedPowerspectrum,
                      AveragedCrossspectrum)

# --- Setup ---
dt   = 1/512
mean = 100
T    = 1024 / dt
rms  = 0.5

# Two independent broadband realizations (one per QPO driver)
sim_a = simulator.Simulator(N=int(T), mean=mean, dt=dt, rms=rms, poisson=True, random_state=1)
w     = np.fft.rfftfreq(sim_a.N, d=sim_a.dt)[1:]

def lorentzian(f, f0, q, norm):
    df = f0 / q
    return norm * (df / np.pi) / ((f - f0)**2 + df**2)

psd_model = lorentzian(w, 0.4, 0.3, 1.0) + lorentzian(w, 2.0, 0.5, 1.2)
lc_a      = sim_a.simulate(psd_model)

broadband = lc_a.counts - lc_a.counts.mean()   # drives primary QPO

# Observable broadband (use realization A so band 1/2 broadband is the shared one)
broadband = broadband

# --- DHO filter ---
def dho_filter(x, dt, f0, zeta, gain=1.0):
    N    = len(x)
    freq = np.fft.rfftfreq(N, d=dt)
    w_   = 2*np.pi*freq
    w0   = 2*np.pi*f0
    H    = w0**2 / (w0**2 - w_**2 + 2j*zeta*w0*w_)
    X    = np.fft.rfft(x - x.mean())
    return np.fft.irfft(gain * H * X, n=N)

# --- Fixed broadband levels ---
mean1, mean2         = 100., 40.
frac_bb_1, frac_bb_2 = 0.25, 0.12

broadband1 = broadband * (frac_bb_1 * mean1 / broadband.std())
broadband2 = broadband * (frac_bb_2 * mean2 / broadband.std())

# --- Two QPOs ---
# Primary: f0 = 1.5 Hz, narrow (Q=10) — driven by broadband A
# Secondary (hidden): f0 = 2.0 Hz, broader (Q=5) — driven by INDEPENDENT broadband B
f0_pri, zeta_pri = 1.5, 0.05    # Q = 10
f0_sec, zeta_sec = 1.6, 0.05    # Q = 10

resp_pri = dho_filter(broadband, dt, f0=f0_pri, zeta=zeta_pri, gain=1.0)
resp_sec = dho_filter(broadband, dt, f0=f0_sec, zeta=zeta_sec, gain=1.0)

qpo_pri_unit = resp_pri / resp_pri.std()
qpo_sec_unit = resp_sec / resp_sec.std()

# Secondary QPO fixed at small fractional rms (the "hidden" one)
frac_qpo_sec = 5e-3

# --- Sweep over PRIMARY QPO fractional rms ---
frac_qpo_vals = np.geomspace(1e-3, 5e-2, 6)
seg = 16.0
t   = np.arange(lc_a.n) * dt
rng = np.random.default_rng(11)

results = []
for frac_qpo in frac_qpo_vals:
    qpo_pri_1 = qpo_pri_unit * (frac_qpo     * mean1)
    qpo_pri_2 = qpo_pri_unit * (frac_qpo     * mean2)
    qpo_sec_1 = qpo_sec_unit * (frac_qpo_sec * mean1)
    qpo_sec_2 = qpo_sec_unit * (frac_qpo_sec * mean2)

    rate1 = mean1 + qpo_pri_1 + qpo_sec_1
    rate2 = mean2 + qpo_pri_2 + qpo_sec_2

    lc1 = Lightcurve(t, rng.poisson(np.clip(rate1, 0, None)), dt=dt, skip_checks=True)
    lc2 = Lightcurve(t, rng.poisson(np.clip(rate2, 0, None)), dt=dt, skip_checks=True)

    ps1 = AveragedPowerspectrum(lc1, segment_size=seg, norm='frac')
    ps2 = AveragedPowerspectrum(lc2, segment_size=seg, norm='frac')
    cs  = AveragedCrossspectrum(lc2, lc1, segment_size=seg, norm='frac')
    coh, _     = cs.coherence()
    lag, lag_e = cs.time_lag()

    results.append(dict(frac_qpo=frac_qpo, ps1=ps1, ps2=ps2, cs=cs,
                        coh=coh, lag=lag, lag_e=lag_e))

# --- Plot ---
font = {'size': 14}
plt.rc('font', **font)

fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
cmap = plt.cm.viridis(np.linspace(0.15, 0.9, len(frac_qpo_vals)))

for ax in axes:
    ax.axvline(f0_pri, ls='--', color='red',    alpha=0.4)
    ax.axvline(f0_sec, ls='--', color='orange', alpha=0.4)

for r, c in zip(results, cmap):
    lbl = f'{r["frac_qpo"]*100:.2f}%'
    axes[0].loglog(r['ps1'].freq, r['ps1'].freq*r['ps1'].power,
                   color=c, alpha=0.85, label=lbl)
    axes[0].loglog(r['ps2'].freq, r['ps2'].freq*r['ps2'].power,
                   color=c, alpha=0.5, ls='--')
    axes[1].errorbar(r['cs'].freq, r['lag']*1e3, yerr=r['lag_e']*1e3,
                     color=c, alpha=0.8, ms=3, lw=1)
    axes[2].semilogx(r['cs'].freq, r['coh'], color=c, alpha=0.85,
                     drawstyle='steps-mid')

axes[0].set(ylabel=r'PSD × ν  [(rms/mean)$^2$]',
            title=(rf'Multiple coherent QPOs — hidden secondary at {f0_sec}Hz (rms = {frac_qpo_sec*100:.2f}%)'
                   '\n'
                   rf'$Q_1$ = {1/(2*zeta_pri):.1f}, $Q_2$ = {1/(2*zeta_sec):.1f}'))
axes[0].legend(title='Primary QPO frac rms', fontsize=8, ncol=2)

axes[1].axhline(0.0, ls=':', color='gray', alpha=0.5)
axes[1].set(ylabel='Time lag [ms]', ylim=(-120, 120))

axes[2].set(xscale='log', xlabel='Frequency [Hz]', ylabel='Coherence',
            ylim=(0.0, 1.02), xlim=(0.1, 20))

plt.tight_layout()
plt.savefig("figures/two_coherent_qpos.pdf")
plt.show()
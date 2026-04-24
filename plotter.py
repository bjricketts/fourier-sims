def plot(t_plot,x1_noise,x2_noise,ps1,ps2,omega1,omega2,avgcs12,x_lim_lower,x_lim_upper,
         file = "figures/decaying_sine_wave.png",omega3=None,theory=False,hil=False):
    """
    Plots results

    Parameters
    -----------
    t_plot : array
        Time array for the plot.
    x1_noise : Lightcurve
        Lightcurve object for the first signal with noise.
    x2_noise : Lightcurve
        Lightcurve object for the second signal with noise.
    ps1 : Powerspectrum
        Power spectrum object for the first signal.
    ps2 : Powerspectrum
        Power spectrum object for the second signal.
    omega1 : float
        Frequency of the first signal.
    omega2 : float
        Frequency of the second signal.
    cs12 : Crossspectrum
        Cross spectrum object between the two signals.
    avgcs12 : AveragedCrossspectrum
        Averaged cross spectrum object between the two signals.
    x_lim_lower : float
        Lower limit for the x-axis.
    x_lim_upper : float
        Upper limit for the x-axis.

    """
    fig = plt.figure(figsize=(24, 12))

    ax1=plt.subplot(2,1,1)
    ax1.plot(t_plot, x1_noise.counts, c="#e66101", label=r"$\phi_1$")
    ax1.plot(t_plot, x2_noise.counts, c="#5e3c99", label=r"$\phi_2$")
    ax1.set_xlabel('Time (s)')
    ax1.legend()
    ax1.set_ylabel('Amplitude')
    ax1.set_xlim(0,5)

    ax2=plt.subplot(2,4,5)
    ax2.errorbar(ps1.freq, ps1.freq*ps1.power, yerr=ps1.freq*ps1.power_err, c="#e66101",label=r"$\phi_1$ PSD")
    ax2.errorbar(ps2.freq, ps2.freq*ps2.power, yerr=ps2.freq*ps2.power_err, c="#5e3c99",label=r"$\phi_2$ PSD")
    ax2.axvline(omega1, color="#b2abd2", linestyle='--')
    ax2.axvline(omega2, color="#b2abd2", linestyle='--')
    if omega3 is not None:
        ax2.axvline(omega3, color="#b2abd2", linestyle='--')
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.legend(loc="lower left")
    ax2.set_xlabel("Frequency (Hz)")
    ax2.set_ylabel(r"Frequency * Power (rms$^2$)")
    ax2.set_ylim(1e-8, 1e-1)
    ax2.set_xlim(x_lim_lower, x_lim_upper)

    ax3=plt.subplot(2,4,6)
    ax3.plot(avgcs12.freq, np.real(avgcs12.power), c="#e66101",label="Real")
    ax3.plot(avgcs12.freq, np.imag(avgcs12.power), c="#5e3c99",label="Imaginary")
    ax3.axvline(omega1, color="#b2abd2", linestyle='--')
    ax3.axvline(omega2, color="#b2abd2", linestyle='--')
    if omega3 is not None:
        ax3.axvline(omega3, color="#b2abd2", linestyle='--')
    ax3.set_yscale('symlog',linthresh=1e-4)
    ax3.set_xscale("log")
    ax3.set_ylabel('Cross-spectrum')
    ax3.legend()
    ax3.set_ylim(-2e-5, 1e-4)
    ax3.set_xlim(x_lim_lower, x_lim_upper)
    ax3.set_xlabel("Frequency (Hz)")

    ax4 = plt.subplot(2,4,7)
    ax4.errorbar(avgcs12.freq, avgcs12.phase_lag()[0], yerr=avgcs12.lag_err, c="#fdb863", label=r"$\phi_1/\phi_2$")
    #ax4.axvline(omega1, color="#b2abd2", linestyle='--')
    #ax4.axvline(omega2, color="#b2abd2", linestyle='--')
    if omega3 is not None:
        ax4.axvline(omega3, color="#b2abd2", linestyle='--')
    if theory == True:
        omegax = omega1
        omegay = omega2
        dampenx = 0.25
        dampeny = 0.5
        numerator = (dampenx)*ps1.freq*(omegay**2 - ps1.freq**2) - (dampeny)*ps2.freq*(omegax**2 - ps2.freq**2)
        denominator = (omega2**2-ps1.freq**2)*(omega1**2-ps1.freq**2) + (dampenx*dampeny)*(ps1.freq)**2
        manual_delta_theta = np.arctan(numerator/denominator)
        theory_corrected = np.angle(avgcs12.power) - manual_delta_theta
        ax4.plot(ps1.freq,manual_delta_theta,label="theory",c="black",ls="--")
        max_phase_pos = 0.5*(np.sqrt(dampenx*dampeny) + np.sqrt(dampenx*dampeny + 4*omegax*omegay))
        max_phase_neg = 0.5*(-np.sqrt(dampenx*dampeny) + np.sqrt(dampenx*dampeny + 4*omegax*omegay))
        #ax4.axvline(max_phase_pos, color="black", linestyle="--")
        #ax4.axvline(max_phase_neg, color="black", linestyle="--")
    if hil == True:
        concat_freq = np.concatenate([-ps1.freq[::-1],ps1.freq])
        subject_band = np.concatenate([ps1.power[::-1],ps1.power])
        ref_band = np.concatenate([ps2.power[::-1],ps2.power])
        #subject_band = ps1.power
        #ref_band = ps2.power
        print(subject_band)
        log_power_ratio = np.log(subject_band/ref_band)
        delta_theta = (np.imag(hilbert(log_power_ratio))/2)[:-(ps2.power.shape[0])]

        #int_delta_theta = -kramers_kronig_onesided_fast(ps1.freq, log_power_ratio)/2
        print(delta_theta)
        corrected_phase = np.angle(avgcs12.power) - delta_theta
        print(corrected_phase)
        ax4.plot(ps1.freq, delta_theta, label="Hilbert", c="black", ls=":")
        ax4.plot(ps1.freq, corrected_phase, label="Corrected", c="red", ls="-")

    ax4.set_xscale("log")
    ax4.set_ylabel('Phase (radians)')
    ax4.set_xlabel("Frequency (Hz)")
    ax4.set_xlim(x_lim_lower, x_lim_upper)
    ax4.set_ylim(np.min(avgcs12.phase_lag()[0][(avgcs12.freq>x_lim_lower)&(avgcs12.freq<x_lim_upper)]), 
                 np.max(avgcs12.phase_lag()[0][(avgcs12.freq>x_lim_lower)&(avgcs12.freq<x_lim_upper)]))
    ax4.legend(loc="upper left")

    ax5 = plt.subplot(2,4,8)
    ax5.plot(avgcs12.freq, avgcs12.coherence()[0], c="#fdb863", label=r"$\phi_1/\phi_2$")
    ax5.axvline(omega1, color="#b2abd2", linestyle='--')
    ax5.axvline(omega2, color="#b2abd2", linestyle='--')
    if omega3 is not None:
        ax5.axvline(omega3, color="#b2abd2", linestyle='--')
    ax5.set_xscale("log")
    ax5.set_ylabel('Coherence')
    ax5.set_xlabel("Frequency (Hz)")
    ax5.set_ylim(0,1.01)
    ax5.set_xlim(x_lim_lower, x_lim_upper)
    #ax5.legend()

    plt.tight_layout()
    plt.savefig(file)
    plt.show()
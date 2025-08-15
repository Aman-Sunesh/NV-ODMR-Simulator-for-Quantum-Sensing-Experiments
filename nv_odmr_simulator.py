import numpy as np
import matplotlib.pyplot as plt
import math

rng = np.random.default_rng(0)  # for shot noise & reproducibility

# ---------- Acquisition / Integration Controls ----------
baseline_slope = 0.0        # unitless per GHz: small linear tilt of baseline vs frequency (0.0 = flat)
read_time = 2e-5            # seconds per read gate (longer -> more counts -> lower noise)
dynamic_steps = 400         # number of frequency points in the sweep (df = (f_end - f_start)/(steps-1))
frames = 10000              # repeats per frequency point (more -> more counts -> lower noise)
averages = 1                # outer averages of the whole sweep (rarely needed; multiplies counts)
cps = 200_000               # counts per second (MW OFF baseline brightness of your setup)
use_noise = True            # True: add Poisson shot noise; False: ideal noiseless traces

# ---------- Physics / Sweep Parameters ----------
f_start = 2.60e9            # sweep start frequency in Hz
f_end   = 3.00e9            # sweep end frequency in Hz
D_mag   = 2.870e9           # zero-field splitting D (Hz) for NV ground state (~2.870 GHz)
gamma_e = 28.0e9            # electron gyromagnetic ratio in Hz/T (~28 GHz/T)
theta   = 0                 # angle (degrees) between B field and NV axis (0° = aligned)
C_max   = 0.04              # max fractional contrast at high power (saturation limit)
A_PARALLEL = 2.16e6         # 14N hyperfine splitting A∥ in Hz (~2.16 MHz spacing)
hyperfine_on = True         # include hyperfine triplets (three lines per transition)

fwhm0   = 8e6               # low-power Lorentzian FWHM in Hz (~1/(pi*T2)); sets base linewidth
power_fac = 1               # saturation parameter proxy s (dimensionless): linewidth → fwhm0*sqrt(1+s)
B_mag = 0.0                 # Tesla. 0 → single dip at D; ~2e-3–3e-3 T (θ=0°) → two dips split by ~2*γe*B∥
B_par = B_mag * math.cos(math.radians(theta))  # Tesla: projection of B along NV axis; sets dip positions

f = np.linspace(f_start, f_end, dynamic_steps)

def lorentzian(f, f0, fwhm, A=1.0):
    hwhm = 0.5 * fwhm
    return ((A * (hwhm**2)) / (hwhm**2 + (f - f0)**2))

f_minus = D_mag - (gamma_e * B_par)
f_plus = D_mag + (gamma_e * B_par)

s = float(power_fac)                       
fwhm = fwhm0 * np.sqrt(1.0 + s)
C_on = C_max * (s / (1.0 + s))   # depth saturates with power

def triplet(f, fc, A, fwhm):
    # peak-normalize the triplet so its max is 1.0
    L0 = lorentzian(np.array([fc]), fc, fwhm)[0]
    Ls = lorentzian(np.array([fc + A]), fc, fwhm)[0]
    norm = L0 + 2*Ls

    return (lorentzian(f, fc - A, fwhm) +
            lorentzian(f, fc, fwhm)     +
            lorentzian(f, fc + A, fwhm)) / norm

if hyperfine_on:
    dips = triplet(f, f_minus, A_PARALLEL, fwhm) + triplet(f, f_plus, A_PARALLEL, fwhm)
else:
    dips = lorentzian(f, f_minus, fwhm) + lorentzian(f, f_plus, fwhm)

baseline = 1.0 + baseline_slope * ((f - f.mean())/1e9)

I0_total = cps * read_time * frames * averages
I_off_ideal = I0_total * baseline
I_on_ideal = I0_total * baseline * np.maximum(1.0 - C_on * dips, 0.02)

if (use_noise == True):
    I_off = rng.poisson(I_off_ideal)
    I_on  = rng.poisson(I_on_ideal)
    ref = (I_off - I_on) / np.clip(I_off, 1, None)
else:
    I_off = I_off_ideal
    I_on  = I_on_ideal
    ref   = (I_off - I_on) / I_off

fig, axes = plt.subplots(3, 1, figsize=(7, 12), sharex=True)

axes[0].plot(f*1e-9, I_on, label="MW ON (counts)")
axes[0].set_ylabel("Counts per Frequency Point")
axes[0].set_title("Simulated NV ODMR — MW ON")
axes[0].legend()

axes[1].plot(f*1e-9, I_off, label="MW OFF (counts)")
axes[1].set_ylabel("Counts per Frequency Point")
axes[1].set_title("Simulated NV ODMR — MW OFF")
axes[1].legend()

axes[2].plot(f*1e-9, ref, label="Referenced ODMR")
axes[2].set_xlabel("Frequency (GHz)")
axes[2].set_ylabel("Referenced Contrast (I_off - I_on)/I_off")
axes[2].set_title("Simulated NV ODMR — Referenced")
axes[2].legend()

# mark centers
for ax in axes:
    for fc in (f_minus, f_plus):
        if hyperfine_on:
            for k in (-1, 0, 1):
                ax.axvline((fc + k*A_PARALLEL)*1e-9, ls='--', alpha=0.25)
        else:
            ax.axvline(fc*1e-9, ls='--', alpha=0.4)

fig.tight_layout()
plt.show()

# readouts
print(f"Transition freqs: {f_minus/1e9:.6f} GHz and {f_plus/1e9:.6f} GHz")
print(f"FWHM used: {fwhm/1e6:.2f} MHz;  C_on={C_on:.3f}")
df = (f_end - f_start) / (dynamic_steps - 1)
print(f"Points per FWHM: {fwhm/df:.1f} (aim ≥ 8–10)")
print(f"Peak referenced contrast (approx): {ref.max():.4f}")


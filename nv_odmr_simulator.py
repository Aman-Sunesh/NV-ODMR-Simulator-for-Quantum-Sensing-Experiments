import numpy as np
import matplotlib.pyplot as plt
import math

rng = np.random.default_rng(0)  # for shot noise & reproducibility

baseline_slope = 0.0        
read_time = 2e-5
dynamic_steps = 400
frames = 10000  
averages = 1
cps = 200_000
use_noise = True


f_start = 2.60e9
f_end = 3.00e9
D_mag = 2.870e9
gamma_e = 28.0e9
theta = 0
C_max = 0.04
A_PARALLEL = 2.16e6   # ^14N hyperfine (Hz)
hyperfine_on = True
fwhm0    = 8e6 
power_fac = 1
B_mag = 0 # 3e-3
B_par = B_mag * math.cos(theta * (math.pi/180))

f = np.linspace(f_start, f_end, dynamic_steps)

def lorentzian(f, f0, fwhm, A=1.0):
    hwhm = 0.5 * fwhm
    return ((A * (hwhm**2)) / (hwhm**2 + (f - f0)**2))

f_minus = D_mag - (gamma_e * B_par)
f_plus = D_mag + (gamma_e * B_par)

s = float(power_fac)                       
fwhm = fwhm0 * np.sqrt(1.0 + s)
C_on = C_max * (s / (1.0 + s)) # depth saturates with power

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

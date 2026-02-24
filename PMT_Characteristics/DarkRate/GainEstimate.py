#gain estimate from dark rate

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
from scipy.optimize import curve_fit
from scipy.special import gamma

allfiles = glob.glob('PMT1_Cleaned/*.csv')
print("Number of files found:", len(allfiles))
e = 1.6e-19

lights_on = {}
lights_off = {}

#need to extract the mean integral from a gaussian fit for each file

for file in allfiles:
    print("\nProcessing file:", file)
    df = pd.read_csv(file)

    if 'integral_pC' not in df.columns:
        print("ERROR: 'integral_pC' column not found in CSV!")
        continue

    all_integrals_picocharge = df['integral_pC'].dropna().to_numpy()
    if all_integrals_picocharge.size == 0:
        print(f"  Skipping file {file} because integral_pC is empty after dropping NaNs")
        continue

    #plt.figure(figsize=(6, 6))
    all_integrals_picocharge = df['integral_pC'].to_numpy()

    def gaussian_pulse(all_integrals_picocharge, amp,  mu, sigma):
        return amp * (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp((-0.5 * (all_integrals_picocharge - mu)**2)/ sigma**2)
    int_number_of_bins = int((np.nanmax(all_integrals_picocharge)-np.nanmin(all_integrals_picocharge))/0.1) #chosen 1 bin = 0.1 pC
    n, bins, patches = plt.hist(all_integrals_picocharge, bins=int_number_of_bins, color = 'skyblue', edgecolor = 'black')
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    fit_mask = (bin_centers >= 0) & (bin_centers <= 6)
    X_fit = bin_centers[fit_mask]
    Y_fit = n[fit_mask]
    p0 = [np.max(Y_fit), np.mean(X_fit), np.std(X_fit)]

    
    # Initialize values to NaN in case of failure
    mu_int = np.nan
    err_mu = np.nan
    err_sigma = np.nan
    file_gain = np.nan
    
    try:
        popt, pcov = curve_fit(gaussian_pulse, X_fit, Y_fit, p0=p0)
        if not np.all(np.isfinite(pcov)):
            raise RuntimeError ("Covariance matrix is non-finite.")
        
        mu_int = popt[1]
        sigma_int = popt[2]
        err_mu = np.sqrt(pcov[1,1])
        err_sigma = np.sqrt(pcov[2,2])
        pulse_fit = gaussian_pulse(bin_centers, *popt)
        plt.plot(bin_centers, pulse_fit, 'r--', linewidth=1)
        plt.xlim([0,10])
        plt.xlabel('Pulse integral (pC)')
        plt.ylabel('Counts per pC')
        plt.legend([r'$\mathrm{Fit:}\ \mu=%.6f,\ \sigma=%.6f$' % (popt[1], popt[2])])
        #plt.show()

        file_gain = (mu_int * 1e-12) / e   
    except Exception as e:
        print(f"  WARNING: Gaussian fit failed for {file}. Data skipped. Error: {e}")

    print(f"  Computed file gain: {file_gain:.3e}")

    match = re.search(r"darkrate_(\d+)_10000_(\d+)_(on|off)", file)
    if not match:
        print(f"No name match for file: {file}")
        continue

    voltage = int(match.group(1))
    trial = int(match.group(2))
    state = match.group(3)

    if state == "on":
        lights_on.setdefault(voltage,[]).append((file_gain, err_mu))
    if state == "off":
        lights_off.setdefault(voltage,[]).append((file_gain, err_mu))

for voltage, rate_list in lights_on.items():
    file_gain, err_mu = zip(*rate_list)
    file_gain = np.array(file_gain)
    err_mu = (np.array(err_mu) * 1e-12) / (e*1e7)
    #x = len(file_gain)
    #mean_err = np.sqrt((1/x)**2 * np.sum((err_mu**2)))  # combine in quadrature
    mean_gain = np.nanmean(file_gain) / 1e7
    valid_errs = err_mu[np.isfinite(err_mu)]
    valid_x = len(valid_errs)
    if valid_x > 0:
        mean_err = np.sqrt((1/valid_x)**2 * np.sum((valid_errs**2))) 
    else:
         mean_err = np.nan
    lights_on[voltage] = (mean_gain, mean_err) #replace the list with a single average rate element



for voltage, rate_list in lights_off.items():
    file_gain, err_mu = zip(*rate_list)
    file_gain = np.array(file_gain)
    err_mu = (np.array(err_mu) * 1e-12) / (e*1e7) 
    #x = len(file_gain)
    #mean_err = np.sqrt((1/x)**2 * np.sum((err_mu**2)))  # combine in quadrature
    mean_gain = np.nanmean(file_gain) / 1e7
    valid_errs = err_mu[np.isfinite(err_mu)]
    valid_x = len(valid_errs)
    if valid_x > 0:
        mean_err = np.sqrt((1/valid_x)**2 * np.sum((valid_errs**2))) 
    else:
         mean_err = np.nan
    lights_off[voltage] = (mean_gain, mean_err)

voltages_on = sorted(lights_on.keys())
average_gain_on = [lights_on[v][0] for v in voltages_on]
average_gain_err_on = [lights_on[v][1] for v in voltages_on]
#voltages_off = sorted(lights_off.keys())
#average_gain_off = [lights_off[v][0] for v in voltages_off]
#average_gain_err_off = [lights_off[v][1] for v in voltages_off]

"""
plt.figure(figsize=(6,6))
plt.errorbar(voltages_on, average_gain_on, average_gain_err_on, ms=4, capsize=3, fmt='o', color='blue', label = 'Lights On')
plt.errorbar(voltages_off, average_gain_off, average_gain_err_off, ms=4, capsize=3, fmt='o', color='magenta', label = 'Lights Off')
plt.xlabel("High Voltage Source (V)")
plt.ylabel("Gain Factor")
plt.title("15σ threshold cut applied")
plt.legend()
plt.show()
"""

def exponential_C(V, A, B, C):
    return A * np.exp(B * V) + C

def exponential(V, A, B):
    return A * np.exp(B * V)

def powerlaw(V, A, B):
    return A*(V**(14*B))

voltages_on_np = np.array(voltages_on)
average_gain_on_np = np.array(average_gain_on)
average_err_on_np = np.array(average_gain_err_on)

finite_mask_on = np.isfinite(average_gain_on_np) & np.isfinite(average_err_on_np)

V_fit_on = voltages_on_np[finite_mask_on]
Gain_fit_on = average_gain_on_np[finite_mask_on]
Err_fit_on = average_err_on_np[finite_mask_on]

try:
    popt_on, pcov_on = curve_fit(
        f=powerlaw,
        xdata=V_fit_on,
        ydata=Gain_fit_on,
        sigma=Err_fit_on,
        absolute_sigma=True,
        p0=[1e-30, 0.8]  # Initial guess for A, B
    )
    gain_fit_on_points = powerlaw(V_fit_on, *popt_on)
    perr_on = np.sqrt(np.diag(pcov_on))
except RuntimeError:
    print("Lights On curve fit failed.")
    popt_on = None

#voltages_off_np = np.array(voltages_off)
#average_gain_off_np = np.array(average_gain_off)
#average_err_off_np = np.array(average_gain_err_off)

#finite_mask_off = np.isfinite(average_gain_off_np) & np.isfinite(average_err_off_np)

#V_fit_off = voltages_off_np[finite_mask_off]
#Gain_fit_off = average_gain_off_np[finite_mask_off]
#Err_fit_off = average_err_off_np[finite_mask_off]


#try:
 #   popt_off, pcov_off = curve_fit(
  #      f=powerlaw,
   #     xdata=V_fit_off,
    #    ydata=Gain_fit_off,
     #   sigma=Err_fit_off,
      #  absolute_sigma=True,
       # p0=[1e-30, 0.8]  # Initial guess for A, B
    #)
    #gain_fit_off_points = powerlaw(V_fit_off, *popt_off)
    #perr_off = np.sqrt(np.diag(pcov_off))
#except RuntimeError:
    #print("Lights Off curve fit failed.")
    #popt_off = None


fig, (ax_main, ax_res) = plt.subplots(
    nrows=2, 
    ncols=1, 
    sharex=True, 
    figsize=(7, 7), 
    gridspec_kw={'height_ratios': [3, 1], 'hspace': 0} # 3:1 ratio and no vertical space
)

ax_main.errorbar(voltages_on, average_gain_on, yerr = average_gain_err_on, ms=2, capsize=3, fmt='o', color='blue', label = 'Lights On')
#ax_main.errorbar(voltages_off, average_gain_off, yerr=average_gain_err_off, ms=2, capsize=3, fmt='o', color='magenta', label = 'Lights Off')
if popt_on is not None:
    ax_main.plot(V_fit_on, gain_fit_on_points, 'b-', label=(r'Gain = $A \cdot V^{14 \cdot B}$' + f': A = {popt_on[0]:.3e} ± {perr_on[0]:.3e} \n $\qquad$ $\qquad$ $\qquad$     B = {popt_on[1]:.3e} ± {perr_on[1]:.3e}'))
# Plot the 'Lights Off' fit line, connecting the fitted points
#if popt_off is not None:
    #ax_main.plot(V_fit_off, gain_fit_off_points, color='magenta', linestyle='-', label=(r'Gain = $A \cdot V^{14 \cdot B}$' + f': A = {popt_on[0]:.3e} ± {perr_on[0]:.3e} \n $\qquad$ $\qquad$ $\qquad$     B = {popt_on[1]:.3e} ± {perr_on[1]:.3e}'))

#plt.xlabel('High Voltage Source (V)')
ax_main.set_ylabel('Gain (x10$^7$)', fontsize=14)
#plt.title('Dark Rate vs. Voltage with Exponential Fit (15-sigma threshold)')
ax_main.legend(loc='upper left', fontsize=14, frameon=False)
ax_main.tick_params(axis='both', labelsize=14)
# Remove the x-axis tick labels for the main plot since the x-axis is shared
ax_main.tick_params(axis='x', labelbottom=False)

# --- Residuals Calculation and Plotting ---
ax_res.tick_params(axis='both', labelsize=14)
ax_res.axhline(0, color='gray', linestyle='--') # Draw a horizontal line at y=0

# 1. Calculate and Plot Lights On Residuals
if popt_on is not None:
    # Residuals = (Measured Rate) - (Fitted Rate at that Voltage)
    residuals_on = Gain_fit_on - gain_fit_on_points
    
    # Optional: Plot standardized residuals (Residuals / Error) if you want to normalize by uncertainty
    # standardized_residuals_on = residuals_on / average_err_on_np
    ax_res.errorbar(
        V_fit_on, 
        residuals_on, 
        yerr=Err_fit_on, 
        ms=3, 
        capsize=4, 
        fmt='o', 
        color='blue', 
        label='Lights On Residuals'
    )

# 2. Calculate and Plot Lights Off Residuals
#if popt_off is not None:
    #residuals_off = Gain_fit_off - gain_fit_off_points
    
    #ax_res.errorbar(
     #   V_fit_off, 
      #  residuals_off, 
       # yerr=Err_fit_off, 
        #ms=3, 
        #capsize=4, 
        #fmt='o', 
        #color='magenta', 
        #label='Lights Off Residuals'
    #)

ax_res.set_xlabel('High Voltage Source (V)', fontsize=14)
ax_res.set_ylabel('Residuals (x10$^7$)', fontsize=14)
ax_res.grid(True, linestyle=':', alpha=0.7)
#ax_res.legend(loc='upper left', fontsize='large', frameon=False)

# 3. Final Display
plt.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust layout to make room for suptitle
plt.show()

if popt_on is not None:
    print("\n##  Lights On Fit Parameters")
    print(f"A: {popt_on[0]:.4e} ± {perr_on[0]:.4e}")
    print(f"B: {popt_on[1]:.4f} ± {perr_on[1]:.4f}")
    #print(f"C: {popt_on[2]:.2f} ± {perr_on[2]:.2f} Hz")
    
#if popt_off is not None:
    #print("\n##  Lights Off Fit Parameters")
    #print(f"A: {popt_off[0]:.4e} ± {perr_off[0]:.4e}")
    #print(f"B: {popt_off[1]:.4f} ± {perr_off[1]:.4f}")
    #print(f"C: {popt_off[2]:.2f} ± {perr_off[2]:.2f} Hz")

    print("\n--- Exponential Fit Results: Rate(V) = A * exp(B * V) + C ---")

chi2_on = None
if popt_on is not None:
    # Calculate the chi-squared statistic for Lights On
    # Residuals = average_rate_on_np - rate_fit_on_points (calculated in previous step)
    # The error is average_err_on_np
    
    # Check for zero errors to avoid division by zero
    non_zero_errors_on = Err_fit_on > 0
    if np.any(non_zero_errors_on):
        chi2_on = np.sum(
            (residuals_on[non_zero_errors_on] / Err_fit_on[non_zero_errors_on])**2
        )
    else:
        print("Cannot calculate Chi-squared for Lights On: all errors are zero.")

#chi2_off = None
#if popt_off is not None:
    # Calculate the chi-squared statistic for Lights Off
    # Residuals = average_rate_off_np - rate_fit_off_points
    # The error is average_err_off_np
    
 #   non_zero_errors_off = Err_fit_off > 0
  #  if np.any(non_zero_errors_off):
   #     chi2_off = np.sum(
    #        (residuals_off[non_zero_errors_off] / Err_fit_off[non_zero_errors_off])**2
     #   )
    #else:
     #   print("Cannot calculate Chi-squared for Lights Off: all errors are zero.")

if popt_on is not None:
    N_on = len(V_fit_on)
    P = len(popt_on)
    nu_on = N_on - P
    
    print("\n## Lights On Fit Parameters")
    # ... (Your existing A, B, C prints) ...
    
    if chi2_on is not None and nu_on > 0:
        chi2_nu_on = chi2_on / nu_on
        print(f"\nChi-Squared ($\chi^2$): {chi2_on:.3f}")
        print(f"Degrees of Freedom ($\nu$): {nu_on}")
        print(f"Reduced Chi-Squared ($\chi^2_\nu$): {chi2_nu_on:.3f}")
    
#if popt_off is not None:
 #   N_off = len(V_fit_off)
  #  P = len(popt_off)
   # nu_off = N_off - P

    #print("\n## Lights Off Fit Parameters")
    # ... (Your existing A, B, C prints) ...
    
    #if chi2_off is not None and nu_off > 0:
     #   chi2_nu_off = chi2_off / nu_off
      #  print(f"\nChi-Squared ($\chi^2$): {chi2_off:.3f}")
       # print(f"Degrees of Freedom ($\nu$): {nu_off}")
        #print(f"Reduced Chi-Squared ($\chi^2_\nu$): {chi2_nu_off:.3f}")



print("\n--- Lights ON (Data Points) ---")
print("voltages_on =", voltages_on)
print("average_gain_on =", average_gain_on)
print("average_gain_err_on =", average_gain_err_on)

if popt_on is not None:
    V_model_on = np.linspace(min(V_fit_on), max(V_fit_on), 400)
    gain_model_on = powerlaw(V_model_on, *popt_on)
    print("V_model_on =", V_model_on)
    print("gain_model_on =", gain_model_on)


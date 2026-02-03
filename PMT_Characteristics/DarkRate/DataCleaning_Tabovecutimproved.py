#import csv
import pandas as pd
import os
from os.path import exists
import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#arguments for terminal
folder_in = sys.argv[1]       #new argument: folder where input CSV is
filename_in = sys.argv[2]     #base name of CSV (without .csv)
folder_out = sys.argv[3]      #output folder

#outputted file
filename_out = f"{filename_in}_cleaned.csv"
file_path = os.path.join(os.getcwd(), folder_out, filename_out)

#read CSV
csv_path = os.path.join(os.getcwd(), folder_in, f"{filename_in}.csv")
df = pd.read_csv(csv_path)

#make data frame for cleaned data
cleaned_df = df.copy()



#15sigma threshold clean
#sigma_thresh = df["time_above_threshold"] > 0 #keep events which cross the threshold
#cleaned_df = cleaned_df[sigma_thresh]

sigma_thresh = df["total_time_above"] > 0 #keep events which cross the threshold
cleaned_df = cleaned_df[sigma_thresh]



#baseline cut
sigma_b_array_mV = df['sd_baseline'].to_numpy() * 1e3
def gaussian_sb(sigma_b_array_mV, amp,  mu, sigma):
    return amp * (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp((-0.5 * (sigma_b_array_mV - mu)**2)/ sigma**2)
sb_number_of_bins = int((np.max(sigma_b_array_mV)-np.min(sigma_b_array_mV))/0.01) #chosen 1 bin = 0.01mV
n, bins, patches = plt.hist(sigma_b_array_mV, bins=sb_number_of_bins, color='skyblue', edgecolor='black')
bin_centers = 0.5 * (bins[1:] + bins[:-1])
MAX_ITERATIONS = 3
SIGMA_THRESHOLD = 3.0  # Remove points > 3 standard deviations away
fit_mask = (bin_centers >=0) & (bin_centers <=2)
bin_centers_filtered = bin_centers[fit_mask]
n_filtered = n[fit_mask]
X_fit = bin_centers_filtered.copy()
Y_fit = n_filtered.copy()
p0 = [np.max(n_filtered), np.mean(bin_centers_filtered), np.std(bin_centers_filtered)]
final_popt = None
final_pcov = None
final_mask = None
for i in range(MAX_ITERATIONS):
    try:
        # 1. Fit the current data
        popt, pcov = curve_fit(gaussian_sb, X_fit, Y_fit, p0=p0)
        # store results
        final_popt = popt
        final_pcov = pcov
        # 2. Calculate the fitted curve and residuals
        Y_fitted = gaussian_sb(X_fit, *popt)
        residuals = Y_fit - Y_fitted
        # 3. Calculate the standard deviation of the residuals (a measure of fit scatter)
        residual_std = np.std(residuals)
        
        # 4. Create a mask to identify outliers
        mask = np.abs(residuals) < (SIGMA_THRESHOLD * residual_std)
        
        # If no points are removed, the process is stable; stop
        if np.all(mask):
            print(f"Converged after {i+1} iterations.")
            break
            
        # 5. Apply the mask to remove outliers for the next iteration
        X_fit = X_fit[mask]
        Y_fit = Y_fit[mask]
        
        # Update the initial guess (p0) for the next iteration with the current fit
        p0 = popt
    except RuntimeError:
        print(f"Fit failed during iteration {i+1}. Stopping.")
        break
else:
    print(f"Fit completed max {MAX_ITERATIONS} iterations.")
# calculate final errors and plot
if final_popt is None:
    print("Error: No successful fit achieved.")
else:
    sb_fit_full = gaussian_sb(bin_centers, *final_popt)
    mu_sb = final_popt[1]
    sigma_sb = final_popt[2]
    err_mu = np.sqrt(final_pcov[1,1])
    err_sigma = np.sqrt(final_pcov[2,2])

good_sigma = cleaned_df['sd_baseline'] < ((mu_sb +15*sigma_sb) * 1e-3)
cleaned_df = cleaned_df[good_sigma]

good_time = cleaned_df['total_time_above'] < 1.75e-8
cleaned_df = cleaned_df[good_time]

#export the cleaned data to another csv
cleaned_df.to_csv(file_path, sep=',', encoding='utf-8-sig', index=True, header=True)
#%% Import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import glob
import os
import matplotlib.cm as cm
from matplotlib.lines import Line2D

#%%
# --- MISSION PARAMETERS ---
mission_duration_years = 15  
mission_duration_seconds = mission_duration_years * 365.25 * 24 * 3600
# --------------------------

# Define the True Weibull function for log-space fitting
def weibull_fit_log(E, ln_k, E_0, b):
    return ln_k + (b - 1) * np.log(E) - (E / E_0)**b

# Target directory and file pattern
folder_path = r"C:\Users\wzt0020\Box\HERT_Box\Energy Resolution"
file_pattern = os.path.join(folder_path, "*_SEP_fluences.txt")
file_list = glob.glob(file_pattern)

if not file_list:
    print("No files matching '*_SEP_fluences.txt' were found.")

# Initialize the single combined figure BEFORE the loop
fig_combined, ax_combined = plt.subplots(figsize=(16, 10))

# Generate a continuous color map using viridis based on the number of files
colors = cm.viridis(np.linspace(0, 0.9, len(file_list)))

# Dictionary to store all calculated arrays, parameters, AND figure objects
compiled_fits = {}

# Loop through every matching file
for idx, file_path in enumerate(file_list):
    filename = os.path.basename(file_path)
    conf_interval = filename.split('_')[0]
    color = colors[idx]
    
    print(f"\n========================================================")
    print(f"Processing Confidence Interval: {conf_interval}%")
    print(f"========================================================")

    # 1. Read and parse the data
    data = {}
    all_x = []
    all_y = []
    current_curve = None
    skip_current = False

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            
            if line.startswith("# Curve"):
                current_curve = line.split(":", 1)[1].strip()
                if "Diff" not in current_curve or "ESP" not in current_curve:
                    skip_current = True
                else:
                    skip_current = False
                    data[current_curve] = {'energy': [], 'flux': []}
                    
            elif not line.startswith("#") and line and not skip_current:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        energy = float(parts[0])
                        raw_fluence = float(parts[1])
                        
                        if energy > 0 and raw_fluence > 0:
                            avg_flux = raw_fluence / mission_duration_seconds
                            
                            data[current_curve]['energy'].append(energy)
                            data[current_curve]['flux'].append(avg_flux)
                            all_x.append(energy)
                            all_y.append(avg_flux)
                    except ValueError:
                        pass

    # 2. Process data for fitting
    x = np.array(all_x)
    y = np.array(all_y)
    
    e_filter = (x > 52.29) & (x < 1000)
    x_filtered = x[e_filter]
    y_filtered = y[e_filter]
    
    if len(x_filtered) < 3:
        print(f"Not enough data points > 10 MeV to perform fit for {filename}. Skipping.")
        continue
        
    ln_y = np.log(y_filtered)

    # 3. Fit True Weibull
    p0 = [np.max(ln_y), 10.0, 0.5]
    popt, pcov = curve_fit(weibull_fit_log, x_filtered, ln_y, p0=p0, bounds=([-np.inf, 1e-5, 1e-5], [np.inf, np.inf, 5]))
    
    ln_k_opt, E_0_opt, b_opt = popt
    k_opt = np.exp(ln_k_opt)

    x_fit = np.logspace(np.log10(np.min(x_filtered)), np.log10(np.max(x_filtered)), 100)
    y_fit_weibull = k_opt * (x_fit**(b_opt - 1)) * np.exp(-(x_fit / E_0_opt)**b_opt)

# 4. Fit Exponential
    # Add cov=True to get the covariance matrix (V)
    p, V = np.polyfit(x_filtered, ln_y, 1, cov=True)
    slope = p[0]
    intercept = p[1]
    
    a = -1.0 / slope
    A = np.exp(intercept)
    y_fit_exp = A * np.exp(-x_fit / a)

    # Calculate 95% CI for the Exponential E_0
    # 1. Degrees of freedom = number of data points - number of fit parameters (2)
    dof = len(x_filtered) - 2
    
    # 2. Critical t-value for 95% CI (two-tailed means alpha/2 = 0.025, so we use 0.975)
    from scipy.stats import t
    t_val = t.ppf(0.975, dof)
    
    # 3. Standard error of the slope is the square root of the variance (diagonal element)
    se_slope = np.sqrt(V[0, 0])
    
    # 4. Margin of error for the slope
    moe_slope = t_val * se_slope
    
    # 5. Confidence interval bounds for the slope
    slope_lower = slope - moe_slope
    slope_upper = slope + moe_slope
    
    # 6. Propagate slope bounds to E_0 bounds
    # Because a = -1/slope, we calculate both and sort them to ensure correct [lower, upper] ordering
    a_bound_1 = -1.0 / slope_lower
    a_bound_2 = -1.0 / slope_upper
    
    a_ci_lower = min(a_bound_1, a_bound_2)
    a_ci_upper = max(a_bound_1, a_bound_2)

    # Print mathematical outputs to console
    print(f"Exponential E_0 (a): {a:.4f} MeV")
    print(f"  --> 95% CI for Exp E_0: [{a_ci_lower:.4f}, {a_ci_upper:.4f}] MeV")
    print(f"Weibull E_0: {E_0_opt:.4f} MeV")
    print(f"Weibull Equation: Avg Flux = {k_opt:.4e} * E^{b_opt-1:.4f} * exp(-(E / {E_0_opt:.4f})^{b_opt:.4f})")

    # ==========================================
    # 5A. INDIVIDUAL PLOT GENERATION
    # ==========================================
    fig_indiv, ax_indiv = plt.subplots(figsize=(12, 6))

    for curve_name, values in data.items():
        # Filtering original data to >10 MeV for the individual plot view
        e_arr = np.array(values['energy'])
        f_arr = np.array(values['flux'])
        mask = e_arr > 10 
        ax_indiv.plot(e_arr[mask], f_arr[mask], marker='o', linestyle='', alpha=0.4, label=curve_name)

    # Plot Fits on individual graph
    label_weibull_indiv = rf'Weibull Fit: $\Phi = {k_opt:.2e} \cdot E^{{{b_opt-1:.2f}}} \cdot \exp(-[E / {E_0_opt:.2f}]^{{{b_opt:.2f}}})$'
    ax_indiv.plot(x_fit, y_fit_weibull, color='red', linewidth=2, label=label_weibull_indiv)

    label_exp_indiv = rf'Exponential Fit: $\Phi = {A:.2e} \cdot e^{{-E / {a:.2f}}}$'
    ax_indiv.plot(x_fit, y_fit_exp, color='blue', linewidth=2, linestyle='--', label=label_exp_indiv)

    min_E = np.max([10, np.min(x_filtered)])

    # Individual Formatting
    ax_indiv.set_xscale('log')
    ax_indiv.set_yscale('log')
    ax_indiv.set_xlim(min_E, max(x_fit))
    ax_indiv.set_ylim(1e-6, 1e2) 
    ax_indiv.set_xlabel('Energy (MeV)')
    ax_indiv.set_ylabel('Average Flux (proton cm$^{-2}$ s$^{-1}$ (MeV/nuc)$^{-1}$)')
    ax_indiv.set_title(f'{conf_interval}% CI Protons Differential Flux Fits')
    ax_indiv.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax_indiv.grid(True, which="both", ls="--", alpha=0.5)
    
    fig_indiv.tight_layout()
    fig_indiv.savefig(f'protons_flux_weibull_{conf_interval}CI.png', dpi=300)
    
    # STORE ALL VARIABLES AND FIGURE OBJECTS FOR LATER MANIPULATION
    compiled_fits[conf_interval] = {
        'x_data': x_filtered,
        'y_data': y_filtered,
        'x_fit': x_fit,
        'weibull_params': {'k': k_opt, 'E_0': E_0_opt, 'b': b_opt},
        'exp_params': {'A': A, 'a': a}
    }

    # ==========================================
    # 5B. ADD TO COMBINED PLOT
    # ==========================================
    ax_combined.plot(x_filtered, y_filtered, marker='o', linestyle='', alpha=0.3, color=color)

    label_weibull_comb = f'{conf_interval}% CI Weibull'
    ax_combined.plot(x_fit, y_fit_weibull, color=color, linewidth=3, linestyle='-', label=label_weibull_comb)

    label_exp_comb = f'{conf_interval}% CI Exponential'
    ax_combined.plot(x_fit, y_fit_exp, color=color, linewidth=3, linestyle='--', label=label_exp_comb)


# ==========================================
# 6. FINALIZE AND SHOW COMBINED PLOT
# ==========================================
ax_combined.set_xscale('log')
ax_combined.set_yscale('log')
ax_combined.set_xlim(min_E, max(x_fit))
ax_combined.set_ylim(1e-6, 1e2) 
ax_combined.set_xlabel('Energy (MeV)', fontsize=16)
ax_combined.set_ylabel('Average Flux (cm$^{-2}$ sr$^{-1}$ s$^{-1}$ MeV$^{-1}$)', fontsize=16)
ax_combined.set_title('Combined Protons High Differential Flux Fits', fontsize=20)
ax_combined.tick_params(axis='both', which='major', labelsize=14)

# Scaled-up legend
leg = ax_combined.legend(
    bbox_to_anchor=(1.05, 1), 
    loc='upper left',
    fontsize=16,          
    markerscale=2.5,      
    labelspacing=1.0,     
    borderpad=1.0,        
    handlelength=3.0      
)

for line in leg.get_lines():
    line.set_linewidth(4.0)

ax_combined.grid(True, which="both", ls="--", alpha=0.5)
fig_combined.tight_layout()
plt.show()

#%% For each CI, find the energy width to reach a cm^2 sr of 0.2
bins = 400
energy_min = 10
energy_max = 1000
energy_edges = np.logspace(np.log10(energy_min), np.log10(energy_max), bins + 1)
energy_midpoints = (energy_edges[:-1] + energy_edges[1:]) / 2

generated_spectra = {}
energy_start = {}
energy_widths = {}  
for conf_interval_str, fit_data in compiled_fits.items():
    conf_interval = int(conf_interval_str)

    k = fit_data['weibull_params']['k']
    E_0 = fit_data['weibull_params']['E_0']
    b = fit_data['weibull_params']['b']
    
    # Apply the Weibull function to the midpoints array
    # Equation: Phi(E) = k * E^(b-1) * exp(-(E/E_0)^b)
    flux_spectrum = k * (energy_midpoints**(b - 1)) * np.exp(-(energy_midpoints / E_0)**b)
    
    # Store the result for later manipulation
    generated_spectra[conf_interval] = flux_spectrum

    G = 1/flux_spectrum * (energy_edges[1:] - energy_edges[:-1])  # Flux * Energy Bin Width

    energy_midpoints_temp = energy_midpoints.copy()  # Create a copy to avoid modifying the original energy_midpoints array
    G_temp = G.copy()  # Create a copy to avoid modifying the original G array
    energy_start[conf_interval_str] = []
    energy_widths[conf_interval_str] = []
    
    sum_G = 0.0
    start_idx = 0  # Anchor index for the current integration window

    for i in range(len(energy_midpoints)):
        sum_G += G[i]
        
        # When the threshold is met, record data and reset for the next cycle
        if sum_G >= 0.2:
            energy_width = energy_midpoints[i] - energy_midpoints[start_idx]
            energy_start[conf_interval_str].append((energy_edges[i] - energy_edges[start_idx])/2 + energy_edges[start_idx])
            energy_widths[conf_interval_str].append(energy_width)
            
            # Reset the accumulator
            sum_G = 0.0
            # Set the anchor index for the next cycle to the immediate next bin
            start_idx = i + 1

plt.figure(figsize=(10, 6))
for conf_interval_str in energy_start.keys():
    plt.scatter(energy_start[conf_interval_str], energy_widths[conf_interval_str], marker='o', label=f'{conf_interval_str}% CI')
plt.xlabel('Energy Start (MeV)')
plt.ylabel('Energy Width (MeV)')
plt.title('Energy Width to Reach 0.2 cm$^2$ sr vs. Confidence Interval')
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()
plt.xlim(0, 70)
plt.tight_layout()
plt.show()

#%% Export Weibull parameters to a text file for later use
import csv

# Define the output path (using the same directory structure you already have)
output_filepath = r"C:\Users\wzt0020\Box\HERT_Box\Energy Resolution\Weibull_Parameters_Export.txt"

# Open the text file and write the data
with open(output_filepath, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    
    # Write the header row so MATLAB knows what each column is
    writer.writerow(['CI', 'k', 'E_0', 'b'])
    
    # Loop through compiled_fits, sort them numerically by CI, and write to file
    for ci_str in sorted(compiled_fits.keys(), key=int):
        params = compiled_fits[ci_str]['weibull_params']
        writer.writerow([ci_str, params['k'], params['E_0'], params['b']])

print(f"Weibull parameters successfully exported to: {output_filepath}")

#%% Determine the energy width to reach 100 counts
# create flux spectra for each CI using the fitted parameters
# Read the comma-separated file directly into a NumPy array
total_geo = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\Proton_FS\proton_FS_14000_v5_range_totalGF.txt")
energy_midpoints_geo = total_geo[:, 0]  # Assuming the first column contains energy midpoints
total_geo_values = total_geo[:, 1]  # Assuming the second column contains the total

total_geo_pen = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\Proton_FS\proton_FS_14000_v5_pen_totalGF.txt")
energy_midpoints_geo_pen = total_geo_pen[:, 0]  # Assuming the first column contains energy midpoints
total_geo_values_pen = total_geo_pen[:, 1]  # Assuming the second column contains the total

generated_spectra = {}
count_rate = {}
for conf_interval, fit_data in compiled_fits.items():
    
    # Extract Weibull parameters
    k = fit_data['weibull_params']['k']
    E_0 = fit_data['weibull_params']['E_0']
    b = fit_data['weibull_params']['b']
    
    # 3. Apply the Weibull function to the midpoints array
    # Equation: Phi(E) = k * E^(b-1) * exp(-(E/E_0)^b)
    flux_spectrum = k * (energy_midpoints_geo**(b - 1)) * np.exp(-(energy_midpoints_geo / E_0)**b)
    
    # Store the result for later manipulation
    generated_spectra[conf_interval] = flux_spectrum

    # Calculate count rate for each confidence interval
    count_rate[conf_interval] = flux_spectrum * (total_geo_values[:] + total_geo_values_pen[:])  # Flux * Geometric Factor

# Plot count rates for each confidence interval
for conf_interval, rates in count_rate.items():
    plt.plot(energy_midpoints_geo, rates, label=f'{conf_interval}% CI')
plt.xscale('log')
plt.yscale('log')
plt.xlim(10, 1000)
plt.ylim(1e-4, 1e1)
plt.xlabel('Energy (MeV)')
plt.ylabel('Count Rate (counts/s/MeV)')
plt.title('Count Rate vs Energy for Each Confidence Interval')
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.show()

#%% Set the value of counts/sec we want to reach and determine energy width to reach it at low end
target_count_rate = 10/60

# Lists to store the numerical x and y data for plotting
numeric_cis = []
energy_widths = []

# Single consolidated loop for calculation, printing, and data storage
for conf_interval_str, rates in count_rate.items():
    # Convert the string key (e.g., '15') to a numerical integer
    conf_interval = int(conf_interval_str)
    
    # Filter energy range 10-70 MeV
    energy_filter = (energy_midpoints >= 10) & (energy_midpoints <= 70)
    max_rate = np.max(rates[energy_filter])
    half_max_rate = max_rate / 2
    
    # Find the starting energy (half of max count rate)
    above_half_max = (rates >= half_max_rate) & energy_filter
    if np.any(above_half_max):
        energy_at_half_max = energy_midpoints[above_half_max][0]
        
        # Define the integration window starting at half max
        energy_filter2 = (energy_midpoints >= energy_at_half_max) & (energy_midpoints <= 70)
        valid_midpoints = energy_midpoints[energy_filter2]
        valid_rates = rates[energy_filter2]
        
        sum_rate = 0
        i = 0
        while sum_rate < target_count_rate and i < len(valid_rates):
            sum_rate += valid_rates[i]
            i += 1

        # If we've reached the target count rate, calculate the energy width
        if sum_rate >= target_count_rate:
            energy_width = valid_midpoints[i] - valid_midpoints[0]
            
            # Output to console
            print(f"{conf_interval}% CI: Energy width to reach {target_count_rate:.4f} counts/s is approximately {energy_width:.4f} MeV")
            
            # Store numerical data for plotting
            numeric_cis.append(conf_interval)
            energy_widths.append(energy_width)

# Plotting the results using the numerical arrays
plt.figure(figsize=(10, 6))
plt.scatter(numeric_cis, energy_widths, marker='o', color='b', s=50)

# Formatting
plt.xlabel('Confidence Interval (%)')   
plt.ylabel('Energy Width (MeV)')
plt.title('Energy Width to Reach Target Count Rate vs. Confidence Interval')
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout()
plt.show()

# %% Plot Delta E vs Energy for each CI and determine the energy width to reach 0.2 cm^2 sr
target_count_rate = 6/60 # counts per second
target_G = 0.1 # cm^2 sr

generated_spectra = {}
delta_E = {}
for conf_interval_str, fit_data in compiled_fits.items():
    conf_interval = int(conf_interval_str)

    k = fit_data['weibull_params']['k']
    E_0 = fit_data['weibull_params']['E_0']
    b = fit_data['weibull_params']['b']
    
    # Apply the Weibull function to the midpoints array
    # Equation: Phi(E) = k * E^(b-1) * exp(-(E/E_0)^b)
    flux_spectrum = k * (energy_midpoints**(b - 1)) * np.exp(-(energy_midpoints / E_0)**b)
    
    # Store the result for later manipulation
    generated_spectra[conf_interval_str] = flux_spectrum

    delta_E[conf_interval_str] = target_count_rate / target_G / flux_spectrum

plt.figure(figsize=(10, 6))
for conf_interval_str in delta_E.keys():
    plt.scatter(energy_edges[:-1], delta_E[conf_interval_str], label=f'{conf_interval_str}% CI')
plt.xscale('log')
plt.yscale('log')
plt.xlim(10, 200)
plt.ylim(1e-2, 1e4)
plt.xlabel('Energy (MeV)')
plt.ylabel('Energy Width (MeV)')
plt.title('Energy Width to Reach 0.2 cm$^2$ sr')
plt.text(15, 1e3, f'Target Count Rate: {target_count_rate*60:.0f} counts/min\nTarget G: {target_G:.2f} cm$^2$ sr', fontsize=12, bbox=dict(facecolor='white', alpha=0.8))
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.show()


#%% Plot modified geometric factor by multiplying geometric factor by flux for each channel to visualize the contribution of each channel to the count rate
loaded_geo = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\Proton_FS\proton_FS_14000_v5_range_GFbyEC.txt")
energy_midpoints_geo = loaded_geo[:, 0]  # Assuming the first column contains energy midpoints
loaded_geo_values = loaded_geo[:, 1:]  # Assuming the rest of the columns contain geometric factor values for each channel

# FIX 1: Add delimiter, transpose the matrix to match (400x24), and share the midpoints array
loaded_geo_pen = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\Proton_FS\proton_FS_14000_v6_pen_GFbyEC.txt", delimiter=',')
loaded_geo_pen_values = loaded_geo_pen[:, 1:].T # Transpose the 24x400 GF data
energy_midpoints_geo_pen = energy_midpoints_geo # Safely reuse the exact same 400 energy bins

# Load in energy channels
energy_channels = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\channel_select\proton_channels_v5.txt",delimiter=',')

select_flux_CI = '50'
num_channels_range = loaded_geo_values.shape[1] # 31
num_channels_pen = loaded_geo_pen_values.shape[1] # 24
color_denom = max(1, num_channels_range - 1)

cmap_range = plt.get_cmap('plasma')
cmap_pen = plt.get_cmap('winter')

# ==========================================
# PLOT 1: RANGE CHANNELS ONLY
# ==========================================
plt.figure(figsize=(14.2, 8)) 

for idx in range(num_channels_range):
    color_val_range = cmap_range(idx / color_denom)
    
    # Format the center-aligned label using Figure Spaces
    e_low = energy_channels[idx, 0]
    e_high = energy_channels[idx, 1]
    label_str = f'{e_low:>5.1f} - {e_high:<5.1f}'.replace(' ', '\u2007')

    # Calculate and Plot Range Data
    geo_factor_temp = loaded_geo_values[:, idx]
    modified_geo = geo_factor_temp * generated_spectra[select_flux_CI]
    
    plt.plot(energy_midpoints_geo, modified_geo, 
             linestyle='-', color=color_val_range, label=label_str)

# Formatting Plot 1
plt.xscale('log')
plt.xticks(fontsize=12)
plt.yscale('log')
plt.yticks(fontsize=12)
plt.xlim(10, 1000)
plt.ylim(1e-6, 10**0.1)
plt.xlabel('Energy (MeV)', fontsize=14)
plt.ylabel('Modified Geometric Factor\n(counts/s/MeV)', fontsize=14) 
plt.title(f'Modified Geometric Factor by Multiplying with Flux ({select_flux_CI}% CI)', fontsize=16, pad=15)
plt.legend(title='   Energy (MeV)   ', loc='center left', bbox_to_anchor=(1.02, 0.5), 
           fontsize=11, title_fontsize=12, frameon=True)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout() 
plt.show()


# ==========================================
# PLOT 2: PENETRATING CHANNELS ONLY
# ==========================================
plt.figure(figsize=(14.2, 8)) 

# FIX 2: Decouple loop and map colors/labels back to the original indices
for back_idx in range(num_channels_pen):
    # Map the index to grab the correct color and text label
    if back_idx == 0:
        idx = 0
        e_low = energy_channels[0, 0]
        e_high = energy_channels[7, 1]
    else:
        idx = back_idx + 7
        e_low = energy_channels[idx, 0]
        e_high = energy_channels[idx, 1]
        
    color_val_pen = cmap_pen(idx / color_denom)
    label_str = f'{e_low:>5.1f} - {e_high:<5.1f}'.replace(' ', '\u2007')

    # Calculate and Plot Penetrating Data
    geo_factor_temp_pen = loaded_geo_pen_values[:, back_idx]
    modified_geo_pen = geo_factor_temp_pen * generated_spectra[select_flux_CI]
    
    plt.plot(energy_midpoints_geo_pen, modified_geo_pen, 
             linestyle='-', color=color_val_pen, label=label_str)

# Formatting Plot 2
plt.xscale('log')
plt.xticks(fontsize=12)
plt.yscale('log')
plt.yticks(fontsize=12)
plt.xlim(10, 1000)
plt.ylim(1e-6, 10**0.1)
plt.xlabel('Energy (MeV)', fontsize=14)
plt.ylabel('Modified Geometric Factor\n(counts/s/MeV)', fontsize=14) 
plt.title(f'Modified Penetrating GF by Multiplying with Flux ({select_flux_CI}% CI)', fontsize=16, pad=15)
plt.legend(title='   Energy (MeV)   ', loc='center left', bbox_to_anchor=(1.02, 0.5), 
           fontsize=11, title_fontsize=12, frameon=True)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout() 
plt.show()


# ==========================================
# PLOT 3: RANGE AND PENETRATING CHANNELS
# ==========================================
plt.figure(figsize=(16, 8)) 

col1_handles = [Line2D([], [], linestyle='none')]
col1_labels = ['    Range   ']  

col2_handles = [Line2D([], [], linestyle='none')]
col2_labels = ['  Energy (MeV)  ']

col3_handles = [Line2D([], [], linestyle='none')]
col3_labels = ['  Penetrating'] 

# Loop based on the full 31 channels so the table spacing is perfect
for idx in range(num_channels_range):
    color_val_range = cmap_range(idx / color_denom)

    # 1. Plot Range Data
    geo_factor_temp = loaded_geo_values[:, idx]
    modified_geo = geo_factor_temp * generated_spectra[select_flux_CI]
    line_r, = plt.plot(energy_midpoints_geo, modified_geo, linestyle='-', color=color_val_range)
    
    col1_handles.append(line_r)
    col1_labels.append('')
    
    # 2. Add Center Text
    e_low = energy_channels[idx, 0]
    e_high = energy_channels[idx, 1]
    label_str = f'{e_low:>5.2f} - {e_high:<5.2f}'.replace(' ', '\u2007')
    
    col2_handles.append(Line2D([], [], linestyle='none'))
    col2_labels.append(label_str)
    
    # 3. Plot Penetrating Data (using mapped index)
    if idx == 0:
        back_idx = 0
    elif idx > 7:
        back_idx = idx - 7
    else:
        back_idx = None # Channels 1-7 don't have their own individual penetrating line

    if back_idx is not None:
        color_val_pen = cmap_pen(idx / color_denom)
        geo_factor_temp_pen = loaded_geo_pen_values[:, back_idx]
        modified_geo_pen = geo_factor_temp_pen * generated_spectra[select_flux_CI]
        line_p, = plt.plot(energy_midpoints_geo_pen, modified_geo_pen, linestyle='-', color=color_val_pen)
        col3_handles.append(line_p)
    else:
        col3_handles.append(Line2D([], [], linestyle='none'))
        
    col3_labels.append('')

custom_handles = col1_handles + col2_handles + col3_handles
custom_labels = col1_labels + col2_labels + col3_labels

plt.xscale('log')
plt.xticks(fontsize=12)
plt.yscale('log')
plt.yticks(fontsize=12)
plt.xlim(10, 1000)
plt.ylim(1e-6, 10**0.1)

plt.xlabel('Energy (MeV)', fontsize=14)
plt.ylabel('Modified Geometric Factor\n(counts/s/MeV)', fontsize=14) 
plt.title(f'Modified Geometric Factor by Multiplying with Flux ({select_flux_CI}% CI)', fontsize=16, pad=15)

plt.legend(custom_handles, custom_labels, 
           ncol=3, 
           loc='center left', 
           bbox_to_anchor=(1.02, 0.5), 
           fontsize=11, 
           columnspacing=3.0,   
           handletextpad=-3.8,  
           frameon=True)

plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout() 
plt.show()

# %% Find count rate from multiplying flux by geometric factor and summing for each channel

# Arrays are already loaded properly in the earlier step, no need to reload them here!
select_flux_CI = '25'  

channel_count_rates = []
for idx in range(num_channels_range):
    geo_factor_temp = loaded_geo_values[:,idx]
    count_rate = np.sum(generated_spectra[select_flux_CI] * geo_factor_temp)
    channel_count_rates.append(count_rate)

channel_count_rates_pen = []
# FIX 3: Iterate based on the correct transposed column index
for back_idx in range(num_channels_pen):
    geo_factor_temp = loaded_geo_pen_values[:, back_idx] 
    count_rate_pen = np.sum(generated_spectra[select_flux_CI] * geo_factor_temp)
    channel_count_rates_pen.append(count_rate_pen)

# Plot the count rates for each channel
plt.figure(figsize=(7, 8))
plt.bar(np.linspace(1, len(channel_count_rates), len(channel_count_rates)), channel_count_rates, alpha=0.5, label='Range GF')

# X-axis mapping for penetrating bar chart so they align with their true channel equivalents
pen_x_axis = [1] + list(range(9, len(channel_count_rates_pen) + 8))
plt.bar(pen_x_axis, channel_count_rates_pen, alpha=0.5, label='Penetrating GF')

plt.text(0.98, 0.98, f'Confidence Interval: {select_flux_CI}%', transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='right', fontsize=12, bbox=dict(facecolor='white', alpha=0.8))
plt.xlabel('Channel Index', fontsize=14)
plt.ylabel('Count Rate (counts/s)', fontsize=14)
plt.title('Count Rates for Each Energy Channel', fontsize=16)
plt.legend()
plt.show()  

# Plot counts after time periods
time_periods = [1, 10, 60, 600, 3600]  # in seconds
time_in_minutes = [t/60 for t in time_periods]

channel_ids = list(range(1, len(channel_count_rates) + 1))
channel_pen_ids = pen_x_axis  # Reuse the mapped array calculated above

plt.figure(figsize=(10, 8)) 

for t_sec, t_min in zip(time_periods, time_in_minutes):
    counts_for_this_time = [rate * t_sec for rate in channel_count_rates]
    counts_for_this_time_pen = [rate * t_sec for rate in channel_count_rates_pen]
    
    if t_sec < 60:
        base_label = f'Time: {t_sec} sec'
    else:
        base_label = f'Time: {t_min:.0f} min'    
    
    line, = plt.plot(channel_ids, counts_for_this_time, marker='o', linestyle='-', label=f'{base_label} (Range)')
    plt.plot(channel_pen_ids, counts_for_this_time_pen, marker='x', linestyle='--', color=line.get_color(), label=f'{base_label} (Pen)')

plt.yscale('log')
plt.xlabel('Channel Index', fontsize=14)
plt.ylabel('Counts', fontsize=14)
plt.yticks(fontsize=12)
plt.ylim(1e-4, 1e4)  
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.title('Counts Across Energy Channels for Specific Time Durations', fontsize=16, pad=10)
plt.text(0.98, 0.98, f'Confidence Interval: {select_flux_CI}%', transform=plt.gca().transAxes, 
         verticalalignment='top', horizontalalignment='right', fontsize=12, bbox=dict(facecolor='white', alpha=0.8))
plt.legend(title='Duration & Type', loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize=12, title_fontsize=14)
plt.tight_layout()
plt.show()
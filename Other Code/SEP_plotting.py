import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import glob
import os
import matplotlib.cm as cm

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
    
    e_filter = (x > 10) & (x <70)
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
    slope, intercept = np.polyfit(x_filtered, ln_y, 1)
    a = -1.0 / slope
    A = np.exp(intercept)
    y_fit_exp = A * np.exp(-x_fit / a)

    # Print mathematical outputs to console
    print(f"Exponential E_0 (a): {a:.4f} MeV | Weibull E_0: {E_0_opt:.4f} MeV")
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
        

#%% Determine the energy width to reach 100 counts
# Set up energy bin system for protons like in matlab
bins = 1000
energy_min = 0.01
energy_max = 1000
# Calculate logarithmic bin edges
energy_edges = np.logspace(np.log10(energy_min), np.log10(energy_max), bins + 1)
energy_midpoints = (energy_edges[:-1] + energy_edges[1:]) / 2

# create flux spectra for each CI using the fitted parameters
file_path = "total_geo_proton_newtest4_nobackveto.txt"
# Read the comma-separated file directly into a NumPy array
total_geo = np.loadtxt(file_path, delimiter=',')

generated_spectra = {}
count_rate = {}
for conf_interval, fit_data in compiled_fits.items():
    
    # Extract Weibull parameters
    k = fit_data['weibull_params']['k']
    E_0 = fit_data['weibull_params']['E_0']
    b = fit_data['weibull_params']['b']
    
    # 3. Apply the Weibull function to the midpoints array
    # Equation: Phi(E) = k * E^(b-1) * exp(-(E/E_0)^b)
    flux_spectrum = k * (energy_midpoints**(b - 1)) * np.exp(-(energy_midpoints / E_0)**b)
    
    # Store the result for later manipulation
    generated_spectra[conf_interval] = flux_spectrum

    # Calculate count rate for each confidence interval
    count_rate[conf_interval] = flux_spectrum * total_geo[:] * (energy_edges[1:] - energy_edges[:-1])  # Flux * Geometric Factor * Energy Bin Width

# Plot count rates for each confidence interval
for conf_interval, rates in count_rate.items():
    plt.plot(energy_midpoints, rates, label=f'{conf_interval}% CI')
plt.xscale('log')
plt.yscale('log')
plt.xlim(10, 70)
plt.ylim(1e-4, 1e1)
plt.xlabel('Energy (MeV)')
plt.ylabel('Count Rate (counts/s)')
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
# %%

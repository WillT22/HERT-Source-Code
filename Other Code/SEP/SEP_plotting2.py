#%% Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.transforms as transforms

#%% Trapped Protons AP9 model
data_dir = r"C:\Users\wzt0020\Documents\OMERE 5.9"
files_ap9 = {
    'Mean': 'trappedProtons_mean.flx',
    '75th Percentile': 'trappedProtons_CORA_75.flx',
    '95th Percentile': 'trappedProtons_CORA_95.flx'
}

# Dictionaries to store raw AP9 data
ae9_raw_flux = {}
ae9_energy = None 

plt.figure(figsize=(12, 5))
for label, filename in files_ap9.items():
    filepath = os.path.join(data_dir, filename)
    try:
        data = np.loadtxt(filepath, comments='#')
        ae9_energy = data[:, 0]
        flux = data[:, 1] / (4 * np.pi) 
        
        ae9_raw_flux[label] = flux # Save to dictionary
        
        plt.loglog(ae9_energy, flux, marker='.', markersize=12, linestyle='-', linewidth=2, label=label)
    except Exception as e:
        print(f"Error processing {filename}: {e}")

plt.xscale('log')
plt.xlim(10, 1000)
custom_ticks = [10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 400, 500, 750, 1000]
plt.xticks(ticks=custom_ticks, labels=custom_ticks, fontsize=16, rotation=45)
plt.yscale('log')
plt.ylim(1e-4, 1e2)
plt.yticks(fontsize=16)
plt.xlabel('Energy (MeV)', fontsize=16)
plt.ylabel('Differential Flux (cm$^{-2}$ sr$^{-1}$ s$^{-1}$ MeV$^{-1}$)', fontsize=16, y=0.45)
plt.grid(True, which="major", ls="--", alpha=0.5)
plt.legend(fontsize=16)
plt.title('AP9 Trapped Protons Differential Flux', fontsize=20)
plt.tight_layout()
plt.show()

#%% SEP ESP model
files_esp = {
    'ESP Mean': 'solarMeanProtons_CORA_mean.flx',
    'ESP 75th Percentile': 'solarMeanProtons_CORA_75.flx', 
    'ESP 95th Percentile': 'solarMeanProtons_CORA_95.flx'
}

# Dictionaries to store raw ESP data
esp_raw_flux = {}
esp_energy = None

plt.figure(figsize=(12, 6))
for label, filename in files_esp.items():
    filepath = os.path.join(data_dir, filename)
    try:
        data = np.loadtxt(filepath, comments='#')
        esp_energy = data[:, 0]
        flux = data[:, 1] / (4 * np.pi) 
        
        esp_raw_flux[label] = flux # Save to dictionary
        
        plt.loglog(esp_energy, flux, marker='.', markersize=12, linestyle='-', linewidth=2, label=label)
    except Exception as e:
        print(f"Error processing {filename}: {e}")

plt.xscale('log')
plt.xlim(10, 300)
custom_ticks = [10, 20, 30, 50, 75, 100, 150, 200, 250, 300]
plt.xticks(ticks=custom_ticks, labels=custom_ticks, fontsize=16, rotation=45)
plt.yscale('log')
plt.ylim(1e-5, 1e1)
plt.yticks(fontsize=16)
plt.xlabel('Energy (MeV)', fontsize=16)
plt.ylabel('Differential Flux (cm$^{-2}$ sr$^{-1}$ s$^{-1}$ MeV$^{-1}$)', fontsize=16)
plt.grid(True, which="major", ls="--", alpha=0.5)
plt.legend(fontsize=16)
plt.title('ESP Solar Protons Differential Flux', fontsize=20)
plt.tight_layout()
plt.show()

#%% SEP SAPPHIRE Model
files_sapphire = {
    '1 per year': 'solarFlareProtons_GEO_50CI_pt5yr.flx',
    '1 per 4 years': 'solarFlareProtons_GEO_75CI_1yr.flx',
    '1 per 20 years': 'solarFlareProtons_GEO_95CI_1yr.flx'
}

# Dictionaries to store raw SAPPHIRE data
sapphire_raw_flux = {}
sapphire_energy = None

plt.figure(figsize=(12, 5))
for label, filename in files_sapphire.items():
    filepath = os.path.join(data_dir, filename)
    try:
        data = np.loadtxt(filepath, comments='#')
        sapphire_energy = data[:, 0]
        flux = data[:, 1] / (4 * np.pi)
        
        sapphire_raw_flux[label] = flux # Save to dictionary
        
        plt.loglog(sapphire_energy, flux, marker='.', markersize=12, linestyle='-', linewidth=2, label=label)
    except Exception as e:
        print(f"Error processing {filename}: {e}")

plt.xscale('log')
plt.xlim(10, 1000)
custom_ticks = [10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 400, 500, 750, 1000]
plt.xticks(ticks=custom_ticks, labels=custom_ticks, fontsize=16, rotation=45)
plt.yscale('log')
plt.ylim(1e-4, 1e4)
plt.yticks(fontsize=16)
plt.xlabel('Energy (MeV)', fontsize=16)
plt.ylabel('Differential Flux (cm$^{-2}$ sr$^{-1}$ s$^{-1}$ MeV$^{-1}$)', fontsize=16, y=0.45)
plt.grid(True, which="major", ls="--", alpha=0.5)
plt.legend(fontsize=16)
plt.title('SAPPHIRE Solar Protons Differential Flux', fontsize=20)
plt.tight_layout()
plt.show()

#%% Plot modified geometric factor
loaded_geo = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\Proton_FS\proton_FS_17000_v6comb_range_GFbyEC.txt")
energy_midpoints_geo = loaded_geo[:, 0]
loaded_geo_values = loaded_geo[:, 1:]

loaded_geo_pen = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\Proton_FS\proton_FS_17000_v6comb_pen_GFbyEC.txt")
loaded_geo_pen_values = loaded_geo_pen[:, 1:] 
energy_midpoints_geo_pen = energy_midpoints_geo 

energy_channels = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\channel_select\proton_channels_v6_combined.txt", delimiter=',')

num_channels_range = loaded_geo_values.shape[1] 
num_channels_pen = loaded_geo_pen_values.shape[1] 
color_denom = max(1, num_channels_range - 1)

cmap_range = plt.get_cmap('plasma')
cmap_pen = plt.get_cmap('winter')

plt.figure(figsize=(12, 8))
for i in range(num_channels_range):
    label_str = 'Range Channels (Solid)' if i == 0 else None
    plt.loglog(energy_midpoints_geo, loaded_geo_values[:, i], color=cmap_range(i / color_denom), linewidth=1.5, linestyle='-', label=label_str)

for i in range(num_channels_pen):
    label_str = 'Penetrating Channels (Dashed)' if i == 0 else None
    plt.loglog(energy_midpoints_geo_pen, loaded_geo_pen_values[:, i], color=cmap_pen(i / color_denom), linewidth=1.5, linestyle='-', alpha=0.8, label=label_str)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy (MeV)', fontsize=16)
plt.ylabel('Geometric Factor (cm$^2$ sr)', fontsize=16)
plt.title('HERT Range vs. Penetrating Geometric Factors', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend(fontsize=14, loc='best')
plt.tight_layout()
plt.show()

#%% Interpolate spectra for ALL models and percentiles
def log_interp(x_target, x_source, y_source):
    valid_mask = y_source > 0
    x_valid = x_source[valid_mask]
    y_valid = y_source[valid_mask]
    log_y_target = np.interp(np.log10(x_target), np.log10(x_valid), np.log10(y_valid), left=-np.inf, right=-np.inf)
    return 10 ** log_y_target

# Dictionaries to hold the final interpolated arrays
ae9_spectra = {}
sapphire_spectra = {}

# Interpolate AE9
for label, raw_flux in ae9_raw_flux.items():
    temp_interp = log_interp(energy_midpoints_geo, ae9_energy, raw_flux)
    temp_interp[energy_midpoints_geo < 10] = 0.0
    ae9_spectra[label] = temp_interp

# Interpolate SAPPHIRE
for label, raw_flux in sapphire_raw_flux.items():
    temp_interp = log_interp(energy_midpoints_geo, sapphire_energy, raw_flux)
    temp_interp[energy_midpoints_geo < 10] = 0.0
    sapphire_spectra[label] = temp_interp


# Plot ONLY the 'Mean' interpolated spectra for verification to keep the plot clean
plt.figure(figsize=(12, 6))

if 'Mean' in ae9_spectra:
    mask = ae9_spectra['Mean'] > 0
    plt.loglog(energy_midpoints_geo[mask], ae9_spectra['Mean'][mask], marker='.', markersize=12, linestyle='-', linewidth=2, label='AE9 Mean (Interp)')

if 'Mean' in sapphire_spectra:
    mask = sapphire_spectra['Mean'] > 0
    plt.loglog(energy_midpoints_geo[mask], sapphire_spectra['Mean'][mask], marker='.', markersize=12, linestyle='-', linewidth=2, label='SAPPHIRE Mean (Interp)')

plt.xscale('log')
plt.xlim(10, 1000)  
custom_ticks = [10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 400, 500, 750, 1000]
plt.xticks(ticks=custom_ticks, labels=custom_ticks, fontsize=16, rotation=45)
plt.yscale('log')
plt.ylim(1e-8, 1e2)
plt.xlabel('Energy (MeV)', fontsize=16)
plt.ylabel('Differential Flux (cm$^{-2}$ sr$^{-1}$ s$^{-1}$ MeV$^{-1}$)', fontsize=16)
plt.title('Log-Interpolated Mean Spectra for AE9, ESP, and SAPPHIRE', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True, which="major", ls="--", alpha=0.5)
plt.legend(fontsize=14, loc='best')
plt.tight_layout()
plt.show()

#%% AE9: Calculate channel count rate by multiplying the interpolated spectra by the geometric factors and summing over all energy bins for each channel
print("--- Generating Plots for AP9 (Trapped Protons) ---")
CI_select_ae9 = '95th Percentile'  # Change this to 'Mean', '75th Percentile', or '95th Percentile' as needed

# 1. Calculate count rates
ae9_channel_count_rates = []
for idx in range(num_channels_range):
    geo_factor_temp = loaded_geo_values[:, idx]
    count_rate = np.sum(ae9_spectra[CI_select_ae9] * geo_factor_temp)
    ae9_channel_count_rates.append(count_rate)

ae9_channel_count_rates_pen = []
for back_idx in range(num_channels_pen):
    geo_factor_temp = loaded_geo_pen_values[:, back_idx] 
    count_rate_pen = np.sum(ae9_spectra[CI_select_ae9] * geo_factor_temp)
    ae9_channel_count_rates_pen.append(count_rate_pen)

# Reverse the penetrating data so it plots from 26 down to 1
ae9_pen_reversed = ae9_channel_count_rates_pen[::-1]

# --- COLORMAP GENERATION ---
# Range: Plasma
cmap_range = plt.colormaps['plasma']
colors_range = cmap_range(np.linspace(0, 1, num_channels_range))

# Penetrating: Winter (Subset matching MATLAB logic)
cmap_pen = plt.colormaps['winter']
total_ec_len = len(energy_channels) # Assumes energy_channels is loaded earlier
full_winter = cmap_pen(np.linspace(0, 1, total_ec_len))
colors_pen = full_winter[-num_channels_pen:] # Take the last 'num_channels_pen' colors

# Reverse the penetrating colors to match the reversed bar placement
colors_pen_reversed = colors_pen[::-1]
# ---------------------------

# Set up the split x-axis with a 1-column gap
# Range is 1 to 26. We leave 27 empty. Penetrating starts at 28 and ends at 53.
range_x_axis = list(range(1, num_channels_range + 1))
pen_x_axis = list(range(num_channels_range + 2, num_channels_range + num_channels_pen + 2))

# Custom ticks correctly mapped to the new layout
tick_locs = [1, 6, 11, 16, 21, 26, 28, 33, 38, 43, 48, 53]
tick_labels = [1, 6, 11, 16, 21, 26, 26, 21, 16, 11, 6, 1]

# ---------------------------------------------------------
# Bar Chart (Count Rates)
# ---------------------------------------------------------
plt.figure(figsize=(12, 4.5))

# Pass the generated color arrays to the color argument
plt.bar(range_x_axis, ae9_channel_count_rates, color=colors_range)
plt.bar(pen_x_axis, ae9_pen_reversed, color=colors_pen_reversed)

# Dashed vertical dividing line placed exactly in the empty gap (x=27)
plt.axvline(x=num_channels_range + 1, color='black', linestyle='--', linewidth=2)

# Set up a blended transform (X is data coordinates, Y is axes coordinates 0 to 1)
trans = transforms.blended_transform_factory(plt.gca().transData, plt.gca().transAxes)

# Place text directly over the true centers of the respective bar groups (X=13.5 and X=40.5)
bbox_props = dict(facecolor='white', alpha=0.8, edgecolor='none')
plt.text(13.5, 0.95, 'Range Channels', transform=trans, ha='center', va='top', 
         fontsize=16, color='black', bbox=bbox_props)
plt.text(40.5, 0.95, 'Penetrating Channels', transform=trans, ha='center', va='top', 
         fontsize=16, color='black', bbox=bbox_props)

plt.text(0.5, 0.95, f'CI: {CI_select_ae9}', transform=plt.gca().transAxes, 
         verticalalignment='top', horizontalalignment='center', fontsize=16, bbox=dict(facecolor='white'))

plt.xlabel('Energy Channel Index', fontsize=16)
# Increased xlim to 54 to ensure the right edge of the 53rd bar isn't clipped
plt.xlim(0, num_channels_range + num_channels_pen + 2)
plt.ylabel('Count Rate (counts/s)', fontsize=16)

plt.yscale('log')
plt.ylim(1e-4, None) 

plt.xticks(ticks=tick_locs, labels=tick_labels, fontsize=16)
plt.yticks(fontsize=16)
plt.title('AP9 (Trapped Protons): Count Rates for Each Energy Channel', fontsize=20)
plt.tight_layout()
plt.show()  


# ---------------------------------------------------------
# Overlapping Bar Chart (1 Minute Duration)
# ---------------------------------------------------------
t_sec_1m = 60  # 1 minute in seconds
counts_1min_range = [rate * t_sec_1m for rate in ae9_channel_count_rates]
counts_1min_pen_reversed = [rate * t_sec_1m for rate in ae9_pen_reversed]

plt.figure(figsize=(12, 4.5))
# Pass the generated color arrays here as well
plt.bar(range_x_axis, counts_1min_range, color=colors_range)
plt.bar(pen_x_axis, counts_1min_pen_reversed, color=colors_pen_reversed)

plt.axvline(x=num_channels_range + 1, color='black', linestyle='--', linewidth=2)

# Use the same blended transform to label this plot
trans2 = transforms.blended_transform_factory(plt.gca().transData, plt.gca().transAxes)
plt.text(13.5, 0.95, 'Range Channels', transform=trans2, ha='center', va='top', 
         fontsize=16, color='black', bbox=bbox_props)
plt.text(40.5, 0.95, 'Penetrating Channels', transform=trans2, ha='center', va='top', 
         fontsize=16, color='black', bbox=bbox_props)

plt.text(0.5, 0.95, f'CI: {CI_select_ae9}\nDuration: 1 min', 
         transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='center', 
         fontsize=16, bbox=dict(facecolor='white'))

plt.xlabel('Energy Channel Index', fontsize=16)
plt.xlim(0, num_channels_range + num_channels_pen + 2)
plt.ylabel('Total Counts', fontsize=16)

plt.yscale('log')
plt.ylim(1e0, 1e5)

plt.xticks(ticks=tick_locs, labels=tick_labels, fontsize=16)
plt.yticks(fontsize=16)
plt.title('AP9 (Trapped Protons): Total Counts per Channel', fontsize=20)
plt.tight_layout()
plt.show()

#%% SAPPHIRE: Calculate channel count rate by multiplying the interpolated spectra by the geometric factors and summing over all energy bins for each channel
print("--- Generating Plots for SAPPHIRE (Solar Protons) ---")
freq_select = '1 per 4 years'  # Change this to '1 per year', '1 per 4 years', or '1 per 20 years' as needed

# 1. Calculate count rates
sapphire_channel_count_rates = []
for idx in range(num_channels_range):
    geo_factor_temp = loaded_geo_values[:, idx]
    count_rate = np.sum(sapphire_spectra[freq_select] * geo_factor_temp)
    sapphire_channel_count_rates.append(count_rate)

sapphire_channel_count_rates_pen = []
for back_idx in range(num_channels_pen):
    geo_factor_temp = loaded_geo_pen_values[:, back_idx] 
    count_rate_pen = np.sum(sapphire_spectra[freq_select] * geo_factor_temp)
    sapphire_channel_count_rates_pen.append(count_rate_pen)

# Reverse the penetrating data so it plots from 26 down to 1
sapphire_pen_reversed = sapphire_channel_count_rates_pen[::-1]

# --- COLORMAP GENERATION ---
# Range: Plasma
cmap_range = plt.colormaps['plasma']
colors_range = cmap_range(np.linspace(0, 1, num_channels_range))

# Penetrating: Winter (Subset matching MATLAB logic)
cmap_pen = plt.colormaps['winter']
total_ec_len = len(energy_channels) # Assumes energy_channels is loaded earlier
full_winter = cmap_pen(np.linspace(0, 1, total_ec_len))
colors_pen = full_winter[-num_channels_pen:] # Take the last 'num_channels_pen' colors

# Reverse the penetrating colors to match the reversed bar placement
colors_pen_reversed = colors_pen[::-1]
# ---------------------------

# Set up the split x-axis with a 1-column gap
# Range is 1 to 26. We leave 27 empty. Penetrating starts at 28 and ends at 53.
range_x_axis = list(range(1, num_channels_range + 1))
pen_x_axis = list(range(num_channels_range + 2, num_channels_range + num_channels_pen + 2))

# Custom ticks correctly mapped to the new layout
tick_locs = [1, 6, 11, 16, 21, 26, 28, 33, 38, 43, 48, 53]
tick_labels = [1, 6, 11, 16, 21, 26, 26, 21, 16, 11, 6, 1]

# ---------------------------------------------------------
# Bar Chart (Count Rates)
# ---------------------------------------------------------
plt.figure(figsize=(12, 4.5))

# Pass the generated color arrays to the color argument
plt.bar(range_x_axis, sapphire_channel_count_rates, color=colors_range)
plt.bar(pen_x_axis, sapphire_pen_reversed, color=colors_pen_reversed)

# Dashed vertical dividing line placed exactly in the empty gap (x=27)
plt.axvline(x=num_channels_range + 1, color='black', linestyle='--', linewidth=2)

# Set up a blended transform (X is data coordinates, Y is axes coordinates 0 to 1)
trans = transforms.blended_transform_factory(plt.gca().transData, plt.gca().transAxes)

# Place text directly over the true centers of the respective bar groups (X=13.5 and X=40.5)
bbox_props = dict(facecolor='white', alpha=0.8, edgecolor='none')
plt.text(13.5, 0.95, 'Range Channels', transform=trans, ha='center', va='top', 
         fontsize=16, color='black', bbox=bbox_props)
plt.text(40.5, 0.95, 'Penetrating Channels', transform=trans, ha='center', va='top', 
         fontsize=16, color='black', bbox=bbox_props)

plt.text(0.5, 0.95, f'Frequency: {freq_select}', transform=plt.gca().transAxes, 
         verticalalignment='top', horizontalalignment='center', fontsize=16, bbox=dict(facecolor='white'))

plt.xlabel('Energy Channel Index', fontsize=16)
# Increased xlim to 54 to ensure the right edge of the 53rd bar isn't clipped
plt.xlim(0, num_channels_range + num_channels_pen + 2)
plt.ylabel('Count Rate (counts/s)', fontsize=16)

plt.yscale('log')
plt.ylim(1e-6, None) # 95th Percentile

plt.xticks(ticks=tick_locs, labels=tick_labels, fontsize=16)
plt.yticks(fontsize=16)
plt.title('SAPPHIRE (Solar Protons): Count Rates for Each Energy Channel', fontsize=20)
plt.tight_layout()
plt.show()  


# ---------------------------------------------------------
# Overlapping Bar Chart (10 Minute Duration)
# ---------------------------------------------------------
t_sec_10m = 600  # 10 minutes in seconds
counts_10min_range = [rate * t_sec_10m for rate in sapphire_channel_count_rates]
counts_10min_pen_reversed = [rate * t_sec_10m for rate in sapphire_pen_reversed]

plt.figure(figsize=(12, 4.5))
# Pass the generated color arrays here as well
plt.bar(range_x_axis, counts_10min_range, color=colors_range)
plt.bar(pen_x_axis, counts_10min_pen_reversed, color=colors_pen_reversed)

plt.axvline(x=num_channels_range + 1, color='black', linestyle='--', linewidth=2)

# Use the same blended transform to label this plot
trans2 = transforms.blended_transform_factory(plt.gca().transData, plt.gca().transAxes)
plt.text(13.5, 0.95, 'Range Channels', transform=trans2, ha='center', va='top', 
         fontsize=16, color='black', bbox=bbox_props)
plt.text(40.5, 0.9, 'Penetrating Channels', transform=trans2, ha='center', va='top', 
         fontsize=16, color='black', bbox=bbox_props)

plt.text(0.5, 0.95, f'Frequency: {freq_select}\nDuration: 10 min', 
         transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='center', 
         fontsize=16, bbox=dict(facecolor='white'))

plt.xlabel('Energy Channel Index', fontsize=16)
plt.xlim(0, num_channels_range + num_channels_pen + 2)
plt.ylabel('Total Counts', fontsize=16)

plt.yscale('log')
plt.ylim(1e0, 1e10)

plt.xticks(ticks=tick_locs, labels=tick_labels, fontsize=16)
plt.yticks(fontsize=16)
plt.title('SAPPHIRE (Solar Protons): Total Counts per Channel', fontsize=20)
plt.tight_layout()
plt.show()




#%% Create-a-spectra using the Weibull fits for events from Laurenza 2015
# Define the True Weibull function for log-space fitting
def weibull_fit_log(E, ln_k, E_0, b):
    return ln_k + (b - 1) * np.log(E) - (E / E_0)**b

def solve_weibull_k(target_pfu, E_0, b):
    """
    Solves for k and ln_k given a target >10 MeV integral flux (pfu), E_0, and b.
    Uses analytical rearrangement of the Weibull integral to solve for k.
    """
    ln_k = np.log(target_pfu) + np.log(b) - b * np.log(E_0) + (10.0 / E_0)**b
    k = np.exp(ln_k)
    return k, ln_k

def classify_s_scale(pfu):
    """
    Classifies the event into the NOAA Solar Radiation Storm Scale based on peak >10 MeV pfu.
    """
    if pfu >= 100000:
        return "S5"
    elif pfu >= 10000:
        return "S4"
    elif pfu >= 1000:
        return "S3"
    elif pfu >= 100:
        return "S2"
    elif pfu >= 10:
        return "S1"
    else:
        return "Below S1"

events = {
    '26 December 2001': {'target_pfu': 779, 'E_0': 6.1, 'b': 0.89}, #S2 
    '21 March 2011':    {'target_pfu': 14,  'E_0': 0.3, 'b': 0.45} #S1
}

for event, params in events.items():
    k, ln_k = solve_weibull_k(params['target_pfu'], params['E_0'], params['b'])
    
    # Determine the NOAA S-class
    s_class = classify_s_scale(params['target_pfu'])
    
    # Include the parameters in the dictionary for use in the Weibull fit
    params['ln_k'] = ln_k
    params['k'] = k
    params['s_class'] = s_class
    
    print(f"{event}: {s_class} (pfu: {params['target_pfu']}), ln_k = {ln_k:.2f}, k = {k:.2e}")

weibull_spectra = {}
for event, params in events.items():
    ln_k = params['ln_k']
    E_0 = params['E_0']
    b = params['b']
    weibull_spectra[event] = np.exp(weibull_fit_log(energy_midpoints_geo, ln_k, E_0, b))

plt.figure(figsize=(12, 6))
for event, spectrum in weibull_spectra.items():
    plt.loglog(energy_midpoints_geo, spectrum, marker='.', markersize=12, linestyle='-', linewidth=2, label=event)
plt.xscale('log')
plt.xlim(10, 1000)
custom_ticks = [10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 400, 500, 750, 1000]
plt.xticks(ticks=custom_ticks, labels=custom_ticks, fontsize=16, rotation=45)
plt.yscale('log')
plt.ylim(1e-12, 1e4)
plt.xlabel('Energy (MeV)', fontsize=16)
plt.ylabel('Differential Flux (cm$^{-2}$ sr$^{-1}$ s$^{-1}$ MeV$^{-1}$)', fontsize=16)
plt.title('Weibull Spectra for Selected SEP Events', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True, which="major", ls="--", alpha=0.5)
plt.legend(fontsize=14, loc='best')


#%% SAPPHIRE: Calculate channel count rate by multiplying the interpolated spectra by the geometric factors and summing over all energy bins for each channel
print("--- Generating Plots for SAPPHIRE (Solar Protons) ---")
event_select = '21 March 2011'  # Change this to the desired event, '26 December 2001' or '21 March 2011'

# 1. Calculate count rates
sapphire_channel_count_rates = []
for idx in range(num_channels_range):
    geo_factor_temp = loaded_geo_values[:, idx]
    count_rate = np.sum(weibull_spectra[event_select] * geo_factor_temp)
    sapphire_channel_count_rates.append(count_rate)

sapphire_channel_count_rates_pen = []
for back_idx in range(num_channels_pen):
    geo_factor_temp = loaded_geo_pen_values[:, back_idx] 
    count_rate_pen = np.sum(weibull_spectra[event_select] * geo_factor_temp)
    sapphire_channel_count_rates_pen.append(count_rate_pen)

# Reverse the penetrating data so it plots from 26 down to 1
sapphire_pen_reversed = sapphire_channel_count_rates_pen[::-1]

# --- COLORMAP GENERATION ---
# Range: Plasma
cmap_range = plt.colormaps['plasma']
colors_range = cmap_range(np.linspace(0, 1, num_channels_range))

# Penetrating: Winter (Subset matching MATLAB logic)
cmap_pen = plt.colormaps['winter']
total_ec_len = len(energy_channels) # Assumes energy_channels is loaded earlier
full_winter = cmap_pen(np.linspace(0, 1, total_ec_len))
colors_pen = full_winter[-num_channels_pen:] # Take the last 'num_channels_pen' colors

# Reverse the penetrating colors to match the reversed bar placement
colors_pen_reversed = colors_pen[::-1]
# ---------------------------

# Set up the split x-axis with a 1-column gap
# Range is 1 to 26. We leave 27 empty. Penetrating starts at 28 and ends at 53.
range_x_axis = list(range(1, num_channels_range + 1))
pen_x_axis = list(range(num_channels_range + 2, num_channels_range + num_channels_pen + 2))

# Custom ticks correctly mapped to the new layout
tick_locs = [1, 6, 11, 16, 21, 26, 28, 33, 38, 43, 48, 53]
tick_labels = [1, 6, 11, 16, 21, 26, 26, 21, 16, 11, 6, 1]

# ---------------------------------------------------------
# Bar Chart (Count Rates)
# ---------------------------------------------------------
plt.figure(figsize=(12, 4.5))

# Pass the generated color arrays to the color argument
plt.bar(range_x_axis, sapphire_channel_count_rates, color=colors_range)
plt.bar(pen_x_axis, sapphire_pen_reversed, color=colors_pen_reversed)

# Dashed vertical dividing line placed exactly in the empty gap (x=27)
plt.axvline(x=num_channels_range + 1, color='black', linestyle='--', linewidth=2)

# Set up a blended transform (X is data coordinates, Y is axes coordinates 0 to 1)
trans = transforms.blended_transform_factory(plt.gca().transData, plt.gca().transAxes)

# Place text directly over the true centers of the respective bar groups (X=13.5 and X=40.5)
bbox_props = dict(facecolor='white', alpha=0.8, edgecolor='none')
plt.text(13.5, 0.95, 'Range Channels', transform=trans, ha='center', va='top', 
         fontsize=16, color='black', bbox=bbox_props)
plt.text(40.5, 0.95, 'Penetrating Channels', transform=trans, ha='center', va='top', 
         fontsize=16, color='black', bbox=bbox_props)

plt.text(0.5, 0.95, f'Event Class: {events[event_select]["s_class"]}\nCadence: 10 min', 
         transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='center', fontsize=16, bbox=dict(facecolor='white'))

plt.xlabel('Energy Channel Index', fontsize=16)
# Increased xlim to 54 to ensure the right edge of the 53rd bar isn't clipped
plt.xlim(0, num_channels_range + num_channels_pen + 2)
plt.ylabel('Count Rate (counts/s)', fontsize=16)

plt.yscale('log')
plt.ylim(1e-6, None) # 95th Percentile

plt.xticks(ticks=tick_locs, labels=tick_labels, fontsize=16)
plt.yticks(fontsize=16)
plt.title(f'SAPPHIRE (Solar Protons): {event_select}', fontsize=20)
plt.tight_layout()
plt.show()  


# ---------------------------------------------------------
# Overlapping Bar Chart (10 Minute Duration)
# ---------------------------------------------------------
t_sec_10m = 600  # 10 minutes in seconds
counts_10min_range = [rate * t_sec_10m for rate in sapphire_channel_count_rates]
counts_10min_pen_reversed = [rate * t_sec_10m for rate in sapphire_pen_reversed]

plt.figure(figsize=(12, 4.5))
# Pass the generated color arrays here as well
plt.bar(range_x_axis, counts_10min_range, color=colors_range)
plt.bar(pen_x_axis, counts_10min_pen_reversed, color=colors_pen_reversed)

plt.axvline(x=num_channels_range + 1, color='black', linestyle='--', linewidth=2)

# Use the same blended transform to label this plot
trans2 = transforms.blended_transform_factory(plt.gca().transData, plt.gca().transAxes)
plt.text(13.5, 0.95, 'Range Channels', transform=trans2, ha='center', va='top', 
         fontsize=16, color='black', bbox=bbox_props)
plt.text(40.5, 0.95, 'Penetrating Channels', transform=trans2, ha='center', va='top', 
         fontsize=16, color='black', bbox=bbox_props)

plt.text(0.5, 0.95, f'Event Class: {events[event_select]["s_class"]}\nCadence: 10 min', 
         transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='center', 
         fontsize=16, bbox=dict(facecolor='white'))

plt.xlabel('Energy Channel Index', fontsize=16)
plt.xlim(0, num_channels_range + num_channels_pen + 2)
plt.ylabel('Total Counts', fontsize=16)

plt.yscale('log')
plt.ylim(1e0, 1e4)

plt.xticks(ticks=tick_locs, labels=tick_labels, fontsize=16)
plt.yticks(fontsize=16)
plt.title(f'SAPPHIRE (Solar Protons): {event_select}', fontsize=20)
plt.tight_layout()
plt.show()

#%% Import necessary libraries
import glob
import os
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

#%% Set up the AE9 flux data based on Skyler's MATLAB script's hardcoded values
# The 17 discrete energy channels defined in the MATLAB script (MeV)
energies_AE9 = np.array([0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0])

# The hardcoded AE9 flux arrays (Counts/cm^2/s/sr/MeV or similar)
MaxFlux50_AE9 = np.array([3043556.52, 506402.542, 211356.954, 62748.4502, 15499.6257, 4167.08963, 1362.59182, 508.053908, 196.611127, 82.051654, 37.154501, 17.4685642, 8.97551526, 5.28536699, 2.81916956, 2.33, 1.84])

MaxFlux95_AE9 = np.array([12559884.1, 3917196.38, 2094895.39, 826481.191, 256040.313, 79982.653, 30444.5749, 13671.416, 6070.77165, 2795.68906, 1292.31374, 654.009611, 366.498171, 235.290299, 139.450819, 119.126969, 98.8031184])

# Create logarithmic interpolation functions (since flux drops exponentially with energy)
# Note: The MATLAB script used linear interpolation (interp1 without 'log'), but log-log interpolation is standard practice for AE9/AP9 data to prevent negative flux values between data points.
interp_25th = interp1d(np.log10(energies_AE9), np.log10(MaxFlux50_AE9), kind='linear', fill_value="extrapolate")
interp_50th = interp1d(np.log10(energies_AE9), np.log10(MaxFlux50_AE9), kind='linear', fill_value="extrapolate")
interp_95th = interp1d(np.log10(energies_AE9), np.log10(MaxFlux95_AE9), kind='linear', fill_value="extrapolate")

# Define the target energies you want to evaluate (e.g., matching your geometric factor arrays)
target_energies = np.linspace(0.5, 8.0, 100)

# Evaluate the interpolation and convert back from log space
flux_25th = 10**interp_25th(np.log10(target_energies))
flux_50th = 10**interp_50th(np.log10(target_energies))
flux_95th = 10**interp_95th(np.log10(target_energies))

# Plot the extracted spectra
plt.figure(figsize=(8, 6))
#plt.plot(target_energies, flux_25th, label='AE9 25th Percentile (Low Flux)') 
plt.plot(target_energies, flux_50th, label='AE9 50th Percentile (Nominal)')
plt.plot(target_energies, flux_95th, label='AE9 95th Percentile (Worst-Case)')
plt.yscale('log')
plt.xlabel('Energy (MeV)')
plt.ylabel('Flux')
plt.title('Extracted AE9 Spectra (Max Flux Condition)')
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.show()

#%% Import electron geometric factor
loaded_geo_range = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\Electron_FS\electron_FS_31000_v2_range_GFbyEC.txt")
energy_midpoints_geo_range = loaded_geo_range[:, 0]  # Assuming the first column contains energy midpoints
loaded_geo_range_values = loaded_geo_range[:, 1:]  # Assuming the rest of the columns contain geometric factor values for each channel

loaded_geo_pen = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\Electron_FS\electron_FS_31000_v2_pen_GFbyEC.txt")
energy_midpoints_geo_pen = loaded_geo_pen[:, 0]  # Assuming the first column contains energy midpoints
loaded_geo_pen_values = loaded_geo_pen[:, 1:]  # Assuming the rest of the columns contain geometric factor values for each channel

# Load in energy channels
energy_channels = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\channel_select\electron_channels_v2.txt",delimiter=',')

generated_spectra = {}
generated_spectra['25'] = (10**interp_25th(np.log10(energy_midpoints_geo_range)))/(4*np.pi)  # Convert to flux units /cm^2/sr/s/MeV
generated_spectra['50'] = (10**interp_50th(np.log10(energy_midpoints_geo_range)))/(4*np.pi)
generated_spectra['95'] = (10**interp_95th(np.log10(energy_midpoints_geo_range)))/(4*np.pi)

select_flux_CI = '50'
num_channels = loaded_geo_range_values.shape[1]
color_denom = max(1, num_channels - 1)

cmap_range = plt.get_cmap('plasma')
cmap_pen = plt.get_cmap('winter')

# ==========================================
# PLOT 1: RANGE CHANNELS ONLY
# ==========================================
plt.figure(figsize=(18.65, 10.5)) 

for idx in range(num_channels):
    color_val_range = cmap_range(idx / color_denom)
    
    # Format the center-aligned label using Figure Spaces
    e_low = energy_channels[idx, 0]
    e_high = energy_channels[idx, 1]
    label_str = f'{e_low:>5.3f} - {e_high:<5.3f}'.replace(' ', '\u2007')

    # Calculate and Plot Range Data
    geo_factor_temp = loaded_geo_range_values[:, idx]
    modified_geo = geo_factor_temp * generated_spectra[select_flux_CI]
    
    plt.plot(energy_midpoints_geo_range, modified_geo, 
             linestyle='-', color=color_val_range, label=label_str)

# Formatting Plot 1
plt.xscale('log')
plt.xticks(fontsize=12)
plt.yscale('log')
plt.yticks(fontsize=12)
plt.xlim(0.1, 10)
plt.ylim(1e-4, 10**6.1)
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
plt.figure(figsize=(18.65, 10.5)) 

for idx in range(num_channels):
    color_val_pen = cmap_pen(idx / color_denom)
    
    # Format the center-aligned label using Figure Spaces
    e_low = energy_channels[idx, 0]
    e_high = energy_channels[idx, 1]
    label_str = f'{e_low:>5.3f} - {e_high:<5.3f}'.replace(' ', '\u2007')

    # Calculate and Plot Penetrating Data
    geo_factor_temp_pen = loaded_geo_pen_values[:, idx]
    modified_geo_pen = geo_factor_temp_pen * generated_spectra[select_flux_CI]
    
    plt.plot(energy_midpoints_geo_pen, modified_geo_pen, 
             linestyle='-', color=color_val_pen, label=label_str)

# Formatting Plot 2
plt.xscale('log')
plt.xticks(fontsize=12)
plt.yscale('log')
plt.yticks(fontsize=12)
plt.xlim(0.1, 10)
plt.ylim(1e-4, 10**6.1)
plt.xlabel('Energy (MeV)', fontsize=14)
plt.ylabel('Modified Geometric Factor\n(counts/s/MeV)', fontsize=14) 
plt.title(f'Modified Geometric Factor by Multiplying with Flux ({select_flux_CI}% CI)', fontsize=16, pad=15)
plt.legend(title='   Energy (MeV)   ', loc='center left', bbox_to_anchor=(1.02, 0.5), 
           fontsize=11, title_fontsize=12, frameon=True)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout() 
plt.show()


# ==========================================
# PLOT 3: RANGE AND PENETRATING CHANNELS
# ==========================================

# 1. Increased figure height from 6 to 8 to better match the tall legend
plt.figure(figsize=(21, 10.5)) 

# 2. Initialize separate lists and add spaces to visually center the headers
col1_handles = [Line2D([], [], linestyle='none')]
col1_labels = ['     Range   ']  # Added spaces to center over the short line

col2_handles = [Line2D([], [], linestyle='none')]
col2_labels = ['   Energy (MeV)   ']

col3_handles = [Line2D([], [], linestyle='none')]
col3_labels = ['   Penetrating'] # Added a space to balance the longer word

for idx in range(num_channels):
    color_val_range = cmap_range(idx / color_denom)
    color_val_pen = cmap_pen(idx / color_denom)

    # Calculate Data
    geo_factor_temp = loaded_geo_range_values[:, idx]
    modified_geo = geo_factor_temp * generated_spectra[select_flux_CI]
    
    geo_factor_temp_pen = loaded_geo_pen_values[:, idx]
    modified_geo_pen = geo_factor_temp_pen * generated_spectra[select_flux_CI]
    
    # Plot Data
    line_r, = plt.plot(energy_midpoints_geo_range, modified_geo, 
                       linestyle='-', color=color_val_range)
    
    line_p, = plt.plot(energy_midpoints_geo_pen, modified_geo_pen, 
                       linestyle='-', color=color_val_pen)

    # Append to lists
    col1_handles.append(line_r)
    col1_labels.append('')
    
    e_low = energy_channels[idx, 0]
    e_high = energy_channels[idx, 1]
    
    # 2. Format both numbers to a fixed 5-character width.
    # The '>' right-aligns the first number. The '<' left-aligns the second number.
    label_str = f'{e_low:>5.3f} - {e_high:<5.3f}'
    
    # 3. Replace the standard spaces with Figure Spaces (\u2007)
    # This forces Matplotlib to treat the invisible padding exactly like visible digits.
    label_str = label_str.replace(' ', '\u2007')
    
    col2_handles.append(Line2D([], [], linestyle='none'))
    col2_labels.append(label_str)
    
    col3_handles.append(line_p)
    col3_labels.append('')

# Combine the columns sequentially
custom_handles = col1_handles + col2_handles + col3_handles
custom_labels = col1_labels + col2_labels + col3_labels

plt.xscale('log')
plt.xticks(fontsize=12)
plt.yscale('log')
plt.yticks(fontsize=12)
plt.xlim(0.1, 10)
plt.ylim(1e-4, 10**6.1)

plt.xlabel('Energy (MeV)', fontsize=14)
plt.ylabel('Modified Geometric Factor\n(counts/s/MeV)', fontsize=14) 
plt.title(f'Modified Geometric Factor by Multiplying with Flux ({select_flux_CI}% CI)', fontsize=16, pad=15)

# 3. Render the perfectly aligned Legend
plt.legend(custom_handles, custom_labels, 
           ncol=3, 
           loc='center left', 
           bbox_to_anchor=(1.02, 0.5), 
           fontsize=11, 
           columnspacing=3.0,   # Increased space between columns for better separation
           handletextpad=-3.8,   # CRITICAL: Pulls the text left to start exactly where the lines start
           frameon=True)

plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout() 
plt.show()

# %% Now do the same thing for electrons without the penetrating condition
loaded_geo = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\Electron_FS\electron_FS_31000_v2_GFbyEC.txt")
energy_midpoints_geo = loaded_geo[:, 0]  # Assuming the first column contains energy midpoints
loaded_geo_values = loaded_geo[:, 1:]  # Assuming the rest of the columns contain geometric factor values for each channel

# Load in energy channels
energy_channels = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\channel_select\electron_channels_v2.txt",delimiter=',')

generated_spectra = {}
generated_spectra['25'] = (10**interp_25th(np.log10(energy_midpoints_geo)))/(4*np.pi)
generated_spectra['50'] = (10**interp_50th(np.log10(energy_midpoints_geo)))/(4*np.pi)
generated_spectra['95'] = (10**interp_95th(np.log10(energy_midpoints_geo)))/(4*np.pi)

select_flux_CI = '50'
num_channels = loaded_geo_values.shape[1]
color_denom = max(1, num_channels - 1)

cmap_range = plt.get_cmap('plasma')
cmap_pen = plt.get_cmap('winter')

# ==========================================
# PLOT 1: RANGE CHANNELS ONLY
# ==========================================
plt.figure(figsize=(18.65, 10.5)) 

for idx in range(num_channels):
    color_val_range = cmap_range(idx / color_denom)
    
    # Format the center-aligned label using Figure Spaces
    e_low = energy_channels[idx, 0]
    e_high = energy_channels[idx, 1]
    label_str = f'{e_low:>5.3f} - {e_high:<5.3f}'.replace(' ', '\u2007')

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
plt.xlim(0.1, 10)
plt.ylim(1e-4, 10**6.1)
plt.xlabel('Energy (MeV)', fontsize=14)
plt.ylabel('Modified Geometric Factor\n(counts/s/MeV)', fontsize=14) 
plt.title(f'Modified Geometric Factor by Multiplying with Flux ({select_flux_CI}% CI)', fontsize=16, pad=15)
plt.legend(title='   Energy (MeV)   ', loc='center left', bbox_to_anchor=(1.02, 0.5), 
           fontsize=11, title_fontsize=12, frameon=True)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout() 
plt.show()

#%% Now using OMERE calculated AE9 fluxes
# Target directory and file pattern
folder_path = r"C:\Users\wzt0020\Box\HERT_Box\Energy Resolution"
file_pattern = os.path.join(folder_path, "AE9_GTO_*.txt")
file_list = glob.glob(file_pattern)

# Initialize a dictionary to store the spectra
omere_spectra = {}

# Loop through each file and extract the spectra
for file_path in file_list:
    # Extract the percentile from the filename (e.g., "AE9_GTO_25.txt" -> "25")
    # Using [-1] indexing is slightly safer in case there are multiple underscores
    percentile = os.path.basename(file_path).split('_')[-1].split('th.')[0]
    
    energies = []
    fluxes = []
    reading_diff = False
    
    with open(file_path, 'r') as file:
        for line in file:
            # Trigger the start of reading when Curve 0 is found
            if "Curve 0: Diff flux" in line:
                reading_diff = True
                continue
            
            # Trigger the end of reading when Curve 1 is found
            if "Curve 1: Integ flux" in line:
                break
            
            # If we are currently in the diff flux section and it's not a comment/empty line
            if reading_diff and not line.startswith('#') and line.strip():
                parts = line.split()
                if len(parts) == 2:
                    energies.append(float(parts[0]))
                    fluxes.append(float(parts[1]))
    fluxes = np.array(fluxes)/(4*np.pi)  # Convert to flux units /cm^2/sr/s/MeV
    
    # Convert lists to numpy arrays and store in the dictionary
    omere_spectra[percentile] = (np.array(energies), fluxes)  # Convert to flux units /cm^2/sr/s/MeV

# Plot the OMERE spectra
plt.figure(figsize=(8, 6))
for percentile, (energies, fluxes) in omere_spectra.items():
    plt.plot(energies, fluxes, label=f'OMERE AE9 {percentile}th Percentile')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy (MeV)', fontsize=12)
plt.ylabel('Flux ($cm^{-2}sr^{-1}s^{-1}MeV^{-1}$)', fontsize=12)
plt.title('OMERE AE9 Flux Spectra', fontsize=14)
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout()
plt.show()

#%% Plot OMERE spectra with the original AE9 spectra for comparison
plt.figure(figsize=(8, 6))
# Plot original AE9 spectra
plt.plot(target_energies, flux_50th, label='Original AE9 50th Percentile', linestyle='--')
# Plot OMERE spectra
for percentile, (energies, fluxes) in omere_spectra.items():
    plt.plot(energies, fluxes, label=f'OMERE AE9 {percentile}th Percentile')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy (MeV)', fontsize=12)
plt.ylabel('Flux ($cm^{-2}sr^{-1}s^{-1}MeV^{-1}$)', fontsize=12)
plt.title('Comparison of AE9 Flux Spectra', fontsize=14)
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout()
plt.show()

#%% Plot modified geometric factore from OMERE AE9 spectra
select_flux_CI = '75'
num_channels = loaded_geo_range_values.shape[1]
color_denom = max(1, num_channels - 1)

cmap_range = plt.get_cmap('plasma')
cmap_pen = plt.get_cmap('winter')

# ==========================================
# PLOT 1: RANGE CHANNELS ONLY
# ==========================================
plt.figure(figsize=(18.65, 10.5)) 

for idx in range(num_channels):
    color_val_range = cmap_range(idx / color_denom)
    
    # Format the center-aligned label using Figure Spaces
    e_low = energy_channels[idx, 0]
    e_high = energy_channels[idx, 1]
    label_str = f'{e_low:>5.3f} - {e_high:<5.3f}'.replace(' ', '\u2007')

    # Calculate and Plot Range Data
    geo_factor_temp = loaded_geo_range_values[:, idx]
    modified_geo = geo_factor_temp * omere_spectra[select_flux_CI][1]
    
    plt.plot(energy_midpoints_geo_range, modified_geo, 
             linestyle='-', color=color_val_range, label=label_str)

# Formatting Plot 1
plt.xscale('log')
plt.xticks(fontsize=12)
plt.yscale('log')
plt.yticks(fontsize=12)
plt.xlim(0.1, 10)
plt.ylim(1e-4, 10**4)
plt.xlabel('Energy (MeV)', fontsize=14)
plt.ylabel('Modified Geometric Factor\n(counts/s/MeV)', fontsize=14) 
plt.title(f'Modified Geometric Factor by Multiplying with Flux ({select_flux_CI}% Percentile)', fontsize=16, pad=15)
plt.legend(title='   Energy (MeV)   ', loc='center left', bbox_to_anchor=(1.02, 0.5), 
           fontsize=11, title_fontsize=12, frameon=True)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout() 
plt.show()


# ==========================================
# PLOT 2: PENETRATING CHANNELS ONLY
# ==========================================
plt.figure(figsize=(18.65, 10.5)) 

for idx in range(num_channels):
    color_val_pen = cmap_pen(idx / color_denom)
    
    # Format the center-aligned label using Figure Spaces
    e_low = energy_channels[idx, 0]
    e_high = energy_channels[idx, 1]
    label_str = f'{e_low:>5.3f} - {e_high:<5.3f}'.replace(' ', '\u2007')

    # Calculate and Plot Penetrating Data
    geo_factor_temp_pen = loaded_geo_pen_values[:, idx]
    modified_geo_pen = geo_factor_temp_pen * omere_spectra[select_flux_CI][1]
    
    plt.plot(energy_midpoints_geo_pen, modified_geo_pen, 
             linestyle='-', color=color_val_pen, label=label_str)

# Formatting Plot 2
plt.xscale('log')
plt.xticks(fontsize=12)
plt.yscale('log')
plt.yticks(fontsize=12)
plt.xlim(0.1, 10)
plt.ylim(1e-4, 10**4)
plt.xlabel('Energy (MeV)', fontsize=14)
plt.ylabel('Modified Geometric Factor\n(counts/s/MeV)', fontsize=14) 
plt.title(f'Modified Geometric Factor by Multiplying with Flux ({select_flux_CI}% Percentile)', fontsize=16, pad=15)
plt.legend(title='   Energy (MeV)   ', loc='center left', bbox_to_anchor=(1.02, 0.5), 
           fontsize=11, title_fontsize=12, frameon=True)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout() 
plt.show()


# ==========================================
# PLOT 3: RANGE AND PENETRATING CHANNELS
# ==========================================

# 1. Increased figure height from 6 to 8 to better match the tall legend
plt.figure(figsize=(21, 10.5)) 

# 2. Initialize separate lists and add spaces to visually center the headers
col1_handles = [Line2D([], [], linestyle='none')]
col1_labels = ['     Range   ']  # Added spaces to center over the short line

col2_handles = [Line2D([], [], linestyle='none')]
col2_labels = ['   Energy (MeV)   ']

col3_handles = [Line2D([], [], linestyle='none')]
col3_labels = ['   Penetrating'] # Added a space to balance the longer word

for idx in range(num_channels):
    color_val_range = cmap_range(idx / color_denom)
    color_val_pen = cmap_pen(idx / color_denom)

    # Calculate Data
    geo_factor_temp = loaded_geo_range_values[:, idx]
    modified_geo = geo_factor_temp * omere_spectra[select_flux_CI][1]
    
    geo_factor_temp_pen = loaded_geo_pen_values[:, idx]
    modified_geo_pen = geo_factor_temp_pen * omere_spectra[select_flux_CI][1]
    
    # Plot Data
    line_r, = plt.plot(energy_midpoints_geo_range, modified_geo, 
                       linestyle='-', color=color_val_range)
    
    line_p, = plt.plot(energy_midpoints_geo_pen, modified_geo_pen, 
                       linestyle='-', color=color_val_pen)

    # Append to lists
    col1_handles.append(line_r)
    col1_labels.append('')
    
    e_low = energy_channels[idx, 0]
    e_high = energy_channels[idx, 1]
    
    # 2. Format both numbers to a fixed 5-character width.
    # The '>' right-aligns the first number. The '<' left-aligns the second number.
    label_str = f'{e_low:>5.3f} - {e_high:<5.3f}'
    
    # 3. Replace the standard spaces with Figure Spaces (\u2007)
    # This forces Matplotlib to treat the invisible padding exactly like visible digits.
    label_str = label_str.replace(' ', '\u2007')
    
    col2_handles.append(Line2D([], [], linestyle='none'))
    col2_labels.append(label_str)
    
    col3_handles.append(line_p)
    col3_labels.append('')

# Combine the columns sequentially
custom_handles = col1_handles + col2_handles + col3_handles
custom_labels = col1_labels + col2_labels + col3_labels

plt.xscale('log')
plt.xticks(fontsize=12)
plt.yscale('log')
plt.yticks(fontsize=12)
plt.xlim(0.1, 10)
plt.ylim(1e-4, 10**4)

plt.xlabel('Energy (MeV)', fontsize=14)
plt.ylabel('Modified Geometric Factor\n(counts/s/MeV)', fontsize=14) 
plt.title(f'Modified Geometric Factor by Multiplying with Flux ({select_flux_CI}% Percentile)', fontsize=16, pad=15)

# 3. Render the perfectly aligned Legend
plt.legend(custom_handles, custom_labels, 
           ncol=3, 
           loc='center left', 
           bbox_to_anchor=(1.02, 0.5), 
           fontsize=11, 
           columnspacing=3.0,   # Increased space between columns for better separation
           handletextpad=-3.8,   # CRITICAL: Pulls the text left to start exactly where the lines start
           frameon=True)

plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout() 
plt.show()

# %% Now do the same thing for electrons without the penetrating condition
loaded_geo = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\Electron_FS\electron_FS_31000_v2_GFbyEC.txt")
energy_midpoints_geo = loaded_geo[:, 0]  # Assuming the first column contains energy midpoints
loaded_geo_values = loaded_geo[:, 1:]  # Assuming the rest of the columns contain geometric factor values for each channel

# Load in energy channels
energy_channels = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\channel_select\electron_channels_v2.txt",delimiter=',')

num_channels = loaded_geo_values.shape[1]
color_denom = max(1, num_channels - 1)

cmap_range = plt.get_cmap('plasma')
cmap_pen = plt.get_cmap('winter')

# ==========================================
# PLOT 1: RANGE CHANNELS ONLY
# ==========================================
plt.figure(figsize=(18.65, 10.5)) 

for idx in range(num_channels):
    color_val_range = cmap_range(idx / color_denom)
    
    # Format the center-aligned label using Figure Spaces
    e_low = energy_channels[idx, 0]
    e_high = energy_channels[idx, 1]
    label_str = f'{e_low:>5.3f} - {e_high:<5.3f}'.replace(' ', '\u2007')

    # Calculate and Plot Range Data
    geo_factor_temp = loaded_geo_values[:, idx]
    modified_geo = geo_factor_temp * omere_spectra[select_flux_CI][1]
    
    plt.plot(energy_midpoints_geo, modified_geo, 
             linestyle='-', color=color_val_range, label=label_str)

# Formatting Plot 1
plt.xscale('log')
plt.xticks(fontsize=12)
plt.yscale('log')
plt.yticks(fontsize=12)
plt.xlim(0.1, 10)
plt.ylim(1e-4, 10**4)
plt.xlabel('Energy (MeV)', fontsize=14)
plt.ylabel('Modified Geometric Factor\n(counts/s/MeV)', fontsize=14) 
plt.title(f'Modified Geometric Factor by Multiplying with Flux ({select_flux_CI}% Percentile)', fontsize=16, pad=15)
plt.legend(title='   Energy (MeV)   ', loc='center left', bbox_to_anchor=(1.02, 0.5), 
           fontsize=11, title_fontsize=12, frameon=True)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout() 
plt.show()
# %%

#%% Import necessary libraries
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

#%% Import electron geometric factor
loaded_geo_range = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\Proton_FS\proton_FS_14000_v5_range_GFbyEC.txt")
energy_midpoints_geo_range = loaded_geo_range[:, 0]  # Assuming the first column contains energy midpoints
loaded_geo_range_values = loaded_geo_range[:, 1:]  # Assuming the rest of the columns contain geometric factor values for each channel

loaded_geo_pen = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\Proton_FS\proton_FS_14000_v5_pen_GFbyEC.txt")
energy_midpoints_geo_pen = loaded_geo_pen[:, 0]  # Assuming the first column contains energy midpoints
loaded_geo_pen_values = loaded_geo_pen[:, 1:]  # Assuming the rest of the columns contain geometric factor values for each channel

# Load in energy channels
energy_channels = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\channel_select\proton_channels_v5.txt",delimiter=',')

select_flux_CI = '50'
num_channels = loaded_geo_range_values.shape[1]
color_denom = max(1, num_channels - 1)

cmap_range = plt.get_cmap('plasma')
cmap_pen = plt.get_cmap('winter')

#%% Now using OMERE calculated AP9 fluxes
# Target directory and file pattern
folder_path = r"C:\Users\wzt0020\Box\HERT_Box\Energy Resolution"
file_pattern = os.path.join(folder_path, "AP9_GTO_*.txt")
file_list = glob.glob(file_pattern)

# Initialize a dictionary to store the spectra
omere_spectra = {}

# Loop through each file and extract the spectra
for file_path in file_list:
    # Extract the percentile from the filename (e.g., "AP9_GTO_25.txt" -> "25")
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
    fluxes = np.array(fluxes)/(4*np.pi)  # Convert from omnidirectional to directional flux

    # Convert lists to numpy arrays and store in the dictionary
    omere_spectra[percentile] = (np.array(energies), fluxes)

# Plot the OMERE spectra
plt.figure(figsize=(8, 6))
for percentile, (energies, fluxes) in omere_spectra.items():
    plt.plot(energies, fluxes, label=f'OMERE AP9 {percentile}th Percentile')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy (MeV)', fontsize=12)
plt.ylabel('Flux ($cm^{-2}s^{-1}MeV^{-1}$)', fontsize=12)
plt.title('OMERE AP9 Flux Spectra', fontsize=14)
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout()
plt.show()

#%% Plot modified geometric factor from OMERE AP9 spectra
select_flux_CI = '25'
num_channels = loaded_geo_range_values.shape[1]
color_denom = max(1, num_channels - 1)

cmap_range = plt.get_cmap('plasma')
cmap_pen = plt.get_cmap('winter')

# ==========================================
# PLOT 1: RANGE CHANNELS ONLY
# ==========================================
plt.figure(figsize=(18.65, 8.5)) 

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
plt.xticks(fontsize=14)
plt.yscale('log')
plt.yticks(fontsize=14)
plt.xlim(10, 1000)
plt.ylim(1e-6, 10**1.1)
plt.xlabel('Energy (MeV)', fontsize=14)
plt.ylabel('Modified Geometric Factor\n(counts/s/MeV)', fontsize=14) 
plt.title(f'Modified Geometric Factor by Multiplying with AP9 Flux ({select_flux_CI}% Percentile)', fontsize=16, pad=15)
plt.legend(title='   Energy (MeV)   ', loc='center left', bbox_to_anchor=(1.02, 0.5), 
           fontsize=11, title_fontsize=12, frameon=True)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout() 
plt.show()


# ==========================================
# PLOT 2: PENETRATING CHANNELS ONLY
# ==========================================
plt.figure(figsize=(18.65, 8.5)) 

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
plt.xticks(fontsize=14)
plt.yscale('log')
plt.yticks(fontsize=14)
plt.xlim(10, 1000)
plt.ylim(1e-6, 10**1.1)
plt.xlabel('Energy (MeV)', fontsize=14)
plt.ylabel('Modified Geometric Factor\n(counts/s/MeV)', fontsize=14) 
plt.title(f'Modified Geometric Factor by Multiplying with AP9 Flux ({select_flux_CI}% Percentile)', fontsize=16, pad=15)
plt.legend(title='   Energy (MeV)   ', loc='center left', bbox_to_anchor=(1.02, 0.5), 
           fontsize=11, title_fontsize=12, frameon=True)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout() 
plt.show()


# ==========================================
# PLOT 3: RANGE AND PENETRATING CHANNELS
# ==========================================

# 1. Increased figure height from 6 to 8 to better match the tall legend
plt.figure(figsize=(21, 8.5)) 

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
plt.xticks(fontsize=14)
plt.yscale('log')
plt.yticks(fontsize=14)
plt.xlim(10, 1000)
plt.ylim(1e-6, 10**1.1)

plt.xlabel('Energy (MeV)', fontsize=14)
plt.ylabel('Modified Geometric Factor\n(counts/s/MeV)', fontsize=14) 
plt.title(f'Modified Geometric Factor by Multiplying with AP9 Flux ({select_flux_CI}% Percentile)', fontsize=16, pad=15)

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
# %%

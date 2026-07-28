#%% Import necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

#%% --- SEP Ion Plotting with OMERE Data ---
# --- 1. Constants & Conversion Factors ---
file_path = r"SEP_Lunar_50.txt"
mission_duration_years = 15

# Calculate total duration in seconds (accounting for leap years)
total_seconds = mission_duration_years * 365.25 * 24 * 3600
print(f"Converting using total mission duration: {total_seconds:,.0f} seconds")

# --- 2. OMERE Text File Parser ---
def load_and_convert_omere_sep(file_path):
    """
    Parses an OMERE SEP text file, separates the curves by Z-number and type 
    (Differential vs Integral), and converts Fluence to Average Flux.
    """
    data = {}
    current_key = None
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            
            # Detect a new curve block
            if line.startswith("# Curve"):
                # Determine if it is Differential or Integral
                is_differential = "Differential" in line
                type_str = "Differential" if is_differential else "Integral"
                
                # Extract the Z number (Ion species). 
                # Note: OMERE often omits Z=1 for protons in the header.
                z_match = re.search(r'Z=(\d+)', line)
                z_val = int(z_match.group(1)) if z_match else 1
                
                # Create a unique dictionary key for this curve
                current_key = f"Z{z_val}_{type_str}"
                data[current_key] = {'Energy_MeV_nuc': [], 'Fluence': []}
                
            # Skip empty lines, headers, or comments
            elif not line or line.startswith("#"):
                continue
                
            # Parse the numeric data lines
            else:
                parts = line.split()
                if len(parts) == 2 and current_key is not None:
                    data[current_key]['Energy_MeV_nuc'].append(float(parts[0]))
                    data[current_key]['Fluence'].append(float(parts[1]))
                    
    # Convert lists to Pandas DataFrames and calculate Flux
    dfs = {}
    for key, dic in data.items():
        df = pd.DataFrame(dic)
        
        # CONVERSION: Average Flux = Total Fluence / Total Seconds
        df['Flux'] = df['Fluence'] / total_seconds
        
        dfs[key] = df
        
    return dfs

# --- 3. Execute Parser ---
# Note: Ensure SEP_Lunar_50.txt is in your working directory or update file_path
sep_data = load_and_convert_omere_sep(file_path)

# Extract just the Differential Fluxes for Z=1 (Protons) and Z=2 (Helium)
df_H_diff = sep_data.get('Z1_Differential')
df_He_diff = sep_data.get('Z2_Differential')

# --- 4. Plot to Verify ---
fig, ax = plt.subplots(figsize=(10, 6))

color_H = 'orange'
color_He = 'purple'

if df_H_diff is not None:
    ax.plot(df_H_diff['Energy_MeV_nuc'], df_H_diff['Flux'], 
             label='Z=1 (Protons)', color=color_H, linewidth=2)
if df_He_diff is not None:
    ax.plot(df_He_diff['Energy_MeV_nuc'], df_He_diff['Flux'], 
             label='Z=2 (Helium)', color=color_He, linewidth=2)

# --- Step C: Add Vertical Lines and Unique Labels ---
vlines_H = {
    14.15: "Trigger D1",
    52.29: "Trigger D9",
    68.09: "Ta Pen",
    107.06: "Rear Pen"
}

vlines_He = {
    56.44/4: "Trigger D1",
    208.07/4: "Trigger D9",
    272.00/4: "Ta Pen",
    407.93/4: "Rear Pen"
}

# Create a transform that uses Data coordinates for X, but Axis percentages (0 to 1) for Y
blended_trans = ax.get_xaxis_transform()

# Plot Hydrogen Lines (Dashed, labels at the TOP left)
for energy, custom_label in vlines_H.items():
    ax.axvline(x=energy, color=color_H, linestyle='--', linewidth=1.5, alpha=0.8)
    
    # Offset X slightly to the right (multiply by 1.05 in log space)
    ax.text(energy * 1.05, 0.95, custom_label, transform=blended_trans, 
            color=color_H, rotation=90, verticalalignment='top', horizontalalignment='left', 
            fontsize=11, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=2))

# Plot Helium Lines (Dotted, labels at the BOTTOM right)
for energy, custom_label in vlines_He.items():
    ax.axvline(x=energy, color=color_He, linestyle=':', linewidth=1.5, alpha=0.8)
    
    # Offset X slightly to the left (multiply by 0.95 in log space)
    ax.text(energy * 0.95, 0.05, custom_label, transform=blended_trans, 
            color=color_He, rotation=90, verticalalignment='bottom', horizontalalignment='right', 
            fontsize=11, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=2))

# Format the plot according to space physics standards
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(10, 1000)
ax.set_xlabel('Energy (MeV/nuc)', fontsize=14)
ax.set_ylabel('Average Flux (ions / cm$^2$ / s / MeV/nuc)', fontsize=14)
ax.set_title('15-Year Time-Averaged SEP Differential Flux', fontsize=16)
ax.grid(True, which="both", ls="--", alpha=0.5)

# Ensure the legend doesn't cover up the new vertical lines
ax.legend(fontsize=12, loc='upper right')
plt.tight_layout()
plt.show()


# %% Calculate the energy at flux of 10^-2 for Z=1 and Z=2
def find_energy_at_flux(df, target_flux):
    """
    Finds the energy corresponding to a given flux level using log-log interpolation.
    """
    if df is None or df.empty:
        return None
    
    # Extract the numpy arrays
    energies = df['Energy_MeV_nuc'].values
    fluxes = df['Flux'].values
    
    if target_flux < fluxes.min() or target_flux > fluxes.max():
        return None  # Target flux is out of bounds
    
    # Step 1: Convert to Log10 space for accurate power-law interpolation
    log_energies = np.log10(energies)
    log_fluxes = np.log10(fluxes)
    log_target_flux = np.log10(target_flux)
    
    # Step 2: Sort the arrays so that the independent variable (Flux) is strictly increasing
    # This is a hard requirement for np.interp to function correctly
    sort_idx = np.argsort(log_fluxes)
    log_fluxes_sorted = log_fluxes[sort_idx]
    log_energies_sorted = log_energies[sort_idx]
    
    # Step 3: Interpolate
    log_energy_at_target = np.interp(log_target_flux, log_fluxes_sorted, log_energies_sorted)
    
    # Step 4: Convert back to linear space
    energy_at_target_flux = 10**(log_energy_at_target)
    
    return energy_at_target_flux

# Apply function
target_flux = 1e-2  # 10^-2 ions/cm^2/s/MeV/nuc

energy_H_at_target_flux = find_energy_at_flux(df_H_diff, target_flux)
energy_He_at_target_flux = find_energy_at_flux(df_He_diff, target_flux)

print(f"Energy at Flux = 10^-2 for Z=1 (Protons): {energy_H_at_target_flux:.2f} MeV/nuc")
print(f"Energy at Flux = 10^-2 for Z=2 (Helium): {energy_He_at_target_flux:.2f} MeV/nuc")  
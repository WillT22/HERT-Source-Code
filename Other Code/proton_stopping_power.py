import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Get PSTAR Data from NIST at https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html

def read_pstar_file_pandas(filepath):
    """
    Reads PSTAR proton stopping power/range table ASCII files into a 
    Pandas DataFrame. Skips initial header lines and uses the defined 
    column structure.
    
    The resulting DataFrame columns are:
    1. Kinetic_Energy_MeV
    2. Electron_Stp_Pow
    3. Nuclear_Stp_Pow
    4. Total_Stp_Pow
    5. CSDA_Range_gcm2
    6. Projected_Range_gcm2
    7. Detour_Factor
    """
    
    # The numerical data starts after 4 lines of header and 3 lines of blank space/titles.
    # The data starts at row index 7 (0-indexed).
    SKIP_HEADER_ROWS = 7 
    
    # Define descriptive column names for the DataFrame
    column_names = [
        'Kinetic_Energy_MeV', 'Electron_Stp_Pow', 'Nuclear_Stp_Pow', 
        'Total_Stp_Pow', 'CSDA_Range_gcm2', 'Projected_Range_gcm2', 'Detour_Factor'
    ]
    
    try:
        # Use pandas.read_csv to handle the fixed-width, space-delimited data.
        # We use delim_whitespace=True because the columns are separated by multiple spaces.
        data_df = pd.read_csv(
            filepath,
            sep='\s+',             # Use one or more spaces as the delimiter
            skiprows=SKIP_HEADER_ROWS,
            names=column_names,    # Assign the descriptive column headers
            engine='python'        # Use python engine for robustness with sep='\s+'
        )
        
        # Ensure the columns are numeric (in case of rogue header rows)
        data_df = data_df.apply(pd.to_numeric, errors='coerce')
        
        # Drop any rows that became entirely NaN after cleaning (e.g., if trailing header rows were missed)
        data_df = data_df.dropna(how='all')
        
        return data_df
        
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
        return None
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None

#%% Plot Stopping Power and Range for Silicon
material_name = "Silicon"
data_silicon = read_pstar_file_pandas("C:\\Users\\wzt0020\\Box\\HERT_Box\\Problem Particles\\proton_range_data_Si.txt")
kinetic_energy = data_silicon['Kinetic_Energy_MeV'].values

# 1. STOPPING POWER PLOT
plt.figure(figsize=(10, 6))

# Plot Total Stopping Power
plt.plot(kinetic_energy, data_silicon['Total_Stp_Pow'], 
            label='Total Stopping Power', linewidth=2, color='C3')

# Plot components for comparison (optional)
plt.plot(kinetic_energy, data_silicon['Electron_Stp_Pow'], 
            label='Electronic Component', linestyle='--', color='C0')

plt.xscale('log')
plt.yscale('log')

plt.title(f'Proton Stopping Power in {material_name} (PSTAR)', fontsize=14)
plt.xlabel('Kinetic Energy (MeV)', fontsize=12)
plt.ylabel(r'Stopping Power ($\text{MeV} \cdot \text{cm}^2/\text{g}$)', fontsize=12)
plt.grid(which='both', linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()


# 2. RANGE PLOT
plt.figure(figsize=(10, 6))

# Plot CSDA Range
plt.plot(kinetic_energy, data_silicon['CSDA_Range_gcm2'], 
            label='CSDA Range', linewidth=3, color='C2')

# Plot Projected Range for comparison
plt.plot(kinetic_energy, data_silicon['Projected_Range_gcm2'], 
            label='Projected Range', linestyle='--', color='C4')

plt.xscale('log')
plt.yscale('log')

plt.title(f'Proton Range in {material_name} (PSTAR)', fontsize=14)
plt.xlabel('Kinetic Energy (MeV)', fontsize=12)
plt.ylabel(r'Range ($\text{g}/\text{cm}^2$)', fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.grid(which='both', linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()

plt.show()

#%% Determine energy to penetrate different thicknesses of Si
density = 2.33  # g/cm^3 for Silicon
detector_thickness = 1.5e-1  # cm
gap_thickness = 1e-1  # cm
gapR_thickness = 1.08e-1  # cm
full_stack_thickness =      detector_thickness * 9 # cm
first_detector_diameter =   10e-1 *2  # cm
detector_diameter =         20e-1 *2  # cm

diag_thicknesses = np.sqrt(detector_diameter**2 + (detector_thickness * 8)**2) #cm
diag_through_first_max = np.sqrt((detector_diameter/2+first_detector_diameter/2)**2 + (full_stack_thickness)**2)  # cm
diag_through_first_min = np.sqrt((detector_diameter/2-first_detector_diameter/2)**2 + (detector_thickness * 7)**2)  # cm
diag_d5_d1_thickness = np.sqrt((detector_diameter/2-first_detector_diameter/2)**2 + (detector_thickness * 3)**2)  # cm

# 1) Full stack
required_range_full_stack = full_stack_thickness * density
energy_full_stack = np.interp(required_range_full_stack, 
                               data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate full stack ({full_stack_thickness:.2f} cm): {energy_full_stack:.2f} MeV")

# 2) Detector Diameter
required_range_detector = detector_diameter * density
energy_detector = np.interp(required_range_detector, 
                               data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate detector diameter ({detector_diameter:.2f} cm): {energy_detector:.2f} MeV")

# 3) First Detector Diameter
required_range_first_detector = first_detector_diameter * density
energy_first_detector = np.interp(required_range_first_detector, 
                               data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate first detector diameter ({first_detector_diameter:.2f} cm): {energy_first_detector:.2f} MeV")

# 4) Diagonal from D9 to D1 (max)
first_max_totalw = 9*detector_thickness + 4*gap_thickness + 4*gapR_thickness
diag_first_max_angle = np.arctan(first_max_totalw / (detector_diameter/2 + first_detector_diameter/2))
diag_through_first_max = 9*(detector_thickness/np.sin(diag_first_max_angle))  # cm
required_range_diag_first = diag_through_first_max * density
energy_diag_first = np.interp(required_range_diag_first, 
                               data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate diagonal through first detector (max) ({diag_through_first_max:.2f} cm): {energy_diag_first:.2f} MeV")

# 5) Diagonal from D9 to D1 (min)
first_min_totalw = 7*detector_thickness + 2*gap_thickness + 2*gapR_thickness
diag_first_min_angle = np.arctan(first_min_totalw / (detector_diameter/2 - first_detector_diameter/2))
diag_through_first_min = 7*(detector_thickness/np.sin(diag_first_min_angle))  # cm
required_range_diag_first_min = (diag_through_first_min * density + 2*data_silicon[data_silicon['Kinetic_Energy_MeV']==0.1]['Projected_Range_gcm2'].values[0]/density)
energy_diag_first_min = np.interp(required_range_diag_first_min, 
                               data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate diagonal through first detector (min) ({diag_through_first_min:.2f} cm): {energy_diag_first_min:.2f} MeV")

# 6) Diagonal from D5 to D1
diag_d5_d1_angle = np.arctan((4*detector_thickness+2*gap_thickness+2*gapR_thickness) / (detector_diameter/2 - first_detector_diameter/2))
diag_d5_d1_thickness = 3*(detector_thickness/np.sin(diag_d5_d1_angle))  # cm
required_range_diag_d5_d1 = (diag_d5_d1_thickness * density + 2*data_silicon[data_silicon['Kinetic_Energy_MeV']==0.1]['Projected_Range_gcm2'].values[0])
energy_diag_d5_d1 = np.interp(required_range_diag_d5_d1,
                               data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate diagonal from D5 to D1 ({diag_d5_d1_thickness:.2f} cm): {energy_diag_d5_d1:.2f} MeV")


#%% Plot Stopping Power and Range for Silicon
material_name = "Tungsten"
data_tungsten = read_pstar_file_pandas("C:\\Users\\wzt0020\\Box\\HERT_Box\\Problem Particles\\proton_range_data_W.txt")
kinetic_energy = data_tungsten['Kinetic_Energy_MeV'].values

# 1. STOPPING POWER PLOT
plt.figure(figsize=(10, 6))

# Plot Total Stopping Power
plt.plot(kinetic_energy, data_tungsten['Total_Stp_Pow'], 
            label='Total Stopping Power', linewidth=2, color='C3')

# Plot components for comparison (optional)
plt.plot(kinetic_energy, data_tungsten['Electron_Stp_Pow'], 
            label='Electronic Component', linestyle='--', color='C0')

plt.xscale('log')
plt.yscale('log')

plt.title(f'Proton Stopping Power in {material_name} (PSTAR)', fontsize=14)
plt.xlabel('Kinetic Energy (MeV)', fontsize=12)
plt.ylabel(r'Stopping Power ($\text{MeV} \cdot \text{cm}^2/\text{g}$)', fontsize=12)
plt.grid(which='both', linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()


# 2. RANGE PLOT
plt.figure(figsize=(10, 6))

# Plot CSDA Range
plt.plot(kinetic_energy, data_tungsten['CSDA_Range_gcm2'], 
            label='CSDA Range', linewidth=3, color='C2')

# Plot Projected Range for comparison
plt.plot(kinetic_energy, data_tungsten['Projected_Range_gcm2'], 
            label='Projected Range', linestyle='--', color='C4')

plt.xscale('log')
plt.yscale('log')

plt.title(f'Proton Range in {material_name} (PSTAR)', fontsize=14)
plt.xlabel('Kinetic Energy (MeV)', fontsize=12)
plt.ylabel(r'Range ($\text{g}/\text{cm}^2$)', fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.grid(which='both', linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()

plt.show()

#%% Determine energy to penetrate different thicknesses of W
density = 16.7  # g/cm^3 for 90% Tungsten, 10% Copper alloy
front_thickness = 1.5e-1 + 3.5e-1  # cm
rear_thickness = 5e-1  # cm

# 1) Front Layer
required_range_front = front_thickness * density
energy_front = np.interp(required_range_front, 
                               data_tungsten['Projected_Range_gcm2'], data_tungsten['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate front layer ({front_thickness:.2f} cm): {energy_front:.2f} MeV")

# 2) Rear Layer
required_range_rear = rear_thickness * density
energy_rear = np.interp(required_range_rear, 
                               data_tungsten['Projected_Range_gcm2'], data_tungsten['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate rear layer ({rear_thickness:.2f} cm): {energy_rear:.2f} MeV")
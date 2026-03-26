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
density_Si = 2.33  # g/cm^3 for Silicon
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
density_Si = 2.33  # g/cm^3 for Silicon
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
required_range_full_stack = full_stack_thickness * density_Si
energy_full_stack = np.interp(required_range_full_stack, 
                               data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate full stack ({full_stack_thickness:.2f} cm): {energy_full_stack:.2f} MeV")

# 2) Trigger D9
required_range_trigger_d9 = detector_thickness * 8 * density_Si
energy_trigger_d9 = np.interp(required_range_trigger_d9, 
                               data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate trigger D9 ({detector_thickness * 8:.2f} cm): {energy_trigger_d9+0.1:.2f} MeV")

# 3) Detector Diameter
required_range_detector = detector_diameter * density_Si
energy_detector = np.interp(required_range_detector, 
                               data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate detector diameter ({detector_diameter:.2f} cm): {energy_detector:.2f} MeV")

# 4) First Detector Diameter
required_range_first_detector = first_detector_diameter * density_Si
energy_first_detector = np.interp(required_range_first_detector, 
                               data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate first detector diameter ({first_detector_diameter:.2f} cm): {energy_first_detector:.2f} MeV")

# 5) Diagonal from D9 to D1 (max)
first_max_totalw = 9*detector_thickness + 4*gap_thickness + 4*gapR_thickness
diag_first_max_angle = np.arctan(first_max_totalw / (detector_diameter/2 + first_detector_diameter/2))
diag_through_first_max = 9*(detector_thickness/np.sin(diag_first_max_angle))  # cm
required_range_diag_first = diag_through_first_max * density_Si
energy_diag_first = np.interp(required_range_diag_first, 
                               data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate diagonal through first detector (max) ({diag_through_first_max:.2f} cm): {energy_diag_first:.2f} MeV")

# 6) Diagonal from D9 to D1 (min)
first_min_totalw = 7*detector_thickness + 2*gap_thickness + 2*gapR_thickness
diag_first_min_angle = np.arctan(first_min_totalw / (detector_diameter/2 - first_detector_diameter/2))
diag_through_first_min = 7*(detector_thickness/np.sin(diag_first_min_angle))  # cm
required_range_diag_first_min = (diag_through_first_min * density_Si + 2*data_silicon[data_silicon['Kinetic_Energy_MeV']==0.1]['Projected_Range_gcm2'].values[0]/density_Si)
energy_diag_first_min = np.interp(required_range_diag_first_min, 
                               data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate diagonal through first detector (min) ({diag_through_first_min:.2f} cm): {energy_diag_first_min:.2f} MeV")

# 7) Diagonal from D5 to D1
diag_d5_d1_angle = np.arctan((4*detector_thickness+2*gap_thickness+2*gapR_thickness) / (detector_diameter/2 - first_detector_diameter/2))
diag_d5_d1_thickness = 3*(detector_thickness/np.sin(diag_d5_d1_angle))  # cm
required_range_diag_d5_d1 = (diag_d5_d1_thickness * density_Si + 2*data_silicon[data_silicon['Kinetic_Energy_MeV']==0.1]['Projected_Range_gcm2'].values[0])
energy_diag_d5_d1 = np.interp(required_range_diag_d5_d1,
                               data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate diagonal from D5 to D1 ({diag_d5_d1_thickness:.2f} cm): {energy_diag_d5_d1:.2f} MeV")

#%% Plot Stopping Power and Range for Silicon (ELECTRON)
material_name = "Silicon"
density_Si = 2.33  # g/cm^3 for Silicon
data_silicon_electron = read_pstar_file_pandas("C:\\Users\\wzt0020\\Box\\HERT_Box\\Problem Particles\\electron_range_data_Si.txt")
kinetic_energy = data_silicon_electron['Kinetic_Energy_MeV'].values

# 1. STOPPING POWER PLOT
plt.figure(figsize=(10, 6))

# Plot Total Stopping Power
plt.plot(kinetic_energy, data_silicon_electron['Total_Stp_Pow'], 
            label='Total Stopping Power', linewidth=2, color='C3')

# Plot components for comparison (optional)
plt.plot(kinetic_energy, data_silicon_electron['Electron_Stp_Pow'], 
            label='Electronic Component', linestyle='--', color='C0')

plt.xscale('log')
plt.yscale('log')

plt.title(f'Electron Stopping Power in {material_name} (ESTAR)', fontsize=14)
plt.xlabel('Kinetic Energy (MeV)', fontsize=12)
plt.ylabel(r'Stopping Power ($\text{MeV} \cdot \text{cm}^2/\text{g}$)', fontsize=12)
plt.grid(which='both', linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()


# 2. RANGE PLOT
plt.figure(figsize=(10, 6))

# Plot CSDA Range
plt.plot(kinetic_energy, data_silicon_electron['CSDA_Range_gcm2'], 
            label='CSDA Range', linewidth=3, color='C2')

plt.xscale('log')
plt.yscale('log')

plt.title(f'Electron Range in {material_name} (ESTAR)', fontsize=14)
plt.xlabel('Kinetic Energy (MeV)', fontsize=12)
plt.ylabel(r'Range ($\text{g}/\text{cm}^2$)', fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.grid(which='both', linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()

plt.show()

#%% Determine energy to penetrate different thicknesses of Si
density_Si = 2.33  # g/cm^3 for Silicon
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

# 0) Full stack
required_range_one = detector_thickness * density_Si
energy_one = np.interp(required_range_one, 
                               data_silicon_electron['CSDA_Range_gcm2'], data_silicon_electron['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate one detector thickness ({detector_thickness:.2f} cm): {energy_one:.2f} MeV")

# 1) Full stack
required_range_full_stack = full_stack_thickness * density_Si
energy_full_stack = np.interp(required_range_full_stack, 
                               data_silicon_electron['CSDA_Range_gcm2'], data_silicon_electron['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate full stack ({full_stack_thickness:.2f} cm): {energy_full_stack:.2f} MeV")

# 2) Detector Diameter
required_range_detector = detector_diameter * density_Si
energy_detector = np.interp(required_range_detector, 
                               data_silicon_electron['CSDA_Range_gcm2'], data_silicon_electron['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate detector diameter ({detector_diameter:.2f} cm): {energy_detector:.2f} MeV")

# 3) First Detector Diameter
required_range_first_detector = first_detector_diameter * density_Si
energy_first_detector = np.interp(required_range_first_detector, 
                               data_silicon_electron['CSDA_Range_gcm2'], data_silicon_electron['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate first detector diameter ({first_detector_diameter:.2f} cm): {energy_first_detector:.2f} MeV")

# 4) Diagonal from D9 to D1 (max)
first_max_totalw = 9*detector_thickness + 4*gap_thickness + 4*gapR_thickness
diag_first_max_angle = np.arctan(first_max_totalw / (detector_diameter/2 + first_detector_diameter/2))
diag_through_first_max = 9*(detector_thickness/np.sin(diag_first_max_angle))  # cm
required_range_diag_first = diag_through_first_max * density_Si
energy_diag_first = np.interp(required_range_diag_first, 
                               data_silicon_electron['CSDA_Range_gcm2'], data_silicon_electron['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate diagonal through first detector (max) ({diag_through_first_max:.2f} cm): {energy_diag_first:.2f} MeV")

# 5) Diagonal from D9 to D1 (min)
first_min_totalw = 7*detector_thickness + 2*gap_thickness + 2*gapR_thickness
diag_first_min_angle = np.arctan(first_min_totalw / (detector_diameter/2 - first_detector_diameter/2))
diag_through_first_min = 7*(detector_thickness/np.sin(diag_first_min_angle))  # cm
required_range_diag_first_min = (diag_through_first_min * density_Si + 2*data_silicon_electron[data_silicon_electron['Kinetic_Energy_MeV']==0.1]['CSDA_Range_gcm2'].values[0]/density_Si)
energy_diag_first_min = np.interp(required_range_diag_first_min, 
                               data_silicon_electron['CSDA_Range_gcm2'], data_silicon_electron['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate diagonal through first detector (min) ({diag_through_first_min:.2f} cm): {energy_diag_first_min:.2f} MeV")

# 6) Diagonal from D5 to D1
diag_d5_d1_angle = np.arctan((4*detector_thickness+2*gap_thickness+2*gapR_thickness) / (detector_diameter/2 - first_detector_diameter/2))
diag_d5_d1_thickness = 3*(detector_thickness/np.sin(diag_d5_d1_angle))  # cm
required_range_diag_d5_d1 = (diag_d5_d1_thickness * density_Si + 2*data_silicon_electron[data_silicon_electron['Kinetic_Energy_MeV']==0.1]['CSDA_Range_gcm2'].values[0]/density_Si)
energy_diag_d5_d1 = np.interp(required_range_diag_d5_d1,
                               data_silicon_electron['CSDA_Range_gcm2'], data_silicon_electron['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate diagonal from D5 to D1 ({diag_d5_d1_thickness:.2f} cm): {energy_diag_d5_d1:.2f} MeV")


#%% Plot Stopping Power and Range for Tungsten
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
density_W = 16.7  # g/cm^3 for 90% Tungsten, 10% Copper alloy
front_thickness = 1.5e-1 + 3.5e-1  # cm
rear_thickness = 5e-1  # cm

# 1) Front Layer
required_range_front = front_thickness * density_W
energy_front = np.interp(required_range_front, 
                               data_tungsten['Projected_Range_gcm2'], data_tungsten['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate front layer ({front_thickness:.2f} cm): {energy_front:.2f} MeV")

# 2) Rear Layer
required_range_rear = rear_thickness * density_W
energy_rear = np.interp(required_range_rear, 
                               data_tungsten['Projected_Range_gcm2'], data_tungsten['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate rear layer ({rear_thickness:.2f} cm): {energy_rear:.2f} MeV")


#%% Plot Stopping Power and Range for Beryllium
material_name = "Beryllium"
data_beryllium = read_pstar_file_pandas("C:\\Users\\wzt0020\\Box\\HERT_Box\\Problem Particles\\proton_range_data_Be.txt")
kinetic_energy = data_beryllium['Kinetic_Energy_MeV'].values

# 1. STOPPING POWER PLOT
plt.figure(figsize=(10, 6))

# Plot Total Stopping Power
plt.plot(kinetic_energy, data_beryllium['Total_Stp_Pow'], 
            label='Total Stopping Power', linewidth=2, color='C3')

# Plot components for comparison (optional)
plt.plot(kinetic_energy, data_beryllium['Electron_Stp_Pow'], 
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

#%% Determine energy to penetrate different thicknesses of Be
density_Be = 1.85  # g/cm^3 for 100% Beryllium
window_thickness = 1.5e-1  # cm

# 1) Window
required_range_window = window_thickness * density_Be
energy_window = np.interp(required_range_window, 
                               data_beryllium['Projected_Range_gcm2'], data_beryllium['Kinetic_Energy_MeV'])
print(f"Energy required to penetrate Be Window ({window_thickness:.2f} cm): {energy_window:.2f} MeV")

#%% Plot Stopping Power and Range for Aluminum
material_name = "Aluminum"
data_aluminum = read_pstar_file_pandas("C:\\Users\\wzt0020\\Box\\HERT_Box\\Problem Particles\\proton_range_data_Al.txt")
kinetic_energy = data_aluminum['Kinetic_Energy_MeV'].values

# 1. STOPPING POWER PLOT
plt.figure(figsize=(10, 6))

# Plot Total Stopping Power
plt.plot(kinetic_energy, data_aluminum['Total_Stp_Pow'], 
            label='Total Stopping Power', linewidth=2, color='C3')

# Plot components for comparison (optional)
plt.plot(kinetic_energy, data_aluminum['Electron_Stp_Pow'], 
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
plt.plot(kinetic_energy, data_aluminum['CSDA_Range_gcm2'], 
            label='CSDA Range', linewidth=3, color='C2')

# Plot Projected Range for comparison
plt.plot(kinetic_energy, data_aluminum['Projected_Range_gcm2'], 
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

density_Al = 2.6989  # g/cm^3 for Aluminum

#%% Backward Calculation to find required incident energy to make it through
def calculate_incident_energy(E_exit, energy_array, stp_pow_array, density, thickness_cm, steps=1000):
    """
    Numerically integrates backward through a material to find the required incident energy.
    
    E_exit: Energy the proton has upon exiting the material (MeV)
    energy_array: Reference array of kinetic energies (MeV)
    stp_pow_array: Reference array of mass stopping power (MeV cm^2 / g)
    density: Material density (g / cm^3)
    thickness_cm: Total thickness of the material layer (cm)
    steps: Number of spatial slices for the integration
    """
    dx = thickness_cm / steps
    E_current = E_exit
    
    for _ in range(steps):
        # 1. Interpolate mass stopping power at the current energy
        mass_stp_pow = np.interp(E_current, energy_array, stp_pow_array)
        
        # 2. Convert to linear stopping power (MeV / cm)
        # S_linear = S_mass * density
        linear_stp_pow = mass_stp_pow * density
        
        # 3. Add the energy lost in this tiny step back to the total
        # dE = S_linear * dx
        E_current += linear_stp_pow * dx
        
    return E_current

#%% Find the amount of energy required to make it through the Be window and trigger D9
# If the proton stops exactly after depositing the trigger energy, 
# the energy entering the detector is simply the trigger energy.
E_entering_D9 = 0.1 # MeV 

# Step 1: Energy required to enter D1 and punch through all 8 prior detectors to reach D9
# 1a. Find the projected range in Si associated with the EXIT energy (E_entering_D9)
range_exit_Si = np.interp(E_entering_D9, data_silicon['Kinetic_Energy_MeV'], data_silicon['Projected_Range_gcm2'])

# 1b. Add the physical mass thickness of the 8 Silicon detectors (g/cm^2)
mass_thickness_D1_D8 = (detector_thickness * 8) * density_Si
range_incident_Si = range_exit_Si + mass_thickness_D1_D8

# 1c. Interpolate back to find the required incident energy entering D1
E_entering_D1 = np.interp(range_incident_Si, data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])

# Step 2: Energy required to hit the outside of the Beryllium window and reach D1
# 2a. Find the projected range in Be associated with the EXIT energy (E_entering_D1)
range_exit_Be = np.interp(E_entering_D1, data_beryllium['Kinetic_Energy_MeV'], data_beryllium['Projected_Range_gcm2'])

# 2b. Add the physical mass thickness of the Beryllium window (g/cm^2)
mass_thickness_Be = window_thickness * density_Be
range_incident_Be = range_exit_Be + mass_thickness_Be

# 2c. Interpolate back to find the required incident energy outside the window
E_incident_D9_trigger = np.interp(range_incident_Be, data_beryllium['Projected_Range_gcm2'], data_beryllium['Kinetic_Energy_MeV'])

print(f"Minimum incident energy required outside Be window to trigger D9: {E_incident_D9_trigger:.2f} MeV")


#%% Calculate energy required to make it through Be Window and trigger D1
# If the proton stops exactly after depositing the trigger energy, 
# the energy entering the detector is simply the trigger energy.
E_entering_D1_only = 0.1 # MeV

# Step 1: Energy required to hit the outside of the Beryllium window and reach D1
# 1a. Find the projected range in Be associated with the EXIT energy (E_entering_D1_only)
range_exit_Be_D1 = np.interp(E_entering_D1_only, data_beryllium['Kinetic_Energy_MeV'], data_beryllium['Projected_Range_gcm2'])

# 1b. Add the physical mass thickness of the Beryllium window (g/cm^2)
# (Reusing mass_thickness_Be from above)
range_incident_Be_D1 = range_exit_Be_D1 + mass_thickness_Be

# 1c. Interpolate back to find the required incident energy outside the window
E_incident_D1_trigger = np.interp(range_incident_Be_D1, data_beryllium['Projected_Range_gcm2'], data_beryllium['Kinetic_Energy_MeV'])

print(f"Minimum incident energy required outside Be window to trigger D1: {E_incident_D1_trigger:.2f} MeV")

#%% Calculate energy required to make it through tungsten later and trigger D1
# The minimum energy required to trigger D1
E_trigger_D1 = 0.1 # MeV

# Step 1: Energy required to enter D1 and deposit the trigger energy
# Using Range Chaining for Silicon
range_exit_Si = 0 # Stopping completely
range_incident_Si = range_exit_Si + (np.interp(E_trigger_D1, data_silicon['Kinetic_Energy_MeV'], data_silicon['Projected_Range_gcm2']))
E_entering_D1 = np.interp(range_incident_Si, data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])

# Step 2: Energy required to punch through the Tungsten layer and reach D1 with E_entering_D1
# 2a. Find the projected range in W associated with the EXIT energy
range_exit_W = np.interp(E_entering_D1, data_tungsten['Kinetic_Energy_MeV'], data_tungsten['Projected_Range_gcm2'])

# 2b. Add the physical mass thickness of the Tungsten layer (in g/cm^2)
mass_thickness_W = front_thickness * density_W
range_incident_W = range_exit_W + mass_thickness_W

# 2c. Interpolate back to find the required incident energy
E_incident_D1_W = np.interp(range_incident_W, data_tungsten['Projected_Range_gcm2'], data_tungsten['Kinetic_Energy_MeV'])

print(f"Minimum incident energy required outside W layer to trigger D1: {E_incident_D1_W:.2f} MeV")
# %% Calculate the energy required to penetrate through W and trigger D9 through D1
# The minimum energy required to trigger D1
E_trigger_D1 = 0.1 # MeV

# Step 1: Base residual range required to deposit the trigger energy in D1
range_incident_Si_D1 = np.interp(E_trigger_D1, data_silicon['Kinetic_Energy_MeV'], data_silicon['Projected_Range_gcm2'])

# Step 2: Energy required to make it through 8 detectors (D9 through D2) and reach D1
range_exit_Si_D9 = range_incident_Si_D1 + (detector_thickness * 8 * density_Si) 
E_entering_D9 = np.interp(range_exit_Si_D9, data_silicon['Projected_Range_gcm2'], data_silicon['Kinetic_Energy_MeV'])

# Step 3: Energy required to punch through the Tungsten layer and reach D9
range_exit_W_D9 = np.interp(E_entering_D9, data_tungsten['Kinetic_Energy_MeV'], data_tungsten['Projected_Range_gcm2'])
mass_thickness_W = rear_thickness * density_W

range_incident_W_D9 = range_exit_W_D9 + mass_thickness_W
E_incident_D9_W = np.interp(range_incident_W_D9, data_tungsten['Projected_Range_gcm2'], data_tungsten['Kinetic_Energy_MeV'])

# Step 4: Energy required to punch through aluminum layer and reach tungsten
Al_back_thickness = 1.5e-1  # cm
mass_thickness_Al_back = Al_back_thickness * density_Al

# FIXED: Convert the boundary energy (entering W) into an Aluminum range first
range_exit_Al_back = np.interp(E_incident_D9_W, data_aluminum['Kinetic_Energy_MeV'], data_aluminum['Projected_Range_gcm2'])

# Add the aluminum mass thickness to find the required incident range in Al
range_incident_Al_back_D9 = range_exit_Al_back + mass_thickness_Al_back

# Convert the incident range back to incident energy
E_incident_D9_Al_back = np.interp(range_incident_Al_back_D9, data_aluminum['Projected_Range_gcm2'], data_aluminum['Kinetic_Energy_MeV'])

print(f"Minimum incident energy required outside Al & W layers to trigger D9-D1: {E_incident_D9_Al_back:.2f} MeV")
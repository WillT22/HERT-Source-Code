#%% Import libraries
import numpy as np
import scipy.constants as sc
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#%% Setting up the problem
# Electron mass in MeV:  Conversion from kg*m^2/s^2 (Joules) to MeV.
e_E0 = sc.electron_mass * sc.c**2 / (sc.electron_volt * 1e6)

# Create electron kinetic energy range
KE = np.linspace(0, 2, num = 201)
KE = KE[1:]

# Relativistic momentum of the electron
p = np.sqrt((KE + e_E0)**2 - e_E0**2) / sc.c 


#%% PARTICLE SOURCE METHOD
# Create electron kinetic energy range
KE = np.linspace(0, 2, num = 201)
KE = KE[1:]

# Relativistic momentum of the electron
p = np.sqrt((KE + e_E0)**2 - e_E0**2) / sc.c 

# Assuming the object and image slits are limiters (not source or detector)
a_O = 5.5e-3       # object slit aperture (meters)
a_I = 8e-3       # image slit aperture (meters)
r = 191e-3 # radius of curverture in mm

d_O = 507.405e-3    # Distance from object slit to magnetic field (meters)
d_B = 299.6e-3      # Effective length in magnetic field (meters)
d_I = 602.693e-3    # Distance from magnetic field to image slit (meters)
d_D = 500e-3        # Distance from image slit to detector (meters)

# Using full slit width as boundary
r_min_abs = (d_B- 1/2 * min(a_O,a_I)) / (np.pi / 2)
r_max_abs = (d_B+ 1/2 * min(a_O,a_I)) / (np.pi / 2)

B_abs = p * (sc.electron_volt * 1e6) / (sc.elementary_charge * r)

# Reverse that: Calculate the max and min momentum that can be seen by the detector
p_min_abs = (B_abs * sc.elementary_charge * r_min_abs) / (sc.electron_volt * 1e6)
p_max_abs = (B_abs * sc.elementary_charge * r_max_abs) / (sc.electron_volt * 1e6)

# Find Kinetic Energy from momentum
KE_min_abs = np.sqrt((p_min_abs * sc.c)**2 + e_E0**2) - e_E0
KE_max_abs = np.sqrt((p_max_abs * sc.c)**2 + e_E0**2) - e_E0
Energy_resolution_abs = (KE_max_abs - KE_min_abs) / KE

'''
plt.figure(figsize=(8, 6)) 
plt.plot(KE, 1/Energy_resolution_abs, label=f"Slit = {min(a_O,a_I) * 10**3:.1f} mm")
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Resolving Power (E/$\Delta$E) ')
plt.grid(True)
plt.legend()
plt.show()
'''

plt.figure(figsize=(8, 6)) 
plt.plot(KE, Energy_resolution_abs, color='C1', label=f"Slit = {min(a_O,a_I) * 10**3:.1f} mm")
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Energy Resolution ($\Delta$E/E) ')
plt.grid(True)
plt.legend()
plt.show()

#%% Bring in Aerospace Data
Aero_data = {} # Dictionary for Aerospace data

Aero_data['csv_data'] = np.genfromtxt('C:/Users/wzt0020/Box/HERT_Box/Aerospace Testing/Resources from Aero/Aerospace Beta-ray Spectrometer 2025-04-30.csv', delimiter=',', filling_values=0, skip_header=1)
Aero_data['channel'] = Aero_data['csv_data'][:, 3]
Aero_data['FWHM'] = Aero_data['csv_data'][:, 4]
Aero_data['KE'] = Aero_data['csv_data'][:, 6]
Aero_data['time'] = Aero_data['csv_data'][:, 7]
Aero_data['max_counts'] = Aero_data['csv_data'][:, 10]
Aero_data['counts_per_day'] = Aero_data['csv_data'][:, 11]

# Fit a linear equation
coefficients = np.polyfit(Aero_data['channel'], Aero_data['KE'], 1)  # 1 for linear fit
slope = coefficients[0]
intercept = coefficients[1]

# Print the linear equation
print(f"Linear equation: KE (keV) = {slope:.4f} * Channel + {intercept:.4f}")

# Calculate the R^2 value
correlation_matrix = np.corrcoef(Aero_data['channel'], Aero_data['KE'])
correlation = correlation_matrix[0, 1]
r_squared = correlation**2

# Print the R^2 value
print(f"R^2 value: {r_squared:.4f}")

Aero_data['energy_resolution'] = Aero_data['FWHM']/Aero_data['channel']

#%% Attempting to find trend in spectrometer manual energy resolution
object_slit = np.array((0.5, 3, 3, 5, 5))
image_slit = np.array((0.4, 3.2, 4.6, 4.8, 5.9))
resolving_power = np.array((600, 83, 68, 58, 45))

# Define the exponential function
def exponential_function(x, a, b):
    return a * np.exp(b * x)

popt, pcov = curve_fit(exponential_function, image_slit, resolving_power)
a, b = popt  # Extract the fitted parameters
print(f"Fitted parameters: a = {a:.3f}, b = {b:.3f}")

# Generate points for the fitted curve
image_slit_curve = np.linspace(min(image_slit), max(image_slit), 100)  # More points for a smooth curve
resolving_power_curve = exponential_function(image_slit_curve, a, b)

plt.figure(figsize=(8, 6)) 
plt.scatter(object_slit, resolving_power, marker='.', s=80, label="Object Slit")
plt.scatter(image_slit, resolving_power, marker='.', s=80, label="Image Slit")
plt.plot(image_slit_curve, resolving_power_curve, label='Fitted Exponential Curve', color='black')
plt.ylabel('Resolving Power (E/$\Delta$E)')
plt.xlabel('Slit Width (mm)')
plt.grid(True)
plt.legend()
plt.show()

#%% Compare theoretical energy resolution with observations
KE_aero = Aero_data['KE']/1000
# Relativistic momentum of the electron
p_aero = np.sqrt((KE_aero + e_E0)**2 - e_E0**2) / sc.c 

# Using full slit width as boundary
r_min_aero = (d_B- 1/2 * min(a_O,a_I)) / (np.pi / 2)
r_max_aero = (d_B+ 1/2 * min(a_O,a_I)) / (np.pi / 2)

B_aero = p_aero * (sc.electron_volt * 1e6) / (sc.elementary_charge * r)

# Reverse that: Calculate the max and min momentum that can be seen by the detector
p_min_aero = (B_aero * sc.elementary_charge * r_min_aero) / (sc.electron_volt * 1e6)
p_max_aero = (B_aero * sc.elementary_charge * r_max_aero) / (sc.electron_volt * 1e6)

# Find Kinetic Energy from momentum
KE_min_aero = np.sqrt((p_min_aero * sc.c)**2 + e_E0**2) - e_E0
KE_max_aero = np.sqrt((p_max_aero * sc.c)**2 + e_E0**2) - e_E0
Energy_resolution_aero = (KE_max_aero - KE_min_aero) / KE_aero
'''
plt.figure(figsize=(8, 6)) 
plt.plot(KE_aero, 1/Energy_resolution_aero, label=f"Slit = {min(a_O,a_I) * 10**3:.1f} mm")
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Resolving Power (E/$\Delta$E)')
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(8, 6)) 
plt.plot(KE_aero, Energy_resolution_aero, label=f"Slit = {min(a_O,a_I) * 10**3:.1f} mm")
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Energy Resolution ($\Delta$E/E) ')
plt.grid(True)
plt.legend()
plt.show()
'''

manual_resolution = 1/exponential_function(min(a_O,a_I)*10**3,a,b)

# Plotting energy resolution v channel
fig, ax1 = plt.subplots(figsize=(8, 6))  # Adjust figure size as needed
ax1.plot(Aero_data['channel'], Aero_data['energy_resolution'], label="Aero")
ax1.plot(Aero_data['channel'], Energy_resolution_aero, label="Theory")
ax1.axhline(y=manual_resolution, color='gray', linestyle='--', label='Manual Resolution')
ax1.set_xlabel('Channel Number')
ax1.set_ylabel(r'Energy Resolution $\Delta E/E$')
ax1.set_title('Energy Resolution vs Channel Number')
ax1.grid(True)
# Get the limits from the first x-axis
ax1.set_xlim(0, 2000)
# Create the second x-axis
ax2 = ax1.twiny()  # Create a twin axis sharing the y-axis
ax2.set_xlim(intercept/1000, (slope * 2000 + intercept)/1000)
ax2.set_xlabel('Kinetic Energy (MeV)')  # Label the top x-axis
ax1.legend(loc='upper right') # Add legend
plt.show()

#%% Plotting FWHM
plt.figure(figsize=(8, 6)) 
plt.scatter(KE_aero*1000, slope * Aero_data['FWHM'], label='Aero')
plt.plot(KE_aero*1000, (KE_max_aero - KE_min_aero)*1000, color='C1', label='Theory')
plt.xlabel('Kinetic Energy (keV)')
plt.ylabel('FWHM (keV)')
plt.grid(True)
plt.legend()
plt.show()

#%% Determine image width required for specified energy resolution
def image_slit(KE_min_test, KE_max_test, KE_test):
    test_resolution = (KE_max_test-KE_min_test)/KE_test

    p_min_test = np.sqrt((KE_min_test + e_E0)**2 - e_E0**2)/sc.c
    p_max_test = np.sqrt((KE_max_test + e_E0)**2 - e_E0**2)/sc.c

    p_ke_test = np.sqrt((KE_test  + e_E0)**2 - e_E0**2) / sc.c 
    B_ke_test = p_ke_test * (sc.electron_volt * 1e6) / (sc.elementary_charge * r)

    # Rearranging to find r_min and r_max
    r_min_test = (p_min_test * sc.electron_volt * 1e6) / (sc.elementary_charge * B_ke_test)
    r_max_test = (p_max_test * sc.electron_volt * 1e6) / (sc.elementary_charge * B_ke_test)

    image_slit1 = (r_min_test * (np.pi / 2) - d_B)*-2
    image_slit2 = (r_max_test * (np.pi / 2) - d_B)*2
    image_slit_required = np.min([image_slit1, image_slit2])
    return image_slit_required, test_resolution


KE_min_channel = 1.0450 # channel 7 = 1.0450 MeV
KE_max_channel = 1.1172 # channel 9 = 1.1172 MeV
KE_test = 1.0805
image_slit_channel, test_resolution_channel = image_slit(KE_min_channel, KE_max_channel, KE_test)
print(f"Image slit required for {test_resolution_channel*100:.2f}% resolution at {KE_test:.2f} MeV: {image_slit_channel*1e3:.2f} mm")

KE_min_between = 1.0590 # channel 7 = 1.0450 MeV
KE_max_between = 1.1024 # channel 9 = 1.1172 MeV
KE_test = 1.0805
image_slit_between, test_resolution_between = image_slit(KE_min_between, KE_max_between, KE_test)
print(f"Image slit required for {test_resolution_between*100:.2f}% resolution at {KE_test:.2f} MeV: {image_slit_between*1e3:.2f} mm")

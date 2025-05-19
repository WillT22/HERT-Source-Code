#%% Import 
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

#%% Calculate B field using theory
KE = np.linspace(0, 3, num=301) # Kinetic Energy in MeV

# Electron rest energy in MeV
e_E0 = sc.electron_mass * sc.c**2 / (sc.electron_volt * 1e6)

# Calculating momentum in units of MeV/c
p_MeV_c = np.sqrt((KE + e_E0)**2 - e_E0**2) / sc.c

# Convert momentum to SI units (kg*m/s)
p_SI = p_MeV_c * (sc.electron_volt * 1e6)

r = 191e-3  # Radius of curvature in meters
q = sc.elementary_charge  # Elementary charge in Coulombs

# Calculate the magnetic field required to bend an electron beam around a set radius of curvature
B = p_SI / (q * r)
'''
# Plotting B vs KE
plt.figure(figsize=(8, 6))  # Adjust figure size as needed
plt.plot(KE, B*10000)
plt.xlim(0, 2)
plt.ylim(0, 450)
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Magnetic Field (G)')  # Label the y-axis with the correct unit
plt.title('Magnetic Field Strength vs Electron Kinetic Energy')
plt.grid(True)
plt.show()
'''

#%% Test theory based on manual suggestions
'''
test_KE = 0.03
test_mass = sc.proton_mass*76 * sc.c**2 / (sc.electron_volt * 1e6)
p_MeV_c_test = np.sqrt((test_KE + test_mass)**2 - test_mass**2) / sc.c
p_SI_test = p_MeV_c_test * 1.602176634e-13
B_test = p_SI_test / (q * r)
print(B_test)
'''

#%% Plot measured data from Aerospace
Aero_Bfield = {}
Aero_Bfield['csv_data'] = np.genfromtxt('C:/Users/Will/Box/HERT_Box/Aerospace Testing/Resources from Aero/Aerospace Beta-ray Spectrometer 2025-04-30.csv', delimiter=',', filling_values=0)
Aero_Bfield['KE'] = Aero_Bfield['csv_data'][:, 6]
Aero_Bfield['Bfield'] = Aero_Bfield['csv_data'][:, 1]
'''
# Plotting B vs KE
plt.figure(figsize=(8, 6)) 
plt.plot(Aero_Bfield['KE']/1000, Aero_Bfield['Bfield'])
plt.xlim(0, 2)
plt.ylim(0, 450)
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Magnetic Field (G)')
plt.title('Magnetic Field Strength vs Electron Kinetic Energy')
plt.grid(True)
plt.show()
'''

#%% Plot comparing theory and Aerospace measurement
# Plotting B vs KE for both datasets on the same figure
plt.figure(figsize=(8, 6))
plt.plot(KE, B*10000, label='Theory')
plt.plot(Aero_Bfield['KE'] / 1000, Aero_Bfield['Bfield'], label='Measured')  # Convert KE to MeV and B to G, label added
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Magnetic Field (G)')
plt.title('Magnetic Field Strength vs Electron Kinetic Energy')
plt.grid(True)
plt.xlim(0, 2)
plt.ylim(0, 450)
plt.legend()
plt.show()

# Plotting B vs KE for both datasets on the same figure
plt.figure(figsize=(8, 6))
plt.plot(B*10000, KE, label='Theory')
plt.plot(Aero_Bfield['Bfield'], Aero_Bfield['KE'] / 1000, label='Measured')  # Convert KE to MeV and B to G, label added
plt.xlabel('Magnetic Field (G)')
plt.ylabel('Kinetic Energy (MeV)')
plt.title('Electron Kinetic Energy vs Magnetic Field Strength')
plt.grid(True)
plt.xlim(0, 450)
plt.ylim(0, 2)
plt.legend()
plt.show()

#%% Find conversion 
# Calculating momentum in units of MeV/c
p_MeV_c_Aero = np.sqrt((Aero_Bfield['KE']/1000 + e_E0)**2 - e_E0**2) / sc.c

# Convert momentum to SI units (kg*m/s)
p_SI_Aero = p_MeV_c_Aero * (sc.electron_volt * 1e6)

# Calculate the magnetic field required to bend an electron beam around a set radius of curvature
B_Aero = p_SI_Aero / (q * r)

conversion = Aero_Bfield['Bfield'] / (B_Aero*10000)

plt.figure(figsize=(8, 6))
plt.plot(Aero_Bfield['KE'], conversion)
plt.xlabel('Kinetic Energy (keV)')
plt.ylabel('Conversion Factor')
plt.title('Magnetic Field Strength vs Electron Kinetic Energy')
plt.grid(True)
plt.show()
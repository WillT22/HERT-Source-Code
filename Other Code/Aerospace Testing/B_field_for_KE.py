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
Aero_Bfield['csv_data'] = np.genfromtxt('C:/Users/wzt0020/Box/HERT_Box/Aerospace Testing/Resources from Aero/Aerospace Beta-ray Spectrometer 2025-04-30.csv', delimiter=',', filling_values=0)
Aero_Bfield['KE'] = Aero_Bfield['csv_data'][:, 6]
Aero_Bfield['Bfield'] = Aero_Bfield['csv_data'][:, 1]
Aero_Bfield['Current'] = Aero_Bfield['csv_data'][:, 0]

KE_Aero = np.sqrt((Aero_Bfield['Bfield']/10000 * q * r / (sc.electron_volt * 1e6) * sc.c)**2 + e_E0**2) -  e_E0

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
# Calculating momentum in units of MeV/c
p_MeV_c_Aero = np.sqrt((Aero_Bfield['KE']/1000 + e_E0)**2 - e_E0**2) / sc.c
# Convert momentum to SI units (kg*m/s)
p_SI_Aero = p_MeV_c_Aero * (sc.electron_volt * 1e6)
# Calculate the magnetic field required to bend an electron beam around a set radius of curvature
B_Aero = p_SI_Aero / (q * r)

plt.figure(figsize=(8, 6))
plt.scatter(Aero_Bfield['KE'], B_Aero*10000, label='Theory')
plt.scatter(Aero_Bfield['KE'], Aero_Bfield['Bfield'], label='Measured')  # Convert KE to MeV and B to G, label added
plt.xlabel('Kinetic Energy (keV)')
plt.ylabel('Magnetic Field (G)')
plt.title('Magnetic Field Strength vs Electron Kinetic Energy')
plt.grid(True)
plt.xlim(0, 2000)
plt.ylim(0, 450)
plt.legend()
plt.show()

# Plotting B vs KE for both datasets on the same figure
plt.figure(figsize=(8, 6))
plt.scatter(Aero_Bfield['Bfield'], KE_Aero*1000, label='Theory')
plt.scatter(Aero_Bfield['Bfield'], Aero_Bfield['KE'], label='Measured')  # Convert KE to MeV and B to G, label added
plt.xlabel('Magnetic Field (G)')
plt.ylabel('Kinetic Energy (keV)')
plt.title('Electron Kinetic Energy vs Magnetic Field Strength')
plt.grid(True)
plt.xlim(0, 450)
plt.ylim(0, 2000)
plt.legend()
plt.show()

#%% Plot Current
plt.figure(figsize=(8, 6))
plt.scatter(Aero_Bfield['Current'], Aero_Bfield['KE'], label='Measured')  # Convert KE to MeV and B to G, label added
plt.xlabel('Current (A)')
plt.ylabel('Kinetic Energy (keV)')
plt.title('Electron Kinetic Energy vs Magnetic Field Strength')
plt.grid(True)
#plt.xlim(0, 450)
plt.ylim(0, 2000)
plt.legend()
plt.show()

#%% Find conversion 
conversion = Aero_Bfield['Bfield'] / (B_Aero*10000)
'''
plt.figure(figsize=(8, 6))
plt.plot(Aero_Bfield['KE'], conversion)
plt.xlabel('Kinetic Energy (keV)')
plt.ylabel('Conversion Factor')
plt.title('Magnetic Field Strength vs Electron Kinetic Energy')
plt.grid(True)
plt.show()
'''
x_data = Aero_Bfield['KE'][1:]
y_data = conversion[1:]
from scipy.optimize import curve_fit
def power_law(x, a, b, c):
    return a * x**-b + c
popt, pcov = curve_fit(power_law, x_data, y_data)
a_fit, b_fit, c_fit = popt

perr = np.sqrt(np.diag(pcov))
a_err, b_err, c_err = perr

print("Fitted parameters and their standard errors:")
print(f"a = {a_fit:.2f} +/- {a_err:.2f}")
print(f"b = {b_fit:.2f} +/- {b_err:.2f}")
print(f"c = {c_fit:.2f} +/- {c_err:.2f}")

plt.figure()
plt.plot(x_data, y_data, 'C1o', label="Original Data")
plt.plot(x_data, power_law(x_data, *popt), 'C1--', label="Fitted Curve")
plt.legend()
plt.xlabel('Kinetic Energy (keV)')
plt.ylabel('Conversion Factor [B/KE]')
plt.grid(True)
plt.show()

def get_fit_error(x_val, popt, pcov):
    a, b, c = popt
    # Calculate the values of the partial derivatives (Jacobian vector) at x_val
    # J = [df/da, df/db, df/dc]
    J = np.array([
        x_val**-b,
        a * -b * x_val**(-b-1),
        np.ones_like(x_val)
    ])
    # Calculate the variance of the fit
    # np.dot(J.T, np.dot(pcov, J)) is equivalent to J * pcov * J.T for a vector J.
    variance = np.dot(J.T, np.dot(pcov, J))
    # The error is the square root of the variance
    return np.sqrt(variance)

x_test = 1000
error_at_x = get_fit_error(x_test, popt, pcov)
fit_value_at_x = power_law(x_test, *popt)

# %% Select nominal energies 
selected_energies = np.array((1.04,1.06,1.08,1.10,1.12,.915,.655,.87,1.25,1.35))*1000  # in keV
b_theory = np.sqrt((selected_energies/1000 + e_E0)**2 - e_E0**2) / sc.c * (sc.electron_volt * 1e6) / (q * r)*1e4
b_predict = b_theory * power_law(selected_energies, *popt)
b_error =  b_theory * np.diag(get_fit_error(selected_energies, popt, pcov))

# Predict currents
current_predict = np.interp(b_predict, Aero_Bfield['Bfield'], Aero_Bfield['Current'])

plt.figure()
plt.plot(Aero_Bfield['KE'], Aero_Bfield['Bfield'], color='C1', label='Measured')
plt.errorbar(selected_energies, b_predict, yerr=b_error, color='black', marker='o',
             linestyle='None', capsize=5, label='Selected Energies', zorder=3)
plt.xlabel('Kinetic Energy (keV)')
plt.ylabel('Magnetic Field (G)')
plt.legend()
plt.grid(True)
plt.show()
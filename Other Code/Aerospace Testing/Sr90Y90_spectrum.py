#%% Import libraries
import numpy as np
import scipy.constants as sc
from scipy.special import gamma
import scipy.integrate as integrate
import matplotlib.pyplot as plt

#%% Fermi's Theory of Beta Decay
def fermi_theory(KE, Emax, Z, A, first_forbidden):
    """
    Calculates the Fermi theory prediction for the beta decay spectrum.

    Args:
        KE (numpy.ndarray): 1D array of kinetic energy range (MeV).
        Emax (float): Endpoint energy (MeV).  Maximum kinetic energy of the emitted beta particle.
        Z (int): Atomic number of the daughter nucleus.
        A (int): Atomic mass number of the daughter nucleus.
        first_forbidden (float): First forbidden shape factor coefficient (dimensionless, experimentally determined).
            Used to correct the spectrum shape for forbidden transitions.

    Returns:
        tuple:
            - N (numpy.ndarray): 1D array of the probability distribution as a function of KE.
              Has the same shape as KE.  Values are NaN where KE > Emax.
            - debug_vars (dict): A dictionary containing all local variables calculated within the function,
              useful for debugging and analysis.
    """
    # Electron mass in MeV:  Conversion from kg*m^2/s^2 (Joules) to MeV.
    e_E0 = sc.electron_mass * sc.c**2 / (sc.electron_volt * 1e6)

    N = np.zeros_like(KE)  # Initialize N with NaN values, same shape as KE.

    # Create a boolean mask for valid KE values.
    valid_ke_mask = (KE > 0) & (KE <= Emax)
    KE_valid = KE[valid_ke_mask]

    # Calculate values only for valid KE
    p = np.sqrt((KE_valid + e_E0)**2 - e_E0**2) / sc.c # Relativistic momentum of the electron.
    E = KE_valid + e_E0**2 # Relativistic energy of the electron.
    
    alpha = 1/137 # Fine-structure constant (dimensionless).
    S = np.sqrt(1 - alpha**2 * Z**2) # Relativistic correction factor.
    R = 1.2E-15 * A**(1/3) # Nuclear radius (meters).  Empirical formula.
    g = R / sc.hbar # dimensionless
    eta = (alpha * Z * E) / (p * sc.c) # Coulomb parameter (dimensionless).

    gamma_top = gamma(S + 1j * eta)
    gamma_bot = gamma(2 * S + 1)
    gamma_ratio = (gamma_top / gamma_bot)**2 # Fermi function related term.

    F = 2 * (1 + S) * (2 * p * g)**(2 * S - 2) * np.exp(sc.pi * eta) * (gamma_ratio * np.conj(gamma_ratio)) # Fermi function (dimensionless).
    #  Describes the Coulomb interaction between the emitted electron and the daughter nucleus.

    shape_factor = 1 + first_forbidden * E # Shape factor (dimensionless).
    #  Corrects the spectrum for first-forbidden transition

    N[valid_ke_mask] = p * E * (Emax - KE_valid)**2 * F * shape_factor
    #  Calculates the Fermi theory prediction for the beta spectrum.
    
    # Normalize N so that the integral over the valid KE range is 1
    integral_N = integrate.trapezoid(N[valid_ke_mask], KE_valid)  # Integrate only over valid KE
    N[valid_ke_mask] /= integral_N

    # Store local variables in a dictionary for debugging
    debug_vars = {'N': N, 'p': p, 'E': E, 'S': S, 'R': R, 'g': g, 'eta': eta, 
                  'gamma_ratio': gamma_ratio, 'Fermi': F, 'shape_factor': shape_factor, 'norm': integral_N}

    return N, debug_vars  # Return the dictionary

#%% Generate spectra
Beta_theory = {} # Dictionary to hold generated spectra data

# Create electron kinetic energy range
Beta_theory['KE'] = np.linspace(0, 3, num = 301)

# Find probability distribution for Sr90 decay
Beta_theory['Sr90'], Beta_theory['sr90_debug_vars'] = fermi_theory(Beta_theory['KE'], 0.5459, 38, 90, -54E-3)
# Find probability distribution for Y90 decay
Beta_theory['Y90'], Beta_theory['y90_debug_vars'] = fermi_theory(Beta_theory['KE'], 2.2785, 39, 90, -7.2E-3)

Beta_theory['combined_spectrum'] = Beta_theory['Sr90'] + Beta_theory['Y90']

# Normalize probabilities according to combined spectrum
Beta_theory['norm_combined'] = np.sum(Beta_theory['combined_spectrum'] * Beta_theory['KE'])

Beta_theory['Sr90_normalized'] = Beta_theory['Sr90'] / Beta_theory['norm_combined'] 
Beta_theory['Y90_normalized'] = Beta_theory['Y90'] / Beta_theory['norm_combined'] 
Beta_theory['combined_normalized'] = Beta_theory['combined_spectrum'] / Beta_theory['norm_combined'] 

Beta_theory['Sr90_mean'] = np.sum(Beta_theory['KE'] * Beta_theory['Sr90']) / np.sum(Beta_theory['Sr90'])
Beta_theory['Sr90_mean_diff'] = abs((Beta_theory['Sr90_mean'] - 0.196) / 0.196) * 100
Beta_theory['Y90_mean'] = np.sum(Beta_theory['KE'] * Beta_theory['Y90']) / np.sum(Beta_theory['Y90'])
Beta_theory['Y90_mean_diff'] = abs((Beta_theory['Y90_mean'] - 0.933) / 0.933) * 100

# Plot Sr90 and Y90
fig, ax = plt.subplots(figsize=(7, 5))

Sr90_scatter = ax.scatter(Beta_theory['KE'], Beta_theory['Sr90'], label="Sr90", color='C2', s=8)
Y90_scatter = ax.scatter(Beta_theory['KE'], Beta_theory['Y90'], label="Y90", color='C0', s=8)

# Plot the sum, making it more prominent
combined_scatter = ax.scatter(Beta_theory['KE'], Beta_theory['combined_spectrum'], label="Combined", s=8, color='C1')

ax.set_xlim(0, 2.3)
ax.set_xlabel("Kinetic Energy (MeV)")
ax.set_ylabel("Probability Distribution")
ax.set_title("Beta Decay Spectra of Sr90 and Y90")
ax.grid(True)

# Add vertical lines for mean energies
plt.axvline(x=0.196, color='C2', linestyle='--', label='Sr90 Mean: 0.196 MeV')
plt.axvline(x=0.933, color='C0', linestyle='--', label='Y90 Mean: 0.933 MeV')
ax.legend()

plt.show()

#%% Import LASP spectra
LASP_spectrum = {} # Dictionary to hold LASP data

LASP_spectrum['csv_data'] = np.genfromtxt('C:/Users/Will/Box/HERT_Box/Sr90 Testing/Sr90Y90.csv', delimiter=',', filling_values=0)
LASP_spectrum['KE'] = LASP_spectrum['csv_data'][:, 0]
LASP_spectrum['Sr90'] = LASP_spectrum['csv_data'][:, 1]/0.02
LASP_spectrum['Y90'] = LASP_spectrum['csv_data'][:, 2]/0.02
LASP_spectrum['combined'] = LASP_spectrum['csv_data'][:, 3]/0.02


# Calculate means and differences and store in the LASP_spectrum dictionary
LASP_spectrum['Sr90_mean'] = np.sum(LASP_spectrum['KE'] * LASP_spectrum['Sr90']) / np.sum(LASP_spectrum['Sr90'])
LASP_spectrum['Sr90_mean_diff'] = abs((LASP_spectrum['Sr90_mean'] - 0.196) / 0.196) * 100
LASP_spectrum['Y90_mean'] = np.sum(LASP_spectrum['KE'] * LASP_spectrum['Y90']) / np.sum(LASP_spectrum['Y90'])
LASP_spectrum['Y90_mean_diff'] = abs((LASP_spectrum['Y90_mean'] - 0.933) / 0.933) * 100
LASP_spectrum['Sr90_total'] = np.sum(LASP_spectrum['Sr90'])
LASP_spectrum['Y90_total'] = np.sum(LASP_spectrum['Y90'])
LASP_spectrum['combined_total'] = np.sum(LASP_spectrum['combined'])

# Plot Sr90 and Y90
fig, ax = plt.subplots(figsize=(7, 5))

Sr90_scatter = ax.scatter(LASP_spectrum['KE'], LASP_spectrum['Sr90'], label="Sr90", color='C2', s=8)
Y90_scatter = ax.scatter(LASP_spectrum['KE'], LASP_spectrum['Y90'], label="Y90", color='C0', s=8)

# Plot the sum, making it more prominent
combined_scatter = ax.scatter(LASP_spectrum['KE'], LASP_spectrum['combined'], label="Combined", s=8, color='C1')

ax.set_xlim(0, 2.3)
ax.set_xlabel("Kinetic Energy (MeV)")
ax.set_ylabel("Counts/s/MeV")
ax.set_title("LASP Beta Decay Spectra of Sr90 and Y90")
ax.grid(True)

'''
# Add vertical lines for mean energies
plt.axvline(x=0.196, color='C2', linestyle='--', label='Sr90 Mean: 0.196 MeV')
plt.axvline(x=0.933, color='C0', linestyle='--', label='Y90 Mean: 0.933 MeV')
ax.legend()
'''

plt.show()

#%% Compare generated and LASP spectra
LASP_comp = {} # Dictionary for comparison data

# Find probability distribution for Sr90 decay
Sr90_LASP, _ = fermi_theory(LASP_spectrum['KE'], 0.5459, 38, 90, -54E-3)
# Find probability distribution for Y90 decay
Y90_LASP, _ = fermi_theory(LASP_spectrum['KE'], 2.2785, 39, 90, -7.2E-3)

combined_LASP = Sr90_LASP + Y90_LASP

LASP_comp['int_gen'] = integrate.trapezoid(combined_LASP, LASP_spectrum['KE'])
LASP_comp['int_LASP'] = integrate.trapezoid(LASP_spectrum['combined'], LASP_spectrum['KE'])

LASP_comp['scale_LASP'] = LASP_comp['int_LASP'] / LASP_comp['int_gen']

LASP_comp['Sr90'] = Sr90_LASP * LASP_comp['scale_LASP']
LASP_comp['Y90'] = Y90_LASP * LASP_comp['scale_LASP']
LASP_comp['combined'] = combined_LASP * LASP_comp['scale_LASP']

LASP_comp['Sr90_mean'] = np.sum(LASP_spectrum['KE'] * LASP_comp['Sr90']) / np.sum(LASP_comp['Sr90'])
LASP_comp['Sr90_mean_diff'] = abs((LASP_comp['Sr90_mean'] - 0.196) / 0.196) * 100
LASP_comp['Y90_mean'] = np.sum(LASP_spectrum['KE'] * LASP_comp['Y90']) / np.sum(LASP_comp['Y90'])
LASP_comp['Y90_mean_diff'] = abs((LASP_comp['Y90_mean'] - 0.933) / 0.933) * 100

# Calculate total count rate for each element
LASP_comp['Sr90_total_scaled'] = sum(LASP_comp['Sr90'])
LASP_comp['Y90_total_scaled'] = sum(LASP_comp['Y90'])
LASP_comp['combined_total_scaled'] = sum(LASP_comp['combined'])

# Plot Sr90 and Y90
fig, ax = plt.subplots(figsize=(7, 5))

Sr90_scatter = ax.scatter(LASP_spectrum['KE'], LASP_comp['Sr90'], label="Sr90", color='C2', s=8)
Y90_scatter = ax.scatter(LASP_spectrum['KE'], LASP_comp['Y90'], label="Y90", color='C0', s=8)

# Plot the sum, making it more prominent
combined_scatter = ax.scatter(LASP_spectrum['KE'], LASP_comp['combined'], label="Combined", s=8, color='C1')

ax.set_xlim(0, 2.3)
ax.set_xlabel("Kinetic Energy (MeV)")
ax.set_ylabel("Counts/s/MeV")
ax.set_title("Scaled Beta Decay Spectra of Sr90 and Y90 to LASP Spectrum")
ax.grid(True)

'''
# Add vertical lines for mean energies
plt.axvline(x=LASP_comp['Sr90_mean'], color='C2', linestyle='--', label=f"Sr90 Mean: {Beta_theory['Sr90_mean']:.3f} MeV")
plt.axvline(x=LASP_comp['Y90_mean'], color='C0', linestyle='--', label=f"Y90 Mean: {Beta_theory['Y90_mean']:.3f} MeV")
ax.legend()
'''

plt.show()

Error = abs(LASP_comp['combined'] - LASP_spectrum['combined'])
SSE = sum((LASP_comp['combined'] - LASP_spectrum['combined'])**2)
MSE = SSE /len(LASP_spectrum['KE'])

#%% Adjust to time since creation
# The rate in 2021 was combined_total=28.818116725999996
# The total count rate of each source when it was generated
R_0 = LASP_spectrum["combined_total"]/2 * 0.02
# The total count rate of each source now
R_now = R_0 * (1/2)**((2025-2010)/28.91)

#
today_Sr90 = R_now * Beta_theory['Sr90_normalized'] / 0.01
today_Y90 = R_now * Beta_theory['Y90_normalized'] / 0.01
today_combined = today_Sr90 + today_Y90

# Plot Sr90 and Y90
fig, ax = plt.subplots(figsize=(7, 5))

Sr90_scatter = ax.scatter(Beta_theory['KE'], today_Sr90, label="Sr90", color='C2', s=8)
Y90_scatter = ax.scatter(Beta_theory['KE'], today_Y90, label="Y90", color='C0', s=8)

# Plot the sum, making it more prominent
combined_scatter = ax.scatter(Beta_theory['KE'], today_combined, label="Combined", s=8, color='C1')

ax.set_xlim(0, 2.3)
ax.set_xlabel("Kinetic Energy (MeV)")
ax.set_ylabel("Counts/s/MeV")
ax.set_title("Beta Decay Spectra of Sr90 and Y90 Today")
ax.grid(True)
'''
# Add vertical lines for mean energies
plt.axvline(x=Beta_theory['Sr90_mean'], color='C2', linestyle='--', label=f"Sr90 Mean: {Beta_theory['Sr90_mean']:.3f} MeV")
plt.axvline(x=Beta_theory['Y90_mean'], color='C0', linestyle='--', label=f"Y90 Mean: {Beta_theory['Y90_mean']:.3f} MeV")
ax.legend()
'''
plt.show()


#%% Aerospace data
Aero_data = {} # Dictionary for Aerospace data

Aero_data['csv_data'] = np.genfromtxt('C:/Users/Will/Box/HERT_Box/Aerospace Testing/Aerospace Beta-ray Spectrometer 2025-04-30.csv', delimiter=',', filling_values=0, skip_header=1)
Aero_data['channel'] = Aero_data['csv_data'][:, 3]
Aero_data['FWHM'] = Aero_data['csv_data'][:, 4]
Aero_data['KE'] = Aero_data['csv_data'][:, 6]
Aero_data['time'] = Aero_data['csv_data'][:, 7]
Aero_data['max_counts'] = Aero_data['csv_data'][:, 10]
Aero_data['counts_per_day'] = Aero_data['csv_data'][:, 11]

#%% Comparing Aero channels and energies
# Plotting E v channel
plt.figure(figsize=(8, 6))  # Adjust figure size as needed
plt.plot(Aero_data['channel'], Aero_data['KE'])
plt.xlabel('Channel Number')
plt.ylabel('Kinetic Energy (keV)')  # Label the y-axis with the correct unit
plt.title('Kinetic Energy vs Channel Number')
plt.grid(True)
plt.show()

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

#%% Calculate Aero energy resolution
Aero_data['energy_resolution'] = Aero_data['FWHM']/Aero_data['channel']

# Plotting energy resolution v channel
fig, ax1 = plt.subplots(figsize=(8, 6))  # Adjust figure size as needed
ax1.plot(Aero_data['channel'], Aero_data['energy_resolution'])
ax1.set_xlabel('Channel Number')
ax1.set_ylabel(r'Energy Resolution $\Delta E/E$')
ax1.set_title('Energy Resolution vs Channel Number')
ax1.grid(True)

# Get the limits from the first x-axis
ax1.set_xlim(0, 2000)
# Create the second x-axis
ax2 = ax1.twiny()  # Create a twin axis sharing the y-axis
ax2.set_xlim(intercept, slope * 2000 + intercept)
ax2.set_xlabel('Kinetic Energy (keV)')  # Label the top x-axis
plt.show()

#%% Look into Aero count rates
Beta_theory['Sr90_per_energy'] = Beta_theory['Sr90_normalized'] / 0.01
Beta_theory['Y90_per_energy'] = Beta_theory['Y90_normalized'] / 0.01
Beta_theory['combined_per_energy'] = Beta_theory['Sr90_per_energy'] + Beta_theory['Y90_per_energy']

Aero_data['countrate'] = Aero_data['max_counts']/Aero_data['time']
Aero_data['1bin_keV'] = slope * 1 + intercept
Aero_data['countrate_per_MeV'] = Aero_data['countrate'] / (Aero_data['1bin_keV']/1000*2)

intmax_theory_per_E = integrate.trapezoid(Beta_theory['combined_per_energy'], Beta_theory['KE'])
intmax_Aero_per_E = integrate.trapezoid(Aero_data['countrate_per_MeV'], Aero_data['KE']/1000)
Aeromax_scale = intmax_Aero_per_E / intmax_theory_per_E

Aeromax_maxscale_theory = max(Aero_data['countrate_per_MeV'])/max(Beta_theory['combined_per_energy'])
Aeromax_maxscale_Sr = max(Aero_data['countrate_per_MeV'])/max(Beta_theory['Sr90_per_energy'])
Aeromax_maxscale_LASP = max(Aero_data['countrate_per_MeV'])/max(LASP_spectrum['combined'])

Sr90_Aeromax_theory = Aeromax_scale * Beta_theory['Sr90_per_energy']
Y90_Aeromax_theory = Aeromax_scale * Beta_theory['Y90_per_energy']
Combined_Aeromax_theory = Aeromax_scale * Beta_theory['combined_per_energy']

# Plotting count rate v KE
plt.figure(figsize=(8, 6))  # Adjust figure size as needed
plt.scatter(Aero_data['KE']/1000, Aero_data['countrate'], color='black', marker='*', label='Aerospace')
plt.xlim(0, 2)
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Count Rate (counts/s)')
plt.title('Count Rate vs KE')
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(8, 6))  # Adjust figure size as needed
plt.scatter(Beta_theory['KE'], Beta_theory['Sr90_per_energy'] * Aeromax_maxscale_Sr, label="Sr90", color='C2', s=8)
#plt.scatter(Beta_theory['KE'], Y90_Aeromax_theory , label="Y90", color='C0', s=8)
plt.scatter(Beta_theory['KE'], Beta_theory['combined_per_energy'] * Aeromax_maxscale_theory, label="Theory", s=10, color='C1')
plt.scatter(LASP_spectrum['KE'], LASP_spectrum['combined'] * Aeromax_maxscale_LASP, label="LASP (scaled)", s=10, color='blue')
plt.scatter(Aero_data['KE']/1000, Aero_data['countrate_per_MeV'], color='black', marker='*', label='Aerospace')
plt.xlim(0, 2)
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Count Rate / MeV (counts/s/MeV)')
plt.title('Count Rate per MeV vs KE')
plt.grid(True)
plt.legend()
plt.show()


max_counts_index_theory = np.argmax(Beta_theory['combined_per_energy'])
ke_at_max_counts_theory = Beta_theory['KE'][max_counts_index_theory]

max_counts_index_LASP = np.argmax(LASP_spectrum['combined'])
ke_at_max_counts_LASP = LASP_spectrum['KE'][max_counts_index_LASP]

max_counts_index_Aero = np.argmax(Aero_data['countrate_per_MeV'])
ke_at_max_counts_Aero = Aero_data['KE'][max_counts_index_Aero]/1000

plt.figure(figsize=(8, 6))  # Adjust figure size as needed
plt.scatter(Beta_theory['KE'] - ke_at_max_counts_theory, Beta_theory['Sr90_per_energy'] * Aeromax_maxscale_Sr, label="Sr90", color='C2', s=8)
plt.scatter(Beta_theory['KE'] - ke_at_max_counts_theory, Beta_theory['combined_per_energy'] * Aeromax_maxscale_theory, label="Theory", s=10, color='C1')
plt.scatter(LASP_spectrum['KE'] - ke_at_max_counts_LASP, LASP_spectrum['combined'] * Aeromax_maxscale_LASP, label="LASP (scaled)", s=10, color='blue')
plt.scatter(Aero_data['KE']/1000 - ke_at_max_counts_Aero, Aero_data['countrate_per_MeV'], color='black', marker='*', label='Aerospace')
plt.xlim(-0.5, 2)
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Count Rate / MeV (counts/s/MeV)')
plt.title('Count Rate per MeV vs KE centerted at Sr90 peak energy')
plt.grid(True)
plt.legend()
plt.show()
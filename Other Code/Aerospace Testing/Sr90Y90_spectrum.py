#%% Import libraries
import numpy as np
import scipy.constants as sc
from scipy.special import gamma
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker

FONT_SIZE = 14
plt.rcParams.update({
    'font.size': FONT_SIZE,          # Base font size
    'axes.titlesize': FONT_SIZE + 2, # Subplot titles
    'axes.labelsize': FONT_SIZE,     # Axis labels (L*, Time, Flux)
    'xtick.labelsize': FONT_SIZE,    # X-axis tick numbers
    'ytick.labelsize': FONT_SIZE,    # Y-axis tick numbers
    'legend.fontsize': FONT_SIZE,    # Legend text
    'figure.titlesize': FONT_SIZE + 4, # Main suptitle
    'legend.markerscale': 2.0
})

#%%
'''

        FERMI'S THEORY OF BETA DECAY SPECTRUM GENERATION

'''

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
Beta_theory['KE'] = np.linspace(0, 3, num = 3001)

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
fig, ax = plt.subplots(figsize=(12, 2.5))

Sr90_scatter = ax.scatter(Beta_theory['KE'], Beta_theory['Sr90'], label="Sr90", color='C9')
Y90_scatter = ax.scatter(Beta_theory['KE'], Beta_theory['Y90'], label="Y90", color='C4')

# Plot the sum, making it more prominent
combined_scatter = ax.scatter(Beta_theory['KE'], Beta_theory['combined_spectrum'], label="Combined", color='C1')

ax.set_xlim(0, 2.3)
ax.set_xlabel("Kinetic Energy (MeV)")
ax.set_ylabel("Probability Distribution", labelpad=35)
ax.set_yticklabels([])
ax.set_title("Beta Decay Spectra of Sr90 and Y90")
plt.grid(axis='x')

# Add vertical lines for mean energies
plt.axvline(x=0.196, color='C9', linestyle='--', linewidth=3, label='Sr90 Mean: 0.196 MeV')
plt.axvline(x=0.933, color='C4', linestyle='--', linewidth=3, label='Y90 Mean: 0.933 MeV')
ax.legend()
plt.show()

#%%
'''

        IMPORT LASP BETA SPECTRUM DATA

'''


#%% Import LASP spectra
LASP_spectrum = {} # Dictionary to hold LASP data

LASP_spectrum['csv_data'] = np.genfromtxt('C:/Users/wzt0020/Box/HERT_Box/Sr90 Testing/Sr90Y90.csv', delimiter=',', filling_values=0)
LASP_spectrum['KE'] = LASP_spectrum['csv_data'][:, 0]
LASP_spectrum['Sr90_CR'] = LASP_spectrum['csv_data'][:, 1]
LASP_spectrum['Y90_CR'] = LASP_spectrum['csv_data'][:, 2]
LASP_spectrum['combined_CR'] = LASP_spectrum['csv_data'][:, 3]
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
fig, ax = plt.subplots(figsize=(14, 6))
Sr90_scatter = ax.scatter(LASP_spectrum['KE'], LASP_spectrum['Sr90_CR'], label="Sr90", color='C2', s=8)
Y90_scatter = ax.scatter(LASP_spectrum['KE'], LASP_spectrum['Y90_CR'], label="Y90", color='C0', s=8)
# Plot the sum, making it more prominent
combined_scatter = ax.scatter(LASP_spectrum['KE'], LASP_spectrum['combined_CR'], label="Combined", s=8, color='C1')

ax.set_xlim(0, 2.3)
ax.set_xlabel("Kinetic Energy (MeV)")
ax.set_ylabel("Counts/s")
ax.set_title("LASP Beta Decay Spectra of Sr90 and Y90")
ax.grid(True)
# Add vertical lines for mean energies
# plt.axvline(x=0.196, color='C2', linestyle='--', label='Sr90 Mean: 0.196 MeV')
# plt.axvline(x=0.933, color='C0', linestyle='--', label='Y90 Mean: 0.933 MeV')
ax.legend()
plt.show()

#%% Normalize LASP spectra
total_counts_LASP = sum(LASP_spectrum['combined_CR'])
LASP_spectrum['Sr90_CR_norm'] = LASP_spectrum['Sr90_CR'] / total_counts_LASP
LASP_spectrum['Y90_CR_norm'] = LASP_spectrum['Y90_CR'] / total_counts_LASP
LASP_spectrum['combined_CR_norm'] = LASP_spectrum['combined_CR'] / total_counts_LASP

# Plot Normalized Sr90 and Y90
# fig, ax = plt.subplots(figsize=(14, 6))
# Sr90_scatter = ax.scatter(LASP_spectrum['KE'], LASP_spectrum['Sr90_CR_norm'], label="Sr90", color='C2', s=8)
# Y90_scatter = ax.scatter(LASP_spectrum['KE'], LASP_spectrum['Y90_CR_norm'], label="Y90", color='C0', s=8)
# # Plot the sum, making it more prominent
# combined_scatter = ax.scatter(LASP_spectrum['KE'], LASP_spectrum['combined_CR_norm'], label="Combined", s=8, color='C1')

# ax.set_xlim(0, 2.3)
# ax.set_xlabel("Kinetic Energy (MeV)")
# ax.set_ylabel("Counts/s")
# ax.set_title("Normalized LASP Beta Decay Spectra of Sr90 and Y90")
# ax.grid(True)
# # Add vertical lines for mean energies
# # plt.axvline(x=0.196, color='C2', linestyle='--', label='Sr90 Mean: 0.196 MeV')
# # plt.axvline(x=0.933, color='C0', linestyle='--', label='Y90 Mean: 0.933 MeV')
# ax.legend()
# plt.show()

#%% Calculate Expected LASP Flux
lda = np.log(2) / 28.91  # Decay constant for Sr90 in years

REPTile2_r = 4.35 # distance between source and detector center, cm
REPTile2_detector = 2.0 # detector diameter, cm
REPTile2_FOV = np.rad2deg(np.atan(REPTile2_detector/2 / REPTile2_r)) # degrees

HERT_r = 6.3 + 0.15736 # distance between source and detector center, cm
HERT_detector = 1.8 # detector diameter, cm
HERT_FOV = np.rad2deg(np.atan(HERT_detector/2 / HERT_r)) # degrees
directional_scaling = (REPTile2_FOV / HERT_FOV)
HERT_eff = 0.5307 # Total instrument efficiency

today_count_rate = total_counts_LASP / 2 * np.exp(-lda * (2025 - 2010)) / HERT_eff / directional_scaling # Initial count rate in 2010
d = 0.15736  # distance between front of collimator and radiation source, cm
s = 0.3/2  # radius of radiation source, cm
LASP_today_scale = 2*today_count_rate/(1-np.cos(np.atan(s/(d+6.3))))/(4*np.pi*(d+6.3)**2*0.02)

LASP_spectrum['Sr90_today'] = LASP_spectrum['Sr90_CR_norm'] * LASP_today_scale
LASP_spectrum['Y90_today'] = LASP_spectrum['Y90_CR_norm'] * LASP_today_scale
LASP_spectrum['combined_today'] = LASP_spectrum['combined_CR_norm'] * LASP_today_scale

Theory_maxscale_LASP = max(LASP_spectrum['combined_today'])/max(Beta_theory['combined_spectrum'])

# Plot Today scaled Sr90 and Y90
fig, ax = plt.subplots(figsize=(12, 2.5)) 
#Sr90_scatter = ax.scatter(LASP_spectrum['KE'], LASP_spectrum['Sr90_today'], label="Sr90", color='C2', s=20)
#Y90_scatter = ax.scatter(LASP_spectrum['KE'], LASP_spectrum['Y90_today'], label="Y90", color='C0', s=20)
# Plot the sum, making it more prominent
ax.scatter(Beta_theory['KE'], Beta_theory['combined_spectrum'] * Theory_maxscale_LASP, label="Theory (scaled)", color='C1')
combined_scatter = ax.scatter(LASP_spectrum['KE'], LASP_spectrum['combined_today'], label="LASP", color='C2')

ax.set_xlim(0, 2.3)
#ax.set_yscale('log')
ax.set_xlabel("Kinetic Energy (MeV)")
ax.set_ylabel(r'Flux ($\# / \text{s } \text{sr } \text{cm}^2 \text{ MeV}$)')
ax.tick_params(axis='both', which='major', labelsize=14)
#ax.set_title("Today's LASP Beta Decay Spectra of Sr90 and Y90")
ax.grid(True)

# create secondary x-axis for channel numbers
E_eff_all = np.genfromtxt(r"C:\Users\wzt0020\Box\HERT_Box\Sr90 Testing\effective_energies_DARTBe.txt")
differences = np.diff(E_eff_all)
try:
    # We add 1 because np.diff shrinks the array by 1.
    break_index = np.where(differences <= 0)[0][0] + 1
except IndexError:
    # If the array is fully sorted, there is no negative difference.
    break_index = len(E_eff_all)
E_eff_reduce = (np.arange(len(E_eff_all)) < break_index) & (E_eff_all < 2.3)
E_eff_DART = E_eff_all[E_eff_reduce]

ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
channel_indices = np.arange(1, len(E_eff_DART) + 1)
ax2.set_xticks(LASP_spectrum['KE'][np.searchsorted(LASP_spectrum['KE'], E_eff_DART[E_eff_DART <= LASP_spectrum['KE'][-1]])])
channel_labels = [''] * len(E_eff_DART[E_eff_DART <= LASP_spectrum['KE'][-1]])
for i in range(len(E_eff_DART[E_eff_DART <= LASP_spectrum['KE'][-1]])):
    if (i + 1) % 5 == 0:  # Label every 5th channel
        channel_labels[i] = str(i + 1)
ax2.set_xticklabels(channel_labels)
ax2.set_xlabel("Channel Number", labelpad=10) # labelpad moves the label up
ax2.tick_params(axis='x', which='major', labelsize=14)

ax.legend()
plt.show()

#%% Calculate Expected LASP count rates
geo_factor_LASP = np.genfromtxt(r"C:\Users\wzt0020\Box\HERT_Box\Sr90 Testing\geofactor_EC_SrTest.txt")
geo_factor_total_LASP = np.sum(geo_factor_LASP, axis=0)
bins_LASP = geo_factor_LASP.shape[-1]
energy_edges_LASP = np.linspace(0, 8, bins_LASP + 1)
bin_width_LASP = np.diff(energy_edges_LASP)
energy_midpoints_LASP = (energy_edges_LASP[1:] + energy_edges_LASP[:-1]) / 2
energy_range_LASP = energy_edges_LASP[:-1] < 2
geo_factor_LASP = geo_factor_LASP[E_eff_reduce,:][:,energy_range_LASP] 
geo_factor_total_LASP = geo_factor_total_LASP[energy_range_LASP]
duration = 10*60

LASP_count_rate_EC = np.zeros(geo_factor_LASP.shape[0])
for channel in range(geo_factor_LASP.shape[0]):
    LASP_count_rate_EC[channel] = np.sum(geo_factor_LASP[channel,:]*LASP_spectrum['combined_today'][LASP_spectrum['KE']<2]*bin_width_LASP[energy_range_LASP])

# Plot expected count rates
fig, ax = plt.subplots(figsize=(16, 2.5))
count_rate_scatter = ax.scatter(E_eff_DART, LASP_count_rate_EC, color='C2', s=60)
ax.set_xlim(0, 2.2)
ax.set_xlabel("Kinetic Energy (MeV)")
ax.set_ylabel("Count Rate (#/s)")
ax.tick_params(axis='both', which='major')
ax.set_title("Today's LASP Predicted Count Rates")
ax.grid(True)
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
channel_indices = np.arange(1, len(E_eff_DART) + 1)
ax2.set_xticks(LASP_spectrum['KE'][np.searchsorted(LASP_spectrum['KE'], E_eff_DART[E_eff_DART <= LASP_spectrum['KE'][-1]])])
channel_labels = [''] * len(E_eff_DART[E_eff_DART <= LASP_spectrum['KE'][-1]])
for i in range(len(E_eff_DART[E_eff_DART <= LASP_spectrum['KE'][-1]])):
    if (i + 1) % 5 == 0:  # Label every 5th channel
        channel_labels[i] = str(i + 1)
ax2.set_xticklabels(channel_labels)
ax2.set_xlabel("Channel Number", labelpad=10) # labelpad moves the label up
ax2.tick_params(axis='x', which='major', labelsize=14)
#plt.yscale('log')
plt.ylim(1e-3, 4)
plt.show()

channels_show = 20
channels = np.linspace(1, E_eff_DART.shape[0], E_eff_DART.shape[0])
fig, ax1 = plt.subplots(figsize=(16, 2.5))
ax1.bar(channels[:channels_show], duration*LASP_count_rate_EC[:channels_show],color='C2')
ax1.errorbar(channels[:channels_show], duration*LASP_count_rate_EC[:channels_show], yerr=np.sqrt(duration*LASP_count_rate_EC[:channels_show]), fmt='none', color='black', ecolor='black', capsize=5)
ax1.set_ylim(0, max(duration*LASP_count_rate_EC)*1.2)
ax1.set_xlabel('Channel Number')
ax1.set_ylabel('Predicted Counts')
ax1.tick_params(axis='x', which='major')
ax1.grid(True, axis='y', linestyle='--', alpha=0.7)
max_channel = channels[:channels_show][-1]
tick_locations = np.arange(1, max_channel + 1)
ax1.set_xticks(tick_locations)
ax1.set_xlim(0, 20.5)
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(tick_locations)
ax2.set_xlabel('Effective Energy (MeV)', labelpad=15)
ax2.set_xticklabels([f'{e:.2f}' if i % 2 == 0 else '' for i, e in enumerate(E_eff_DART[:channels_show])], rotation=45, ha='left')
plt.title('Predicted Counts at {:.4f} MeV Per Channel'.format(specific_energy))
plt.show()

# Plot expected counts at 10 minutes for each detector
detector_number = np.linspace(1,9,9)
detector_efficiency = np.genfromtxt(r"C:\Users\wzt0020\Box\HERT_Box\Sr90 Testing\efficiency_Detectors_LASPTest.txt")
fig, ax = plt.subplots(figsize=(16, 2.5))
count_rate_plot = ax.bar(
    detector_number, 
    detector_efficiency * today_count_rate * (10 * 60),               # Y-data (Counts)
    yerr=np.sqrt(detector_efficiency * today_count_rate * (10 * 60)), # Y-Error (square root of counts for Poisson statistics)
    # Appearance Settings for the Error Bar:
    width=0.8,                                      
    color='C2',                                   # Color of both marker and error bars
    capsize=20,                                    # Makes the horizontal end-caps larger
)
ax.set_xlabel("Detector Number", fontsize=16)
ax.set_ylabel("Counts", fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(1, 10, 1))
ax.set_title("Predicted Counts at 10 minutes", fontsize=18)
ax.grid(True)
plt.ylim(0,6000)
#plt.yscale('log')
plt.show()

#%% Calculate Expected LASP time to counts and save
LASP_time_to_10 = 10 / LASP_count_rate_EC /60 #convert to minutes
LASP_time_to_100 = 100 / LASP_count_rate_EC /60 #convert to minutes

fig, ax = plt.subplots(figsize=(16, 4))
time_to_100_LASP_scatter = ax.scatter(E_eff_DART, LASP_time_to_10, color='C0', label='10 Counts', s=60)
time_to_1000_LASP_scatter = ax.scatter(E_eff_DART, LASP_time_to_100, color='C1', label='100 Counts', s=60)
ax.set_xlim(0, 2)
ax.set_ylim(0, 20)
ax.set_xlabel("Kinetic Energy (MeV)", fontsize=16)
ax.set_ylabel('Time (minutes)', fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_title("Expected Exposure Time", fontsize=18)
ax.grid(True)
ax.legend(loc =  'upper left',fontsize=14)
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
channel_indices = np.arange(1, len(E_eff_DART) + 1)
ax2.set_xticks(LASP_spectrum['KE'][np.searchsorted(LASP_spectrum['KE'], E_eff_DART[E_eff_DART <= LASP_spectrum['KE'][-1]])])
channel_labels = [''] * len(E_eff_DART[E_eff_DART <= LASP_spectrum['KE'][-1]])
for i in range(len(E_eff_DART[E_eff_DART <= LASP_spectrum['KE'][-1]])):
    if (i + 1) % 5 == 0:  # Label every 5th channel
        channel_labels[i] = str(i + 1)
ax2.set_xticklabels(channel_labels)
ax2.set_xlabel("Channel Number", fontsize=16, labelpad=10) # labelpad moves the label up
ax2.tick_params(axis='x', which='major', labelsize=14)
plt.show()

channel_numbers = np.arange(1, len(E_eff_DART) + 1)
time_to_counts = np.column_stack((channel_numbers, E_eff_DART, LASP_time_to_10, LASP_time_to_100))
np.savetxt('LASP_time_to_counts.txt', time_to_counts,
    fmt=['%6d', '%30.6f', '%23.6f', '%24.6f'],
    delimiter='\t',
    header="Channel\tEffective_Energy_(MeV)\tTime_to_10_Counts_(min)\tTime_to_100_Counts_(min)",
)

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
# fig, ax = plt.subplots(figsize=(14, 6))

# Sr90_scatter = ax.scatter(LASP_spectrum['KE'], LASP_comp['Sr90'], label="Sr90", color='C2', s=8)
# Y90_scatter = ax.scatter(LASP_spectrum['KE'], LASP_comp['Y90'], label="Y90", color='C0', s=8)

# # Plot the sum, making it more prominent
# combined_scatter = ax.scatter(LASP_spectrum['KE'], LASP_comp['combined'], label="Combined", s=8, color='C1')

# ax.set_xlim(0, 2.3)
# ax.set_xlabel("Kinetic Energy (MeV)")
# ax.set_ylabel("Counts/s/MeV")
# ax.set_title("Scaled Beta Decay Spectra of Sr90 and Y90 to LASP Spectrum")
# ax.grid(True)

# # Add vertical lines for mean energies
# plt.axvline(x=LASP_comp['Sr90_mean'], color='C2', linestyle='--', label=f"Sr90 Mean: {Beta_theory['Sr90_mean']:.3f} MeV")
# plt.axvline(x=LASP_comp['Y90_mean'], color='C0', linestyle='--', label=f"Y90 Mean: {Beta_theory['Y90_mean']:.3f} MeV")

# ax.legend()
# plt.show()


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
fig, ax = plt.subplots(figsize=(14, 6))

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
'''
ax.legend()
plt.show()


#%%
'''

        COMPARE TO AEROSPACE SPECTRUM DATA

'''

#%% Aerospace data
Aero_data = {} # Dictionary for Aerospace data

Aero_data['csv_data'] = np.genfromtxt('C:/Users/wzt0020/Box/HERT_Box/Aerospace Testing/Resources from Aero/Aerospace Beta-ray Spectrometer 2025-04-30.csv', delimiter=',', filling_values=0, skip_header=1)
Aero_data['channel'] = Aero_data['csv_data'][:, 3]
Aero_data['FWHM'] = Aero_data['csv_data'][:, 4]
Aero_data['KE'] = Aero_data['csv_data'][:, 6]
Aero_data['time'] = Aero_data['csv_data'][:, 7]
Aero_data['max_counts'] = Aero_data['csv_data'][:, 10]
Aero_data['counts_per_day'] = Aero_data['csv_data'][:, 11]

#%% Comparing Aero channels and energies
# Plotting E v channel
# plt.figure(figsize=(14, 6))
# plt.plot(Aero_data['channel'], Aero_data['KE'])
# plt.xlabel('Channel Number')
# plt.ylabel('Kinetic Energy (keV)')  # Label the y-axis with the correct unit
# plt.title('Kinetic Energy vs Channel Number')
# plt.grid(True)
# plt.show()

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

#%% Look into Aero count rates
Beta_theory['Sr90_per_energy'] = Beta_theory['Sr90_normalized'] / 0.01
Beta_theory['Y90_per_energy'] = Beta_theory['Y90_normalized'] / 0.01
Beta_theory['combined_per_energy'] = Beta_theory['Sr90_per_energy'] + Beta_theory['Y90_per_energy']

Aero_data['countrate'] = Aero_data['max_counts']/(Aero_data['time']/60/60)
Aero_data['countrate_per_MeV'] = Aero_data['countrate'] / (slope/1000)

intmax_theory_per_E = integrate.trapezoid(Beta_theory['combined_per_energy'], Beta_theory['KE'])
intmax_Aero_per_E = integrate.trapezoid(Aero_data['countrate_per_MeV'], Aero_data['KE']/1000)
Aeromax_scale = intmax_Aero_per_E / intmax_theory_per_E

Aeromax_maxscale_theory = max(Aero_data['countrate'])/(60*60)/max(Beta_theory['combined_spectrum'])
Aeromax_maxscale_Sr = max(Aero_data['countrate'])/max(Beta_theory['Sr90'])
Aeromax_maxscale_LASP = max(Aero_data['countrate'])/max(LASP_spectrum['combined'])

Sr90_Aeromax_theory = Aeromax_scale * Beta_theory['Sr90_per_energy']
Y90_Aeromax_theory = Aeromax_scale * Beta_theory['Y90_per_energy']
Combined_Aeromax_theory = Aeromax_scale * Beta_theory['combined_per_energy']

max_counts_index_theory = np.argmax(Beta_theory['combined_per_energy'])
ke_at_max_counts_theory = Beta_theory['KE'][max_counts_index_theory]

max_counts_index_LASP = np.argmax(LASP_spectrum['combined'])
ke_at_max_counts_LASP = LASP_spectrum['KE'][max_counts_index_LASP]

max_counts_index_Aero = np.argmax(Aero_data['countrate_per_MeV'])
ke_at_max_counts_Aero = Aero_data['KE'][max_counts_index_Aero]/1000

fig, ax = plt.subplots(figsize=(8, 4)) 
#ax.scatter(Beta_theory['KE'], Beta_theory['Sr90_per_energy'] * Aeromax_maxscale_Sr, label="Sr90", color='C2', s=8)
ax.scatter(Beta_theory['KE'], Beta_theory['combined_spectrum'] * Aeromax_maxscale_theory, label="Theory (scaled)", color='C1')
#ax.scatter(LASP_spectrum['KE'], LASP_spectrum['combined_today'], label="LASP", color='blue')
ax.scatter(Aero_data['KE']/1000, Aero_data['countrate']/(60*60), color='black', marker='*', label='Aerospace')
ax.set_xlim(0, 2.3)
ax.set_xlabel('Kinetic Energy (MeV)')
ax.set_ylabel('Count Rate (#/second)')
ax.set_title('Theoretical v. Aerospace Measured Beta Decay of Sr90 and Y90')
#ax.set_yticklabels([])
ax.grid()
ax.legend()
plt.show()

#%% Convert from Aero detector to HERT
aero_detector = 8e-3 /2
hert_detector = 18e-3 /2

#assuming that aero and HERT see only one source (the central source)
scale_factor = 1
hert_countrate = Aero_data['countrate'] * scale_factor

from scipy.optimize import curve_fit
from scipy.stats import beta
def beta_func(x, A, a, b, loc, scale):
    return A * beta.pdf(x, a, b, loc=loc, scale=scale)

A_guess = np.max(Aero_data['countrate'])
a_guess = 5
b_guess = 2
loc_guess = np.min(Aero_data['KE']/1000)
scale_guess = np.max(Aero_data['KE']/1000) - np.min(Aero_data['KE']/1000)
p0_guess = [A_guess, a_guess, b_guess, loc_guess, scale_guess]

popt_aero, pcov_aero = curve_fit(beta_func, Aero_data['KE']/1000, Aero_data['countrate'], p0=p0_guess)
A_fit_aero, a_fit_aero, b_fit_aero, loc_fit_aero, scale_fit_aero = popt_aero
print("Fitted Parameters (Aero):")
print(f"Amplitude (A): {A_fit_aero:.2f}")
print(f"Shape (a): {a_fit_aero:.2f}")
print(f"Shape (b): {b_fit_aero:.2f}")
print(f"Location (loc): {loc_fit_aero:.2f}")
print(f"Scale (scale): {scale_fit_aero:.2f}")
aero_fit = beta_func(Beta_theory['KE'], *popt_aero)

popt_hert, pcov_hert = curve_fit(beta_func, Aero_data['KE']/1000, hert_countrate, p0=p0_guess)
A_fit_hert, a_fit_hert, b_fit_hert, loc_fit_hert, scale_fit_hert = popt_hert
print("Fitted Parameters (HERT):")
print(f"Amplitude (A): {A_fit_hert:.2f}")
print(f"Shape (a): {a_fit_hert:.2f}")
print(f"Shape (b): {b_fit_hert:.2f}")
print(f"Location (loc): {loc_fit_hert:.2f}")
print(f"Scale (scale): {scale_fit_hert:.2f}")
#hert_fit = beta_func(Beta_theory['KE'], *popt_hert)

# import geometric factor
geo_factor = np.genfromtxt(r"C:\Users\wzt0020\Box\HERT_Box\Aerospace Testing\efficiency_EC_AeroTest.txt")
geo_factor_total = np.nansum(geo_factor, axis=0)
bins = geo_factor.shape[-1]
energy_edges = np.logspace(np.log10(0.01), np.log10(8), bins + 1)
bin_width = np.diff(energy_edges)
energy_midpoints = (energy_edges[1:] + energy_edges[:-1]) / 2
energy_range = energy_edges[:-1] < 2
geo_factor_total = geo_factor_total[energy_range]
hert_fit = beta_func(energy_midpoints[energy_range], *popt_hert) * geo_factor_total

#%% Plot instrument efficiency with points on selected energies
# Optional Energies: 0.723, 0.7531, 0.957, 0.983, 1.010, 1.155, 1.202, 1.294, 1.402, 1.459
selected_energies = [1.038, 1.059, 1.080, 1.102, 1.117, .9265, .6950, .8782, 1.25, 1.347, 1.50]
selected_energies =np.sort(selected_energies)
selected_indices = np.searchsorted(energy_midpoints[energy_range], selected_energies)
duration = 10*60

plasma_cmap = cm.get_cmap('plasma')
colors = plasma_cmap(np.linspace(0, 1, geo_factor.shape[0]))
# plt.figure(figsize=(14, 6))
# for channel in range(geo_factor.shape[0]):
#     current_color = colors[channel]
#     plt.plot(energy_midpoints[energy_range], geo_factor[channel, energy_range], color=current_color, linewidth=3, alpha=0.7)
#     plt.scatter(energy_midpoints[selected_indices], geo_factor[channel, energy_range][selected_indices], color=current_color, marker='o', zorder=3)
# for idx in selected_indices:
#     plt.axvline(x=energy_midpoints[idx],color='black',linestyle='--',alpha=0.5,linewidth=2)
# plt.tick_params(axis='both', which='major')
# plt.xlim(0.5, 2)
# plt.ylim(1e-4, 1)
# plt.yscale('log')
# plt.xlabel('Kinetic Energy (MeV)')
# plt.ylabel('Insturment Efficiency')
# plt.grid(True)
# plt.show()

plt.figure(figsize=(16, 3.5))
for channel in range(geo_factor.shape[0]):
    current_color = colors[channel]
    plt.plot(energy_midpoints[energy_range], duration*geo_factor[channel, energy_range], color=current_color, linewidth=3, alpha=0.7)
    plt.scatter(energy_midpoints[selected_indices], duration*geo_factor[channel, energy_range][selected_indices], color=current_color, marker='o', zorder=3)
for idx in selected_indices:
    plt.axvline(x=energy_midpoints[idx],color='black',linestyle='--',alpha=0.5,linewidth=2)
plt.tick_params(axis='both', which='major')
plt.xlim(0.5, 2)
plt.ylim(0, 280)
#plt.yscale('log')
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Counts')
plt.title('Predicted Counts at 10 Minutes')
plt.grid(True)
plt.show()

#%% Show expected counts at a specific energy
specific_index = 5 #5, 6, 9
specific_energy = selected_energies[specific_index]
specific_energy_all = np.searchsorted(energy_midpoints, specific_energy)
specific_geo_factor = geo_factor[:, specific_energy_all]

E_eff = np.genfromtxt(r"C:\Users\wzt0020\Box\HERT_Box\Aerospace Testing\E_eff_Aero.txt")
differences = np.diff(E_eff)

channels = np.linspace(1, E_eff.shape[0], E_eff.shape[0])
duration = 10*60

# fig, ax1 = plt.subplots(figsize=(10, 6))
# ax1.bar(channels[:E_eff.shape[0]], specific_geo_factor[:E_eff.shape[0]])
# ax1.set_xlim(0, E_eff.shape[0]+1)
# ax1.set_xlabel('Channel Number', fontsize=14)
# ax1.set_ylabel('Instrument Efficiency', fontsize=14)
# ax1.tick_params(axis='x', which='major', labelsize=12)
# ax1.grid(True, axis='y', linestyle='--', alpha=0.7)
# ax2 = ax1.twiny()
# ax2.set_xlim(ax1.get_xlim())
# ax2.set_xticks(channels)
# ax2.set_xticklabels([f'{e:.2f}' if i % 5 == 0 else '' for i, e in enumerate(E_eff)], rotation=45, ha='left')
# ax2.set_xlabel('Effective Energy (MeV)', fontsize=14, labelpad=15)
# plt.title('Instrument Efficiency at {:.4f} MeV Per Channel'.format(specific_energy), fontsize=16)
# plt.show()

channels_show = 20
fig, ax1 = plt.subplots(figsize=(8, 5))
ax1.bar(channels[:channels_show], duration*specific_geo_factor[:channels_show],color='C3')
ax1.errorbar(channels[:channels_show], duration*specific_geo_factor[:channels_show], yerr=np.sqrt(duration*specific_geo_factor[:channels_show]), fmt='none', color='black', ecolor='black', capsize=5)
ax1.set_ylim(0, max(duration*specific_geo_factor)*1.2)
ax1.set_xlabel('Channel Number')
ax1.set_ylabel('Predicted Counts')
ax1.tick_params(axis='x', which='major')
ax1.grid(True, axis='y', linestyle='--', alpha=0.7)
max_channel = channels[:channels_show][-1]
tick_locations = np.arange(1, max_channel + 1)
ax1.set_xticks(tick_locations)
ax1.set_xlim(0, 20.5)
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(tick_locations)
ax2.set_xlabel('Effective Energy (MeV)', labelpad=15)
ax2.set_xticklabels([f'{e:.2f}' if i % 2 == 0 else '' for i, e in enumerate(E_eff[:channels_show])], rotation=45, ha='left')
plt.title('Predicted Counts at {:.4f} MeV Per Channel'.format(specific_energy))
plt.show()

#%% Plotting count rate v KE
plt.figure(figsize=(8, 6))
plt.scatter(Aero_data['KE']/1000, Aero_data['countrate'], color='C0', marker='*', label='Aerospace')
plt.plot(Beta_theory['KE'], aero_fit, color='C0', label='Aerospace Fit')
plt.plot(energy_midpoints[energy_range], hert_fit, color='C1', label='HERT Fit')
plt.axvline(x=0.6,color='black',linestyle='--',label='HERT Threshold')
plt.scatter(energy_midpoints[selected_indices],hert_fit[selected_indices],color='black',marker='o',label='selected energies',zorder=3)
plt.xlim(0, 2)
#plt.ylim(0, 20000)
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Count Rate (counts/hour)')
plt.title('Count Rate vs KE')
plt.grid(True)
plt.legend()
plt.show()

#%% Plot counts as function of energy for a given time
# Calculate when Poisson uncertainty is less than 10%
# 10**2 counts is 10%, 4*10**2 is 5%, 10**3 counts is 3.2%, 10**4 counts is 1%
poisson_error = np.array((0.1,0.05,0.032,0.01)) 
poisson_count = (1/poisson_error)**2
plt.figure(figsize=(8, 6))  # Adjust figure size as needed
plt.plot(energy_midpoints[energy_range], hert_fit, label='1 hour')
plt.plot(energy_midpoints[energy_range], hert_fit*10, label='10 hours')
plt.plot(energy_midpoints[energy_range], hert_fit*24, label='24 hours')
plt.errorbar(energy_midpoints[selected_indices], hert_fit[selected_indices], yerr=np.sqrt(hert_fit[selected_indices]), color='C0', marker='o',
             linestyle='None', capsize=5, zorder=3)
plt.errorbar(energy_midpoints[selected_indices], hert_fit[selected_indices]*10, yerr=np.sqrt(hert_fit[selected_indices]*10), color='C1', marker='o',
             linestyle='None', capsize=5, zorder=3)
plt.errorbar(energy_midpoints[selected_indices], hert_fit[selected_indices]*24, yerr=np.sqrt(hert_fit[selected_indices]*24), color='C2', marker='o',
             linestyle='None', capsize=5, zorder=3)
line_styles = [':', '-.','--','-']
for i in range(len(poisson_error)):
    label_text = f'{poisson_error[i]*100:.0f}% Error'
    plt.axhline(y=poisson_count[i], color='black',linestyle=line_styles[i],label=label_text)
plt.axvline(x=0.6,color='black',linestyle='--',label='HERT Threshold')
plt.xlim(0.5, 2)
plt.ylim(10**0, 10**6)
#plt.yticks(np.arange(0, 25, 3))
plt.yscale('log')
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('# Counts')
plt.grid(True)
plt.legend()
plt.show()

#%% Figure out the time it takes to reach the Poisson error
time_to_error = np.zeros((len(hert_fit),len(poisson_count)))
time_to_error_error = np.zeros((len(hert_fit),len(poisson_count)))
time_error_indices = np.zeros((len(hert_fit),len(poisson_count)))
energy_error = np.zeros((len(hert_fit),len(poisson_count)))

selected_times = [1,10,24]
time_indices = np.zeros((len(selected_times),len(poisson_count)),int)
max_energies = np.zeros((len(selected_times),len(poisson_count)))
for i in range(len(poisson_count)):
    time_indices[:,i] = np.searchsorted(time_to_error[energy_midpoints[energy_range]>=0.8,i], selected_times)
    max_energies[:,i] = energy_midpoints[energy_midpoints>=0.8][time_indices[:,i]]

#%% Replot Time to Counts as a function of KE
time_to_100_hert = 100/(hert_fit)
poisson_100_hert = np.sqrt(100)/time_to_100_hert/hert_fit
time_to_400_hert = 400/(hert_fit)
poisson_400_hert = np.sqrt(400)/time_to_400_hert/hert_fit
time_to_1000_hert = 1000/(hert_fit)
poisson_1000_hert = np.sqrt(1000)/time_to_1000_hert/hert_fit
time_to_10000_hert = 10000/(hert_fit)
poisson_10000_hert = np.sqrt(10000)/time_to_10000_hert/hert_fit
plt.figure(figsize=(8, 6))  # Adjust figure size as needed
plt.plot(energy_midpoints[energy_range], time_to_100_hert, color='C0', label='100 counts (HERT)')
plt.errorbar(energy_midpoints[selected_indices], time_to_100_hert[selected_indices], yerr=poisson_100_hert[selected_indices]*time_to_100_hert[selected_indices], color='C0', marker='o',
             linestyle='None', capsize=5, label='Selected Energies', zorder=3)
plt.plot(energy_midpoints[energy_range], time_to_400_hert, color='C1', label='400 counts (HERT)')
plt.errorbar(energy_midpoints[selected_indices], time_to_400_hert[selected_indices], yerr=poisson_400_hert[selected_indices]*time_to_400_hert[selected_indices], color='C1', marker='o',
             linestyle='None', capsize=5, zorder=3)
plt.plot(energy_midpoints[energy_range], time_to_1000_hert, color='C2', label='1000 counts (HERT)')
plt.errorbar(energy_midpoints[selected_indices], time_to_1000_hert[selected_indices], yerr=poisson_1000_hert[selected_indices]*time_to_1000_hert[selected_indices], color='C2', marker='o',
             linestyle='None', capsize=5, zorder=3)
plt.plot(energy_midpoints[energy_range], time_to_10000_hert, color='C3', label='10000 counts (HERT)')
plt.errorbar(energy_midpoints[selected_indices], time_to_10000_hert[selected_indices], yerr=poisson_10000_hert[selected_indices]*time_to_10000_hert[selected_indices], color='C3', marker='o',
             linestyle='None', capsize=5, zorder=3)
plt.axvline(x=0.6,color='black',linestyle='--',label='HERT Threshold')
plt.xlim(0.5, 2)
plt.ylim(0, 25)
plt.yticks(np.arange(0, 26, 3))
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Time to counts (hours)')
plt.grid(True)
plt.legend()
plt.show()


# plt.figure(figsize=(8, 6))  # Adjust figure size as needed
# plt.plot(energy_midpoints[energy_range], time_to_100_hert, color='C0', label='100 counts (HERT)')
# plt.errorbar(energy_midpoints[selected_indices], time_to_100_hert[selected_indices], yerr=poisson_100_hert[selected_indices]*time_to_100_hert[selected_indices], color='C0', marker='o',
#              linestyle='None', capsize=5, label='Selected Energies', zorder=3)
# plt.plot(energy_midpoints[energy_range], time_to_1000_hert, color='C1', label='1000 counts (HERT)')
# plt.errorbar(energy_midpoints[selected_indices], time_to_1000_hert[selected_indices], yerr=poisson_1000_hert[selected_indices]*time_to_1000_hert[selected_indices], color='C1', marker='o',
#              linestyle='None', capsize=5, zorder=3)
# plt.plot(energy_midpoints[energy_range], time_to_10000_hert, color='C2', label='10000 counts (HERT)')
# plt.errorbar(energy_midpoints[selected_indices], time_to_10000_hert[selected_indices], yerr=poisson_10000_hert[selected_indices]*time_to_10000_hert[selected_indices], color='C2', marker='o',
#              linestyle='None', capsize=5, zorder=3)
# plt.axvline(x=0.6,color='black',linestyle='--',label='HERT Threshold')
# plt.xlim(0.5, 2)
# plt.ylim(0, 6)
# plt.xlabel('Kinetic Energy (MeV)')
# plt.ylabel('Time to counts (hours)')
# plt.grid(True)
# plt.legend()
# plt.show()

#%% Save time to counts at selected energies
time_to_counts = np.column_stack((energy_midpoints[selected_indices], time_to_100_hert[selected_indices], time_to_400_hert[selected_indices], time_to_1000_hert[selected_indices], time_to_10000_hert[selected_indices]))
np.savetxt('Test Theory and Data/Aero_time_to_counts.csv', time_to_counts,
    fmt='%.4f, %.4f, %.4f, %.4f, %.4f',
    delimiter=',',
    header="Energy (MeV), 100 Counts (hrs), 400 Counts (hrs), 1000 Counts (hrs), 10000 Counts (hrs)",
)

total_time = np.sum(time_to_1000_hert[selected_indices][:-2]) + time_to_400_hert[selected_indices][-2] + time_to_100_hert[selected_indices][-1]

#%% FOV rotation calculation
# this needs to be corrected based on the actual distance between the front of the instrument and the axis of rotation
d_rotation = 220 / 10 / 2 # cm, distance between axis of rotation and front of collimator
d_FOV = (63 - 32) / 10 # cm, distance between field of view and front of collimator
aperture_hw = 18 / 2 / 10 # cm, half width of the aperture of the instrument
FOV_half_angle = np.rad2deg(np.atan(aperture_hw/d_FOV))
FOV_angle = FOV_half_angle * 2

max_rotation = np.rad2deg(np.atan(aperture_hw/d_rotation))
print(f"Max Angle of Rotation for FOV Angle: {FOV_angle:.2f} degrees is {max_rotation:.2f} degrees")

# Find how angles correspond between 0 and FOV
rotation_angle_steps = np.linspace(0, max_rotation, 1000)
FOV_angle_steps = np.linspace(0, FOV_half_angle, 1000)

plt.plot(rotation_angle_steps, FOV_angle_steps)
plt.xlabel('Rotation Angle (degrees)')
plt.ylabel('FOV Half Angle (degrees)')
plt.title('Rotation Angle vs FOV Angle')
plt.grid(True)
plt.show()

set_rotation_type = 'FOV' # 'axis' or 'FOV'
if set_rotation_type == 'axis':
    set_rotation = 1 # degrees
    FOV_half_angle_rotation = (set_rotation / max_rotation) * FOV_half_angle
    print(f"For a rotation of {set_rotation:.2f} degrees, the FOV half angle is {FOV_half_angle_rotation:.2f} degrees")
elif set_rotation_type == 'FOV':
    set_FOV_half_angle = 5 # degrees
    rotation_angle = (set_FOV_half_angle / FOV_half_angle) * max_rotation
    print(f"For a FOV half angle of {set_FOV_half_angle:.2f} degrees, the rotation angle is {rotation_angle:.2f} degrees")
#%% Import Libraries
from matplotlib import ticker
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import pandas as pd
from scipy.io import loadmat

FONT_SIZE = 20
plt.rcParams.update({
    'font.size': FONT_SIZE,          # Base font size
    'axes.titlesize': FONT_SIZE + 2, # Subplot titles
    'axes.labelsize': FONT_SIZE,     # Axis labels (L*, Time, Flux)
    'xtick.labelsize': FONT_SIZE,    # X-axis tick numbers
    'ytick.labelsize': FONT_SIZE,    # Y-axis tick numbers
    'legend.fontsize': FONT_SIZE,    # Legend text
    'figure.titlesize': FONT_SIZE + 4, # Main suptitle
    'legend.markerscale': 1.0
})

#%% Import effective energies for DART Be window
# Load MAT file
mat_path = r"C:\Users\Will\Box\HERT_Box\Matlab Main\Result\DARTmod_ISO\DARTmod_ISO_new.mat"
mat_data = loadmat(mat_path)

# Extract 1D vectors
energy_edges = mat_data['energy_edges'].squeeze()
energy_midpoints = mat_data['energy_midpoints'].squeeze()
energy_channels = mat_data['energy_channels'].squeeze()

E_eff_DART = mat_data['E_eff'].squeeze()
differences = np.diff(E_eff_DART)
try:
    # We add 1 because np.diff shrinks the array by 1.
    break_index = np.where(differences <= 0)[0][0] + 1
except IndexError:
    # If the array is fully sorted, there is no negative difference.
    break_index = len(E_eff_DART)


# Extract geometric factor as guaranteed 2D matrix
geo_factor_LASP = np.atleast_2d(mat_data['geo_EC'])
geo_factor_LASP[np.isnan(geo_factor_LASP)] = 0
geo_factor_total_LASP = np.sum(geo_factor_LASP, axis=0)

bins_LASP = geo_factor_LASP.shape[-1]
bin_width = np.diff(energy_edges)

#%% Import LASP test2 data and plot
# Load the file
mat_ctrl = scipy.io.loadmat(r"C:\Users\Will\Box\HERT_Box\Sr90 Testing\HERT SR90 Test2_LASP\HERTSR90Test2_CTRL.mat")
mat = scipy.io.loadmat(r"C:\Users\Will\Box\HERT_Box\Sr90 Testing\HERT SR90 Test2_LASP\HERTSR90Test2.mat")

# Flatten the (130, 1) matrix into a 1D 130-element array
dt = 10
counts_ctrl = mat_ctrl['SummedBinCounts'].flatten()
counts_raw = mat['SummedBinCounts'].flatten()

# Extract requested variables
test_EC_counts_raw = counts_raw[:41]            # First 41 bins are energy channel counts per second
test_detector_counts_raw = counts_raw[-(2+8):-2]  # Bins 121–128 are counts per detector
duration_raw = counts_raw[-1]                       # Last index is test time in seconds

test_EC_counts_raw_error = np.sqrt(test_EC_counts_raw)
test_EC_countrate_raw = test_EC_counts_raw / dt
test_EC_countrate_raw_error = test_EC_counts_raw_error / dt

test_detector_counts_raw_error = np.sqrt(test_detector_counts_raw)
test_detector_countrate_raw = test_detector_counts_raw / dt
test_detector_countrate_raw_error = test_detector_counts_raw_error / dt

test_EC_countrate_ctrl = counts_ctrl[:41]            # First 41 bins are energy channel counts per second
test_detector_countrate_ctrl = counts_ctrl[-(2+8):-2]  # Bins 121–128 are counts per detector
duration_ctrl = counts_ctrl[-1]                       # Last index is test time in seconds

#test_EC_countrate = test_EC_countrate_raw - test_EC_countrate_ctrl
test_EC_countrate = test_EC_countrate_raw
test_EC_countrate[test_EC_countrate<0] = 0
test_EC_countrate_error = test_EC_countrate_raw_error

#test_detector_countrate = test_detector_counts_raw - test_detector_counts_ctrl
test_detector_countrate = test_detector_countrate_raw
test_detector_countrate[test_detector_countrate<0] = 0
test_detector_countrate_error = test_detector_countrate_raw_error

total_countrate_test = np.sum(test_detector_countrate)
print('Total count rate for test:', total_countrate_test)


# Plot count rate for detectors
x_locs = np.arange(1, 9)
detector_labels = ['1', '2', '3', '4', '5', '6', '7 & 8', '9']

fig, ax = plt.subplots(figsize=(10, 10))
test_detector_counts_plot = ax.bar(
    x_locs,
    test_detector_countrate,
    yerr=test_detector_countrate_error,
    width=0.8,
    color='C3',
    capsize=10)
ax.set_xlabel("Detector Number")
ax.set_ylabel("Count Rate (#/second)")
ax.tick_params(axis='both', which='major')
ax.set_xticks(x_locs)
ax.set_xticklabels(detector_labels)
ax.set_title(fr"Detector Count Rate from LASP SR90/Y90 Source Test 2",pad=16)
ax.grid(True)
#plt.yscale('log')
plt.ylim(0,100/dt)
#plt.ylim(10^-4,None)
#ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
plt.show()

# Plot count rate for channels
channels = np.arange(1, 25 + 1, dtype=int)
channel_countrate = test_EC_countrate

capsize = 6   # Proportional capsize for side-by-side bars

fig, ax1 = plt.subplots(figsize=(12, 10))
# 2. Plot Measured Counts (Shifted Right)
ax1.bar(
    channels, 
    test_EC_countrate[channels - 1], 
    yerr=test_EC_countrate_error[channels - 1], 
    color='C3', 
    edgecolor='black',
    linewidth=1.0,
    capsize=capsize,
    error_kw={
        'ecolor': 'black',
        'elinewidth': 1.5,
        'capthick': 1.5
    },
    label='Measured'
)
# Configure Bottom Axis & Log Scaling
ax1.set_xlabel('Channel Number')
ax1.set_ylabel('Count Rate (#/second)')
#ax1.set_yscale('log')
#ax1.set_ylim(10**-2, None)
ax1.set_ylim(0, 100/dt)
#ax1.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax1.set_xticks(channels)
ax1.set_xlim(0.3, len(channels) + 0.7)
ax1.grid(True, which="both", ls="--", alpha=0.5)
# 3. Create Twinned Axis for Effective Energies (Top)
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(channels)
ax2.set_xlabel('Effective Energy (MeV)', labelpad=10)
ax2.set_xticklabels(
    [f'{e:.2f}' if i % 2 == 0 else '' for i, e in enumerate(E_eff_DART[channels - 1])], 
    rotation=45, 
    ha='left'
)
plt.tight_layout()
plt.show()

#%%
'''

        IMPORT THEORETICAL LASP BETA SPECTRUM DATA

'''

#%% Import LASP spectra
LASP_spectrum = {} # Dictionary to hold LASP data
LASP_spectrum['csv_data'] = np.genfromtxt('C:/Users/Will/Box/HERT_Box/Sr90 Testing/Sr90Y90.csv', delimiter=',', filling_values=0)

# Original CSV Data (0.02 MeV bins)
LASP_spectrum['KE_orig'] = LASP_spectrum['csv_data'][:, 0]
LASP_spectrum['Sr90_CR_orig'] = LASP_spectrum['csv_data'][:, 1]
LASP_spectrum['Y90_CR_orig'] = LASP_spectrum['csv_data'][:, 2]
LASP_spectrum['combined_CR_orig'] = LASP_spectrum['csv_data'][:, 3]

# Divide by original energy bin width (0.02) to get differential flux
LASP_spectrum['Sr90_flux_orig'] = LASP_spectrum['Sr90_CR_orig'] / 0.02
LASP_spectrum['Y90_flux_orig'] = LASP_spectrum['Y90_CR_orig'] / 0.02
LASP_spectrum['combined_flux_orig'] = LASP_spectrum['combined_CR_orig'] / 0.02

#%% Interpolate onto new logarithmic midpoints
LASP_spectrum['combined_flux_interp'] = np.interp(
    energy_midpoints,
    LASP_spectrum['KE_orig'],
    LASP_spectrum['combined_flux_orig'],
    left=0.0,
    right=0.0
)

LASP_spectrum['Sr90_flux_interp'] = np.interp(
    energy_midpoints,
    LASP_spectrum['KE_orig'],
    LASP_spectrum['Sr90_flux_orig'],
    left=0.0,
    right=0.0
)

LASP_spectrum['Y90_flux_interp'] = np.interp(
    energy_midpoints,
    LASP_spectrum['KE_orig'],
    LASP_spectrum['Y90_flux_orig'],
    left=0.0,
    right=0.0
)

# Multiply by new logarithmic bin_width to get re-binned count rates
LASP_spectrum['Sr90_CR_interp'] = LASP_spectrum['Sr90_flux_interp'] * bin_width
LASP_spectrum['Y90_CR_interp'] = LASP_spectrum['Y90_flux_interp'] * bin_width
LASP_spectrum['combined_CR_interp'] = LASP_spectrum['combined_flux_interp'] * bin_width

#%% Normalize LASP spectra
total_counts_LASP = np.sum(LASP_spectrum['combined_CR_interp'])

LASP_spectrum['Sr90_CR_norm'] = LASP_spectrum['Sr90_CR_interp'] / total_counts_LASP
LASP_spectrum['Y90_CR_norm'] = LASP_spectrum['Y90_CR_interp'] / total_counts_LASP
LASP_spectrum['combined_CR_norm'] = LASP_spectrum['combined_CR_interp'] / total_counts_LASP

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

today_count_rate = total_counts_LASP / 2 * np.exp(-lda * (2026 - 2010)) / HERT_eff / directional_scaling # Initial count rate in 2010
d = 0.15736  # distance between front of collimator and radiation source, cm
s = 0.3/2  # radius of radiation source, cm
LASP_today_scale = 2*today_count_rate/(1-np.cos(np.atan(s/(d+6.3))))/(4*np.pi*(d+6.3)**2)

LASP_spectrum['Sr90_today'] = LASP_spectrum['Sr90_CR_norm'] * LASP_today_scale / bin_width
LASP_spectrum['Y90_today'] = LASP_spectrum['Y90_CR_norm'] * LASP_today_scale / bin_width
LASP_spectrum['combined_today'] = LASP_spectrum['combined_CR_norm'] * LASP_today_scale / bin_width

# Plot Today scaled Sr90 and Y90
fig, ax = plt.subplots(figsize=(14, 4)) 
combined_scatter = ax.scatter(energy_midpoints, LASP_spectrum['combined_today'], label="LASP", color='C2')

ax.set_xlim(0, 2.5)
ax.set_xlabel("Kinetic Energy (MeV)")
ax.set_ylabel(r'Flux ($\# / \text{s } \text{sr } \text{cm}^2 \text{ MeV}$)')
ax.tick_params(axis='both', which='both', labelsize=14)
#ax.set_title("Today's LASP Beta Decay Spectra of Sr90 and Y90")
ax.grid(True)

ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())

# 1. Isolate the effective energies that fit within your data bounds
valid_E_eff = E_eff_DART[E_eff_DART <= energy_midpoints[-1]]

# 2. SLICE the array to start at Channel 3 (Index 2 in a 0-indexed array)
valid_E_eff_ch3_up = valid_E_eff[2:]

# 3. Set the tick positions using only the sliced array
tick_positions = energy_midpoints[np.searchsorted(energy_midpoints, valid_E_eff_ch3_up)]
ax2.set_xticks(tick_positions)

# 4. Generate the labels, keeping track of the offset channel numbers
channel_numbers = np.arange(3, 3 + len(valid_E_eff_ch3_up))
channel_labels = []

for ch in channel_numbers:
    if ch % 5 == 0:  # Preserving your logic to label every 5th channel
        channel_labels.append(str(ch))
    else:
        channel_labels.append('')

ax2.set_xticklabels(channel_labels)
ax2.set_xlabel("Channel Number", labelpad=10)
ax2.tick_params(axis='x', which='major', labelsize=14)

plt.ylim(10**-2,10**3)
plt.yscale('log')
major_ticks = 10**np.arange(np.log10(10**-2), np.log10(10**3) + 1)
ax.set_yticks(major_ticks)
ax.legend()
ax.set_title(r"Beta Decay Spectra of Sr$^{90}$ and Y$^{90}$", pad=12)
plt.show()

#%% Calculate Expected LASP count rates
LASP_count_rate_EC = np.zeros(geo_factor_LASP.shape[0])
for channel in range(geo_factor_LASP.shape[0]):
    LASP_count_rate_EC[channel] = np.sum(geo_factor_LASP[channel,:]*LASP_spectrum['combined_today']*bin_width)

# Plot expected count rates
fig, ax = plt.subplots(figsize=(16, 3))
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
ax2.set_xticks(energy_midpoints[np.searchsorted(energy_midpoints, E_eff_DART[E_eff_DART <= energy_midpoints[-1]])])
channel_labels = [''] * len(E_eff_DART[E_eff_DART <= energy_midpoints[-1]])
for i in range(len(E_eff_DART[E_eff_DART <= energy_midpoints[-1]])):
    if (i + 1) % 5 == 0:  # Label every 5th channel
        channel_labels[i] = str(i + 1)
ax2.set_xticklabels(channel_labels)
ax2.set_xlabel("Channel Number", labelpad=10) # labelpad moves the label up
ax2.tick_params(axis='x', which='major')
plt.yscale('log')
plt.ylim(1e-3, None)
plt.show()

#%% Compare theoretical and measured counts for each energy channel
channels = np.arange(1, 41 + 1, dtype=int)
last_channel = 26
channel_counts = LASP_count_rate_EC * dt
channel_counts_error = np.sqrt(channel_counts)
channel_countrate = channel_counts / dt
channel_countrate_error = channel_counts_error / dt

width = 0.38  # Width of each bar
capsize = 6   # Proportional capsize for side-by-side bars

fig, ax1 = plt.subplots(figsize=(12, 5))

# 1. Plot Theoretical Counts (Shifted Left)
ax1.bar(
    channels - width / 2, 
    channel_countrate, 
    width=width,
    yerr=channel_countrate_error, 
    color='C2', 
    edgecolor='black',
    linewidth=1.0,
    capsize=capsize,
    error_kw={
        'ecolor': 'black',
        'elinewidth': 1.5,
        'capthick': 1.5
    },
    label='Theoretical'
)

# 2. Plot Measured Counts (Shifted Right)
ax1.bar(
    channels + width / 2, 
    test_EC_countrate, 
    width=width,
    yerr=test_EC_countrate_error, 
    color='C3', 
    edgecolor='black',
    linewidth=1.0,
    capsize=capsize,
    error_kw={
        'ecolor': 'black',
        'elinewidth': 1.5,
        'capthick': 1.5
    },
    label='Measured'
)

# Configure Bottom Axis & Log Scaling
ax1.set_xlabel('Channel Number')
ax1.set_ylabel('Count Rate (#/second)')
ax1.set_yscale('log')
ax1.set_ylim(10**-2, None)
ax1.set_xticks(channels)
ax1.set_xlim(0.3, last_channel + 0.7)
ax1.grid(True, which="major", ls="--", alpha=0.5)

# 3. Create Twinned Axis for Effective Energies (Top)
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())

# Slice the arrays so we only place ticks/labels for the visible channels
visible_channels = channels[:last_channel]
visible_energies = E_eff_DART[:last_channel]
ax2.set_xticks(visible_channels)
ax2.set_xlabel('Effective Energy (MeV)', labelpad=10)
ax2.set_xticklabels(
    [f'{e:.2f}' if i % 2 == 0 else '' for i, e in enumerate(visible_energies)], 
    rotation=45, 
    ha='left'
)
ax1.legend(loc='upper right')
plt.tight_layout()
plt.show()

# Mask: True ONLY if BOTH theoretical and measured counts are > 0
valid_channel_mask = (channel_countrate > 0.1) & (test_EC_countrate > 0.1)

with np.errstate(divide='ignore', invalid='ignore'):
    # 1-sigma Poisson % uncertainties (requires count > 0)
    theo_stat_err = np.where(
        channel_countrate > 0, 
        (np.sqrt(channel_countrate) / channel_countrate) * 100, 
        np.nan
    )
    meas_stat_err = np.where(
        test_EC_countrate > 0, 
        (np.sqrt(test_EC_countrate) / test_EC_countrate) * 100, 
        np.nan
    )
    
    # % Relative Difference (returns NaN if EITHER count is 0)
    rel_diff_pct = np.where(
        valid_channel_mask, 
        ((test_EC_countrate - channel_countrate) / channel_countrate) * 100, 
        np.nan
    )

# Formatters for displaying clean 'N/A' instead of NaN
def fmt_pct(val):
    return 'N/A' if np.isnan(val) or np.isinf(val) else f"{val:+.2f}%"

def fmt_stat_err(val):
    return 'N/A' if np.isnan(val) or np.isinf(val) else f"{val:.2f}%"

# Construct DataFrame
df_summary = pd.DataFrame({
    'Channel': channels[valid_channel_mask],
    'Theoretical Count Rate': channel_countrate[valid_channel_mask],
    'Measured Count Rate': test_EC_countrate[valid_channel_mask],
    '% Rel Difference': rel_diff_pct[valid_channel_mask]
})

print(df_summary.to_string(index=False, formatters={
    'Theoretical Count Rate': '{:,.0f}'.format,
    'Measured Count Rate': '{:,.0f}'.format,
    '% Rel Difference': fmt_pct
}))

#%% Compare theoretical and measured counts for each detectors
detector_number = np.linspace(1,9,9)
detector_efficiency = np.genfromtxt(r"C:\Users\Will\Box\HERT_Box\Sr90 Testing\efficiency_Detectors_LASPTest.txt")
detector_efficiency = detector_efficiency/detector_efficiency[0]
raw_detector_countrate = detector_efficiency * today_count_rate
detector_countrate = np.array([
    *raw_detector_countrate[:6],                        # Detectors 1-6
    raw_detector_countrate[6] + raw_detector_countrate[7], # Detectors 7 & 8 combined
    raw_detector_countrate[8]                           # Detector 9
])
detector_counts = detector_countrate * dt
detector_counts_error = np.sqrt(detector_counts)
detector_countrate_error = detector_counts_error / dt

x_locs = np.arange(1, 9)
detector_labels = ['1', '2', '3', '4', '5', '6', '7 & 8', '9']

width = 0.38  # Width of each bar
capsize = 10   # Proportional capsize for side-by-side bars

fig, ax = plt.subplots(figsize=(12, 4.5))

# Plot Theory (Left)
theory_bar = ax.bar(
    x_locs - width / 2, 
    detector_countrate, 
    width=width,
    yerr=detector_countrate_error, 
    color='C2', 
    error_kw={
        'elinewidth': 1.5,  # Thicker vertical error lines
        'capthick': 1.5     # Thicker horizontal cap lines
    },
    ecolor='black',
    capsize=capsize,
    label='Theoretical'
)

# Plot Measurement (Right)
measure_bar = ax.bar(
    x_locs + width / 2, 
    test_detector_countrate, 
    width=width,
    yerr=test_detector_countrate_error, 
    color='C3', 
    error_kw={
        'elinewidth': 1.5,  # Thicker vertical error lines
        'capthick': 1.5     # Thicker horizontal cap lines
    },
    ecolor='black',
    capsize=capsize,
    label='Measured'
)

ax.set_xlabel("Detector Number")
ax.set_ylabel("Count Rate (#/second)")
ax.set_xticks(x_locs)
ax.set_xticklabels(detector_labels)
ax.set_yscale('log')
plt.ylim(10**-4,None)
ax.grid(True, which="major", ls="--", alpha=0.5)
ax.legend()
plt.show()

# Mask: True ONLY if BOTH theoretical and measured counts are > 0
valid_diff_mask = (detector_countrate > 0) & (test_detector_countrate > 0)

with np.errstate(divide='ignore', invalid='ignore'):
    # 1-sigma Poisson % uncertainties (requires count > 0)
    theo_stat_err = np.where(
        detector_countrate > 0, 
        (np.sqrt(detector_countrate) / detector_countrate) * 100, 
        np.nan
    )
    meas_stat_err = np.where(
        test_detector_countrate > 0, 
        (np.sqrt(test_detector_countrate) / test_detector_countrate) * 100, 
        np.nan
    )
    
    # % Relative Difference (returns NaN if EITHER count is 0)
    rel_diff_pct = np.where(
        valid_diff_mask, 
        ((test_detector_countrate - detector_countrate) / detector_countrate) * 100, 
        np.nan
    )

# Formatters for displaying clean 'N/A' instead of NaN
def fmt_pct(val):
    return 'N/A' if np.isnan(val) or np.isinf(val) else f"{val:+.2f}%"

def fmt_stat_err(val):
    return 'N/A' if np.isnan(val) or np.isinf(val) else f"{val:.2f}%"

# Construct DataFrame
df_summary = pd.DataFrame({
    'Detector': detector_labels,
    'Theoretical Count Rate': detector_countrate,
    'Measured Count Rate': test_detector_countrate,
    '% Rel Difference': rel_diff_pct
})

print(df_summary.to_string(index=False, formatters={
    'Theoretical Count Rate': '{:,.0f}'.format,
    'Measured Count Rate': '{:,.0f}'.format,
    '% Rel Difference': fmt_pct
}))
#%% Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.transforms as transforms
from scipy.optimize import curve_fit

#%% Core Utility & Classification Functions
def classify_s_scale(pfu):
    """
    Classifies the event into the NOAA Solar Radiation Storm Scale based on peak >10 MeV pfu.
    """
    if pfu >= 100000:
        return "S5"
    elif pfu >= 10000:
        return "S4"
    elif pfu >= 1000:
        return "S3"
    elif pfu >= 100:
        return "S2"
    elif pfu >= 10:
        return "S1"
    return "Below S1"

def log_interp(x_target, x_source, y_source):
    """
    Log-log interpolation that safely handles zero/negative values.
    """
    valid_mask = y_source > 0
    x_valid = x_source[valid_mask]
    y_valid = y_source[valid_mask]
    log_y_target = np.interp(np.log10(x_target), np.log10(x_valid), np.log10(y_valid), left=-np.inf, right=-np.inf)
    return 10 ** log_y_target

def plot_channel_bar_chart(range_data, pen_reversed_data, title, ylabel, ylim, center_text, 
                           colors_range, colors_pen_rev, range_x, pen_x, tick_locs, tick_labels,
                           num_range, pen_label_y=0.95):
    """
    Unified plotting helper for count-rate and total-count bar charts.
    Preserves exact font sizes, labels, transforms, and bounding boxes.
    """
    plt.figure(figsize=(12, 4.5))
    plt.bar(range_x, range_data, color=colors_range)
    plt.bar(pen_x, pen_reversed_data, color=colors_pen_rev)
    
    plt.axvline(x=num_range + 1, color='black', linestyle='--', linewidth=2)
    
    trans = transforms.blended_transform_factory(plt.gca().transData, plt.gca().transAxes)
    bbox_props = dict(facecolor='white', alpha=0.8, edgecolor='none')
    
    plt.text(13.5, 0.95, 'Range Channels', transform=trans, ha='center', va='top', 
             fontsize=16, color='black', bbox=bbox_props)
    plt.text(40.5, pen_label_y, 'Penetrating Channels', transform=trans, ha='center', va='top', 
             fontsize=16, color='black', bbox=bbox_props)
    
    plt.text(0.5, 0.95, center_text, transform=plt.gca().transAxes, 
             verticalalignment='top', horizontalalignment='center', fontsize=16, bbox=dict(facecolor='white'))
    
    plt.xlabel('Energy Channel Index', fontsize=16)
    plt.xlim(0, len(range_x) + len(pen_x) + 2)
    plt.ylabel(ylabel, fontsize=16)
    plt.yscale('log')
    plt.ylim(ylim)
    
    plt.xticks(ticks=tick_locs, labels=tick_labels, fontsize=16)
    plt.yticks(fontsize=16)
    plt.title(title, fontsize=20)
    plt.tight_layout()
    plt.show()

#%% Trapped Protons AP9 model
data_dir = r"C:\Users\wzt0020\Documents\OMERE 5.9"
files_ap9 = {
    'Mean': 'trappedProtons_mean.flx',
    '75th Percentile': 'trappedProtons_CORA_75.flx',
    '95th Percentile': 'trappedProtons_CORA_95.flx'
}

ae9_raw_flux = {}
ae9_energy = None 

plt.figure(figsize=(12, 5))
for label, filename in files_ap9.items():
    filepath = os.path.join(data_dir, filename)
    try:
        data = np.loadtxt(filepath, comments='#')
        ae9_energy = data[:, 0]
        flux = data[:, 1] / (4 * np.pi) 
        ae9_raw_flux[label] = flux
        plt.loglog(ae9_energy, flux, marker='.', markersize=12, linestyle='-', linewidth=2, label=label)
    except Exception as e:
        print(f"Error processing {filename}: {e}")

plt.xscale('log')
plt.xlim(10, 1000)
custom_ticks = [10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 400, 500, 750, 1000]
plt.xticks(ticks=custom_ticks, labels=custom_ticks, fontsize=16, rotation=45)
plt.yscale('log')
plt.ylim(1e-4, 1e2)
plt.yticks(fontsize=16)
plt.xlabel('Energy (MeV)', fontsize=16)
plt.ylabel('Differential Flux (cm$^{-2}$ sr$^{-1}$ s$^{-1}$ MeV$^{-1}$)', fontsize=16, y=0.45)
plt.grid(True, which="major", ls="--", alpha=0.5)
plt.legend(fontsize=16)
plt.title('AP9 Trapped Protons Differential Flux', fontsize=20)
plt.tight_layout()
plt.show()

#%% SEP ESP model
files_esp = {
    'ESP Mean': 'solarMeanProtons_CORA_mean.flx',
    'ESP 75th Percentile': 'solarMeanProtons_CORA_75.flx', 
    'ESP 95th Percentile': 'solarMeanProtons_CORA_95.flx'
}

esp_raw_flux = {}
esp_energy = None

plt.figure(figsize=(12, 6))
for label, filename in files_esp.items():
    filepath = os.path.join(data_dir, filename)
    try:
        data = np.loadtxt(filepath, comments='#')
        esp_energy = data[:, 0]
        flux = data[:, 1] / (4 * np.pi) 
        esp_raw_flux[label] = flux
        plt.loglog(esp_energy, flux, marker='.', markersize=12, linestyle='-', linewidth=2, label=label)
    except Exception as e:
        print(f"Error processing {filename}: {e}")

plt.xscale('log')
plt.xlim(10, 300)
custom_ticks_esp = [10, 20, 30, 50, 75, 100, 150, 200, 250, 300]
plt.xticks(ticks=custom_ticks_esp, labels=custom_ticks_esp, fontsize=16, rotation=45)
plt.yscale('log')
plt.ylim(1e-5, 1e1)
plt.yticks(fontsize=16)
plt.xlabel('Energy (MeV)', fontsize=16)
plt.ylabel('Differential Flux (cm$^{-2}$ sr$^{-1}$ s$^{-1}$ MeV$^{-1}$)', fontsize=16)
plt.grid(True, which="major", ls="--", alpha=0.5)
plt.legend(fontsize=16)
plt.title('ESP Solar Protons Differential Flux', fontsize=20)
plt.tight_layout()
plt.show()

#%% SEP SAPPHIRE Model
files_sapphire = {
    '1 per year': 'solarFlareProtons_GEO_50CI_pt5yr.flx',
    '1 per 4 years': 'solarFlareProtons_GEO_75CI_1yr.flx',
    '1 per 20 years': 'solarFlareProtons_GEO_95CI_1yr.flx'
}

sapphire_raw_flux = {}
sapphire_energy = None
sapphire_pfu_results = {}

plt.figure(figsize=(12, 5))
for label, filename in files_sapphire.items():
    filepath = os.path.join(data_dir, filename)
    try:
        data = np.loadtxt(filepath, comments='#')
        sapphire_energy = data[:, 0]
        flux = data[:, 1] / (4 * np.pi)
        sapphire_raw_flux[label] = flux

        # Numerical Integration for >10 MeV pfu
        idx_10mev = np.where(sapphire_energy >= 10.0)[0]
        if len(idx_10mev) > 0:
            e_slice = sapphire_energy[idx_10mev]
            flux_slice = flux[idx_10mev]
            calculated_pfu = np.trapezoid(flux_slice, e_slice)
            s_class = classify_s_scale(calculated_pfu)
            sapphire_pfu_results[label] = {'pfu': calculated_pfu, 'class': s_class}
            print(f"{label:15s} | Integrated Flux: {calculated_pfu:8.1f} pfu | Classification: {s_class}")
        else:
            print(f"{label:15s} | No data points >= 10 MeV found.")
            
        plt.loglog(sapphire_energy, flux, marker='.', markersize=12, linestyle='-', linewidth=2, label=label)
    except Exception as e:
        print(f"Error processing {filename}: {e}")

plt.xscale('log')
plt.xlim(10, 1000)
plt.xticks(ticks=custom_ticks, labels=custom_ticks, fontsize=16, rotation=45)
plt.yscale('log')
plt.ylim(1e-4, 1e4)
plt.yticks(fontsize=16)
plt.xlabel('Energy (MeV)', fontsize=16)
plt.ylabel('Differential Flux (cm$^{-2}$ sr$^{-1}$ s$^{-1}$ MeV$^{-1}$)', fontsize=16, y=0.45)
plt.grid(True, which="major", ls="--", alpha=0.5)
plt.legend(fontsize=16)
plt.title('SAPPHIRE Solar Protons Differential Flux', fontsize=20)
plt.tight_layout()
plt.show()

#%% Test Band Fit
def band_flux(E, A, gamma1, gamma2, E0):
    """
    Standard Band function for SEP differential flux.
    """
    E_break = (gamma2 - gamma1) * E0
    E = np.asarray(E, dtype=float)
    flux = np.zeros_like(E)
    
    mask_low = (E < E_break)
    flux[mask_low] = A * (E[mask_low] ** -gamma1) * np.exp(-E[mask_low] / E0)
    
    mask_high = ~mask_low
    continuity_factor = A * ((gamma2 - gamma1) * E0) ** (gamma2 - gamma1) * np.exp(gamma1 - gamma2)
    flux[mask_high] = continuity_factor * (E[mask_high] ** -gamma2)
    return flux

def log_band_flux(E, log10_A, gamma1, gamma2, E0):
    """
    Log10 wrapper around the Band function for stable least-squares fitting.
    """
    A = 10.0 ** log10_A
    y = band_flux(E, A, gamma1, gamma2, E0)
    return np.log10(np.clip(y, 1e-30, np.inf))

sapphire_flux_select = sapphire_raw_flux['1 per 4 years'] 
valid_mask = (sapphire_energy > 0) & (sapphire_flux_select > 0)
E_fit = sapphire_energy[valid_mask]
y_fit = np.log10(sapphire_flux_select[valid_mask])

p0 = [4.0, 1.5, 4.0, 30.0]
bounds = ([-10.0, 0.1, 2.0, 1.0], [15.0, 3.0, 10.0, 500.0])

try:
    popt, pcov = curve_fit(
        log_band_flux, 
        E_fit, 
        y_fit, 
        p0=p0, 
        bounds=bounds,
        maxfev=10000
    )
    A_best = 10.0 ** popt[0]
    gamma1_best = popt[1]
    gamma2_best = popt[2]
    E0_best = popt[3]
    E_break_best = (gamma2_best - gamma1_best) * E0_best
    
    print("="*50)
    print("BEST BAND FUNCTION FIT PARAMETERS ('1 per year')")
    print("="*50)
    print(f"Amplitude (A)        : {A_best:.3e}")
    print(f"Low-E index (γ1)     : {gamma1_best:.3f}")
    print(f"High-E index (γ2)    : {gamma2_best:.3f}")
    print(f"Roll-off energy (E0) : {E0_best:.2f} MeV")
    print(f"Transition (E_break) : {E_break_best:.2f} MeV")
    print("="*50)
except RuntimeError as e:
    print(f"Band fit failed to converge: {e}")

E_smooth = np.logspace(np.log10(E_fit.min()), np.log10(E_fit.max()), 200)
flux_band_smooth = band_flux(E_smooth, A_best, gamma1_best, gamma2_best, E0_best)

plt.figure(figsize=(10, 6))
plt.loglog(E_fit, sapphire_flux_select[valid_mask], 
           marker='o', markersize=8, linestyle='none', color='black', 
           label="OMERE SAPPHIRE ('1 per year' Data)")
plt.loglog(E_smooth, flux_band_smooth, 
           linestyle='-', linewidth=3, color='crimson', 
           label=rf"Band Fit ($\gamma_1$={gamma1_best:.2f}, $\gamma_2$={gamma2_best:.2f}, $E_0$={E0_best:.1f} MeV)")
plt.axvline(E_break_best, color='blue', linestyle='--', alpha=0.7, 
            label=rf"Transition Break ($E_{{break}}$ = {E_break_best:.1f} MeV)")

plt.xlim(1, 1000)
plt.xlabel('Energy (MeV)', fontsize=14)
plt.ylabel('Differential Flux (cm$^{-2}$ sr$^{-1}$ s$^{-1}$ MeV$^{-1}$)', fontsize=14)
plt.title("Band Function Fit — SAPPHIRE '1 per year' Spectrum", fontsize=16)
plt.grid(True, which="both", ls="--", alpha=0.4)
plt.legend(fontsize=12, loc="lower left")
plt.tight_layout()
plt.show()

#%% Plot modified geometric factor
loaded_geo = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\Proton_FS\proton_FS_17000_v6comb_range_GFbyEC.txt")
energy_midpoints_geo = loaded_geo[:, 0]
loaded_geo_values = loaded_geo[:, 1:]

loaded_geo_pen = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\Proton_FS\proton_FS_17000_v6comb_pen_GFbyEC.txt")
loaded_geo_pen_values = loaded_geo_pen[:, 1:] 
energy_midpoints_geo_pen = energy_midpoints_geo 

energy_channels = np.loadtxt(r"D:\HERT_Drive\Matlab Main\Result\channel_select\proton_channels_v6_combined.txt", delimiter=',')

num_channels_range = loaded_geo_values.shape[1] 
num_channels_pen = loaded_geo_pen_values.shape[1] 
color_denom = max(1, num_channels_range - 1)

cmap_range_fn = plt.get_cmap('plasma')
cmap_pen_fn = plt.get_cmap('winter')

plt.figure(figsize=(12, 8))
for i in range(num_channels_range):
    label_str = 'Range Channels (Solid)' if i == 0 else None
    plt.loglog(energy_midpoints_geo, loaded_geo_values[:, i], color=cmap_range_fn(i / color_denom), linewidth=1.5, linestyle='-', label=label_str)

for i in range(num_channels_pen):
    label_str = 'Penetrating Channels (Dashed)' if i == 0 else None
    plt.loglog(energy_midpoints_geo_pen, loaded_geo_pen_values[:, i], color=cmap_pen_fn(i / color_denom), linewidth=1.5, linestyle='-', alpha=0.8, label=label_str)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy (MeV)', fontsize=16)
plt.ylabel('Geometric Factor (cm$^2$ sr)', fontsize=16)
plt.title('HERT Range vs. Penetrating Geometric Factors', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend(fontsize=14, loc='best')
plt.tight_layout()
plt.show()

# --- Pre-compute Shared Bar-Chart Parameters & Colormaps ---
colors_range = plt.colormaps['plasma'](np.linspace(0, 1, num_channels_range))
full_winter = plt.colormaps['winter'](np.linspace(0, 1, len(energy_channels)))
colors_pen = full_winter[-num_channels_pen:]
colors_pen_reversed = colors_pen[::-1]

range_x_axis = list(range(1, num_channels_range + 1))
pen_x_axis = list(range(num_channels_range + 2, num_channels_range + num_channels_pen + 2))
tick_locs = [1, 6, 11, 16, 21, 26, 28, 33, 38, 43, 48, 53]
tick_labels = [1, 6, 11, 16, 21, 26, 26, 21, 16, 11, 6, 1]
# -----------------------------------------------------------

#%% Interpolate spectra for ALL models and percentiles
ae9_spectra = {}
sapphire_spectra = {}

for label, raw_flux in ae9_raw_flux.items():
    temp_interp = log_interp(energy_midpoints_geo, ae9_energy, raw_flux)
    temp_interp[energy_midpoints_geo < 10] = 0.0
    ae9_spectra[label] = temp_interp

for label, raw_flux in sapphire_raw_flux.items():
    temp_interp = log_interp(energy_midpoints_geo, sapphire_energy, raw_flux)
    temp_interp[energy_midpoints_geo < 10] = 0.0
    sapphire_spectra[label] = temp_interp

plt.figure(figsize=(12, 6))
if 'Mean' in ae9_spectra:
    mask = ae9_spectra['Mean'] > 0
    plt.loglog(energy_midpoints_geo[mask], ae9_spectra['Mean'][mask], marker='.', markersize=12, linestyle='-', linewidth=2, label='AE9 Mean (Interp)')

if 'Mean' in sapphire_spectra:
    mask = sapphire_spectra['Mean'] > 0
    plt.loglog(energy_midpoints_geo[mask], sapphire_spectra['Mean'][mask], marker='.', markersize=12, linestyle='-', linewidth=2, label='SAPPHIRE Mean (Interp)')

plt.xscale('log')
plt.xlim(10, 1000)  
plt.xticks(ticks=custom_ticks, labels=custom_ticks, fontsize=16, rotation=45)
plt.yscale('log')
plt.ylim(1e-8, 1e2)
plt.xlabel('Energy (MeV)', fontsize=16)
plt.ylabel('Differential Flux (cm$^{-2}$ sr$^{-1}$ s$^{-1}$ MeV$^{-1}$)', fontsize=16)
plt.title('Log-Interpolated Mean Spectra for AE9, ESP, and SAPPHIRE', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True, which="major", ls="--", alpha=0.5)
plt.legend(fontsize=14, loc='best')
plt.tight_layout()
plt.show()

#%% AE9: Calculate channel count rate & total counts
print("--- Generating Plots for AP9 (Trapped Protons) ---")
CI_select_ae9 = '95th Percentile'

# Fast vectorized dot products
ae9_channel_count_rates = ae9_spectra[CI_select_ae9] @ loaded_geo_values
ae9_channel_count_rates_pen = ae9_spectra[CI_select_ae9] @ loaded_geo_pen_values
ae9_pen_reversed = ae9_channel_count_rates_pen[::-1]

# Plot Count Rates
plot_channel_bar_chart(
    ae9_channel_count_rates, ae9_pen_reversed,
    title='AP9 (Trapped Protons): Count Rates for Each Energy Channel',
    ylabel='Count Rate (counts/s)', ylim=(1e-4, None),
    center_text=f'CI: {CI_select_ae9}',
    colors_range=colors_range, colors_pen_rev=colors_pen_reversed,
    range_x=range_x_axis, pen_x=pen_x_axis, tick_locs=tick_locs, tick_labels=tick_labels,
    num_range=num_channels_range
)

# Plot Overlapping Bar Chart (1 Minute Duration)
t_sec_1m = 60
counts_1min_range = ae9_channel_count_rates * t_sec_1m
counts_1min_pen_reversed = ae9_pen_reversed * t_sec_1m

plot_channel_bar_chart(
    counts_1min_range, counts_1min_pen_reversed,
    title='AP9 (Trapped Protons): Total Counts per Channel',
    ylabel='Total Counts', ylim=(1e0, 1e5),
    center_text=f'CI: {CI_select_ae9}\nDuration: 1 min',
    colors_range=colors_range, colors_pen_rev=colors_pen_reversed,
    range_x=range_x_axis, pen_x=pen_x_axis, tick_locs=tick_locs, tick_labels=tick_labels,
    num_range=num_channels_range
)

#%% SAPPHIRE: Calculate channel count rate & total counts
print("--- Generating Plots for SAPPHIRE (Solar Protons) ---")
freq_select = '1 per 4 years'

# Fast vectorized dot products
sapphire_channel_count_rates = sapphire_spectra[freq_select] @ loaded_geo_values
sapphire_channel_count_rates_pen = sapphire_spectra[freq_select] @ loaded_geo_pen_values
sapphire_pen_reversed = sapphire_channel_count_rates_pen[::-1]

# Plot Count Rates
plot_channel_bar_chart(
    sapphire_channel_count_rates, sapphire_pen_reversed,
    title='SAPPHIRE (Solar Protons): Count Rates for Each Energy Channel',
    ylabel='Count Rate (counts/s)', ylim=(1e-6, None),
    center_text=f'Frequency: {freq_select}',
    colors_range=colors_range, colors_pen_rev=colors_pen_reversed,
    range_x=range_x_axis, pen_x=pen_x_axis, tick_locs=tick_locs, tick_labels=tick_labels,
    num_range=num_channels_range
)

# Plot Overlapping Bar Chart (10 Minute Duration)
t_sec_10m = 600
counts_10min_range = sapphire_channel_count_rates * t_sec_10m
counts_10min_pen_reversed = sapphire_pen_reversed * t_sec_10m

plot_channel_bar_chart(
    counts_10min_range, counts_10min_pen_reversed,
    title='SAPPHIRE (Solar Protons): Total Counts per Channel',
    ylabel='Total Counts', ylim=(1e0, 1e10),
    center_text=f'Frequency: {freq_select}\nDuration: 10 min',
    colors_range=colors_range, colors_pen_rev=colors_pen_reversed,
    range_x=range_x_axis, pen_x=pen_x_axis, tick_locs=tick_locs, tick_labels=tick_labels,
    num_range=num_channels_range, pen_label_y=0.9
)

#%% Create-a-spectra using the Weibull fits for events from Laurenza 2015
def weibull_fit_log(E, ln_k, E_0, b):
    return ln_k + (b - 1) * np.log(E) - (E / E_0)**b

def solve_weibull_k(target_pfu, E_0, b):
    """
    Solves for k and ln_k given a target >10 MeV integral flux (pfu), E_0, and b.
    """
    ln_k = np.log(target_pfu) + np.log(b) - b * np.log(E_0) + (10.0 / E_0)**b
    k = np.exp(ln_k)
    return k, ln_k

events = {
    '26 December 2001': {'target_pfu': 2080, 'E_0': 6.1, 'b': 0.89}, #S3, from ACE >10MeV at 2001 12 26 0830
    '21 March 2011':    {'target_pfu': 21.9,  'E_0': 0.3, 'b': 0.45} #S1, from ACE >10MeV at 2011 03 21 2355
}

weibull_spectra = {}
for event, params in events.items():
    k, ln_k = solve_weibull_k(params['target_pfu'], params['E_0'], params['b'])
    s_class = classify_s_scale(params['target_pfu'])
    
    params['ln_k'] = ln_k
    params['k'] = k
    params['s_class'] = s_class
    print(f"{event}: {s_class} (pfu: {params['target_pfu']}), ln_k = {ln_k:.2f}, k = {k:.2e}")
    weibull_spectra[event] = np.exp(weibull_fit_log(energy_midpoints_geo, ln_k, params['E_0'], params['b']))

plt.figure(figsize=(12, 6))
for event, spectrum in weibull_spectra.items():
    plt.loglog(energy_midpoints_geo, spectrum, marker='.', markersize=12, linestyle='-', linewidth=2, label=event)
plt.xscale('log')
plt.xlim(10, 1000)
plt.xticks(ticks=custom_ticks, labels=custom_ticks, fontsize=16, rotation=45)
plt.yscale('log')
plt.ylim(1e-12, 1e4)
plt.xlabel('Energy (MeV)', fontsize=16)
plt.ylabel('Differential Flux (cm$^{-2}$ sr$^{-1}$ s$^{-1}$ MeV$^{-1}$)', fontsize=16)
plt.title('Weibull Spectra for Selected SEP Events', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True, which="major", ls="--", alpha=0.5)
plt.legend(fontsize=14, loc='best')

#%% Calculate channel count rate for Weibull Fits (Solar Protons)
print("--- Generating Plots for Weibull Fits (Solar Protons) ---")
event_select = '21 March 2011'

# Fast vectorized dot products
weibull_count_rates = weibull_spectra[event_select] @ loaded_geo_values
weibull_count_rates_pen = weibull_spectra[event_select] @ loaded_geo_pen_values
weibull_pen_reversed = weibull_count_rates_pen[::-1]

# Plot Count Rates
plot_channel_bar_chart(
    weibull_count_rates, weibull_pen_reversed,
    title=f'SEP Solar Protons: {event_select}',
    ylabel='Count Rate (counts/s)', ylim=(1e-6, None),
    center_text=f'Event Class: {events[event_select]["s_class"]}\nCadence: 10 min',
    colors_range=colors_range, colors_pen_rev=colors_pen_reversed,
    range_x=range_x_axis, pen_x=pen_x_axis, tick_locs=tick_locs, tick_labels=tick_labels,
    num_range=num_channels_range
)

# Plot Overlapping Bar Chart (10 Minute Duration)
counts_10min_range = weibull_count_rates * t_sec_10m
counts_10min_pen_reversed = weibull_pen_reversed * t_sec_10m

plot_channel_bar_chart(
    counts_10min_range, counts_10min_pen_reversed,
    title=f'SEP Solar Protons: {event_select}',
    ylabel='Total Counts', ylim=(1e0, 1e7),
    center_text=f'Event Class: {events[event_select]["s_class"]}\nCadence: 10 min',
    colors_range=colors_range, colors_pen_rev=colors_pen_reversed,
    range_x=range_x_axis, pen_x=pen_x_axis, tick_locs=tick_locs, tick_labels=tick_labels,
    num_range=num_channels_range
)

#%% Redo, now with a simple power law fit given by the SEPEM database of events
events_powerlaw = {
    '17 May 2012': {'target_pfu': 255, 'k':2.22e3, 'alpha': -1.81},
    '16 June 2012': {'target_pfu': 14, 'k':1.14e3, 'alpha': -2.45}, 
    '19 February 2014': {'target_pfu': 22, 'k':9.07e2, 'alpha': -2.23},
    '18 April 2014': {'target_pfu': 58,  'k':1.3e4, 'alpha': -2.78}
}

def power_law(E, k, alpha):
    return k * E**alpha

powerlaw_spectra = {}
for event, params in events_powerlaw.items():
    s_class = classify_s_scale(params['target_pfu'])
    params['s_class'] = s_class
    print(f"{event}: {s_class} (pfu: {params['target_pfu']}), k = {params['k']:.2e}, alpha = {params['alpha']:.2f}")
    powerlaw_spectra[event] = power_law(energy_midpoints_geo, params['k'], params['alpha'])

plt.figure(figsize=(12, 6))
for event, spectrum in powerlaw_spectra.items():
    plt.loglog(energy_midpoints_geo, spectrum, marker='.', markersize=12, linestyle='-', linewidth=2, label=event)
plt.xscale('log')
plt.xlim(10, 1000)
plt.xticks(ticks=custom_ticks, labels=custom_ticks, fontsize=16, rotation=45)
plt.yscale('log')
plt.ylim(1e-5, 1e2)
plt.xlabel('Energy (MeV)', fontsize=16)
plt.ylabel('Differential Flux (cm$^{-2}$ sr$^{-1}$ s$^{-1}$ MeV$^{-1}$)', fontsize=16)
plt.title('Power Law Spectra for Selected SEP Events', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True, which="major", ls="--", alpha=0.5)
plt.legend(fontsize=14, loc='best')

#%% Calculate channel count rate for Power Law Fits (Solar Protons)
print("--- Generating Plots for Power Law Fits (Solar Protons) ---")
event_select = '18 April 2014'

# Fast vectorized dot products
powerlaw_count_rates = powerlaw_spectra[event_select] @ loaded_geo_values
powerlaw_count_rates_pen = powerlaw_spectra[event_select] @ loaded_geo_pen_values
powerlaw_pen_reversed = powerlaw_count_rates_pen[::-1]

# Plot Count Rates
plot_channel_bar_chart(
    powerlaw_count_rates, powerlaw_pen_reversed,
    title=f'SEP Solar Protons: {event_select}',
    ylabel='Count Rate (counts/s)', ylim=(1e-6, None),
    center_text=f'Event Class: {events_powerlaw[event_select]["s_class"]}',
    colors_range=colors_range, colors_pen_rev=colors_pen_reversed,
    range_x=range_x_axis, pen_x=pen_x_axis, tick_locs=tick_locs, tick_labels=tick_labels,
    num_range=num_channels_range
)

# Plot Overlapping Bar Chart (10 Minute Duration)
counts_10min_range = powerlaw_count_rates * t_sec_10m
counts_10min_pen_reversed = powerlaw_pen_reversed * t_sec_10m

plot_channel_bar_chart(
    counts_10min_range, counts_10min_pen_reversed,
    title=f'SEP Solar Protons: {event_select}',
    ylabel='Total Counts', ylim=(1e0, 1e7),
    center_text=f'Event Class: {events_powerlaw[event_select]["s_class"]}\nCadence: 10 min',
    colors_range=colors_range, colors_pen_rev=colors_pen_reversed,
    range_x=range_x_axis, pen_x=pen_x_axis, tick_locs=tick_locs, tick_labels=tick_labels,
    num_range=num_channels_range
)
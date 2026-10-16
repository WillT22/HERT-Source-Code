#%% Import Libraries 
import os
import sys
import importlib
import pandas as pd
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
import time
import joblib
import types
import pickle

current_script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_script_dir)
import spectra_selection_functions, logic_equations, plot_functions
importlib.reload(spectra_selection_functions)
importlib.reload(logic_equations)
importlib.reload(plot_functions)
from spectra_selection_functions import assign_Normalized_Flat_Weights, get_available_counts, generate_flat_counts, generate_exponential_counts, generate_powerlaw_counts, sample_spectrum
from logic_equations import evaluate_rept_logic, evaluate_reptile2_logic
from plot_functions import create_and_save_pairgrid, plot_species_error_distribution

#%% Optional: Load in workspace
load = True
if load == True:    
    # Load the dictionary back into memory
    load_path = r"C:\Users\wzt0020\Box\HERT_Box\Particle Classification\entire_workspace_hadronic.joblib"
    workspace = joblib.load(load_path)
    # Unpack the dictionary directly into your global environment
    globals().update(workspace)
    print(f"Successfully loaded {len(workspace)} variables into the workspace.")

#%% Start Processing
# Process Data
# Start Clock
start_time = time.perf_counter()

# ==========================================
# CONFIGURATION
# ==========================================
os.chdir(r"C:\Users\wzt0020\Box\HERT_Box\Particle Classification")
PROCESS_DATA = False  # Set to True to process raw data, False to load from pre-processed CSVs

electron_file = r"D:\Matlab Main\Results\Electron_FS\Aggregate Data\Aggregate_Electron_FS_Data_reduced.txt"
proton_file   = r"D:\Matlab Main\Results\Proton_FS_Hadronic\Aggregate Data\Aggregate_Proton_FS_Hadronic_Data_reduced.txt"

column_names = ["E_Inc", "Detector1", "Detector2", "Detector3",
                "Detector4", "Detector5", "Detector6",
                "Detector7", "Detector8", "Detector9"]

# ==========================================
# MAIN EXECUTION
# ==========================================
if PROCESS_DATA:
    # ==========================================
    # 1. Load Data
    # ==========================================
    print("Processing raw data...")

    imported_electron_data = pd.read_csv(electron_file, delim_whitespace=True, skiprows=1, names=column_names)
    imported_proton_data   = pd.read_csv(proton_file, delim_whitespace=True, skiprows=1, names=column_names)

    imported_electron_data['Particle_Type'] = "Electron"
    imported_proton_data['Particle_Type']   = "Proton"
    
    # ==========================================
    # 2. Filter & Clean
    # ==========================================
    # Apply baseline threshold
    electron_data_filtered = imported_electron_data[imported_electron_data['Detector1'] >= 0.1].copy()
    proton_data_filtered   = imported_proton_data[imported_proton_data['Detector1'] >= 0.1].copy()

    print(f"Total Electrons after D1 threshold: {len(electron_data_filtered):,}")
    print(f"Total Protons after D1 threshold:   {len(proton_data_filtered):,}")
    print("-" * 40) # Adds a clean visual divider in the console

    del imported_electron_data
    del imported_proton_data

    # Sum Detectors 7 and 8
    electron_data_filtered['Detector7_8_sum'] = electron_data_filtered['Detector7'] + electron_data_filtered['Detector8']
    proton_data_filtered['Detector7_8_sum']   = proton_data_filtered['Detector7'] + proton_data_filtered['Detector8']

    # Zero out low noise values (< 0.1 MeV) across specific detectors
    cols_to_modify = [col for col in electron_data_filtered.columns if col not in ["E_Inc", "Detector7", "Detector8", "Particle_Type"]]
    
    electron_data_filtered[cols_to_modify] = electron_data_filtered[cols_to_modify].where(electron_data_filtered[cols_to_modify] >= 0.1, 0)
    proton_data_filtered[cols_to_modify]   = proton_data_filtered[cols_to_modify].where(proton_data_filtered[cols_to_modify] >= 0.1, 0)

    # ==========================================
    # 3. Stratified Random Split
    # ==========================================
    num_training_rows = 80000
    num_validation_rows = 10000
    num_test_rows = 10000
    num_points = num_training_rows + num_validation_rows + num_test_rows

    print("\nExtracting random samples for Train, Validation, and Test sets...")

    # --- ELECTRONS ---
    e_train = electron_data_filtered.sample(n=num_training_rows, random_state=42)
    e_remaining = electron_data_filtered[~electron_data_filtered.index.isin(e_train.index)]
    e_val = e_remaining.sample(n=num_validation_rows, random_state=44)
    e_remaining = e_remaining[~e_remaining.index.isin(e_val.index)]
    e_test = e_remaining.sample(n=num_test_rows, random_state=46)

    # --- PROTONS ---
    p_train = proton_data_filtered.sample(n=num_training_rows, random_state=43)
    p_remaining = proton_data_filtered[~proton_data_filtered.index.isin(p_train.index)]
    p_val = p_remaining.sample(n=num_validation_rows, random_state=45)
    p_remaining = p_remaining[~p_remaining.index.isin(p_val.index)]
    p_test = p_remaining.sample(n=num_test_rows, random_state=47)

    # --- COMBINE AND SHUFFLE ---
    # Concatenate the classes and shuffle (frac=1) so the SVM doesn't read all electrons before protons
    training_data   = pd.concat([e_train, p_train], ignore_index=True).sample(frac=1, random_state=48).reset_index(drop=True)
    validation_data = pd.concat([e_val, p_val], ignore_index=True).sample(frac=1, random_state=49).reset_index(drop=True)
    test_data       = pd.concat([e_test, p_test], ignore_index=True).sample(frac=1, random_state=50).reset_index(drop=True)

    print(f"Combined Training Set:   {len(training_data):,} rows")
    print(f"Combined Validation Set: {len(validation_data):,} rows")
    print(f"Combined Test Set:       {len(test_data):,} rows")

    # ==========================================
    # 4. Generate Flat Spectrum Weights
    # ==========================================
    e_bins = np.logspace(np.log10(0.1), np.log10(10), 301)
    p_bins = np.logspace(np.log10(10), np.log10(2000), 401)

    training_data['Flat_Weight'] = assign_Normalized_Flat_Weights(training_data, e_bins, p_bins)
    validation_data['Flat_Weight'] = assign_Normalized_Flat_Weights(validation_data, e_bins, p_bins)
    w_train = training_data['Flat_Weight'].values

    # ==========================================
    # 5. Create Test Sets
    # ==========================================
    # Check what Geant4 actually produced in the remainder sets
    e_avail = get_available_counts(e_remaining, e_bins)
    p_avail = get_available_counts(p_remaining, p_bins)

    print("Extracting Physically Flat Test Set...")
    e_flat_targets = generate_flat_counts(e_avail, num_test_rows)
    p_flat_targets = generate_flat_counts(p_avail, num_test_rows)

    e_flat_test = sample_spectrum(e_remaining, e_flat_targets, e_bins, random_state=42)
    p_flat_test = sample_spectrum(p_remaining, p_flat_targets, p_bins, random_state=42)
    flat_test_data = pd.concat([e_flat_test, p_flat_test], ignore_index=True).sample(frac=1, random_state=50).reset_index(drop=True)

    print("Extracting Scaled Space Spectrum Test Set...")
    e_scaled_targets = generate_exponential_counts(e_bins, e_avail, num_test_rows, E0=2)
    p_scaled_targets = generate_powerlaw_counts(p_bins, p_avail, num_test_rows, gamma=2)

    e_scaled_test = sample_spectrum(e_remaining, e_scaled_targets, e_bins, random_state=42)
    p_scaled_test = sample_spectrum(p_remaining, p_scaled_targets, p_bins, random_state=42)
    scaled_test_data = pd.concat([e_scaled_test, p_scaled_test], ignore_index=True).sample(frac=1, random_state=51).reset_index(drop=True)

    del e_remaining, p_remaining

    # Plot the incident energy distributions for both particles in the flat and scaled test set
    e_ticks = [0.1, 0.2, 0.5, 1, 2, 5, 10]
    p_ticks = [10, 20, 50, 100, 200, 500, 1000, 2000]

    # Flat spectra
    plt.figure(figsize=(14, 5))

    # Left Panel: Electrons
    ax1 = plt.subplot(1, 2, 1)
    plt.hist(
        flat_test_data.loc[flat_test_data['Particle_Type'] == 'Electron', 'E_Inc'], 
        bins=e_bins, color='red', alpha=0.7, label='Electrons'
    )
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(e_bins[0], e_bins[-1])
    plt.xticks(e_ticks, labels=[str(v) for v in e_ticks])
    ax1.xaxis.set_minor_formatter(plt.NullFormatter()) # Prevents duplicate minor labels at 2, 5, etc.
    plt.title('Flat Test Set: Electron Incident Energy')
    plt.xlabel('Incident Energy (MeV)')
    plt.ylabel('Counts')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)

    # Right Panel: Protons
    ax2 = plt.subplot(1, 2, 2)
    plt.hist(
        flat_test_data.loc[flat_test_data['Particle_Type'] == 'Proton', 'E_Inc'], 
        bins=p_bins, color='blue', alpha=0.7, label='Protons'
    )
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(p_bins[0], p_bins[-1])
    plt.xticks(p_ticks, labels=[str(v) for v in p_ticks])
    ax2.xaxis.set_minor_formatter(plt.NullFormatter()) # Prevents duplicate minor labels
    plt.title('Flat Test Set: Proton Incident Energy')
    plt.xlabel('Incident Energy (MeV)')
    plt.ylabel('Counts')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.show()

    # Scaled spectra
    plt.figure(figsize=(14, 5))

    # Left Panel: Electrons
    ax1 = plt.subplot(1, 2, 1)
    plt.hist(
        scaled_test_data.loc[scaled_test_data['Particle_Type'] == 'Electron', 'E_Inc'], 
        bins=e_bins, color='red', alpha=0.7, label='Electrons'
    )
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(e_bins[0], e_bins[-1])
    plt.xticks(e_ticks, labels=[str(v) for v in e_ticks])
    ax1.xaxis.set_minor_formatter(plt.NullFormatter()) # Prevents duplicate minor labels at 2, 5, etc.
    plt.title('Scaled Test Set: Electron Incident Energy')
    plt.xlabel('Incident Energy (MeV)')
    plt.ylabel('Counts')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)

    # Right Panel: Protons
    ax2 = plt.subplot(1, 2, 2)
    plt.hist(
        scaled_test_data.loc[scaled_test_data['Particle_Type'] == 'Proton', 'E_Inc'], 
        bins=p_bins, color='blue', alpha=0.7, label='Protons'
    )
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(p_bins[0], p_bins[-1])
    plt.xticks(p_ticks, labels=[str(v) for v in p_ticks])
    ax2.xaxis.set_minor_formatter(plt.NullFormatter()) # Prevents duplicate minor labels
    plt.title('Scaled Test Set: Proton Incident Energy')
    plt.xlabel('Incident Energy (MeV)')
    plt.ylabel('Counts')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.show()


    # Save each data frame to a separate text file (.csv)
    training_data.to_csv("training_data_hadronic.csv", index=False)
    validation_data.to_csv("validation_data_hadronic.csv", index=False)
    flat_test_data.to_csv("flat_test_data_hadronic.csv", index=False)
    scaled_test_data.to_csv("scaled_test_data_hadronic.csv", index=False)

    print("Data processed and saved to CSVs.")

else:
    # Load directly from CSVs to save time
    print("Loading data from pre-processed CSVs...")
    training_data = pd.read_csv("training_data_hadronic.csv")
    w_train = training_data['Flat_Weight'].values
    validation_data = pd.read_csv("validation_data_hadronic.csv")
    flat_test_data = pd.read_csv("flat_test_data_hadronic.csv")
    scaled_test_data = pd.read_csv("scaled_test_data_hadronic.csv")
    print("Data loaded.")

# ==========================================
# COMMON POST-PROCESSING
# ==========================================
# Convert Particle_Type to an ordered Categorical data type
cat_type = pd.CategoricalDtype(categories=["Electron", "Proton"], ordered=True)
training_data['Particle_Type'] = training_data['Particle_Type'].astype(cat_type)
validation_data['Particle_Type'] = validation_data['Particle_Type'].astype(cat_type)
flat_test_data['Particle_Type'] = flat_test_data['Particle_Type'].astype(cat_type)
scaled_test_data['Particle_Type'] = scaled_test_data['Particle_Type'].astype(cat_type)

#%% Create SVM Model
# Create SVM
# ==========================================
# TOGGLE SVM EXECUTION MODE
# ==========================================
# "tune"   : Run GridSearch on validation data to find best C values, then create SVMs
# "create" : Train the SVMs directly using the predefined optimized costs (C=10)
# "skip"   : Bypass SVM training entirely
SVM_MODE = "skip"  # Options: "tune", "create", "skip"

if SVM_MODE != "skip":
    print(f"--- Starting SVM Processing (Mode: {SVM_MODE}) ---")
    
    # Define features and target for the main models
    features = ["Detector1", "Detector2", "Detector3", "Detector4", 
                "Detector5", "Detector6", "Detector7_8_sum", "Detector9"]
    
    X_train = training_data[features]
    y_train = training_data['Particle_Type']
    
    X_val = validation_data[features]
    y_val = validation_data['Particle_Type']
    
    X_test_flat = flat_test_data[features]
    y_test_flat = flat_test_data['Particle_Type']

    X_test_scaled = scaled_test_data[features]
    y_test_scaled = scaled_test_data['Particle_Type']

    # --- TUNE OPTION ---
    if SVM_MODE == "tune":
        print("Tuning SVMs on validation data...")
        param_grid = {'C': [0.001, 0.01, 0.1, 1, 5, 10, 100]}
        
        # Linear Tune
        tuner_lin = GridSearchCV(SVC(kernel='linear'), param_grid, cv=3, n_jobs=-1)
        tuner_lin.fit(X_val, y_val)
        print(f"Best Linear Cost: {tuner_lin.best_params_['C']}")
        
        # Poly Tune (adding 1000 to param grid based on R script)
        param_grid_poly = {'C': [0.001, 0.01, 0.1, 1, 5, 10, 100, 1000]}
        tuner_poly = GridSearchCV(SVC(kernel='poly'), param_grid_poly, cv=3, n_jobs=-1)
        tuner_poly.fit(X_val, y_val)
        print(f"Best Poly Cost: {tuner_poly.best_params_['C']}")
        
        # Radial Tune
        tuner_rad = GridSearchCV(SVC(kernel='rbf'), param_grid_poly, cv=3, n_jobs=-1)
        tuner_rad.fit(X_val, y_val)
        print(f"Best Radial Cost: {tuner_rad.best_params_['C']}")

    # --- CREATE OPTION (or proceeding after tuning) ---
    start_linearSVM_time = time.perf_counter()
    print("\n### Linear SVM ###")
    # In scikit-learn, setting class_weight works identically to class.weights in R
    svm_linear = SVC(kernel='linear', C=10, class_weight={'Electron': 1, 'Proton': 1})
    svm_linear.fit(X_train, y_train, sample_weight=w_train)
    
    # Get the support vectors and coefficients of the hyperplane (Intercept is prepended to match R output)
    linear_suppvec = svm_linear.n_support_
    print("Support Vectors (Electrons, Protons):", linear_suppvec)
    linear_hp_coefs = np.insert(svm_linear.coef_[0], 0, svm_linear.intercept_[0])
    print("Linear Hyperplane Coefficients:\n", linear_hp_coefs)
    
    # Calculate weights relative to the intercept (matching linear_hp_coefs[[1]] in R)
    weights = np.abs(linear_hp_coefs / np.sum(linear_hp_coefs[1:]))  # Exclude intercept for normalization
    print("Weights:\n", weights)
    
    # Predict to validate model (Flat test data)
    linear_predictions_flat = svm_linear.predict(X_test_flat)
    # Print predictions (prop.table equivalent)
    linear_ct = pd.crosstab(index=linear_predictions_flat, columns=y_test_flat, rownames=['Predict'], normalize='columns') * 100
    print("\nLinear SVM Predictions Flat Spectra (%):\n", linear_ct)

    # Predict to validate model (Scaled test data)
    linear_predictions_scaled = svm_linear.predict(X_test_scaled)
    # Print predictions (prop.table equivalent)
    linear_ct = pd.crosstab(index=linear_predictions_scaled, columns=y_test_scaled, rownames=['Predict'], normalize='columns') * 100
    print("\nLinear SVM Predictions Scaled Spectra (%):\n", linear_ct)

    end_linearSVM_time = time.perf_counter()
    elapsed = end_linearSVM_time - start_linearSVM_time
    minutes, seconds = divmod(elapsed, 60)
    print(f"Execution time: {int(minutes)} minutes and {int(seconds)} seconds")


    start_polySVM_time = time.perf_counter()
    print("\n### Polynomial SVM ###")
    svm_poly = SVC(kernel='poly', C=10) # default degree is 3, same as R
    svm_poly.fit(X_train, y_train)

    # Get the support vectors
    poly_suppvec = svm_poly.n_support_
    print("Support Vectors (Electrons, Protons):", poly_suppvec)

    # Predict to validate model
    poly_predictions_flat = svm_poly.predict(X_test_flat)
    poly_ct = pd.crosstab(index=poly_predictions_flat, columns=y_test_flat, rownames=['Predict'], normalize='columns') * 100
    print("\nPolynomial SVM Predictions Flat Spectra (%):\n", poly_ct)

    poly_predictions_scaled = svm_poly.predict(X_test_scaled)
    poly_ct = pd.crosstab(index=poly_predictions_scaled, columns=y_test_scaled, rownames=['Predict'], normalize='columns') * 100
    print("\nPolynomial SVM Predictions Scaled Spectra (%):\n", poly_ct)

    end_polySVM_time = time.perf_counter()
    elapsed = end_polySVM_time - start_polySVM_time
    minutes, seconds = divmod(elapsed, 60)
    print(f"Execution time: {int(minutes)} minutes and {int(seconds)} seconds")


    start_radialSVM_time = time.perf_counter()
    print("\n### Radial SVM ###")
    svm_radial = SVC(kernel='rbf', C=10, class_weight={'Electron': 1, 'Proton': 1})
    svm_radial.fit(X_train, y_train)

    # Get the support vectors
    radial_suppvec = svm_radial.n_support_
    print("Support Vectors (Electrons, Protons):", radial_suppvec)

    radial_predictions_flat = svm_radial.predict(X_test_flat)
    radial_ct = pd.crosstab(index=radial_predictions_flat, columns=y_test_flat, rownames=['Predict'], normalize='columns') * 100
    print("\nRadial SVM Predictions Flat Spectra (%):\n", radial_ct)

    radial_predictions_scaled = svm_radial.predict(X_test_scaled)
    radial_ct = pd.crosstab(index=radial_predictions_scaled, columns=y_test_scaled, rownames=['Predict'], normalize='columns') * 100
    print("\nRadial SVM Predictions Scaled Spectra (%):\n", radial_ct)

    end_radialSVM_time = time.perf_counter()
    elapsed = end_radialSVM_time - start_radialSVM_time
    minutes, seconds = divmod(elapsed, 60)
    print(f"Execution time: {int(minutes)} minutes and {int(seconds)} seconds")

else:
    print("--- Skipping SVM Processing ---")

#%% Evaluate Instrument Performance Comparison
# Run legacy logic equations
print("\n--- Evaluating Legacy Logic Equations on Flat Spectra ---")
evaluate_rept_logic(flat_test_data)
evaluate_reptile2_logic(flat_test_data)

print("\n--- Evaluating Legacy Logic Equations on Scaled Spectra ---")
evaluate_rept_logic(scaled_test_data)
evaluate_reptile2_logic(scaled_test_data)

#%% Find Max Energies
# Find the max energies deposited onto each detector for both particle types
# Create a list of just the detector columns (if not already defined)
features_only = ["Detector1", "Detector2", "Detector3", "Detector4", 
                 "Detector5", "Detector6", "Detector7_8_sum", "Detector9"]

# Subset and calculate the max for each particle
e_max_energies = flat_test_data[flat_test_data['Particle_Type'] == 'Electron'][features_only].max()
p_max_energies = flat_test_data[flat_test_data['Particle_Type'] == 'Proton'][features_only].max()

# Combine them into a single summary dataframe for a clean printout
max_energy_summary = pd.DataFrame({
    'Electron Max (MeV)': e_max_energies,
    'Proton Max (MeV)': p_max_energies
})

print("\n--- Maximum Deposited Energy by Detector ---")
print(max_energy_summary.to_string())

# Define the output file name
output_filename = "max_energies_summary.txt"

# Export to text file
with open(output_filename, "w") as file:
    file.write("--- Maximum Deposited Energy by Detector ---\n")
    file.write(max_energy_summary.to_string())
    file.write("\n")

print(f"\nResults successfully exported to: {os.path.abspath(output_filename)}")

#%% Create and Save PairGrid with Density Plots
print("\n--- Initiating Visualization Pipeline ---")
# Define exactly what columns go into the PairGrid matrix
columns_to_keep = ["Detector1", "Detector2", "Detector3", "Detector4", 
                   "Detector5", "Detector6", "Detector7_8_sum", "Detector9", "Particle_Type"]
flat_test_data_plot = flat_test_data[columns_to_keep].copy()

# Features for tracking weights/indices
features_only = ["Detector1", "Detector2", "Detector3", "Detector4", 
                 "Detector5", "Detector6", "Detector7_8_sum", "Detector9"]

create_and_save_pairgrid(flat_test_data_plot, features_only, 
                         output="square_pairs_plot_flatspectra.png")

features_only_reduced = ["Detector1", "Detector2", "Detector7_8_sum", "Detector9"]
create_and_save_pairgrid(flat_test_data_plot, features_only_reduced,
                         output="square_pairs_plot_reduced_flatspectra.png")  


# Re-run with Scaled Test Data
scaled_test_data_plot = scaled_test_data[columns_to_keep].copy()

# Features for tracking weights/indices
features_only = ["Detector1", "Detector2", "Detector3", "Detector4", 
                 "Detector5", "Detector6", "Detector7_8_sum", "Detector9"]

create_and_save_pairgrid(scaled_test_data_plot, features_only, 
                         output="square_pairs_plot_spectrascale.png")

features_only_reduced = ["Detector1", "Detector2", "Detector7_8_sum", "Detector9"]
create_and_save_pairgrid(scaled_test_data_plot, features_only_reduced,
                         output="square_pairs_plot_reduced_spectrascale.png")   

plt.rcdefaults()               

#%% Save workspace
save = False  # Set to True to save the workspace, False to skip
if save:
    workspace = {}
    
    # Loop through all variables currently in memory
    # Force a static list copy of the globals dictionary before iterating
    for name, value in list(globals().items()):
        
        # Skip Jupyter/IPython hidden variables (which start with '_') and modules/functions
        if not name.startswith('_') and not isinstance(value, (types.ModuleType, types.FunctionType, type)):
            
            # Test if the object is safe to save
            try:
                pickle.dumps(value)
                workspace[name] = value
            except Exception:
                # Skip objects that cannot be saved (like open file handles or active plots)
                continue

    # Actually save the compiled dictionary
    save_path = os.path.join(os.getcwd(), 'entire_workspace_hadronic.joblib')
    print(f"Compressing {len(workspace)} variables...")
    joblib.dump(workspace, save_path)
    print(f"Workspace successfully saved to: {save_path}")
else:
    print("Workspace saving skipped.")

end_time = time.perf_counter()
elapsed = end_time - start_time
minutes, seconds = divmod(elapsed, 60)
print(f"Execution time: {int(minutes)} minutes and {int(seconds)} seconds")


#%% Exploratory: Find the data that are inaccurately classified by the Linear SVM model
# Extract missclassified particle data
misclassified_mask_flat = linear_predictions_flat != y_test_flat
misclassified_data_flat = X_test_flat[misclassified_mask_flat].copy()
misclassified_data_flat.insert(0, 'E_Inc', flat_test_data.loc[misclassified_mask_flat, 'E_Inc'])
misclassified_data_flat['Particle_Type'] = y_test_flat[misclassified_mask_flat]
misclassified_data_flat['Predicted_Label'] = linear_predictions_flat[misclassified_mask_flat]

features_mc_plot = ["Detector1", "Detector2", "Detector3", "Detector4", 
                    "Detector5", "Detector6", "Detector7_8_sum", "Detector9"]

# ---------------------------------------------------------
# Missclassification Pair Plot
# ---------------------------------------------------------
create_and_save_pairgrid(
    test_data_plot=misclassified_data_flat, 
    features_only=features_mc_plot, 
    output="misclassified_particles_8x8_grid_flatspectra.png"
)

# Run for both species in two lines:
plot_species_error_distribution(misclassified_data_flat, 'Electron', e_bins, (0.5, 10), [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10], features_mc_plot, 'red', 'darkred', '(Flat Spectra)')
plot_species_error_distribution(misclassified_data_flat, 'Proton', p_bins, (10, 2000), [10, 20, 50, 100, 200, 500, 1000, 2000], features_mc_plot, 'blue', 'darkblue', '(Flat Spectra)')

plt.rcdefaults()

#%% Redo misclassification analysis for Scaled Spectra
# Extract missclassified particle data
misclassified_mask_scaled = linear_predictions_scaled != y_test_scaled
misclassified_data_scaled = X_test_scaled[misclassified_mask_scaled].copy()
misclassified_data_scaled.insert(0, 'E_Inc', scaled_test_data.loc[misclassified_mask_scaled, 'E_Inc'])
misclassified_data_scaled['Particle_Type'] = y_test_scaled[misclassified_mask_scaled]
misclassified_data_scaled['Predicted_Label'] = linear_predictions_scaled[misclassified_mask_scaled]

features_mc_plot = ["Detector1", "Detector2", "Detector3", "Detector4", 
                    "Detector5", "Detector6", "Detector7_8_sum", "Detector9"]

# ---------------------------------------------------------
# Missclassification Pair Plot
# ---------------------------------------------------------
create_and_save_pairgrid(
    test_data_plot=misclassified_data_scaled, 
    features_only=features_mc_plot, 
    output="misclassified_particles_8x8_grid_scaledspectra.png"
)

# Run for both species in two lines:
plot_species_error_distribution(misclassified_data_scaled, 'Electron', e_bins, (0.5, 10), [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10], features_mc_plot, 'red', 'darkred', '(Scaled Spectra)')
plot_species_error_distribution(misclassified_data_scaled, 'Proton', p_bins, (10, 2000), [10, 20, 50, 100, 200, 500, 1000, 2000], features_mc_plot, 'blue', 'darkblue', '(Scaled Spectra)')

plt.rcdefaults()
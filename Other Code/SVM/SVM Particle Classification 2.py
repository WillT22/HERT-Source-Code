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

current_script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_script_dir)
import logic_equations, plot_functions
importlib.reload(logic_equations)
importlib.reload(plot_functions)
from logic_equations import evaluate_rept_logic, evaluate_reptile2_logic
from plot_functions import check_densities, create_and_save_pairgrid, create_and_save_pairgrid_reduced

#%% Optional: Load in workspace
load = False
if load == True:    
    import joblib
    # Load the dictionary back into memory
    load_path = 'entire_workspace.joblib'
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
PROCESS_DATA = False

electron_file = r"D:\Matlab Main\Results\Electron_FS\Aggregate Data\Aggregate_Electron_FS_Data_new.txt"
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
    num_training_rows = 8000
    num_validation_rows = 1000
    num_test_rows = 1000
    num_points = num_training_rows + num_validation_rows + num_test_rows

    print("\nExtracting random samples for Train, Validation, and Test sets...")

    # Safety Check
    for name, df in [("Electrons", electron_data_filtered), ("Protons", proton_data_filtered)]:
        if len(df) < num_points:
            raise ValueError(f"Not enough {name}! You need {num_points} total, but only have {len(df)}.")

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


    # Optional: Free up RAM
    del electron_data_filtered
    del proton_data_filtered
    del e_remaining
    del p_remaining

    # Define energy bins for histogramming and flat spectrum sampling
    e_bins = np.logspace(np.log10(0.1), np.log10(10), 301)
    p_bins = np.logspace(np.log10(10), np.log10(2000), 401)
    
    e_train_counts, _ = np.histogram(training_data[training_data['Particle_Type'] == 'Electron']['E_Inc'], bins=e_bins)
    e_val_counts, _   = np.histogram(validation_data[validation_data['Particle_Type'] == 'Electron']['E_Inc'], bins=e_bins)
    e_test_counts, _  = np.histogram(test_data[test_data['Particle_Type'] == 'Electron']['E_Inc'], bins=e_bins)

    p_train_counts, _ = np.histogram(training_data[training_data['Particle_Type'] == 'Proton']['E_Inc'], bins=p_bins)
    p_val_counts, _   = np.histogram(validation_data[validation_data['Particle_Type'] == 'Proton']['E_Inc'], bins=p_bins)
    p_test_counts, _  = np.histogram(test_data[test_data['Particle_Type'] == 'Proton']['E_Inc'], bins=p_bins)

    def assign_flat_weights(df, particle_type, bins, actual_counts):
        # Locate the particles
        mask = df['Particle_Type'] == particle_type
        
        # Map each particle's incident energy to its respective bin index
        bin_indices = np.digitize(df.loc[mask, 'E_Inc'], bins) - 1
        
        # Define the mathematically flat target: total particles spread perfectly evenly across all active bins
        active_bins = (actual_counts > 0).sum()
        target_per_bin = actual_counts.sum() / active_bins if active_bins > 0 else 0
        
        # Create the lookup array: Weight = Target Environment / Actual Simulated Counts
        weight_lookup = np.zeros(len(actual_counts))
        for i, count in enumerate(actual_counts):
            if count > 0:
                weight_lookup[i] = target_per_bin / count
                
        # Return the mapped weights for the specific rows
        return weight_lookup[bin_indices]

    # Apply the weights to the Training Data
    training_data['Flat_Weight'] = 0.0

    e_train_mask = training_data['Particle_Type'] == 'Electron'
    p_train_mask = training_data['Particle_Type'] == 'Proton'

    training_data.loc[e_train_mask, 'Flat_Weight'] = assign_flat_weights(training_data, 'Electron', e_bins, e_train_counts)
    training_data.loc[p_train_mask, 'Flat_Weight'] = assign_flat_weights(training_data, 'Proton', p_bins, p_train_counts)

    # Apply the weights to the Validation Data
    validation_data['Flat_Weight'] = 0.0

    e_val_mask = validation_data['Particle_Type'] == 'Electron'
    p_val_mask = validation_data['Particle_Type'] == 'Proton'

    validation_data.loc[e_val_mask, 'Flat_Weight'] = assign_flat_weights(validation_data, 'Electron', e_bins, e_val_counts)
    validation_data.loc[p_val_mask, 'Flat_Weight'] = assign_flat_weights(validation_data, 'Proton', p_bins, p_val_counts)

    # Class-Level Normalization (Critical for SVM Balance)
    # This forces the sum of all electron weights to equal the sum of all proton weights (1.0), 
    # preventing the SVM from simply guessing the class that happens to have a higher total weight sum.
    def normalize_50_50(df, weight_col):
        e_sum = df.loc[df['Particle_Type'] == 'Electron', weight_col].sum()
        p_sum = df.loc[df['Particle_Type'] == 'Proton', weight_col].sum()
        
        if e_sum > 0:
            df.loc[df['Particle_Type'] == 'Electron', weight_col] /= e_sum
        if p_sum > 0:
            df.loc[df['Particle_Type'] == 'Proton', weight_col] /= p_sum
            
        return df

    training_data = normalize_50_50(training_data, 'Flat_Weight')
    validation_data = normalize_50_50(validation_data, 'Flat_Weight')

    # 5. Extract the array for the SVM
    w_train = training_data['Flat_Weight'].values

    # Save each data frame to a separate text file (.csv)
    training_data.to_csv("training_data_hadronic_red.csv", index=False)
    validation_data.to_csv("validation_data_hadronic_red.csv", index=False)
    test_data.to_csv("test_data_hadronic_red.csv", index=False)

    print("Data processed and saved to CSVs.")

else:
    # Load directly from CSVs to save time
    print("Loading data from pre-processed CSVs...")
    training_data = pd.read_csv("training_data_hadronic_red.csv")
    w_train = training_data['Flat_Weight'].values
    validation_data = pd.read_csv("validation_data_hadronic_red.csv")
    test_data = pd.read_csv("test_data_hadronic_red.csv")
    print("Data loaded.")

# ==========================================
# COMMON POST-PROCESSING
# ==========================================
# Convert Particle_Type to an ordered Categorical data type
cat_type = pd.CategoricalDtype(categories=["Electron", "Proton"], ordered=True)
training_data['Particle_Type'] = training_data['Particle_Type'].astype(cat_type)
validation_data['Particle_Type'] = validation_data['Particle_Type'].astype(cat_type)
test_data['Particle_Type'] = test_data['Particle_Type'].astype(cat_type)

#%% Create SVM Model
# Create SVM
# ==========================================
# TOGGLE SVM EXECUTION MODE
# ==========================================
# "tune"   : Run GridSearch on validation data to find best C values, then create SVMs
# "create" : Train the SVMs directly using the predefined optimized costs (C=10)
# "skip"   : Bypass SVM training entirely
SVM_MODE = "create"  # Options: "tune", "create", "skip"

if SVM_MODE != "skip":
    print(f"--- Starting SVM Processing (Mode: {SVM_MODE}) ---")
    
    # Define features and target for the main models
    features = ["Detector1", "Detector2", "Detector3", "Detector4", 
                "Detector5", "Detector6", "Detector7_8_sum", "Detector9"]
    
    X_train = training_data[features]
    y_train = training_data['Particle_Type']
    
    X_val = validation_data[features]
    y_val = validation_data['Particle_Type']
    
    X_test = test_data[features]
    y_test = test_data['Particle_Type']

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
    
    # Predict to validate model
    linear_predictions = svm_linear.predict(X_test)
    
    # Print predictions (prop.table equivalent)
    linear_ct = pd.crosstab(index=linear_predictions, columns=y_test, rownames=['Predict'], normalize='columns') * 100
    print("\nLinear SVM Predictions (%):\n", linear_ct)
    
    # Using sklearn's built-in decision function is safer than manually multiplying,
    # but achieves the exact same math as your manual test block. 
    # > 0 assigns it to the positive class (Protons, based on sklearn's alphabetical class sorting).
    manual_decisions = svm_linear.decision_function(X_test)
    # sklearn sorts classes alphabetically: 0='Electron', 1='Proton'. So < 0 is Electron.
    linear_hp_test = np.where(manual_decisions < 0, "Electron", "Proton")
    
    hp_ct = pd.crosstab(index=linear_hp_test, columns=test_data['Particle_Type'], rownames=['Predict'], normalize='columns') * 100
    print("\nManual Hyperplane Test (%):\n", hp_ct)

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
    poly_predictions = svm_poly.predict(X_test)
    poly_ct = pd.crosstab(index=poly_predictions, columns=y_test, rownames=['Predict'], normalize='columns') * 100
    print("Polynomial SVM Predictions (%):\n", poly_ct)

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

    radial_predictions = svm_radial.predict(X_test)
    radial_ct = pd.crosstab(index=radial_predictions, columns=y_test, rownames=['Predict'], normalize='columns') * 100
    print("Radial SVM Predictions (%):\n", radial_ct)

    end_radialSVM_time = time.perf_counter()
    elapsed = end_radialSVM_time - start_radialSVM_time
    minutes, seconds = divmod(elapsed, 60)
    print(f"Execution time: {int(minutes)} minutes and {int(seconds)} seconds")

    print("\n### Simplified Linear SVM ###")
    features_si = ["Detector1", "Detector2", "Detector9"]
    X_train_si = training_data[features_si]
    y_train_si = training_data['Particle_Type']
    
    svm_linearsi = SVC(kernel='linear', C=10)
    svm_linearsi.fit(X_train_si, y_train_si)
    
    linearsi_hp_coefs = np.insert(svm_linearsi.coef_[0], 0, svm_linearsi.intercept_[0])
    print("Simplified Linear Coefficients:\n", linearsi_hp_coefs)
    
    linearsi_test = svm_linearsi.predict(X_test[features_si])
    linearsi_ct = pd.crosstab(index=linearsi_test, columns=y_test, rownames=['Predict'], normalize='columns') * 100
    print("\nSimplified Linear Predictions (%):\n", linearsi_ct)

    # Manual test for simplified hyperplane
    si_manual_decisions = svm_linearsi.decision_function(test_data[features_si])
    linearsi_hp_test = np.where(si_manual_decisions < 0, "Electron", "Proton")
    
    hp_si_ct = pd.crosstab(index=linearsi_hp_test, columns=test_data['Particle_Type'], rownames=['Predict'], normalize='columns') * 100
    print("\nManual Simplified Hyperplane Test (%):\n", hp_si_ct)

else:
    print("--- Skipping SVM Processing ---")

#%% Evaluate Instrument Performance Comparison
# Run legacy logic equations
evaluate_rept_logic(test_data)
evaluate_reptile2_logic(test_data)


#%% Find Max Energies
# Find the max energies deposited onto each detector for both particle types
# Create a list of just the detector columns (if not already defined)
features_only = ["Detector1", "Detector2", "Detector3", "Detector4", 
                 "Detector5", "Detector6", "Detector7_8_sum", "Detector9"]

# Subset and calculate the max for each particle
e_max_energies = test_data[test_data['Particle_Type'] == 'Electron'][features_only].max()
p_max_energies = test_data[test_data['Particle_Type'] == 'Proton'][features_only].max()

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
current_script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_script_dir)
import logic_equations, plot_functions
importlib.reload(logic_equations)
importlib.reload(plot_functions)
from logic_equations import evaluate_rept_logic, evaluate_reptile2_logic
from plot_functions import check_densities, create_and_save_pairgrid, create_and_save_pairgrid_reduced

print("\n--- Initiating Visualization Pipeline ---")
# Define exactly what columns go into the PairGrid matrix
columns_to_keep = ["Detector1", "Detector2", "Detector3", "Detector4", 
                   "Detector5", "Detector6", "Detector7_8_sum", "Detector9", "Particle_Type"]
test_data_plot = test_data[columns_to_keep].copy()

# Features for tracking weights/indices
features_only = ["Detector1", "Detector2", "Detector3", "Detector4", 
                 "Detector5", "Detector6", "Detector7_8_sum", "Detector9"]
create_and_save_pairgrid(test_data_plot, features_only, 
                         output="square_pairs_plot_with_horizontal_legends.png")

features_only_reduced = ["Detector1", "Detector2", "Detector3", "Detector9"]
create_and_save_pairgrid_reduced(test_data_plot, features_only_reduced,
                         output="square_pairs_plot_reduced.png")                    

#%% Save workspace
import joblib
import types
import pickle
import os

save = True  # Set to True to save the workspace, False to skip
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
    save_path = os.path.join(os.getcwd(), 'entire_workspace_spectrascale.joblib')
    print(f"Compressing {len(workspace)} variables...")
    joblib.dump(workspace, save_path)
    print(f"Workspace successfully saved to: {save_path}")
else:
    print("Workspace saving skipped.")

end_time = time.perf_counter()
elapsed = end_time - start_time
minutes, seconds = divmod(elapsed, 60)
print(f"Execution time: {int(minutes)} minutes and {int(seconds)} seconds")

#%% Load workspace
# import joblib

# # Load the dictionary back into memory
# load_path = 'entire_workspace.joblib'
# workspace = joblib.load(load_path)

# # Unpack the dictionary directly into your global environment
# globals().update(workspace)

# print(f"Successfully loaded {len(workspace)} variables into the workspace.")


#%% Exploratory: Find the data that are inaccurately classified by the Linear SVM model
# Extract missclassified particle data
misclassified_mask = linear_predictions != y_test
misclassified_data = X_test[misclassified_mask].copy()
misclassified_data.insert(0, 'E_Inc', test_data.loc[misclassified_mask, 'E_Inc'])
misclassified_data['Particle_Type'] = y_test[misclassified_mask]
misclassified_data['Predicted_Label'] = linear_predictions[misclassified_mask]

features_mc_plot = ["Detector1", "Detector2", "Detector3", "Detector4", 
                    "Detector5", "Detector6", "Detector7_8_sum", "Detector9"]

# ---------------------------------------------------------
# Missclassification Pair Plot
# ---------------------------------------------------------
create_and_save_pairgrid(
    test_data=misclassified_data, 
    features_only=features, 
    output="misclassified_particles_8x8_grid.png"
)

# =============================================================================
# 1. ELECTRONS (Misclassified as Protons)
# =============================================================================
electrons_as_protons = misclassified_data[misclassified_data['Particle_Type'] == 'Electron']
total_electrons = test_data[test_data['Particle_Type'] == 'Electron']

max_edep_e = electrons_as_protons[features_mc_plot].max().max()
d_bins_e = np.linspace(0, max_edep_e, 101)

e_totals_df = pd.read_csv('D:/Matlab Main/Results/Electron_FS/PostProcess_Histograms/Electron_Incident_Energy_Histogram.txt', sep='\t')
total_per_MeV_e = e_totals_df['Counts_per_MeV'].values
e_ticks = [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10]

fig_e, axes_e = plt.subplots(3, 3, figsize=(20, 18))
fig_e.suptitle('Distribution of Misclassified Electrons', fontsize=40, fontweight='bold', y=0.97)

for i, feature in enumerate(features_mc_plot + ['E_Inc']):
    ax = axes_e.flat[i]
    
    if feature != 'E_Inc':
        mc_counts, _ = np.histogram(electrons_as_protons[feature], bins=d_bins_e)
        prob = (mc_counts / np.sum(mc_counts)) * 100
        
        ax.hist(d_bins_e[:-1], bins=d_bins_e, weights=prob, color='red', alpha=0.7)
        ax.set_xscale('linear')
        ax.set_xlim(0, max_edep_e)
        ax.set_xlabel('Deposited Energy (MeV)', fontsize=24)
        
    else:
        mc_counts, _ = np.histogram(electrons_as_protons['E_Inc'], bins=e_bins)
        tot_counts_e = total_per_MeV_e * np.diff(e_bins)
        
        corrected_mc = np.divide(mc_counts, tot_counts_e, out=np.zeros_like(mc_counts, dtype=float), where=tot_counts_e!=0)
        prob = (corrected_mc / np.sum(corrected_mc)) * 100
        
        ax.hist(e_bins[:-1], bins=e_bins, weights=prob, color='darkred', alpha=0.9)
        ax.set_xscale('log')
        ax.set_xlim(0.5, 10)
        ax.set_xticks(e_ticks, labels=[str(v) for v in e_ticks])
        ax.set_xlabel('Incident Energy (MeV)', fontsize=24)
        
    ax.set_title(feature, fontsize=28, fontweight='bold')
    ax.set_ylabel('% of Total Errors', fontsize=24)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.grid(True, which='major', linestyle='--', alpha=0.7)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

# =============================================================================
# 2. PROTONS (Misclassified as Electrons)
# =============================================================================
protons_as_electrons = misclassified_data[misclassified_data['Particle_Type'] == 'Proton']
total_protons = test_data[test_data['Particle_Type'] == 'Proton']

max_edep_p = protons_as_electrons[features_mc_plot].max().max()
d_bins_p = np.linspace(0, max_edep_p, 101)

p_totals_df = pd.read_csv('D:/Matlab Main/Results/Proton_FS_Hadronic/PostProcess_Results/Proton_Incident_Energy_Histogram.txt', sep='\t')
total_per_MeV_p = p_totals_df['Counts_per_MeV'].values
p_ticks = [10, 100, 1000, 2000]

fig_p, axes_p = plt.subplots(3, 3, figsize=(20, 18))
fig_p.suptitle('Distribution of Misclassified Protons', fontsize=40, fontweight='bold', y=0.97)

for i, feature in enumerate(features_mc_plot + ['E_Inc']):
    ax = axes_p.flat[i]
    
    if feature != 'E_Inc':
        mc_counts, _ = np.histogram(protons_as_electrons[feature], bins=d_bins_p)
        prob = (mc_counts / np.sum(mc_counts)) * 100
        
        ax.hist(d_bins_p[:-1], bins=d_bins_p, weights=prob, color='blue', alpha=0.7)
        ax.set_xscale('linear')
        ax.set_xlim(0, max_edep_p)
        ax.set_xlabel('Deposited Energy (MeV)', fontsize=24)
        
    else:
        mc_counts, _ = np.histogram(protons_as_electrons['E_Inc'], bins=p_bins)
        tot_counts_p = total_per_MeV_p * np.diff(p_bins)
        
        corrected_mc = np.divide(mc_counts, tot_counts_p, out=np.zeros_like(mc_counts, dtype=float), where=tot_counts_p!=0)
        prob = (corrected_mc / np.sum(corrected_mc)) * 100
        
        ax.hist(p_bins[:-1], bins=p_bins, weights=prob, color='darkblue', alpha=0.9)
        ax.set_xscale('log')
        ax.set_xlim(10, 2000)
        ax.set_xticks(p_ticks, labels=[str(v) for v in p_ticks])
        ax.set_xlabel('Incident Energy (MeV)', fontsize=24)
        
    ax.set_title(feature, fontsize=28, fontweight='bold')
    ax.set_ylabel('% of Total Errors', fontsize=24)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.grid(True, which='both', linestyle='--', alpha=0.7)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()
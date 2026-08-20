#%% Import Libraries 
import os
import sys
import importlib
import pandas as pd
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

current_script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_script_dir)
import logic_equations, plot_functions
importlib.reload(logic_equations)
importlib.reload(plot_functions)
from logic_equations import evaluate_rept_logic, evaluate_reptile2_logic
from plot_functions import check_densities, create_and_save_pairgrid, get_bin_edges

#%% Optional: Load in workspace
import joblib
# Load the dictionary back into memory
load_path = 'entire_workspace.joblib'
workspace = joblib.load(load_path)
# Unpack the dictionary directly into your global environment
globals().update(workspace)
print(f"Successfully loaded {len(workspace)} variables into the workspace.")

#%% Load Data
# Set working directory
os.chdir(r"C:\Users\wzt0020\Box\HERT_Box\Particle Classification")

# ==========================================
# TOGGLE DATA PIPELINE MODE
# ==========================================
# True: Read raw .txt files, filter, sample, and save to .csv
# False: Skip processing and load directly from saved .csv files
PROCESS_DATA = False

# Set the number of rows you want for each data set
num_training_rows = 400000
num_validation_rows = 50000
num_test_rows = 50000
num_points = (num_training_rows + num_validation_rows + num_test_rows) * 2

if PROCESS_DATA:
    print("Processing raw data...")
    # Set the file paths
    electron_file = r"D:\HERT_Drive\Matlab Main\Result\Electron_FS\Aggregate Data\Aggregate_Electron_FS_Data_new.txt"
    proton_file   = r"D:\HERT_Drive\Matlab Main\Result\Proton_FS\Aggregate Data\Aggregate_Proton_FS_Data_new.txt"

    # Define column names
    column_names = ["E_Inc", "Detector1", "Detector2", "Detector3",
                    "Detector4", "Detector5", "Detector6",
                    "Detector7", "Detector8", "Detector9"]

    # Read the data
    imported_electron_data = pd.read_csv(electron_file, delim_whitespace=True, skiprows=1, names=column_names)
    imported_proton_data   = pd.read_csv(proton_file, delim_whitespace=True, skiprows=1, names=column_names)

    # Create a new column named "Particle_Type"
    imported_electron_data['Particle_Type'] = "Electron"
    imported_proton_data['Particle_Type'] = "Proton"

    # Apply the condition to the full data sets first
    electron_data_filtered = imported_electron_data[imported_electron_data['Detector1'] >= 0.1].copy()
    proton_data_filtered   = imported_proton_data[imported_proton_data['Detector1'] >= 0.1].copy()

    del imported_electron_data
    del imported_proton_data

    # Sum detectors 7 and 8 after the data is filtered
    electron_data_filtered['Detector7_8_sum'] = electron_data_filtered['Detector7'] + electron_data_filtered['Detector8']
    proton_data_filtered['Detector7_8_sum'] = proton_data_filtered['Detector7'] + proton_data_filtered['Detector8']

    # Filter values < 0.1 to 0 for specific columns
    cols_to_modify = [col for col in electron_data_filtered.columns if col not in ["E_Inc", "Detector7", "Detector8", "Particle_Type"]]

    electron_data_filtered[cols_to_modify] = electron_data_filtered[cols_to_modify].where(electron_data_filtered[cols_to_modify] >= 0.1, 0)
    proton_data_filtered[cols_to_modify] = proton_data_filtered[cols_to_modify].where(proton_data_filtered[cols_to_modify] >= 0.1, 0)

    # Randomly sample row indices from the filtered data
    num_electron_rows = len(electron_data_filtered)
    num_proton_rows = len(proton_data_filtered)

    # ==========================================
    # Generate Mutually Exclusive Indices
    # ==========================================
    
    # For Electrons: Shuffle all indices, then slice
    np.random.seed(42)
    all_e_indices = np.random.permutation(num_electron_rows)
    
    e_training_indices   = all_e_indices[:num_training_rows]
    e_validation_indices = all_e_indices[num_training_rows : num_training_rows + num_validation_rows]
    e_test_indices       = all_e_indices[num_training_rows + num_validation_rows : num_training_rows + num_validation_rows + num_test_rows]

    # For Protons: Shuffle all indices, then slice
    np.random.seed(43)
    all_p_indices = np.random.permutation(num_proton_rows)
    
    p_training_indices   = all_p_indices[:num_training_rows]
    p_validation_indices = all_p_indices[num_training_rows : num_training_rows + num_validation_rows]
    p_test_indices       = all_p_indices[num_training_rows + num_validation_rows : num_training_rows + num_validation_rows + num_test_rows]

    # Select the random rows and combine datasets
    training_data = pd.concat([
        electron_data_filtered.iloc[e_training_indices], 
        proton_data_filtered.iloc[p_training_indices]
    ], ignore_index=True)

    validation_data = pd.concat([
        electron_data_filtered.iloc[e_validation_indices], 
        proton_data_filtered.iloc[p_validation_indices]
    ], ignore_index=True)

    test_data = pd.concat([
        electron_data_filtered.iloc[e_test_indices], 
        proton_data_filtered.iloc[p_test_indices]
    ], ignore_index=True)

    # Save each data frame to a separate text file (.csv)
    training_data.to_csv("training_data2.csv", index=False)
    validation_data.to_csv("validation_data2.csv", index=False)
    test_data.to_csv("test_data2.csv", index=False)

    del electron_data_filtered
    del proton_data_filtered
    print("Data processed and saved to CSVs.")

else:
    # Load directly from CSVs to save time
    print("Loading data from pre-processed CSVs...")
    training_data = pd.read_csv("training_data2.csv")
    validation_data = pd.read_csv("validation_data2.csv")
    test_data = pd.read_csv("test_data2.csv")
    print("Data loaded.")

# ==========================================
# COMMON POST-PROCESSING
# ==========================================
# Convert Particle_Type to an ordered Categorical data type
cat_type = pd.CategoricalDtype(categories=["Electron", "Proton"], ordered=True)
training_data['Particle_Type'] = training_data['Particle_Type'].astype(cat_type)
validation_data['Particle_Type'] = validation_data['Particle_Type'].astype(cat_type)
test_data['Particle_Type'] = test_data['Particle_Type'].astype(cat_type)

# Identify ideal HERT test data
HERT_data_electrons = test_data[(test_data['Particle_Type'] == 'Electron') & 
                                (test_data['E_Inc'] >= 0.6) & 
                                (test_data['E_Inc'] <= 7.5)]

HERT_data_protons = test_data[(test_data['Particle_Type'] == 'Proton') & 
                              (test_data['E_Inc'] >= 14) & 
                              (test_data['E_Inc'] <= 70)]

HERT_data = pd.concat([HERT_data_electrons, HERT_data_protons], ignore_index=True)
HERT_data = test_data.copy() # Overwriting as specified in original script

# Combine 10% of data
training_data_sample = training_data.sample(n=int(num_training_rows * 0.1), random_state=1, replace=False)

# Combine 1% of data
training_data_sample2 = training_data.sample(n=int(num_training_rows * 0.01), random_state=2, replace=False)

#%% Create SVM Model
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
    
    X_test = HERT_data[features]
    y_test = HERT_data['Particle_Type']

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
    print("\n### Linear SVM ###")
    # In scikit-learn, setting class_weight works identically to class.weights in R
    svm_linear = SVC(kernel='linear', C=10, class_weight={'Electron': 1, 'Proton': 1})
    svm_linear.fit(X_train, y_train)
    
    # Get the support vectors and coefficients of the hyperplane (Intercept is prepended to match R output)
    linear_suppvec = svm_linear.n_support_
    print("Support Vectors (Electrons, Protons):", linear_suppvec)
    linear_hp_coefs = np.insert(svm_linear.coef_[0], 0, svm_linear.intercept_[0])
    print("Linear Hyperplane Coefficients:\n", linear_hp_coefs)
    
    # Calculate weights relative to the intercept (matching linear_hp_coefs[[1]] in R)
    weights = np.abs(linear_hp_coefs / linear_hp_coefs[0])
    print("Weights:\n", weights)
    
    # Predict to validate model
    linear_predictions = svm_linear.predict(X_test)
    
    # Print predictions (prop.table equivalent)
    linear_ct = pd.crosstab(index=linear_predictions, columns=y_test, rownames=['Predict'], normalize='columns') * 100
    print("\nLinear SVM Predictions (%):\n", linear_ct)

    # Test the hyperplane manually (filtered for Detector1 > 0.1 as in R script)
    hert_filtered = HERT_data[HERT_data['Detector1'] > 0.1]
    X_test_filtered = hert_filtered[features]
    
    # Using sklearn's built-in decision function is safer than manually multiplying,
    # but achieves the exact same math as your manual test block. 
    # > 0 assigns it to the positive class (Protons, based on sklearn's alphabetical class sorting).
    manual_decisions = svm_linear.decision_function(X_test_filtered)
    # sklearn sorts classes alphabetically: 0='Electron', 1='Proton'. So < 0 is Electron.
    linear_hp_test = np.where(manual_decisions < 0, "Electron", "Proton")
    
    hp_ct = pd.crosstab(index=linear_hp_test, columns=hert_filtered['Particle_Type'], rownames=['Predict'], normalize='columns') * 100
    print("\nManual Hyperplane Test (%):\n", hp_ct)

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

    print("\n### Radial SVM ###")
    svm_radial = SVC(kernel='rbf', C=10, class_weight={'Electron': 1, 'Proton': 1})
    svm_radial.fit(X_train, y_train)

    # Get the support vectors
    radial_suppvec = svm_radial.n_support_
    print("Support Vectors (Electrons, Protons):", radial_suppvec)

    radial_predictions = svm_radial.predict(X_test)
    radial_ct = pd.crosstab(index=radial_predictions, columns=y_test, rownames=['Predict'], normalize='columns') * 100
    print("Radial SVM Predictions (%):\n", radial_ct)

    print("\n### Simplified Linear SVM ###")
    features_si = ["Detector1", "Detector2"]
    X_train_si = training_data_sample2[features_si]
    y_train_si = training_data_sample2['Particle_Type']
    
    svm_linearsi = SVC(kernel='linear', C=10)
    svm_linearsi.fit(X_train_si, y_train_si)
    
    linearsi_hp_coefs = np.insert(svm_linearsi.coef_[0], 0, svm_linearsi.intercept_[0])
    print("Simplified Linear Coefficients:\n", linearsi_hp_coefs)
    
    linearsi_test = svm_linearsi.predict(X_test[features_si])
    linearsi_ct = pd.crosstab(index=linearsi_test, columns=y_test, rownames=['Predict'], normalize='columns') * 100
    print("\nSimplified Linear Predictions (%):\n", linearsi_ct)

    # Manual test for simplified hyperplane
    si_manual_decisions = svm_linearsi.decision_function(hert_filtered[features_si])
    linearsi_hp_test = np.where(si_manual_decisions < 0, "Electron", "Proton")
    
    hp_si_ct = pd.crosstab(index=linearsi_hp_test, columns=hert_filtered['Particle_Type'], rownames=['Predict'], normalize='columns') * 100
    print("\nManual Simplified Hyperplane Test (%):\n", hp_si_ct)

else:
    print("--- Skipping SVM Processing ---")

#%% Evaluate Instrument Performance Comparison
# Run legacy logic equations
evaluate_rept_logic(test_data)
evaluate_reptile2_logic(test_data)


#%% Explore Means and Medians
# ==========================================
# Conditional Means and Medians Matrices
# ==========================================

columns_to_keep = ["Detector1", "Detector2", "Detector3", "Detector4", 
                   "Detector5", "Detector6", "Detector7_8_sum", "Detector9", "Particle_Type"]
test_data_plot = test_data[columns_to_keep].copy()

# Extract just the detector names for the grid
features_only = columns_to_keep[:-1]

output_col_names = features_only
n_cols = len(output_col_names)

# Initialize empty DataFrames with 0.0
test_data_means_electrons = pd.DataFrame(0.0, index=output_col_names, columns=output_col_names)
test_data_means_protons   = pd.DataFrame(0.0, index=output_col_names, columns=output_col_names)
test_data_medians_electrons = pd.DataFrame(0.0, index=output_col_names, columns=output_col_names)
test_data_medians_protons   = pd.DataFrame(0.0, index=output_col_names, columns=output_col_names)

electrons = test_data_plot[test_data_plot['Particle_Type'] == 'Electron']
protons   = test_data_plot[test_data_plot['Particle_Type'] == 'Proton']

# Loop through each detector column index
for i in range(n_cols):
    for j in range(i + 1):  # Runs from 0 to i inclusive
        col_j = output_col_names[j]
        
        # Identify columns physically "higher" in the detector stack
        higher_cols = output_col_names[i+1:]
        
        if len(higher_cols) > 0:
            # Mask: Current col != 0 AND sum of all higher columns == 0
            mask_electrons = (electrons[col_j] != 0) & (electrons[higher_cols].sum(axis=1) == 0)
            mask_protons   = (protons[col_j] != 0) & (protons[higher_cols].sum(axis=1) == 0)
        else:
            # If checking the final column, there are no "higher" columns
            mask_electrons = (electrons[col_j] != 0)
            mask_protons   = (protons[col_j] != 0)
            
        # Calculate Mean
        test_data_means_electrons.iloc[i, j] = electrons.loc[mask_electrons, col_j].mean()
        test_data_means_protons.iloc[i, j]   = protons.loc[mask_protons, col_j].mean()
        
        # Calculate Median
        test_data_medians_electrons.iloc[i, j] = electrons.loc[mask_electrons, col_j].median()
        test_data_medians_protons.iloc[i, j]   = protons.loc[mask_protons, col_j].median()

# Fill any NaN values (caused by trying to find the mean/median of an empty slice) with 0
test_data_means_electrons.fillna(0, inplace=True)
test_data_means_protons.fillna(0, inplace=True)
test_data_medians_electrons.fillna(0, inplace=True)
test_data_medians_protons.fillna(0, inplace=True)

# Calculate Mean of Medians
mean_of_medians = (test_data_medians_electrons + test_data_medians_protons) / 2

# Display Results
print("\n--- Test Data Means (Electrons) ---")
print(test_data_means_electrons.round(3))

print("\n--- Test Data Means (Protons) ---")
print(test_data_means_protons.round(3))

print("\n--- Mean of Medians ---")
print(mean_of_medians.round(3))


#%% Simple Plot
# Find the max energies deposited onto each detector for both particle types
# Create a list of just the detector columns (if not already defined)
features_only = ["Detector1", "Detector2", "Detector3", "Detector4", 
                 "Detector5", "Detector6", "Detector7_8_sum", "Detector9"]

# Subset and calculate the max for each particle
e_max_energies = test_data_plot[test_data_plot['Particle_Type'] == 'Electron'][features_only].max()
p_max_energies = test_data_plot[test_data_plot['Particle_Type'] == 'Proton'][features_only].max()

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

# ==========================================
# 1. Plot Edep v Einc (Single Scatter)
# ==========================================
plt.figure(figsize=(8, 6))

# Map colors based on Particle_Type
colors = np.where(training_data['Particle_Type'] == 'Electron', 'red', 'blue')

plt.scatter(training_data['Detector4'], training_data['Detector3'], 
            c=colors, s=1, alpha=0.5)

plt.xlim(0, 25)
plt.ylim(0, 25)
plt.xlabel('E_dep on Detector 4 (MeV)', fontsize=14, fontweight='bold')
plt.ylabel('E_dep on Detector 3 (MeV)', fontsize=14, fontweight='bold')
plt.tick_params(axis='both', labelsize=12)

# Custom Legend
red_patch = mpatches.Patch(color='red', label='Electron')
blue_patch = mpatches.Patch(color='blue', label='Proton')
plt.legend(handles=[red_patch, blue_patch], loc='upper right', frameon=False, fontsize=12)

plt.tight_layout()
plt.show()

# ==========================================
# 2. Pairs Plot with Hyperplanes
# ==========================================
def plot_hyperplane_scatter(x, y, **kwargs):
    """Custom function to plot scatter points and overlay the SVM hyperplane."""
    ax = plt.gca()
    
    # Scatter plot (hue and palette are now automatically inherited via **kwargs)
    sns.scatterplot(x=x, y=y, s=5, ax=ax, legend=False, linewidth=0, alpha=0.5, **kwargs)
    
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 20)
    
    # Identify indices of the current variables to grab the right coefficients
    i = features_only.index(x.name)
    j = features_only.index(y.name)
    
    # linear_hp_coefs mapping: [Intercept, D1, D2, D3, D4, D5, D6, D7_8_sum, D9]
    coef_x = linear_hp_coefs[i + 1]
    coef_y = linear_hp_coefs[j + 1]
    intercept = linear_hp_coefs[0]
    
    # Avoid division by zero if a coefficient is 0
    if coef_y != 0:
        slope = -coef_x / coef_y
        int_y = -intercept / coef_y
        
        # Plot the line across the 0-20 limits
        x_vals = np.array([0, 20])
        y_vals = int_y + slope * x_vals
        ax.plot(x_vals, y_vals, color='black', linestyle='--', linewidth=2)


print("Generating PairGrid. This may take a moment...")

# 1. DEFINE HUE AND PALETTE GLOBALLY HERE
g = sns.PairGrid(test_data_plot, vars=features_only, height=2.5, 
                 hue='Particle_Type', palette={'Electron': 'red', 'Proton': 'blue'})

# 2. Map the custom scatter+line function to the off-diagonal plots
g.map_offdiag(plot_hyperplane_scatter)

# 3. Map standard histograms to the diagonal (it automatically uses the global hue)
g.map_diag(sns.histplot, multiple="stack")

plt.show()

#%% Create and Save PairGrid with Density Plots
print("\n--- Initiating Visualization Pipeline ---")
# Define exactly what columns go into the PairGrid matrix
columns_to_keep = ["Detector1", "Detector2", "Detector3", "Detector4", 
                   "Detector5", "Detector6", "Detector7_8_sum", "Detector9", "Particle_Type"]
test_data_plot = test_data[columns_to_keep].copy()

# Features for tracking weights/indices
features_only = ["Detector1", "Detector2", "Detector3", "Detector4", 
                 "Detector5", "Detector6", "Detector7_8_sum", "Detector9"]

# ## Check Densities for Detector1 and Detector2
# # Generate the edges array using your main dataset
# edges = get_bin_edges(test_data_plot, features_only)
# # Now run the diagnostic function using the edges you just calculated
# d1_array = check_densities(test_data_plot, 'Detector1', edges)
# d9_array = check_densities(test_data_plot, 'Detector9', edges)

# def print_particle_max_densities(density_df, detector_name):
#     """Prints the maximum density and corresponding row for electrons and protons."""
#     print(f"\n--- Max Density: {detector_name} ---")
    
#     # Electron max
#     e_idx = density_df['Electron_Pct'].idxmax()
#     print("Electron Peak:")
#     print(density_df.loc[e_idx].to_string())
#     print(f"Max Electron Pct: {density_df['Electron_Pct'].max():.4f}%\n")
    
#     # Proton max
#     p_idx = density_df['Proton_Pct'].idxmax()
#     print("Proton Peak:")
#     print(density_df.loc[p_idx].to_string())
#     print(f"Max Proton Pct: {density_df['Proton_Pct'].max():.4f}%")

# # --- Usage ---
# print_particle_max_densities(d1_array, "Detector 1")
# print_particle_max_densities(d9_array, "Detector 9")

# Call the external function. (Requires the 'weights' variable from your Linear SVM block)
create_and_save_pairgrid(test_data_plot, weights, features_only, 
                         output="square_pairs_plot_with_vertical_legends.png")
                         

#%% Save workspace
import joblib
import types
import pickle
import os

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
    save_path = os.path.join(os.getcwd(), 'entire_workspace.joblib')
    print(f"Compressing {len(workspace)} variables...")
    joblib.dump(workspace, save_path)
    print(f"Workspace successfully saved to: {save_path}")
else:
    print("Workspace saving skipped.")

#%% Load workspace
import joblib

# Load the dictionary back into memory
load_path = 'entire_workspace.joblib'
workspace = joblib.load(load_path)

# Unpack the dictionary directly into your global environment
globals().update(workspace)

print(f"Successfully loaded {len(workspace)} variables into the workspace.")
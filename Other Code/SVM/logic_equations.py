import pandas as pd
import numpy as np

def Rxy(df, x, y):
    """Sums energy across detectors x to y (inclusive)"""
    cols = [f"Detector{i}" for i in range(x, y+1)]
    return df[cols].sum(axis=1)

def Rbarexy(df, x, y):
    """Checks if all detectors from x to y are < 0.4"""
    cols = [f"Detector{i}" for i in range(x, y+1)]
    return (df[cols] < 0.4).all(axis=1)

def Rbarpxy(df, x, y):
    """Checks if all detectors from x to y are < 0.5"""
    cols = [f"Detector{i}" for i in range(x, y+1)]
    return (df[cols] < 0.5).all(axis=1)

def evaluate_reptile2_logic(test_data):
    print("\n--- Evaluating Khoo 2022 (REPTile-2) Logic Equations ---")
    
    detector_cols = [f"Detector{i}" for i in range(1, 10)]
    total_energy = test_data[detector_cols].sum(axis=1)
    
    D1 = test_data['Detector1']
    D2 = test_data['Detector2']
    D3 = test_data['Detector3']
    D4 = test_data['Detector4']
    truth = test_data['Particle_Type']

    slant_eq_D12 = (D1 / 2.8 + D2 / 4.2) > 1
    slant_eq_D34 = (D3 / 13.5 + D2 / 30) > 1 
    
    rng_p = (D1 > 0.1) & slant_eq_D12 & (~(D4 > 0.1) | slant_eq_D12) & (total_energy <= 35)
    pen_p = (D1 > 0.1) & slant_eq_D12 & ((D4 > 0.1) | ~slant_eq_D12) & (total_energy <= 35)
    rng_e = (D1 > 0.1) & ~slant_eq_D12 & ~(D4 > 0.1) & (total_energy <= 4)
    pen_e = (D1 > 0.1) & ~slant_eq_D12 & (D4 > 0.1) & (total_energy <= 4)

    d1_mask = D1 > 0.1
    
    # FIXED: Added .values to all truth references to prevent index misalignment
    predict_rnge = np.where(rng_e[d1_mask], "Electron", "Rejected RNG_E")
    khoo_rngetab = (pd.crosstab(predict_rnge, truth[d1_mask].values, rownames=['Predict'], normalize='columns') * 100).round(3)
    print("\nRNG_E Logic (%):\n", khoo_rngetab)

    predict_pene = np.where(pen_e[d1_mask], "Electron", "Rejected PEN_E")
    khoo_penetab = (pd.crosstab(predict_pene, truth[d1_mask].values, rownames=['Predict'], normalize='columns') * 100).round(3)
    print("\nPEN_E Logic (%):\n", khoo_penetab)
    
    predict_e = np.where(rng_e | pen_e, "Electron", "Rejected Electron")
    khoo_etab = (pd.crosstab(predict_e, truth.values, rownames=['Predict'], normalize='columns') * 100).round(3)
    print("\nTotal Electron Logic (%):\n", khoo_etab)

    predict_rngp = np.where(rng_p, "Proton", "Rejected RNG_P")
    khoo_rngptab = (pd.crosstab(predict_rngp, truth.values, rownames=['Predict'], normalize='columns') * 100).round(3)
    print("\nRNG_P Logic (%):\n", khoo_rngptab)

    predict_penp = np.where(pen_p, "Proton", "Rejected PEN_P")
    khoo_penptab = (pd.crosstab(predict_penp, truth.values, rownames=['Predict'], normalize='columns') * 100).round(3)
    print("\nPEN_P Logic (%):\n", khoo_penptab)
    
    predict_p = np.where(rng_p | pen_p, "Proton", "Rejected Proton")
    khoo_ptab = (pd.crosstab(predict_p, truth.values, rownames=['Predict'], normalize='columns') * 100).round(3)
    print("\nTotal Proton Logic (%):\n", khoo_ptab)


def evaluate_rept_logic(test_data):
    print("\n--- Evaluating Baker 2013 (REPT) Logic Equations ---")
    
    E_Inc = test_data['E_Inc']
    truth = test_data['Particle_Type']
    D1 = test_data['Detector1']
    D2 = test_data['Detector2']
    D9 = test_data['Detector9']

    # 1. Create the Incident Energy Masks (1.6 - 8 MeV for e-, 18 - 75 MeV for p+)
    valid_electrons = (truth == 'Electron') & (E_Inc >= 1.6) & (E_Inc <= 8.0)
    valid_protons   = (truth == 'Proton') & (E_Inc >= 18.0) & (E_Inc <= 75.0)
    energy_mask     = valid_electrons | valid_protons

    # Electron Logic Equations
    EL1 = (D1 >= 1.0) & (D1 <= 1.2) & (D2 <= 1.5) & (Rxy(test_data, 1, 2) >= 1.1) & (Rxy(test_data, 1, 2) <= 1.2) & Rbarexy(test_data, 3, 9)
    EL2 = (D1 >= 0.4) & (D2 >= 0.4) & (Rxy(test_data, 1, 2) >= 1.3) & (Rxy(test_data, 1, 2) <= 1.7) & Rbarexy(test_data, 3, 9)
    EL3 = (D1 >= 0.4) & (D2 >= 0.4) & (Rxy(test_data, 1, 4) >= 1.85) & (Rxy(test_data, 1, 4) <= 2.25) & Rbarexy(test_data, 5, 9)
    EL4 = (D1 >= 0.4) & (D2 >= 0.4) & (Rxy(test_data, 1, 4) >= 2.65) & (Rxy(test_data, 1, 4) <= 2.95) & Rbarexy(test_data, 5, 9)
    EL5 = (D1 >= 0.4) & (Rxy(test_data, 2, 4) >= 0.4) & (Rxy(test_data, 1, 6) >= 3.35) & (Rxy(test_data, 1, 6) <= 3.95) & Rbarexy(test_data, 7, 9)
    EL6 = (D1 >= 0.4) & (Rxy(test_data, 2, 6) >= 0.4) & (Rxy(test_data, 1, 8) >= 4.4) & (Rxy(test_data, 1, 8) <= 5.0) & (D9 < 0.4)
    EL7 = (D1 >= 0.4) & (D1 <= 2.0) & (D2 >= 0.4) & (D2 <= 2.0) & (Rxy(test_data, 3, 6) >= 0.4) & (Rxy(test_data, 1, 8) >= 5.5) & (Rxy(test_data, 1, 8) <= 6.25) & (D9 < 0.4)
    EL8 = (D1 >= 0.4) & (D2 >= 0.4) & (D2 <= 1.0) & (Rxy(test_data, 3, 6) >= 2.4) & (Rxy(test_data, 3, 9) >= 5.75) & (Rxy(test_data, 3, 9) <= 6.6)
    EL9 = (D1 >= 0.4) & (D2 >= 0.4) & (D2 <= 1.0) & (Rxy(test_data, 3, 4) >= 0.4) & (Rxy(test_data, 3, 4) <= 2.0) & (Rxy(test_data, 5, 6) >= 0.4) & (Rxy(test_data, 7, 9) >= 0.4) & (Rxy(test_data, 3, 9) >= 8.0) & (Rxy(test_data, 3, 9) <= 9.0)
    EL10 = (D1 >= 0.4) & (D2 >= 0.4) & (Rxy(test_data, 3, 4) >= 0.4) & (Rxy(test_data, 5, 6) >= 0.4) & (Rxy(test_data, 7, 8) >= 0.4) & (D9 >= 0.1) & (Rxy(test_data, 3, 9) >= 10.3) & (Rxy(test_data, 3, 9) <= 12.5)
    EL11 = (D1 >= 0.4) & (D1 <= 1.0) & (D2 >= 0.4) & (Rxy(test_data, 3, 4) >= 0.4) & (Rxy(test_data, 5, 9) >= 0.4) & (Rxy(test_data, 7, 9) >= 11)
    EL12 = (D1 >= 0.4) & (D1 <= 1.0) & (D2 >= 0.4) & (D2 <= 1.0) & (Rxy(test_data, 3, 4) >= 0.4) & (Rxy(test_data, 3, 4) <= 1.5) & (Rxy(test_data, 5, 9) >= 0.4) & (Rxy(test_data, 7, 9) >= 15)

    ELOGIC = pd.DataFrame({'EL1':EL1, 'EL2':EL2, 'EL3':EL3, 'EL4':EL4, 'EL5':EL5, 'EL6':EL6, 
                           'EL7':EL7, 'EL8':EL8, 'EL9':EL9, 'EL10':EL10, 'EL11':EL11, 'EL12':EL12})
    
    # 2. Combine the D1 threshold mask with the Energy Range mask
    eval_mask = (D1 > 0.4) & energy_mask
    
    E_problem_rows = (ELOGIC[eval_mask].sum(axis=1) > 1).sum()
    ratio_Eproblem = E_problem_rows / eval_mask.sum()
    print(f"\nRatio of Filtered Particles triggering >1 Electron equation: {ratio_Eproblem * 100:.3f}%")

    # FIXED: Added .values to truth[eval_mask]
    Erow_predict = np.where(ELOGIC[eval_mask].sum(axis=1) > 0, "Electron", "Rejected Electron")
    REPT_etab = (pd.crosstab(Erow_predict, truth[eval_mask].values, rownames=['Predict'], normalize='columns') * 100).round(3)
    print("\nREPT Electron Logic (%):\n", REPT_etab)

    # Proton Logic Equations
    PL1 = (D1 > 8.2) & (D1 < 16) & (D2 < 6.5) & (Rxy(test_data, 1, 2) > 8.2) & (Rxy(test_data, 1, 2) < 18) & (Rxy(test_data, 3, 9) < 0.5) & Rbarpxy(test_data, 3, 9)
    PL2 = (D1 > 5.4) & (D1 < 12.2) & (D2 > 5.0) & (D2 < 16.9) & (Rxy(test_data, 3, 4) > 0.1) & (Rxy(test_data, 3, 4) < 11) & (Rxy(test_data, 1, 4) > 15.9) & (Rxy(test_data, 1, 4) < 25.7) & (Rxy(test_data, 5, 9) < 0.5) & Rbarpxy(test_data, 5, 9)
    PL3 = (D1 > 4) & (D1 < 7) & (D2 > 4) & (D2 < 9.5) & (Rxy(test_data, 3, 4) > 10) & (Rxy(test_data, 5, 6) < 12.5) & (Rxy(test_data, 1, 6) > 24) & (Rxy(test_data, 1, 6) < 35.5) & (Rxy(test_data, 7, 9) < 0.5) & Rbarpxy(test_data, 7, 9)
    PL4 = (D1 > 3.1) & (D1 < 4.9) & (D2 > 3.2) & (D2 < 5.7) & (Rxy(test_data, 3, 4) > 7.6) & (Rxy(test_data, 3, 4) < 16.8) & (Rxy(test_data, 5, 6) > 9.2) & (Rxy(test_data, 5, 6) < 24) & (Rxy(test_data, 7, 8) < 23) & (D9 < 4.1) & (Rxy(test_data, 5, 9) > 11.5) & (Rxy(test_data, 5, 9) < 33.0)
    PL5 = (D1 > 2.2) & (D1 < 4) & (D2 > 1.9) & (D2 < 4.2) & (Rxy(test_data, 3, 4) > 5.5) & (Rxy(test_data, 3, 4) < 12.5) & (Rxy(test_data, 5, 6) > 5.8) & (Rxy(test_data, 5, 6) < 12.5) & (Rxy(test_data, 7, 8) > 7) & (Rxy(test_data, 7, 8) < 22.7) & (D9 > 1) & (D9 < 13) & (Rxy(test_data, 7, 9) > 5) & (Rxy(test_data, 7, 9) < 45)
    PL6 = (D1 > 1.5) & (D1 < 3.3) & (D2 > 1.0) & (D2 < 3.3) & (Rxy(test_data, 3, 4) > 4.1) & (Rxy(test_data, 3, 4) < 6.5) & (Rxy(test_data, 5, 6) > 4.5) & (Rxy(test_data, 5, 6) < 7.2) & (Rxy(test_data, 7, 8) > 4.8) & (Rxy(test_data, 7, 8) < 8.0) & (D9 > 2.0) & (D9 < 8.5) & (Rxy(test_data, 1, 6) > 11) & (Rxy(test_data, 1, 6) < 22) & (Rxy(test_data, 1, 9) < 65)
    PL7 = (D1 > 1.4) & (D1 < 2.5) & (D2 > 1.4) & (D2 < 2.8) & (Rxy(test_data, 3, 4) > 3.4) & (Rxy(test_data, 3, 4) < 5.4) & (Rxy(test_data, 5, 6) > 3.4) & (Rxy(test_data, 5, 6) < 5.9) & (Rxy(test_data, 7, 8) > 3.5) & (Rxy(test_data, 7, 8) < 6) & (Rxy(test_data, 1, 9) > 10) & (Rxy(test_data, 1, 9) < 45)
    PL8 = (D1 > 0.8) & (D1 < 3.0) & (D2 > 0.8) & (D2 < 3.0) & (Rxy(test_data, 3, 4) > 2.5) & (Rxy(test_data, 3, 4) < 5) & (Rxy(test_data, 5, 6) > 2.5) & (Rxy(test_data, 5, 6) < 5.5) & (Rxy(test_data, 7, 8) > 2.5) & (Rxy(test_data, 7, 8) < 5.5) & (D9 > 1) & (D9 < 6) & (Rxy(test_data, 1, 9) < 8) & (Rxy(test_data, 1, 9) < 32.0)

    PLOGIC = pd.DataFrame({'PL1':PL1, 'PL2':PL2, 'PL3':PL3, 'PL4':PL4, 'PL5':PL5, 'PL6':PL6, 'PL7':PL7, 'PL8':PL8})

    P_problem_rows = (PLOGIC[eval_mask].sum(axis=1) > 1).sum()
    ratio_Pproblem = P_problem_rows / eval_mask.sum()
    print(f"\nRatio of Filtered Particles triggering >1 Proton equation: {ratio_Pproblem * 100:.3f}%")

    # FIXED: Added .values to truth[eval_mask]
    Prow_predict = np.where(PLOGIC[eval_mask].sum(axis=1) > 0, "Proton", "Rejected Proton")
    REPT_ptab = (pd.crosstab(Prow_predict, truth[eval_mask].values, rownames=['Predict'], normalize='columns') * 100).round(3)
    print("\nREPT Proton Logic (%):\n", REPT_ptab)
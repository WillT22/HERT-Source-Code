# ==========================================
# HELPER FUNCTIONS
# ==========================================
import numpy as np
import pandas as pd
from scipy.integrate import quad

def get_bin_indices(energies, bins):
    """Vectorized helper to map energies to clipped bin indices."""
    return np.clip(np.digitize(energies, bins) - 1, 0, len(bins) - 2)

def get_available_counts(df, bins):
    """Finds exactly how many physical particles exist in each bin."""
    bin_indices = get_bin_indices(df['E_Inc'].values, bins)
    return np.bincount(bin_indices, minlength=len(bins) - 1)

def assign_Normalized_Flat_Weights(df, e_bins, p_bins):
    """
    Computes flat bin weights for both species and normalizes each class 
    to sum to 1.0 (50/50 SVM class balance) in a single pass.
    """
    weights = np.zeros(len(df), dtype=float)
    
    for p_type, bins in [('Electron', e_bins), ('Proton', p_bins)]:
        mask = (df['Particle_Type'] == p_type).values
        if not np.any(mask):
            continue
            
        bin_indices = get_bin_indices(df.loc[mask, 'E_Inc'].values, bins)
        counts = np.bincount(bin_indices, minlength=len(bins) - 1)
        
        # Each occupied bin gets equal total weight; divide by count per bin
        w_lookup = np.zeros_like(counts, dtype=float)
        nonzero = counts > 0
        w_lookup[nonzero] = 1.0 / counts[nonzero]
        
        particle_weights = w_lookup[bin_indices]
        # Normalize so the sum of weights for this species equals 1.0
        weights[mask] = particle_weights / np.sum(particle_weights)
        
    return weights

def allocate_weighted_counts(raw_probs, avail_counts, total_points):
    """
    Universal water-filling allocator. Redistributes slack from physically 
    capped bins proportionally to raw_probs (handles flat, exp, and power-law).
    """
    sample_counts = np.zeros_like(avail_counts, dtype=int)
    remaining = int(total_points)
    
    while remaining > 0:
        open_mask = (avail_counts > sample_counts) & (raw_probs > 0)
        open_bins = np.where(open_mask)[0]
        
        if len(open_bins) == 0:
            print(f"Warning: Exhausted all physical particles in valid bins! Short by {remaining}.")
            break
            
        open_probs = raw_probs[open_bins]
        norm_probs = open_probs / np.sum(open_probs)
        
        exact_add = norm_probs * remaining
        int_add = np.floor(exact_add).astype(int)
        
        if np.sum(int_add) == 0:
            remainders = exact_add - int_add
            priority_order = np.argsort(remainders)[::-1]
            for idx in priority_order[:remaining]:
                sample_counts[open_bins[idx]] += 1
                remaining -= 1
            break
            
        can_take = avail_counts[open_bins] - sample_counts[open_bins]
        actual_add = np.minimum(int_add, can_take)
        
        if np.sum(actual_add) == 0:
            priority_order = np.argsort(exact_add)[::-1]
            for idx in priority_order:
                target_bin = open_bins[idx]
                if avail_counts[target_bin] > sample_counts[target_bin]:
                    sample_counts[target_bin] += 1
                    remaining -= 1
                    break
        else:
            sample_counts[open_bins] += actual_add
            remaining -= np.sum(actual_add)

    return sample_counts

def generate_flat_counts(avail_counts, total_points):
    """Flat distribution is simply equal raw_probs across all bins."""
    raw_probs = np.ones_like(avail_counts, dtype=float)
    return allocate_weighted_counts(raw_probs, avail_counts, total_points)

def generate_exponential_counts(e_bins, avail_counts, total_points, E0=0.5):
    e_lower, e_upper = e_bins[:-1], e_bins[1:]
    raw_probs = np.exp(-e_lower / E0) - np.exp(-e_upper / E0)
    return allocate_weighted_counts(raw_probs, avail_counts, total_points)

def generate_powerlaw_counts(p_bins, avail_counts, total_points, j0=1.32e6, V_sw=460.9, kappa0=3.46e18, alpha=1.94, gamma=3.75, r=1.0):
    conv_diff = (4.0 / 9.0) * ((V_sw * 1e5) / kappa0) * (r * 1.495978707e13)
    flux_fn = lambda E: j0 * (E**alpha + conv_diff)**(-gamma / alpha)

    p_lower, p_upper = p_bins[:-1], p_bins[1:]
    raw_probs = np.array([quad(flux_fn, lo, hi)[0] for lo, hi in zip(p_lower, p_upper)])
    return allocate_weighted_counts(raw_probs, avail_counts, total_points)

def sample_spectrum(df, target_counts, bins, random_state=42):
    """Fast, single-pass stratified sampler without replacement."""
    rng = np.random.default_rng(random_state)
    bin_indices = get_bin_indices(df['E_Inc'].values, bins)
    
    selected_positions = []
    # Group positional indices by bin in a single pass via pandas groupby
    for bin_idx, pos_indices in pd.Series(np.arange(len(df))).groupby(bin_indices):
        target_n = target_counts[bin_idx]
        if target_n <= 0:
            continue
        avail_n = len(pos_indices)
        sample_n = min(target_n, avail_n)
        if sample_n < target_n:
            print(f"Warning: Bin {bin_idx} requested {target_n} but only has {avail_n}.")
        selected_positions.append(rng.choice(pos_indices.values, size=sample_n, replace=False))
        
    all_positions = np.concatenate(selected_positions)
    return df.iloc[all_positions].reset_index(drop=True)
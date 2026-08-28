import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
from PIL import Image

def get_bin_edges(data, features, nbins=100):
    """Generates logspace edges from exactly 0.1 to 100 MeV."""
    # Force the bins to perfectly cover the 0.1 to 100 visual space
    edges = np.logspace(np.log10(0.1), np.log10(100.0), nbins + 1)
    return edges

def get_density_limits(data, particle_type, edges, features):
    """Finds the true physical minimum (1 particle) and global maximum density."""
    subset = data[data['Particle_Type'] == particle_type][features]
    max_density = 0.0
    
    # True physical minimum: 1 particle out of the total particles
    total_particles = len(subset)
    min_density = 1.0 / total_particles if total_particles > 0 else 1e-4

    for i in range(len(features) - 1):
        for j in range(i + 1, len(features)):
            H, _, _ = np.histogram2d(subset.iloc[:, i], subset.iloc[:, j], bins=(edges, edges))
            total = H.sum()
            if total > 0:
                frac = H.max() / total
                if frac > max_density:
                    max_density = frac

    # Round max up for a clean upper limit on the colorbar
    max_density = np.ceil(max_density * 2000) / 2000
    
    return min_density, max_density

# ==========================================
# Seaborn PairGrid Mapping Functions
# ==========================================

def plot_density_lower(x, y, **kwargs):
    """Lower triangle: Electron 2D density (Red)"""
    ax = plt.gca()
    df = kwargs['full_data']
    edges = kwargs['edges']
    # features = kwargs['features'] # Unused, can be removed

    # Subset electrons for this specific x/y pair
    subset = df[df['Particle_Type'] == 'Electron']
    x_data = subset[x.name]
    y_data = subset[y.name]

    H, xedges, yedges = np.histogram2d(x_data, y_data, bins=(edges, edges))
    
    # Normalizes by the total number of electrons in the entire dataset
    total_electrons = len(subset) 
    H_frac = H.T / max(1, total_electrons) 

    # Mask zeros to render as white
    H_masked = np.ma.masked_where(H_frac == 0, H_frac)

    orig_cmap = plt.cm.Reds
    # Extract colors from the 0.1 to 1.0 range
    truncated_colors = orig_cmap(np.linspace(0.1, 1.0, 256)) 
    cmap = mcolors.ListedColormap(truncated_colors)
    norm = mcolors.LogNorm(vmin=min(kwargs['e_min'], kwargs['p_min']), vmax=max(kwargs['e_max'], kwargs['p_max']))

    ax.pcolormesh(xedges, yedges, H_masked, cmap=cmap, norm=norm)

def plot_density_upper(x, y, **kwargs):
    """Upper triangle: Proton 2D density (Blue)"""
    ax = plt.gca()
    df = kwargs['full_data']
    edges = kwargs['edges']
    # features = kwargs['features'] # Unused, can be removed

    subset = df[df['Particle_Type'] == 'Proton']
    x_data = subset[x.name]
    y_data = subset[y.name]

    H, xedges, yedges = np.histogram2d(x_data, y_data, bins=(edges, edges))
    
    # Normalizes by the total number of protons in the entire dataset
    total_protons = len(subset)
    H_frac = H.T / max(1, total_protons)

    H_masked = np.ma.masked_where(H_frac == 0, H_frac)

    orig_cmap = plt.cm.Blues
    # Extract colors from the 0.1 to 1.0 range
    truncated_colors = orig_cmap(np.linspace(0.1, 1.0, 256)) 
    cmap = mcolors.ListedColormap(truncated_colors)
    norm = mcolors.LogNorm(vmin=min(kwargs['p_min'], kwargs['e_min']), vmax=max(kwargs['p_max'], kwargs['e_max']))

    ax.pcolormesh(xedges, yedges, H_masked, cmap=cmap, norm=norm)


def plot_diag_1d(x, **kwargs):
    ax = plt.gca()
    df = kwargs['full_data']
    edges = kwargs['edges']
    # features = kwargs['features'] # Unused, can be removed

    col_name = x.name

    e_data = df[df['Particle_Type'] == 'Electron'][col_name]
    p_data = df[df['Particle_Type'] == 'Proton'][col_name]

    e_hist, _ = np.histogram(e_data, bins=edges)
    p_hist, _ = np.histogram(p_data, bins=edges)

    # --- CORRECTION: Divide by total particles, not just in-bounds particles ---
    total_e = len(e_data)
    total_p = len(p_data)
    e_ratio = e_hist / max(1, total_e)
    p_ratio = p_hist / max(1, total_p)
    # -------------------------------------------------------------------------

    # Calculate the width of every individual bin
    bin_widths = np.diff(edges)

    # Divide the fraction by the bin width to get density
    e_density = e_ratio / bin_widths
    p_density = p_ratio / bin_widths

    midpoints = np.sqrt(edges[:-1] * edges[1:])

    # --- 1. BASE AXIS (Log Scale) ---
    ax.set_xscale('log')
    ax.set_xlim(0.1, 100)
    ax.set_yscale('log')
    ax.set_ylim(0.1, 100)
    
    ax.yaxis.set_visible(False)

    # --- 2. TWIN AXIS (Linear Scale) ---
    ax2 = ax.twinx()
    
    density_red = plt.cm.Reds(0.8)
    density_blue = plt.cm.Blues(0.8)
    
    ax2.plot(midpoints, e_density, color=density_red, lw=2)
    ax2.plot(midpoints, p_density, color=density_blue, lw=2)

    ax2.set_yscale('linear')
    
    max_density = max(np.max(e_density), np.max(p_density))
    y_max = max_density * 1.1 if max_density > 0 else 1.0 
    ax2.set_ylim(0, y_max)
    
    # The 'steps' array restricts the algorithm to 1-decimal friendly intervals
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator(n=4))
    
    # Format labels to a max of 1 decimal place, stripping trailing zeros
    ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{x:.1f}".rstrip('0').rstrip('.')))

    ax2.yaxis.tick_left()
    ax2.yaxis.set_label_position("left")

    # --- 3. DYNAMIC LABELS ---
    ax2.tick_params(axis='y', which='both', left=True, labelleft=True)
    
    ax2.text(0.95, 1.25, "PDF (MeV$^{-1}$)", transform=ax2.transAxes,
             ha='right', va='top', fontsize=24, style='italic', color='gray')


def check_densities(df, col_name, edges):
    """Calculates and returns the density arrays for a specific detector."""
    e_data = df[df['Particle_Type'] == 'Electron'][col_name]
    p_data = df[df['Particle_Type'] == 'Proton'][col_name]

    e_hist, _ = np.histogram(e_data, bins=edges)
    p_hist, _ = np.histogram(p_data, bins=edges)

    # --- CORRECTION: Divide by total particles ---
    total_e = len(e_data)
    total_p = len(p_data)
    e_ratio = e_hist / max(1, total_e)
    p_ratio = p_hist / max(1, total_p)
    # ---------------------------------------------

    midpoints = np.sqrt(edges[:-1] * edges[1:])

    density_df = pd.DataFrame({
        'Bin_Center_MeV': midpoints,
        'Electron_Ratio': e_ratio,
        'Proton_Ratio': p_ratio
    })

    active_bins = density_df[(density_df['Electron_Ratio'] > 0) | (density_df['Proton_Ratio'] > 0)]

    return active_bins

# ==========================================
# Main Generation Function
# ==========================================

def create_and_save_pairgrid(test_data_plot, features_only, output="square_pairs_plot.png"):
    label_mapping = {
        "Detector1": "D1", "Detector2": "D2", "Detector3": "D3",
        "Detector4": "D4", "Detector5": "D5", "Detector6": "D6",
        "Detector7_8_sum": "D7+8", "Detector9": "D9"
    }

    plot_df = test_data_plot.rename(columns=label_mapping)
    display_features = [label_mapping.get(f, f) for f in features_only]

    print("Calculating bin edges and max densities...")
    edges = get_bin_edges(plot_df, display_features)
    e_min, e_max = get_density_limits(plot_df, "Electron", edges, display_features)
    p_min, p_max = get_density_limits(plot_df, "Proton", edges, display_features)

    print("Generating PairGrid mapping... (This may take a few minutes)")

    sns.set_theme(style="ticks")

    plt.rcParams.update({
        'font.size': 36,
        'axes.labelsize': 48,
        'axes.labelweight': 'bold',
        'xtick.labelsize': 28,
        'ytick.labelsize': 28,

        'xtick.major.size': 12,
        'xtick.major.width': 3,
        'xtick.minor.size': 6,
        'xtick.minor.width': 2,

        'ytick.major.size': 12,
        'ytick.major.width': 3,
        'ytick.minor.size': 6,
        'ytick.minor.width': 2
    })

    g = sns.PairGrid(plot_df, vars=display_features, height=3, aspect=1, diag_sharey=False)

    map_kwargs = {
        'full_data': plot_df,
        'edges': edges,
        'e_min': e_min,
        'e_max': e_max,
        'p_min': p_min,
        'p_max': p_max,
        'features': display_features
    }

    g.map_lower(plot_density_lower, **map_kwargs)
    g.map_upper(plot_density_upper, **map_kwargs)
    g.map_diag(plot_diag_1d, **map_kwargs)

    for ax in g.axes.flat:
        if ax is not None:
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(1.5)
                spine.set_edgecolor('black')
            ax.set_box_aspect(1) 

    g.fig.supxlabel("Deposited Energy (MeV)", fontsize=40, fontweight='bold', y=-0.01)
    g.fig.supylabel("Deposited Energy (MeV)", fontsize=40, fontweight='bold', x=-0.01)

    major_ticks = [0.1, 1, 10, 100]
    major_tick_labels = ["0.1", "1", "10", "100"]

    for i in range(len(g.axes)):
        # 1. Left Column (Y-axis labels: D1, D2, D3...)
        ax_left = g.axes[i, 0]
        ax_left.yaxis.label.set_size(40)
        ax_left.yaxis.label.set_weight('bold')
        ax_left.yaxis.set_label_coords(-0.45, 0.45) 
        
        # 2. Bottom Row (X-axis labels: D1, D2, D3...)
        ax_bottom = g.axes[-1, i]
        ax_bottom.xaxis.label.set_size(40)
        ax_bottom.xaxis.label.set_weight('bold')
        ax_bottom.xaxis.labelpad = 20

        for j in range(len(g.axes)):
            ax = g.axes[i, j]
            if ax is None:
                continue

            # ==========================================
            # X-AXIS (Log Scale)
            # ==========================================
            ax.set_xscale('log')
            ax.set_xlim(0.1, 100)
            ax.set_xticks(major_ticks)
            ax.set_xticklabels(major_tick_labels) # Set labels globally
            ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)))
            ax.xaxis.set_minor_formatter(ticker.NullFormatter())

            # Toggle text visibility based on row position (bottom row only)
            if i == len(g.axes) - 1:
                ax.tick_params(axis='x', which='both', bottom=True, labelbottom=True)
            else:
                ax.tick_params(axis='x', which='both', bottom=True, labelbottom=False)

            # ==========================================
            # Y-AXIS (Log Scale)
            # ==========================================
            if i != j:
                ax.set_yscale('log')
                ax.set_ylim(0.1, 100)
                ax.set_yticks(major_ticks)
                ax.set_yticklabels(major_tick_labels) # Set labels globally
                ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)))
                ax.yaxis.set_minor_formatter(ticker.NullFormatter())

                # Toggle text visibility based on column position (left column only)
                if j == 0 or j == i+1:
                    ax.tick_params(axis='y', which='both', left=True, labelleft=True)
                else:
                    ax.tick_params(axis='y', which='both', left=True, labelleft=False)
            else:
                # Diagonal: Hide axis ticks and labels
                ax.tick_params(axis='y', which='both', left=False, labelleft=False)

# --- Stitch Horizontal Colorbar Legends ---
    # Remove the right margin adjustment; rely on savefig's bbox_inches="tight" 
    # to capture the axes below the plot.
    g.fig.subplots_adjust(wspace=0.5)

    # 1. ELECTRON COLORBAR (Bottom Left)
    # [left, bottom, width, height]
    cbar_ax_elec = g.fig.add_axes([0.1, -0.05, 0.35, 0.03])
    orig_cmap = plt.cm.Reds
    # Extract colors from the 0.1 to 1.0 range
    truncated_colors = orig_cmap(np.linspace(0.1, 1.0, 256)) 
    cmap_elec = mcolors.ListedColormap(truncated_colors)
    norm_elec = mcolors.LogNorm(vmin=min(e_min, p_min), vmax=max(e_max, p_max))
    
    cbar_elec = g.fig.colorbar(
        plt.cm.ScalarMappable(norm=norm_elec, cmap=cmap_elec), 
        cax=cbar_ax_elec, 
        orientation='horizontal'
    )
    
    # Text will automatically center below or above the horizontal bar
    cbar_elec.set_label("Electron Fraction", labelpad=10, fontweight='bold', fontsize=40)
    cbar_elec.ax.tick_params(labelsize=32)
    
    cbar_elec.locator = ticker.LogLocator(base=10.0, numticks=5)
    cbar_elec.formatter = ticker.LogFormatterMathtext()
    cbar_elec.update_ticks()

    # 2. PROTON COLORBAR (Bottom Right)
    cbar_ax_prot = g.fig.add_axes([0.55, -0.05, 0.35, 0.03])
    orig_cmap = plt.cm.Blues
    # Extract colors from the 0.1 to 1.0 range
    truncated_colors = orig_cmap(np.linspace(0.1, 1.0, 256)) 
    cmap_prot = mcolors.ListedColormap(truncated_colors)
    norm_prot = mcolors.LogNorm(vmin=min(e_min, p_min), vmax=max(e_max, p_max))
    
    cbar_prot = g.fig.colorbar(
        plt.cm.ScalarMappable(norm=norm_prot, cmap=cmap_prot), 
        cax=cbar_ax_prot, 
        orientation='horizontal'
    )
    
    cbar_prot.set_label("Proton Fraction", labelpad=10, fontweight='bold', fontsize=40)
    cbar_prot.ax.tick_params(labelsize=32)
    
    cbar_prot.locator = ticker.LogLocator(base=10.0, numticks=5)
    cbar_prot.formatter = ticker.LogFormatterMathtext()
    cbar_prot.update_ticks()

    print(f"Saving high-res image to {output}...")
    g.fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(g.fig)

    # --- RESET GLOBAL PLOTTING DEFAULTS ---
    plt.rcdefaults()
    sns.reset_defaults()

    print("Export complete.")


#%% For reduced number of detectors
def plot_diag_1d_reduced(x, **kwargs):
    ax = plt.gca()
    df = kwargs['full_data']
    edges = kwargs['edges']
    features = kwargs['features']

    col_name = x.name

    e_data = df[df['Particle_Type'] == 'Electron'][col_name]
    p_data = df[df['Particle_Type'] == 'Proton'][col_name]

    e_hist, _ = np.histogram(e_data, bins=edges)
    p_hist, _ = np.histogram(p_data, bins=edges)

    e_ratio = (e_hist / max(1, e_hist.sum()))
    p_ratio = (p_hist / max(1, p_hist.sum()))

    # Calculate the width of every individual bin
    bin_widths = np.diff(edges)

    # Divide the fraction by the bin width to get density
    e_density = e_ratio / bin_widths
    p_density = p_ratio / bin_widths

    midpoints = np.sqrt(edges[:-1] * edges[1:])

    # --- 1. BASE AXIS (Log Scale) ---
    # We maintain the log scale to satisfy Seaborn's row-sharing constraint
    ax.set_xscale('log')
    ax.set_xlim(0.1, 100)
    ax.set_yscale('log')
    ax.set_ylim(0.1, 100)
    
    # Disable the base Y-axis so it doesn't conflict
    ax.yaxis.set_visible(False)

    # --- 2. TWIN AXIS (Linear Scale) ---
    # Create the independent axis and plot the data here
    ax2 = ax.twinx()
    
    # Sample the exact colormaps used in the 2D plots (0.8 gives a strong, visible shade)
    density_red = plt.cm.Reds(0.8)
    density_blue = plt.cm.Blues(0.8)
    
    ax2.plot(midpoints, e_density, color=density_red, lw=2)
    ax2.plot(midpoints, p_density, color=density_blue, lw=2)

    ax2.set_yscale('linear')
    
    max_density = max(np.max(e_density), np.max(p_density))
    y_max = max_density * 1.1 if max_density > 0 else 1.0 
    ax2.set_ylim(0, y_max)
    
    # The 'steps' array restricts the algorithm to 1-decimal friendly intervals
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator(n=4))
    
    # Format labels to a max of 1 decimal place, stripping trailing zeros
    ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{x:.1f}".rstrip('0').rstrip('.')))

    # Move the twin axis from the right side to the left side
    ax2.yaxis.tick_left()
    ax2.yaxis.set_label_position("left")

    # --- 3. DYNAMIC LABELS ---
    ax2.tick_params(axis='y', which='both', left=True, labelleft=True)
    
    # Place the PDF label INSIDE the plot area (top-right corner) to save margin space
    ax2.text(0.8, 1.18, "PDF (MeV$^{-1}$)", transform=ax2.transAxes,
             ha='right', va='top', fontsize=18, style='italic', color='gray')

# ==========================================
# Main Generation Function
# ==========================================

def create_and_save_pairgrid_reduced(test_data_plot, features_only, output="square_pairs_plot.png"):
    label_mapping = {
        "Detector1": "D1", "Detector2": "D2", "Detector3": "D3",
        "Detector4": "D4", "Detector5": "D5", "Detector6": "D6",
        "Detector7_8_sum": "D7+8", "Detector9": "D9"
    }

    plot_df = test_data_plot.rename(columns=label_mapping)
    display_features = [label_mapping.get(f, f) for f in features_only]

    print("Calculating bin edges and max densities...")
    edges = get_bin_edges(plot_df, display_features)
    e_min, e_max = get_density_limits(plot_df, "Electron", edges, display_features)
    p_min, p_max = get_density_limits(plot_df, "Proton", edges, display_features)

    print("Generating PairGrid mapping... (This may take a few minutes)")

    sns.set_theme(style="ticks")

    plt.rcParams.update({
        'font.size': 22,
        'axes.labelsize': 22,
        'axes.labelweight': 'bold',
        'xtick.labelsize': 20,
        'ytick.labelsize': 20,

        'xtick.major.size': 12,
        'xtick.major.width': 3,
        'xtick.minor.size': 6,
        'xtick.minor.width': 2,

        'ytick.major.size': 12,
        'ytick.major.width': 3,
        'ytick.minor.size': 6,
        'ytick.minor.width': 2
    })

    g = sns.PairGrid(plot_df, vars=display_features, height=3, aspect=1, diag_sharey=False)

    map_kwargs = {
        'full_data': plot_df,
        'edges': edges,
        'e_min': e_min,
        'e_max': e_max,
        'p_min': p_min,
        'p_max': p_max,
        'features': display_features
    }

    g.map_lower(plot_density_lower, **map_kwargs)
    g.map_upper(plot_density_upper, **map_kwargs)
    g.map_diag(plot_diag_1d_reduced, **map_kwargs)

    for ax in g.axes.flat:
        if ax is not None:
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(1.5)
                spine.set_edgecolor('black')
            ax.set_box_aspect(1) 

    g.fig.supxlabel("Deposited Energy (MeV)", fontsize=22, fontweight='bold', y=0.005)
    g.fig.supylabel("Deposited Energy (MeV)", fontsize=22, fontweight='bold', x=-0.006)

    major_ticks = [0.1, 1, 10, 100]
    major_tick_labels = ["0.1", "1", "10", "100"]

    for i in range(len(g.axes)):
        # 1. Left Column (Y-axis labels: D1, D2, D3...)
        ax_left = g.axes[i, 0]
        ax_left.yaxis.label.set_size(20)
        ax_left.yaxis.label.set_weight('bold')
        ax_left.yaxis.set_label_coords(-0.35, 0.45) 
        
        # 2. Bottom Row (X-axis labels: D1, D2, D3...)
        ax_bottom = g.axes[-1, i]
        ax_bottom.xaxis.label.set_size(20)
        ax_bottom.xaxis.label.set_weight('bold')
        ax_bottom.xaxis.labelpad = 8

        for j in range(len(g.axes)):
            ax = g.axes[i, j]
            if ax is None:
                continue

            # ==========================================
            # X-AXIS (Log Scale)
            # ==========================================
            ax.set_xscale('log')
            ax.set_xlim(0.1, 100)
            ax.set_xticks(major_ticks)
            ax.set_xticklabels(major_tick_labels) # Set labels globally
            ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)))
            ax.xaxis.set_minor_formatter(ticker.NullFormatter())

            # Toggle text visibility based on row position (bottom row only)
            if i == len(g.axes) - 1:
                ax.tick_params(axis='x', which='both', bottom=True, labelbottom=True)
            else:
                ax.tick_params(axis='x', which='both', bottom=True, labelbottom=False)

            # ==========================================
            # Y-AXIS (Log Scale)
            # ==========================================
            if i != j:
                ax.set_yscale('log')
                ax.set_ylim(0.1, 100)
                ax.set_yticks(major_ticks)
                ax.set_yticklabels(major_tick_labels) # Set labels globally
                ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)))
                ax.yaxis.set_minor_formatter(ticker.NullFormatter())

                # Toggle text visibility based on column position (left column only)
                if j == 0 or j == i+1:
                    ax.tick_params(axis='y', which='both', left=True, labelleft=True)
                else:
                    ax.tick_params(axis='y', which='both', left=True, labelleft=False)
            else:
                # Diagonal: Hide axis ticks and labels
                ax.tick_params(axis='y', which='both', left=False, labelleft=False)

# --- Stitch Horizontal Colorbar Legends ---
    # Remove the right margin adjustment; rely on savefig's bbox_inches="tight" 
    # to capture the axes below the plot.
    g.fig.subplots_adjust(wspace=0.4)

    # 1. ELECTRON COLORBAR (Bottom Left)
    # [left, bottom, width, height]
    cbar_ax_elec = g.fig.add_axes([0.1, -0.05, 0.35, 0.03])
    orig_cmap = plt.cm.Reds
    # Extract colors from the 0.1 to 1.0 range
    truncated_colors = orig_cmap(np.linspace(0.1, 1.0, 256)) 
    cmap_elec = mcolors.ListedColormap(truncated_colors)
    norm_elec = mcolors.LogNorm(vmin=min(e_min, p_min), vmax=max(e_max, p_max))
    
    cbar_elec = g.fig.colorbar(
        plt.cm.ScalarMappable(norm=norm_elec, cmap=cmap_elec), 
        cax=cbar_ax_elec, 
        orientation='horizontal'
    )
    
    # Text will automatically center below or above the horizontal bar
    cbar_elec.set_label("Electron Fraction", labelpad=6, fontweight='bold', fontsize=22)
    cbar_elec.ax.tick_params(labelsize=20)
    
    cbar_elec.locator = ticker.LogLocator(base=10.0, numticks=5)
    cbar_elec.formatter = ticker.LogFormatterMathtext()
    cbar_elec.update_ticks()

    # 2. PROTON COLORBAR (Bottom Right)
    cbar_ax_prot = g.fig.add_axes([0.55, -0.05, 0.35, 0.03])
    orig_cmap = plt.cm.Blues
    # Extract colors from the 0.1 to 1.0 range
    truncated_colors = orig_cmap(np.linspace(0.1, 1.0, 256)) 
    cmap_prot = mcolors.ListedColormap(truncated_colors)
    norm_prot = mcolors.LogNorm(vmin=min(e_min, p_min), vmax=max(e_max, p_max))
    
    cbar_prot = g.fig.colorbar(
        plt.cm.ScalarMappable(norm=norm_prot, cmap=cmap_prot), 
        cax=cbar_ax_prot, 
        orientation='horizontal'
    )
    
    cbar_prot.set_label("Proton Fraction", labelpad=6, fontweight='bold', fontsize=22)
    cbar_prot.ax.tick_params(labelsize=20)
    
    cbar_prot.locator = ticker.LogLocator(base=10.0, numticks=5)
    cbar_prot.formatter = ticker.LogFormatterMathtext()
    cbar_prot.update_ticks()

    print(f"Saving high-res image to {output}...")
    g.fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(g.fig)

    # --- RESET GLOBAL PLOTTING DEFAULTS ---
    plt.rcdefaults()
    sns.reset_defaults()

    print("Export complete.")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker


def get_bin_edges(nbins=100):
    """Generates logspace edges from exactly 0.1 to 100 MeV."""
    return np.logspace(np.log10(0.1), np.log10(100.0), nbins + 1)


def get_density_limits(data, particle_type, edges, features):
    """Finds the true physical minimum (1 particle) and global maximum fraction."""
    subset = data[data['Particle_Type'] == particle_type][features]
    total_particles = len(subset)
    
    if total_particles == 0:
        return 1e-4, 1.0

    min_density = 1.0 / total_particles
    max_density = 0.0

    for i in range(len(features) - 1):
        for j in range(i + 1, len(features)):
            H, _, _ = np.histogram2d(subset.iloc[:, i], subset.iloc[:, j], bins=(edges, edges))
            # Divide by total_particles of this species (matches 2D panels)
            frac = H.max() / total_particles
            if frac > max_density:
                max_density = frac

    max_density = np.ceil(max_density * 2000) / 2000 if max_density > 0 else 1.0
    if max_density <= min_density:
        max_density = min_density * 10.0

    return min_density, max_density


def _make_truncated_cmap(base_cmap):
    """Extracts the 0.1 to 1.0 range of a colormap so low counts remain visible."""
    return mcolors.ListedColormap(base_cmap(np.linspace(0.1, 1.0, 256)))


# ==========================================
# Seaborn PairGrid Mapping Functions
# ==========================================

def plot_density_2d(x, y, particle_type, base_cmap, **kwargs):
    """Unified 2D histogram mapper for both upper (Proton) and lower (Electron) triangles."""
    ax = plt.gca()
    df = kwargs['full_data']
    edges = kwargs['edges']

    subset = df[df['Particle_Type'] == particle_type]
    total_species = len(subset)

    H, xedges, yedges = np.histogram2d(subset[x.name], subset[y.name], bins=(edges, edges))
    H_frac = H.T / max(1, total_species)
    H_masked = np.ma.masked_where(H_frac == 0, H_frac)

    cmap = _make_truncated_cmap(base_cmap)
    norm = mcolors.LogNorm(vmin=kwargs['global_vmin'], vmax=kwargs['global_vmax'])

    ax.pcolormesh(xedges, yedges, H_masked, cmap=cmap, norm=norm)


def plot_density_lower(x, y, **kwargs):
    """Lower triangle: Electron 2D density (Red)"""
    plot_density_2d(x, y, particle_type='Electron', base_cmap=plt.cm.Reds, **kwargs)


def plot_density_upper(x, y, **kwargs):
    """Upper triangle: Proton 2D density (Blue)"""
    plot_density_2d(x, y, particle_type='Proton', base_cmap=plt.cm.Blues, **kwargs)


def plot_diag_1d(x, **kwargs):
    """Unified 1D diagonal PDF plotter for both full and reduced grids."""
    ax = plt.gca()
    df = kwargs['full_data']
    edges = kwargs['edges']
    style_cfg = kwargs['style_cfg']

    col_name = x.name
    e_data = df[df['Particle_Type'] == 'Electron'][col_name]
    p_data = df[df['Particle_Type'] == 'Proton'][col_name]

    e_hist, _ = np.histogram(e_data, bins=edges)
    p_hist, _ = np.histogram(p_data, bins=edges)

    # Normalize by total particles of each species in the dataset
    e_ratio = e_hist / max(1, len(e_data))
    p_ratio = p_hist / max(1, len(p_data))

    bin_widths = np.diff(edges)
    e_density = e_ratio / bin_widths
    p_density = p_ratio / bin_widths
    midpoints = np.sqrt(edges[:-1] * edges[1:])

    # 1. Base Axis (Log Scale to satisfy Seaborn's grid sharing)
    ax.set_xscale('log')
    ax.set_xlim(0.1, 100)
    ax.set_yscale('log')
    ax.set_ylim(0.1, 100)
    ax.yaxis.set_visible(False)

    # 2. Twin Axis (Linear Scale for PDF)
    ax2 = ax.twinx()
    ax2.plot(midpoints, e_density, color=plt.cm.Reds(0.8), lw=2)
    ax2.plot(midpoints, p_density, color=plt.cm.Blues(0.8), lw=2)
    ax2.set_yscale('linear')

    max_density = max(np.max(e_density), np.max(p_density))
    ax2.set_ylim(0, max_density * 1.1 if max_density > 0 else 1.0)

    ax2.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator(n=4))
    ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f"{val:.1f}".rstrip('0').rstrip('.')))

    ax2.yaxis.tick_left()
    ax2.yaxis.set_label_position("left")
    ax2.tick_params(axis='y', which='both', left=True, labelleft=True)

    # 3. Dynamic PDF Label based on full vs. reduced preset
    tx, ty = style_cfg['pdf_text_pos']
    ax2.text(tx, ty, "PDF (MeV$^{-1}$)", transform=ax2.transAxes,
             ha='right', va='top', fontsize=style_cfg['pdf_fontsize'], style='italic', color='gray')


# ==========================================
# Unified Main Generation Function
# ==========================================

def create_and_save_pairgrid(test_data_plot, features_only, output="square_pairs_plot.png", reduced=None):
    """
    Generates and saves the 2D/1D detector PairGrid.
    
    Parameters:
    - reduced (bool or None): If True, uses compact font/spacing presets for smaller grids.
                              If None, automatically sets reduced=True when len(features_only) < 8.
    """
    if reduced is None:
        reduced = len(features_only) < 8

    # Style presets encapsulating all differences between full (8x8) and reduced grids
    if reduced:
        style_cfg = {
            'font_size': 22, 'axes_labelsize': 22, 'tick_labelsize': 20,
            'sup_fontsize': 22, 'supx_y': 0.005, 'supy_x': -0.006,
            'row_col_labelsize': 20, 'ylabel_coords': (-0.35, 0.45), 'xlabel_pad': 8,
            'wspace': 0.4, 'cbar_labelsize': 22, 'cbar_ticksize': 20, 'cbar_pad': 6,
            'pdf_text_pos': (0.80, 1.18), 'pdf_fontsize': 18
        }
    else:
        style_cfg = {
            'font_size': 36, 'axes_labelsize': 48, 'tick_labelsize': 28,
            'sup_fontsize': 40, 'supx_y': -0.01, 'supy_x': -0.01,
            'row_col_labelsize': 40, 'ylabel_coords': (-0.45, 0.45), 'xlabel_pad': 20,
            'wspace': 0.5, 'cbar_labelsize': 40, 'cbar_ticksize': 32, 'cbar_pad': 10,
            'pdf_text_pos': (0.95, 1.25), 'pdf_fontsize': 24
        }

    label_mapping = {
        "Detector1": "D1", "Detector2": "D2", "Detector3": "D3",
        "Detector4": "D4", "Detector5": "D5", "Detector6": "D6",
        "Detector7_8_sum": "D7+8", "Detector9": "D9"
    }

    plot_df = test_data_plot.rename(columns=label_mapping)
    display_features = [label_mapping.get(f, f) for f in features_only]

    print("Calculating bin edges and max densities...")
    edges = get_bin_edges()
    e_min, e_max = get_density_limits(plot_df, "Electron", edges, display_features)
    p_min, p_max = get_density_limits(plot_df, "Proton", edges, display_features)
    global_vmin = min(e_min, p_min)
    global_vmax = max(e_max, p_max)

    print(f"Generating PairGrid mapping (reduced={reduced})...")
    sns.set_theme(style="ticks")

    plt.rcParams.update({
        'font.size': style_cfg['font_size'],
        'axes.labelsize': style_cfg['axes_labelsize'],
        'axes.labelweight': 'bold',
        'xtick.labelsize': style_cfg['tick_labelsize'],
        'ytick.labelsize': style_cfg['tick_labelsize'],
        'xtick.major.size': 12, 'xtick.major.width': 3,
        'xtick.minor.size': 6,  'xtick.minor.width': 2,
        'ytick.major.size': 12, 'ytick.major.width': 3,
        'ytick.minor.size': 6,  'ytick.minor.width': 2
    })

    g = sns.PairGrid(plot_df, vars=display_features, height=3, aspect=1, diag_sharey=False)

    map_kwargs = {
        'full_data': plot_df,
        'edges': edges,
        'global_vmin': global_vmin,
        'global_vmax': global_vmax,
        'style_cfg': style_cfg
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

    g.fig.supxlabel("Deposited Energy (MeV)", fontsize=style_cfg['sup_fontsize'], fontweight='bold', y=style_cfg['supx_y'])
    g.fig.supylabel("Deposited Energy (MeV)", fontsize=style_cfg['sup_fontsize'], fontweight='bold', x=style_cfg['supy_x'])

    major_ticks = [0.1, 1, 10, 100]
    major_tick_labels = ["0.1", "1", "10", "100"]
    n_axes = len(g.axes)

    for i in range(n_axes):
        ax_left = g.axes[i, 0]
        ax_left.yaxis.label.set_size(style_cfg['row_col_labelsize'])
        ax_left.yaxis.label.set_weight('bold')
        ax_left.yaxis.set_label_coords(*style_cfg['ylabel_coords'])

        ax_bottom = g.axes[-1, i]
        ax_bottom.xaxis.label.set_size(style_cfg['row_col_labelsize'])
        ax_bottom.xaxis.label.set_weight('bold')
        ax_bottom.xaxis.labelpad = style_cfg['xlabel_pad']

        for j in range(n_axes):
            ax = g.axes[i, j]
            if ax is None:
                continue

            # X-Axis (Log Scale)
            ax.set_xscale('log')
            ax.set_xlim(0.1, 100)
            ax.set_xticks(major_ticks)
            ax.set_xticklabels(major_tick_labels)
            ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)))
            ax.xaxis.set_minor_formatter(ticker.NullFormatter())
            ax.tick_params(axis='x', which='both', bottom=True, labelbottom=(i == n_axes - 1))

            # Y-Axis (Log Scale on off-diagonals)
            if i != j:
                ax.set_yscale('log')
                ax.set_ylim(0.1, 100)
                ax.set_yticks(major_ticks)
                ax.set_yticklabels(major_tick_labels)
                ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)))
                ax.yaxis.set_minor_formatter(ticker.NullFormatter())
                ax.tick_params(axis='y', which='both', left=True, labelleft=(j == 0 or j == i + 1))
            else:
                ax.tick_params(axis='y', which='both', left=False, labelleft=False)

    # --- Stitch Horizontal Colorbar Legends ---
    g.fig.subplots_adjust(wspace=style_cfg['wspace'])
    shared_norm = mcolors.LogNorm(vmin=global_vmin, vmax=global_vmax)

    cbar_specs = [
        ([0.10, -0.05, 0.35, 0.03], plt.cm.Reds,  "Electron Fraction"),
        ([0.55, -0.05, 0.35, 0.03], plt.cm.Blues, "Proton Fraction")
    ]

    for rect, base_cmap, label in cbar_specs:
        cbar_ax = g.fig.add_axes(rect)
        cmap = _make_truncated_cmap(base_cmap)
        cbar = g.fig.colorbar(plt.cm.ScalarMappable(norm=shared_norm, cmap=cmap), cax=cbar_ax, orientation='horizontal')
        cbar.set_label(label, labelpad=style_cfg['cbar_pad'], fontweight='bold', fontsize=style_cfg['cbar_labelsize'])
        cbar.ax.tick_params(labelsize=style_cfg['cbar_ticksize'])
        cbar.locator = ticker.LogLocator(base=10.0, numticks=5)
        cbar.formatter = ticker.LogFormatterMathtext()
        cbar.update_ticks()

    print(f"Saving high-res image to {output}...")
    g.fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(g.fig)

    plt.rcdefaults()
    sns.reset_defaults()
    print("Export complete.")

def plot_species_error_distribution(misclassified_df, species, e_inc_bins, e_inc_xlim, e_inc_ticks,
                                    features, base_color, dark_color, title_suffix=""):

    subset = misclassified_df[misclassified_df['Particle_Type'] == species]
    total_mc = len(subset)

    max_edep = subset[features].max().max() if total_mc > 0 else 1.0
    d_bins = np.linspace(0, max_edep, 101)

    fig, axes = plt.subplots(3, 3, figsize=(20, 18))
    fig.suptitle(f'Distribution of Misclassified {species}s {title_suffix}', fontsize=40, fontweight='bold', y=0.97)

    for i, feature in enumerate(features + ['E_Inc']):
        ax = axes.flat[i]

        if feature != 'E_Inc':
            mc_counts, _ = np.histogram(subset[feature], bins=d_bins)
            prob = (mc_counts / total_mc) * 100 if total_mc > 0 else np.zeros_like(mc_counts, dtype=float)

            ax.hist(d_bins[:-1], bins=d_bins, weights=prob, color=base_color, alpha=0.7)
            ax.set_xscale('linear')
            ax.set_xlim(0, max_edep)
            ax.set_xlabel('Deposited Energy (MeV)', fontsize=24)

        else:
            mc_counts, _ = np.histogram(subset['E_Inc'], bins=e_inc_bins)
            prob = (mc_counts / total_mc) * 100 if total_mc > 0 else np.zeros_like(mc_counts, dtype=float)

            ax.hist(e_inc_bins[:-1], bins=e_inc_bins, weights=prob, color=dark_color, alpha=0.9)
            ax.set_xscale('log')
            ax.set_xlim(*e_inc_xlim)
            ax.set_xticks(e_inc_ticks, labels=[str(v) for v in e_inc_ticks])
            ax.set_xlabel('Incident Energy (MeV)', fontsize=24)

        ax.set_title(feature, fontsize=28, fontweight='bold')
        ax.set_ylabel(f'% of Total {species} Errors', fontsize=24)
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.grid(True, which='major', linestyle='--', alpha=0.7)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    plt.show() 
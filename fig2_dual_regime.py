"""
Figure 2: Dual-Regime Spatial Analysis of Three BCI Tree Species
With Bootstrap Confidence Intervals

This script compares the Structure Factor S(k) profiles of three species
with contrasting ecological strategies:
- Heisteria concinna (balanced/liquid-like)
- Guarea sebifera (random-like)
- Gustavia superba (clustered/habitat specialist)

Author: Talha Aksoy
Repository: https://github.com/tlhksy/bci-spatial-analysis
Data Source: ForestGEO (https://doi.org/10.15146/5xcp-0d46)
"""

import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.special import j0
import matplotlib.pyplot as plt


def calculate_structure_factor(points, k_values):
    """
    Calculate the radially-averaged Structure Factor S(k) for a 2D point pattern.
    
    Parameters
    ----------
    points : ndarray
        Nx2 array of (x, y) coordinates
    k_values : ndarray
        Array of wavenumbers to evaluate
        
    Returns
    -------
    ndarray
        S(k) values for each wavenumber
    """
    n = len(points)
    
    dist_matrix = squareform(pdist(points))
    np.fill_diagonal(dist_matrix, 0)
    
    s_values = []
    for k in k_values:
        bessel_sum = np.sum(j0(k * dist_matrix)) - n
        s_k = 1 + (1/n) * bessel_sum
        s_values.append(s_k)
    
    return np.array(s_values)


def bootstrap_full_profile(points, k_values, n_bootstrap=50):
    """
    Calculate bootstrap confidence interval for full S(k) profile.
    
    Parameters
    ----------
    points : ndarray
        Nx2 array of coordinates
    k_values : ndarray
        Array of wavenumbers
    n_bootstrap : int
        Number of bootstrap samples
        
    Returns
    -------
    tuple
        (s_observed, lower_ci, upper_ci)
    """
    n = len(points)
    s_observed = calculate_structure_factor(points, k_values)
    
    bootstrap_results = []
    for i in range(n_bootstrap):
        indices = np.random.choice(n, size=n, replace=True)
        resampled_points = points[indices]
        resampled_points = np.unique(resampled_points, axis=0)
        
        if len(resampled_points) < 50:
            continue
            
        s_boot = calculate_structure_factor(resampled_points, k_values)
        bootstrap_results.append(s_boot)
    
    bootstrap_results = np.array(bootstrap_results)
    lower_ci = np.percentile(bootstrap_results, 2.5, axis=0)
    upper_ci = np.percentile(bootstrap_results, 97.5, axis=0)
    
    return s_observed, lower_ci, upper_ci


def load_species_data(filepath, species_list, census=8, dbh_min=100, n_bootstrap=50):
    """
    Load and process data for multiple species with bootstrap CI.
    
    Parameters
    ----------
    filepath : str
        Path to the TSV data file
    species_list : list
        List of species names to analyze
    census : int
        Census number (default: 8 for 2015)
    dbh_min : float
        Minimum DBH in mm
    n_bootstrap : int
        Number of bootstrap samples
        
    Returns
    -------
    dict
        Dictionary with species data and S(k) profiles
    """
    df = pd.read_csv(filepath, sep='\t', low_memory=False)
    
    # Define wavenumber range
    k_values = np.linspace(0.02, 4.0, 150)
    
    species_data = {}
    
    for species in species_list:
        print(f"\nProcessing {species}...")
        
        # Filter data
        mask = (
            (df['SpeciesName'] == species) &
            (df['PlotCensusNumber'] == census) &
            (df['Status'].str.lower() == 'alive')
        )
        sub = df[mask].copy()
        
        # Filter adults
        sub['DBH'] = pd.to_numeric(sub['DBH'], errors='coerce')
        sub = sub[sub['DBH'] >= dbh_min].dropna(subset=['PX', 'PY'])
        
        if len(sub) < 50:
            print(f"  Warning: {species} has insufficient data (N={len(sub)})")
            continue
        
        # Extract and clean coordinates
        points = np.column_stack((sub['PX'].values, sub['PY'].values))
        points = np.unique(np.round(points * 2) / 2, axis=0)
        n = len(points)
        
        print(f"  N = {n}, calculating bootstrap CI...")
        
        # Calculate S(k) with bootstrap CI
        s_obs, s_lower, s_upper = bootstrap_full_profile(points, k_values, n_bootstrap)
        
        # Calculate key metrics
        min_sk = np.min(s_obs)
        low_k_max = np.max(s_obs[:50])  # k < ~0.5 region
        
        species_data[species] = {
            'n': n,
            'points': points,
            's_values': s_obs,
            's_lower': s_lower,
            's_upper': s_upper,
            'k_values': k_values,
            'min_sk': min_sk,
            'low_k_max': low_k_max
        }
        
        print(f"  Min S(k)={min_sk:.4f} | Low-k Max={low_k_max:.2f}")
    
    return species_data


def plot_dual_regime(species_data, output_path='fig2_dual_regime.png'):
    """
    Create Figure 2: Dual-regime comparison plot with bootstrap CI.
    
    Parameters
    ----------
    species_data : dict
        Output from load_species_data()
    output_path : str
        Output file path
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Style definitions
    styles = {
        'concinna': {
            'color': 'blue', 
            'ls': '-', 
            'lw': 2.5,
            'label': 'H. concinna'
        },
        'sebifera': {
            'color': 'gray', 
            'ls': '--', 
            'lw': 2.0,
            'label': 'G. sebifera'
        },
        'superba': {
            'color': 'red', 
            'ls': '-.', 
            'lw': 2.5,
            'label': 'G. superba'
        }
    }
    
    # Plot each species on both panels
    for sp, data in species_data.items():
        if sp not in styles:
            continue
        
        k_vals = data['k_values']
        s_vals = data['s_values']
        s_lower = data['s_lower']
        s_upper = data['s_upper']
        n = data['n']
        
        label = f"{styles[sp]['label']} (N={n})"
        
        # Panel A with CI band
        ax1.fill_between(k_vals, s_lower, s_upper,
                         alpha=0.2, color=styles[sp]['color'])
        ax1.plot(k_vals, s_vals, 
                 color=styles[sp]['color'],
                 linestyle=styles[sp]['ls'], 
                 linewidth=styles[sp]['lw'],
                 label=label)
        
        # Panel B with CI band
        ax2.fill_between(k_vals, s_lower, s_upper,
                         alpha=0.2, color=styles[sp]['color'])
        ax2.plot(k_vals, s_vals, 
                 color=styles[sp]['color'],
                 linestyle=styles[sp]['ls'], 
                 linewidth=styles[sp]['lw'],
                 label=label)
    
    # =========================================================================
    # Panel A: Large-Scale Ecological Domain (Log Scale)
    # =========================================================================
    ax1.axhline(1.0, color='black', linewidth=1, linestyle='-', alpha=0.7)
    ax1.set_xlim(0, 1.0)
    ax1.set_yscale('log')
    ax1.set_ylim(0.5, 500)
    ax1.set_xlabel('Wavevector k (m⁻¹)', fontsize=12)
    ax1.set_ylabel('S(k) [Log Scale]', fontsize=12)
    ax1.set_title('(a) Large-Scale Ecological Domain (k < 1.0)\n'
                  'with 95% Bootstrap CI', fontsize=11)
    ax1.legend(loc='upper right', fontsize=9)
    ax1.grid(True, which='both', alpha=0.3)
    
    # Annotation for extreme clustering
    if 'superba' in species_data:
        ax1.annotate('Extreme\nClustering', xy=(0.05, 150), xytext=(0.3, 100),
                     fontsize=9, color='red',
                     arrowprops=dict(arrowstyle='->', color='red', lw=1.5))
    
    # =========================================================================
    # Panel B: Short-Range Physical Domain (Linear Scale)
    # =========================================================================
    ax2.axhline(1.0, color='black', linewidth=1.5, linestyle='-', 
                label='Poisson Baseline')
    ax2.set_xlim(0.5, 4.0)
    ax2.set_ylim(0, 2.0)
    ax2.set_xlabel('Wavevector k (m⁻¹)', fontsize=12)
    ax2.set_ylabel('S(k) [Linear Scale]', fontsize=12)
    ax2.set_title('(b) Short-Range Physical Domain (k > 0.5)\n'
                  'with 95% Bootstrap CI', fontsize=11)
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, alpha=0.3)
    
    # Universal convergence annotation
    ax2.annotate('Universal\nConvergence\n(S(k)≈0.8-1.0)', 
                 xy=(2.5, 0.9), fontsize=9,
                 bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    # Main title
    plt.suptitle('Dual-Regime Spatial Analysis - BCI Census 8 (2015)', 
                 fontsize=13, y=1.02, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✅ Figure saved: {output_path}")
    
    return fig, (ax1, ax2)


def print_summary_table(species_data):
    """Print a summary comparison table."""
    print("\n" + "=" * 70)
    print("SPECIES COMPARISON SUMMARY")
    print("=" * 70)
    print(f"{'Species':<15} {'N':<8} {'Min S(k)':<12} {'Low-k Max':<12} {'Strategy':<15}")
    print("-" * 70)
    
    strategies = {
        'concinna': 'Liquid-like',
        'sebifera': 'Random-like', 
        'superba': 'Clustered'
    }
    
    for sp in ['concinna', 'sebifera', 'superba']:
        if sp in species_data:
            data = species_data[sp]
            print(f"{sp:<15} {data['n']:<8} {data['min_sk']:<12.4f} "
                  f"{data['low_k_max']:<12.2f} {strategies[sp]:<15}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    import os
    
    # Data file options (in order of preference)
    DATA_FILES = [
        "BCI_3species.tsv",           # Filtered data (recommended)
        "FullMeasurementBCI.tsv",     # Full dataset
    ]
    
    # Find available data file
    DATA_FILE = None
    for f in DATA_FILES:
        if os.path.exists(f):
            DATA_FILE = f
            break
    
    if DATA_FILE is None:
        print("ERROR: Data file not found!")
        print("\nPlease ensure one of these files is in the current directory:")
        print("  1. BCI_3species.tsv (filtered, ~4 MB)")
        print("  2. FullMeasurementBCI.tsv (full dataset, ~300 MB)")
        print("\nDownload from: https://doi.org/10.15146/5xcp-0d46")
        print("Or get filtered data from: https://github.com/tlhksy/bci-spatial-analysis")
        exit(1)
    
    print("=" * 60)
    print("Figure 2: Dual-Regime Spatial Analysis")
    print(f"Data file: {DATA_FILE}")
    print("=" * 60)
    
    # Species to analyze
    SPECIES_LIST = ['concinna', 'sebifera', 'superba']
    
    # Load and analyze data with bootstrap
    species_data = load_species_data(DATA_FILE, SPECIES_LIST, census=8, n_bootstrap=50)
    
    # Create figure
    plot_dual_regime(species_data, 'fig2_dual_regime.png')
    
    # Print summary
    print_summary_table(species_data)
    
    print("\n" + "=" * 60)
    print("KEY FINDINGS")
    print("=" * 60)
    print("""
    1. PHYSICAL REGIME (k > 0.5): All species converge to S(k) ≈ 0.8-1.0
       → Universal constraint from crown/root exclusion
       
    2. ECOLOGICAL REGIME (k < 0.5): Species-specific strategies diverge
       → G. superba: Extreme clustering (habitat specialist)
       → G. sebifera: Near-random (generalist)
       → H. concinna: Balanced intermediate (liquid-like)
    """)

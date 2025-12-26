"""
Figure 1: Temporal Evolution of Spatial Order in H. concinna (1982-2015)
With Bootstrap Confidence Intervals

This script calculates the Structure Factor S(k) for Heisteria concinna
across 8 census intervals (1982-2015) from the BCI 50-ha plot dataset.

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
    
    The Structure Factor quantifies density fluctuations at different spatial scales.
    S(k) = 1 for complete spatial randomness (Poisson process)
    S(k) < 1 indicates repulsive order ("liquid-like")
    S(k) > 1 indicates clustering
    
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
    
    # Calculate pairwise distance matrix
    dist_matrix = squareform(pdist(points))
    np.fill_diagonal(dist_matrix, 0)
    
    s_values = []
    for k in k_values:
        # S(k) = 1 + (2/N) * sum_{i<j} J_0(k * r_ij)
        # Using full matrix: S(k) = 1 + (1/N) * (sum_all J_0(k*r) - N)
        bessel_sum = np.sum(j0(k * dist_matrix)) - n
        s_k = 1 + (1/n) * bessel_sum
        s_values.append(s_k)
    
    return np.array(s_values)


def bootstrap_min_sk(points, k_values, n_bootstrap=100):
    """
    Calculate bootstrap confidence interval for minimum S(k).
    
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
        (observed_min_sk, lower_ci, upper_ci)
    """
    n = len(points)
    
    # Observed S(k) and its minimum
    s_observed = calculate_structure_factor(points, k_values)
    min_sk_observed = np.min(s_observed)
    
    # Bootstrap resampling for min S(k)
    min_sk_boots = []
    for i in range(n_bootstrap):
        # Resample with replacement
        indices = np.random.choice(n, size=n, replace=True)
        resampled_points = points[indices]
        
        # Remove duplicate points
        resampled_points = np.unique(resampled_points, axis=0)
        
        if len(resampled_points) < 50:
            continue
            
        s_boot = calculate_structure_factor(resampled_points, k_values)
        min_sk_boots.append(np.min(s_boot))
    
    min_sk_boots = np.array(min_sk_boots)
    
    # 95% CI
    lower_ci = np.percentile(min_sk_boots, 2.5)
    upper_ci = np.percentile(min_sk_boots, 97.5)
    
    return min_sk_observed, lower_ci, upper_ci


def load_and_filter_data(filepath, species, census, dbh_min=100):
    """
    Load BCI data and filter for adult individuals of a specific species.
    
    Parameters
    ----------
    filepath : str
        Path to the TSV data file
    species : str
        Species name to filter (partial match)
    census : int
        Census number (1-8)
    dbh_min : float
        Minimum DBH in mm (default: 100mm = 10cm for adults)
        
    Returns
    -------
    ndarray
        Nx2 array of (x, y) coordinates
    int
        Number of individuals
    """
    df = pd.read_csv(filepath, sep='\t', low_memory=False)
    
    # Filter by species, census, and alive status
    mask = (
        df['SpeciesName'].str.contains(species, case=False, na=False) &
        (df['PlotCensusNumber'] == census) &
        (df['Status'].str.lower() == 'alive')
    )
    sub = df[mask].copy()
    
    # Filter by DBH (adults only)
    sub['DBH'] = pd.to_numeric(sub['DBH'], errors='coerce')
    sub = sub[sub['DBH'] >= dbh_min].dropna(subset=['PX', 'PY'])
    
    if len(sub) < 10:
        return None, 0
    
    # Extract coordinates
    points = np.column_stack((sub['PX'].values, sub['PY'].values))
    
    # Spatial consolidation: merge points within 0.5m (measurement artifacts)
    points_clean = np.unique(np.round(points * 2) / 2, axis=0)
    
    return points_clean, len(points_clean)


def analyze_temporal_evolution(filepath, species='concinna', n_bootstrap=50):
    """
    Analyze the temporal evolution of spatial order across all censuses.
    
    Parameters
    ----------
    filepath : str
        Path to the BCI data file
    species : str
        Target species name
    n_bootstrap : int
        Number of bootstrap samples for CI
        
    Returns
    -------
    list
        List of dictionaries with census results
    """
    # Define wavenumber range
    k_values = np.linspace(0.02, 4.0, 150)
    years = [1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015]
    
    results = []
    
    for census in range(1, 9):
        points, n = load_and_filter_data(filepath, species, census)
        
        if points is None:
            print(f"Census {census}: Insufficient data")
            continue
        
        # Calculate S(k) with bootstrap CI
        min_sk, lower, upper = bootstrap_min_sk(points, k_values, n_bootstrap)
        
        results.append({
            'census': census,
            'year': years[census-1],
            'n': n,
            'min_sk': min_sk,
            'lower_ci': lower,
            'upper_ci': upper
        })
        
        print(f"Census {census} ({years[census-1]}): N={n:3d}, "
              f"Min S(k)={min_sk:.4f} [{lower:.4f}, {upper:.4f}]")
    
    return results


def plot_temporal_evolution(results, output_path='fig1_temporal_evolution.png'):
    """
    Create Figure 1: Temporal evolution plot with bootstrap CI.
    
    Parameters
    ----------
    results : list
        Results from analyze_temporal_evolution()
    output_path : str
        Output file path for the figure
    """
    fig, ax = plt.subplots(figsize=(11, 6))
    
    # Extract data
    census_nums = [r['census'] for r in results]
    min_sks = [r['min_sk'] for r in results]
    lowers = [r['lower_ci'] for r in results]
    uppers = [r['upper_ci'] for r in results]
    ns = [r['n'] for r in results]
    years = [r['year'] for r in results]
    
    # Error bars - ensure non-negative
    err_lower = [max(0, m - l) for m, l in zip(min_sks, lowers)]
    err_upper = [max(0, u - m) for m, u in zip(min_sks, uppers)]
    
    # Random baseline (Poisson)
    ax.axhline(1.0, color='red', linestyle='--', linewidth=2, 
               label='Random (Poisson)', zorder=1)
    
    # Reference line at S(k) = 0.85
    ax.axhline(0.85, color='gray', linestyle=':', alpha=0.7, 
               label='S(k)=0.85 Reference')
    
    # Observed values with error bars
    ax.errorbar(census_nums, min_sks, yerr=[err_lower, err_upper],
                fmt='o-', color='blue', linewidth=2.5, markersize=12,
                capsize=5, capthick=2, elinewidth=2,
                label='H. concinna (95% Bootstrap CI)', zorder=5)
    
    # Annotate with sample sizes
    for c, s, n, u in zip(census_nums, min_sks, ns, uppers):
        ax.annotate(f'N={n}', (c, u + 0.03), ha='center', fontsize=8, color='blue')
    
    # Labels and formatting
    ax.set_xlabel('Census Number', fontsize=12)
    ax.set_ylabel('Minimum Structure Factor S(k)', fontsize=12)
    ax.set_title('Temporal Evolution of Spatial Order in H. concinna (1982-2015)\n'
                 'with 95% Bootstrap Confidence Intervals', fontsize=13)
    ax.set_xlim(0.5, 8.5)
    ax.set_ylim(0.45, 1.25)
    ax.set_xticks(census_nums)
    ax.set_xticklabels([f'{c}\n({y})' for c, y in zip(census_nums, years)])
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    
    # Summary statistics box
    textstr = f'Range: {min(min_sks):.2f} - {max(min_sks):.2f}\nMean: {np.mean(min_sks):.2f}'
    props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✅ Figure saved: {output_path}")
    
    return fig, ax


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
    print("Figure 1: Temporal Evolution of Spatial Order")
    print(f"Data file: {DATA_FILE}")
    print("=" * 60)
    
    # Run analysis with bootstrap
    results = analyze_temporal_evolution(DATA_FILE, species='concinna', n_bootstrap=50)
    
    # Create figure
    plot_temporal_evolution(results, 'fig1_temporal_evolution.png')
    
    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    min_sks = [r['min_sk'] for r in results]
    print(f"S(k) Range: {min(min_sks):.2f} - {max(min_sks):.2f}")
    print(f"Mean S(k): {np.mean(min_sks):.2f}")
    print(f"Population: {results[0]['n']} → {results[-1]['n']}")

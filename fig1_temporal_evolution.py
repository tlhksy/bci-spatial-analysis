"""
Figure 1: Temporal Evolution of Spatial Order in H. concinna (1982-2015)

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


def analyze_temporal_evolution(filepath, species='concinna'):
    """
    Analyze the temporal evolution of spatial order across all censuses.
    
    Parameters
    ----------
    filepath : str
        Path to the BCI data file
    species : str
        Target species name
        
    Returns
    -------
    list
        List of dictionaries with census results
    """
    # Define wavenumber range
    # k_min based on plot size (~1000m), k_max for short-range structure
    k_values = np.linspace(0.02, 4.0, 250)
    
    results = []
    
    for census in range(1, 9):
        points, n = load_and_filter_data(filepath, species, census)
        
        if points is None:
            print(f"Census {census}: Insufficient data")
            continue
        
        # Calculate S(k)
        s_values = calculate_structure_factor(points, k_values)
        min_sk = np.min(s_values)
        min_k = k_values[np.argmin(s_values)]
        
        results.append({
            'census': census,
            'n': n,
            'min_sk': min_sk,
            'min_k': min_k,
            's_values': s_values,
            'k_values': k_values
        })
        
        print(f"Census {census}: N={n:3d}, Min S(k)={min_sk:.4f} at k={min_k:.2f}")
    
    return results


def plot_temporal_evolution(results, output_path='fig1_temporal_evolution.png'):
    """
    Create Figure 1: Temporal evolution plot.
    
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
    ns = [r['n'] for r in results]
    years = [1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015]
    
    # Random baseline (Poisson)
    ax.axhline(1.0, color='red', linestyle='--', linewidth=2, 
               label='Random (Poisson)', zorder=1)
    
    # Reference line at S(k) = 0.85
    ax.axhline(0.85, color='gray', linestyle=':', alpha=0.7, 
               label='S(k)=0.85 Reference')
    
    # Observed values
    ax.plot(census_nums, min_sks, 'o-', color='blue', linewidth=2.5, 
            markersize=12, label='Observed H. concinna', zorder=5)
    
    # Annotate with sample sizes
    for c, s, n in zip(census_nums, min_sks, ns):
        ax.annotate(f'N={n}', (c, s - 0.03), ha='center', fontsize=8, color='blue')
    
    # Labels and formatting
    ax.set_xlabel('Census Number', fontsize=12)
    ax.set_ylabel('Minimum Structure Factor S(k)', fontsize=12)
    ax.set_title('Temporal Evolution of Spatial Order in H. concinna (1982-2015)', 
                 fontsize=13)
    ax.set_xlim(0.5, 8.5)
    ax.set_ylim(0.55, 1.10)
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
        print("  1. BCI_3species.tsv (filtered, ~3 MB)")
        print("  2. FullMeasurementBCI.tsv (full dataset, ~300 MB)")
        print("\nDownload from: https://doi.org/10.15146/5xcp-0d46")
        print("Or get filtered data from: https://github.com/tlhksy/bci-spatial-analysis")
        exit(1)
    
    print("=" * 60)
    print("Figure 1: Temporal Evolution of Spatial Order")
    print(f"Data file: {DATA_FILE}")
    print("=" * 60)
    
    # Run analysis
    results = analyze_temporal_evolution(DATA_FILE, species='concinna')
    
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

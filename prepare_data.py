"""
Data Preparation Script for BCI Spatial Analysis

This script filters the full BCI dataset to extract three focal species:
- Heisteria concinna (Olacaceae) - balanced/liquid-like strategy
- Guarea sebifera (Meliaceae) - random-like/generalist strategy
- Gustavia superba (Lecythidaceae) - clustered/habitat specialist

Author: Talha Aksoy
Repository: https://github.com/tlhksy/bci-spatial-analysis
"""

import pandas as pd
import os


def prepare_filtered_dataset(input_path, output_path='BCI_3species.tsv'):
    """
    Filter the full BCI dataset to include only three focal species.
    
    Parameters
    ----------
    input_path : str
        Path to the full BCI dataset (FullMeasurementBCI.tsv)
    output_path : str
        Path for the filtered output file
        
    Returns
    -------
    pd.DataFrame
        Filtered dataset
    """
    print(f"Loading full dataset: {input_path}")
    df = pd.read_csv(input_path, sep='\t', low_memory=False)
    
    print(f"Full dataset: {len(df):,} rows")
    
    # Target species (using exact SpeciesName matches)
    target_species = ['concinna', 'sebifera', 'superba']
    
    # Filter for target species
    mask = df['SpeciesName'].isin(target_species)
    df_filtered = df[mask].copy()
    
    print(f"\nFiltered dataset: {len(df_filtered):,} rows")
    
    # Summary by species
    print("\nSpecies summary:")
    print("-" * 50)
    for sp in target_species:
        sp_data = df_filtered[df_filtered['SpeciesName'] == sp]
        n_total = len(sp_data)
        n_alive = len(sp_data[sp_data['Status'].str.lower() == 'alive'])
        print(f"  {sp}: {n_total:,} records ({n_alive:,} alive)")
    
    # Save filtered dataset
    df_filtered.to_csv(output_path, sep='\t', index=False)
    print(f"\n✅ Saved filtered dataset: {output_path}")
    print(f"   File size: {os.path.getsize(output_path) / 1e6:.1f} MB")
    
    return df_filtered


def verify_dataset(filepath):
    """
    Verify the filtered dataset contents.
    
    Parameters
    ----------
    filepath : str
        Path to the filtered dataset
    """
    print(f"\nVerifying dataset: {filepath}")
    print("=" * 60)
    
    df = pd.read_csv(filepath, sep='\t', low_memory=False)
    
    print(f"Total rows: {len(df):,}")
    print(f"Columns: {list(df.columns)[:10]}...")
    
    # Check species
    print("\nSpecies in dataset:")
    print(df['SpeciesName'].value_counts())
    
    # Check censuses
    print("\nCensuses available:")
    print(sorted(df['PlotCensusNumber'].unique()))
    
    # Check coordinates
    print("\nCoordinate ranges:")
    print(f"  PX: {df['PX'].min():.1f} - {df['PX'].max():.1f} m")
    print(f"  PY: {df['PY'].min():.1f} - {df['PY'].max():.1f} m")
    
    # Check DBH for adults
    df['DBH'] = pd.to_numeric(df['DBH'], errors='coerce')
    adults = df[df['DBH'] >= 100]
    print(f"\nAdult trees (DBH ≥ 10cm): {len(adults):,}")
    
    for sp in df['SpeciesName'].unique():
        sp_adults = adults[(adults['SpeciesName'] == sp) & 
                          (adults['Status'].str.lower() == 'alive') &
                          (adults['PlotCensusNumber'] == 8)]
        print(f"  {sp} (Census 8): {len(sp_adults)} adults")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # Check for full dataset
    FULL_DATASET = "FullMeasurementBCI.tsv"
    OUTPUT_FILE = "BCI_3species.tsv"
    
    if os.path.exists(FULL_DATASET):
        # Create filtered dataset from full data
        prepare_filtered_dataset(FULL_DATASET, OUTPUT_FILE)
        verify_dataset(OUTPUT_FILE)
    elif os.path.exists(OUTPUT_FILE):
        # Just verify existing filtered dataset
        print("Full dataset not found. Verifying existing filtered dataset...")
        verify_dataset(OUTPUT_FILE)
    else:
        print("ERROR: No data files found!")
        print("\nTo prepare the filtered dataset:")
        print("1. Download FullMeasurementBCI.tsv from:")
        print("   https://doi.org/10.15146/5xcp-0d46")
        print("2. Place it in this directory")
        print("3. Run this script again")
        print("\nAlternatively, use the pre-filtered BCI_3species.tsv from:")
        print("   https://github.com/tlhksy/bci-spatial-analysis")

"""
Data Preparation Script: Extract Target Species from Full BCI Dataset

This script filters the complete BCI measurement file to extract only
the three focal species used in the analysis:
- Heisteria concinna
- Guarea sebifera  
- Gustavia superba

Run this script first if you downloaded the full dataset from ForestGEO.

Author: Talha Aksoy
Repository: https://github.com/tlhksy/bci-spatial-analysis
Data Source: ForestGEO (https://doi.org/10.15146/5xcp-0d46)
"""

import pandas as pd
import os


def prepare_data(input_file="FullMeasurementBCI.tsv", output_file="BCI_3species.tsv"):
    """
    Filter the full BCI dataset to extract target species.
    
    Parameters
    ----------
    input_file : str
        Path to the full BCI measurement file (FullMeasurementBCI.tsv)
        Download from: https://datadryad.org/dataset/doi:10.15146/5xcp-0d46
    output_file : str
        Path for the filtered output file
        
    Returns
    -------
    pd.DataFrame
        Filtered dataframe with only target species
    """
    
    print("=" * 60)
    print("BCI Data Preparation")
    print("=" * 60)
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"\n❌ ERROR: Input file not found: {input_file}")
        print("\nPlease download the BCI dataset from:")
        print("https://datadryad.org/dataset/doi:10.15146/5xcp-0d46")
        print("\nDownload 'FullMeasurementBCI.zip' and extract it to this folder.")
        return None
    
    print(f"\n📂 Loading: {input_file}")
    print("   This may take a moment...")
    
    # Load data
    df = pd.read_csv(input_file, sep='\t', low_memory=False)
    print(f"   Total records: {len(df):,}")
    
    # Target species
    target_species = ['concinna', 'sebifera', 'superba']
    
    print(f"\n🔍 Filtering for target species: {target_species}")
    
    # Filter
    mask = df['SpeciesName'].isin(target_species)
    filtered = df[mask].copy()
    
    print(f"   Filtered records: {len(filtered):,}")
    
    # Show species breakdown
    print("\n📊 Species breakdown:")
    for sp in target_species:
        count = len(filtered[filtered['SpeciesName'] == sp])
        print(f"   {sp}: {count:,} records")
    
    # Save
    print(f"\n💾 Saving to: {output_file}")
    filtered.to_csv(output_file, sep='\t', index=False)
    
    # File size
    size_mb = os.path.getsize(output_file) / (1024 * 1024)
    print(f"   File size: {size_mb:.1f} MB")
    
    print("\n✅ Data preparation complete!")
    print(f"   You can now run the analysis scripts.")
    
    return filtered


def verify_data(filepath="BCI_3species.tsv"):
    """
    Verify the prepared data file and show summary statistics.
    
    Parameters
    ----------
    filepath : str
        Path to the filtered data file
    """
    
    if not os.path.exists(filepath):
        print(f"❌ File not found: {filepath}")
        print("   Run prepare_data() first.")
        return
    
    print("=" * 60)
    print("Data Verification")
    print("=" * 60)
    
    df = pd.read_csv(filepath, sep='\t', low_memory=False)
    
    print(f"\n📂 File: {filepath}")
    print(f"   Total records: {len(df):,}")
    print(f"   Columns: {len(df.columns)}")
    
    print("\n📊 Species summary:")
    print(df['SpeciesName'].value_counts())
    
    print("\n📊 Census 8 adult counts (DBH >= 100mm):")
    for sp in df['SpeciesName'].unique():
        sub = df[(df['SpeciesName'] == sp) & 
                 (df['PlotCensusNumber'] == 8) &
                 (df['Status'].str.lower() == 'alive')]
        sub['DBH'] = pd.to_numeric(sub['DBH'], errors='coerce')
        adults = sub[sub['DBH'] >= 100]
        print(f"   {sp}: N = {len(adults)}")
    
    print("\n✅ Data verification complete!")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    import sys
    
    # Default files
    INPUT_FILE = "FullMeasurementBCI.tsv"
    OUTPUT_FILE = "BCI_3species.tsv"
    
    # Check command line arguments
    if len(sys.argv) > 1:
        INPUT_FILE = sys.argv[1]
    if len(sys.argv) > 2:
        OUTPUT_FILE = sys.argv[2]
    
    # Check if output already exists
    if os.path.exists(OUTPUT_FILE):
        print(f"📂 Found existing data file: {OUTPUT_FILE}")
        verify_data(OUTPUT_FILE)
    else:
        # Prepare data
        prepare_data(INPUT_FILE, OUTPUT_FILE)
        print("\n" + "-" * 60)
        verify_data(OUTPUT_FILE)

# bci-spatial-analysis
Scale-Dependent Spatial Order in Tropical Forests
# Scale-Dependent Spatial Order in Tropical Forests

[![DOI](https://img.shields.io/badge/DOI-10.xxxx/xxxxx-blue)](https://doi.org/10.xxxx/xxxxx)
[![Data](https://img.shields.io/badge/Data-ForestGEO-green)](https://doi.org/10.15146/5xcp-0d46)
[![Python](https://img.shields.io/badge/Python-3.8+-yellow)](https://www.python.org/)

This repository contains the analysis code and reproducible figures for the manuscript:

> **Scale-Dependent Spatial Order in a Tropical Forest: A Dual-Regime Framework with Implications for Urban Planting Design**
> 
> Talha Aksoy, Department of Landscape Architecture, Kırklareli University, Turkey

## Abstract

We employed the Structure Factor S(k) from condensed matter physics to analyze spatial patterns of three tree species in the Barro Colorado Island (BCI) 50-ha plot. Our analysis reveals a **dual-regime spatial architecture**:

- **Physical Regime (k > 0.5 m⁻¹):** Universal convergence to S(k) ≈ 0.8-1.0 driven by crown/root exclusion
- **Ecological Regime (k < 0.5 m⁻¹):** Species-specific strategies from extreme clustering to near-randomness

These findings provide quantifiable biomimetic benchmarks for urban forestry planting design.

## Key Findings

| Species | N | Strategy | S(k) Profile |
|---------|---|----------|--------------|
| *H. concinna* | 312 | Liquid-like (balanced) | Min S(k) ≈ 0.85 |
| *G. sebifera* | 470 | Random-like (generalist) | S(k) ≈ 1.0 |
| *G. superba* | 579 | Clustered (specialist) | Low-k S(k) > 100 |

## Repository Structure

```
├── fig1_temporal_evolution.py    # Figure 1: 40-year temporal analysis
├── fig2_dual_regime.py           # Figure 2: Three-species comparison
├── README.md                     # This file
├── requirements.txt              # Python dependencies
└── figures/
    ├── fig1_temporal_evolution.png
    └── fig2_dual_regime.png
```

## Installation

```bash
# Clone repository
git clone https://github.com/[your-username]/bci-spatial-analysis.git
cd bci-spatial-analysis

# Install dependencies
pip install -r requirements.txt
```

## Requirements

```
numpy>=1.20.0
pandas>=1.3.0
scipy>=1.7.0
matplotlib>=3.4.0
```

## Data

The analysis uses the BCI 50-ha plot census data from ForestGEO:

> Condit R., Perez, R., Aguilar, S., Lao, S., Foster, R., Hubbell, S.P. (2019). 
> Complete data from the Barro Colorado 50-ha plot: 423617 trees, 35 years.
> Dryad. https://doi.org/10.15146/5xcp-0d46

**Download:** [FullMeasurementBCI.zip](https://datadryad.org/dataset/doi:10.15146/5xcp-0d46) (33.7 MB)

### Data Preparation

Extract species of interest from the full dataset:

```python
import pandas as pd

df = pd.read_csv("FullMeasurementBCI.tsv", sep='\t', low_memory=False)

# Filter target species
species = ['concinna', 'sebifera', 'superba']
filtered = df[df['SpeciesName'].isin(species)]
filtered.to_csv("BCI_3species.tsv", sep='\t', index=False)
```

## Usage

### Figure 1: Temporal Evolution

```bash
python fig1_temporal_evolution.py
```

Analyzes *H. concinna* spatial order across 8 censuses (1982-2015), demonstrating persistent "liquid-like" structure (S(k) ≈ 0.80-0.89).

### Figure 2: Dual-Regime Analysis

```bash
python fig2_dual_regime.py
```

Compares three species with contrasting ecological strategies, revealing scale-dependent divergence in spatial organization.

## Methods

### Structure Factor Calculation

The Structure Factor S(k) for a 2D point pattern is calculated as:

$$S(k) = 1 + \frac{2}{N} \sum_{i < j} J_0(k |\mathbf{r}_i - \mathbf{r}_j|)$$

where:
- $J_0$ is the zero-order Bessel function of the first kind
- $|\mathbf{r}_i - \mathbf{r}_j|$ is the Euclidean distance between points
- $k$ is the wavenumber (spatial frequency)

### Interpretation

| S(k) Value | Interpretation |
|------------|----------------|
| S(k) = 1 | Complete Spatial Randomness (Poisson) |
| S(k) < 1 | Repulsive order ("liquid-like") |
| S(k) > 1 | Clustering |
| S(k) → 0 as k → 0 | Hyperuniformity |

## Results

### Temporal Stability (Figure 1)

*H. concinna* maintains stable "liquid-like" order over 40 years:
- S(k) range: 0.80 - 0.89
- Population increase: 245 → 312 (27%)
- No drift toward "jammed" state

### Dual-Regime Structure (Figure 2)

**Physical Regime (short-range, k > 0.5 m⁻¹):**
- All species converge to S(k) ≈ 0.8-1.0
- Driven by physical exclusion (crowns cannot overlap)

**Ecological Regime (large-scale, k < 0.5 m⁻¹):**
- *G. superba*: Extreme clustering (S(k) > 100) - habitat specialist
- *G. sebifera*: Near-random (S(k) ≈ 2-3) - generalist
- *H. concinna*: Intermediate (S(k) ≈ 5-10) - balanced strategy

## Applications

### Urban Forestry Design Protocol

Based on our findings, we propose:

1. **Enforce minimum exclusion zones:** ~2.0-2.5 m for canopy trees
2. **Target S(k) ≈ 0.85** at short ranges using algorithmic placement
3. **Allow large-scale flexibility** depending on design goals

## Citation

If you use this code or methodology, please cite:

```bibtex
@article{aksoy2025spatial,
  title={Scale-Dependent Spatial Order in a Tropical Forest: A Dual-Regime Framework with Implications for Urban Planting Design},
  author={Aksoy, Talha},
  journal={Urban Forestry \& Urban Greening},
  year={2025},
  doi={10.xxxx/xxxxx}
}
```

## License

This code is released under the MIT License. See [LICENSE](LICENSE) for details.

The BCI data is subject to ForestGEO data use policies. Please cite the data source appropriately.

## Contact

**Talha Aksoy**  
Department of Landscape Architecture  
Kırklareli University, Turkey  
📧 talha.aksoy@klu.edu.tr

## Acknowledgments

- BCI forest dynamics plot data from ForestGEO/STRI
- AI assistance from Google Gemini and Anthropic Claude for code development and manuscript preparation

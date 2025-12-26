# BCI Forest Spatial Analysis

**Structure Factor Analysis of Tropical Forest Spatial Patterns**

This repository contains Python scripts and data for analyzing spatial patterns in the Barro Colorado Island (BCI) 50-hectare forest dynamics plot using the Structure Factor S(k) from condensed matter physics.

## Publication

**Scale-Dependent Spatial Order in a Tropical Forest: A Dual-Regime Framework for Urban Planting Design**

*Talha Aksoy*  
Department of Landscape Architecture, Kırklareli University, Turkey

Target Journal: Urban Forestry & Urban Greening

## Key Findings

### The Dual-Regime Hypothesis

We demonstrate that tropical forest spatial structure is fundamentally **scale-dependent**:

| Regime | Scale | Wavenumber | Behavior |
|--------|-------|------------|----------|
| **Physical** | < 12 m | k > 0.5 m⁻¹ | Universal S(k) ≈ 0.8-1.0 (crown exclusion) |
| **Ecological** | > 12 m | k < 0.5 m⁻¹ | Species-specific (clustering to random) |

### Species Comparison (Census 8, 2015)

| Species | N | Min S(k) | Strategy |
|---------|---|----------|----------|
| *H. concinna* | 312 | 0.84 | Liquid-like (balanced) |
| *G. sebifera* | 470 | 0.83 | Random-like (generalist) |
| *G. superba* | 579 | 0.74 | Clustered (specialist) |

## Repository Contents

```
bci-spatial-analysis/
├── README.md
├── requirements.txt
├── BCI_3species.tsv          # Filtered dataset (~4 MB)
├── fig1_temporal_evolution.py # H. concinna 40-year analysis
├── fig2_dual_regime.py        # Three-species comparison
└── prepare_data.py            # Data filtering script
```

## Quick Start

### 1. Clone Repository
```bash
git clone https://github.com/tlhksy/bci-spatial-analysis.git
cd bci-spatial-analysis
```

### 2. Install Dependencies
```bash
pip install -r requirements.txt
```

### 3. Run Analysis
```bash
# Figure 1: Temporal evolution (1982-2015)
python fig1_temporal_evolution.py

# Figure 2: Dual-regime comparison
python fig2_dual_regime.py
```

## Data Source

The BCI dataset is publicly available from the ForestGEO data portal:
- **DOI:** https://doi.org/10.15146/5xcp-0d46
- **Portal:** https://forestgeo.si.edu/

### Included Dataset

`BCI_3species.tsv` contains filtered data for three focal species:
- *Heisteria concinna* (Olacaceae)
- *Guarea sebifera* (Meliaceae)  
- *Gustavia superba* (Lecythidaceae)

## Methods

### Structure Factor S(k)

The radially-averaged Structure Factor for 2D point patterns:

$$S(k) = 1 + \frac{2}{N} \sum_{i<j} J_0(k \cdot r_{ij})$$

Where:
- $k$ = wavenumber (m⁻¹)
- $r_{ij}$ = distance between trees $i$ and $j$
- $J_0$ = zero-order Bessel function
- $N$ = number of trees

### Interpretation

| S(k) Value | Meaning |
|------------|---------|
| S(k) = 1 | Complete Spatial Randomness (Poisson) |
| S(k) < 1 | Repulsive order ("liquid-like") |
| S(k) > 1 | Clustering |

### Wavenumber to Distance Conversion

$$\lambda = \frac{2\pi}{k}$$

| k (m⁻¹) | λ (m) | Scale |
|---------|-------|-------|
| 0.1 | 62.8 | Landscape |
| 0.3 | 20.9 | Grove/patch |
| **0.5** | **12.6** | **Regime boundary** |
| 1.0 | 6.3 | Inter-tree |
| 2.0 | 3.1 | Crown zone |

### Statistical Testing

- **Monte Carlo Null Envelope:** 100 CSR simulations for significance testing
- **Bootstrap Confidence Intervals:** 50 resamples for S(k) uncertainty

## Applications

### Biomimetic Urban Planting Design

1. **Enforce minimum exclusion zones:** 2.0-2.5 m for canopy trees
2. **Target S(k) ≈ 0.85 at short ranges:** Liquid-like order
3. **Allow large-scale flexibility:** Uniform to clustered arrangements

## Citation

If you use this code or methodology, please cite.

## Acknowledgments

- BCI forest dynamics plot: ForestGEO network, Smithsonian Tropical Research Institute
- Data collection: STRI staff and collaborators (1981-2015)

## License

MIT License

## Contact

Talha Aksoy  
talha.aksoy@klu.edu.tr  
Kırklareli University, Turkey

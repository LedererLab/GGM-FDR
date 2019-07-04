# False Discovery Rates in Biological Networks

This repository provides the data and implementation of the methods described in the paper *False Discovery Rates in Biological Networks* by Lu Yu, Tobias Kaufmann and Johannes Lederer.

## Usage
The file `Knockoff.R` contains a function `GraphEstimation` to estimate a connectivity graph using the KO or KO+ method. Our paper contains detailed descriptions of these methods.

## Simulation
We provide an example code `Simulation.R` for the comparison of the estimators' accuracy in terms of power and FDR among KO, graphical lasso, neighborhood selection with the “and-rule” and the “or-rule”, thresholding the correlation matrix, and thresholding the partial correlation matrix.

## Real Data Analyses

###### Brain Connectivity Analysis
The program for Brain Connectivity Analysis can be found in `RealDataAnalysis/BrainConnectivity.R`. The corresponding datasets can be found under `RealDataAnalysis/Data/BrainConnectivity/AAL_YAN.csv`, `RealDataAnalysis/Data/BrainConnectivity/BrainRegions.csv` and `RealDataAnalysis/Data/BrainConnectivity/detrendedData`

###### Human Microbiome Analysis
The program for Human Microbiome Analysis is included, and can be found in `RealDataAnalysis/AmericanGut.R`. The corresponding datasets can be found under `RealDataAnalysis/Data/AmericanGut/ag-cleaned_L2.txt`.

###### Atlantic Amphibians Abundance Analysis
The program for Atlantic Amphibians Abundance Analysis can be found under `RealDataAnalysis/AtlanticAmphibians.R`. The corresponding datasets can be found in `RealDataAnalysis/Data/AtlanticAmphibians/ATLANTIC_AMPHIBIANS_sites.csv`, and `RealDataAnalysis/Data/AtlanticAmphibians/ATLANTIC_AMPHIBIANS_species.csv`


## Repository Authors
- Lu Yu — Ph.D. candidate in Statistics, University of Toronto 

- Tobias Kaufmann — 

- Johannes Lederer — Professor in Mathematical Statistics, Ruhr-University Bochum

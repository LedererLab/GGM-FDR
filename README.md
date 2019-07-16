# False Discovery Rates in Biological Networks

This repository provides the data and implementation of the methods described in [(Yu et al., 2018)](https://arxiv.org/pdf/1907.03808.pdf).

## Usage
The file `Knockoff.R` contains a function `GraphEstimation` to estimate a connectivity graph using the KO or KO+ method. Our paper contains detailed descriptions of these methods.

## Simulation
We provide an example code `Simulation.R` for the comparison of the estimators' accuracy in terms of power and FDR among KO, graphical lasso, neighborhood selection with the “and-rule” and the “or-rule”, thresholding the correlation matrix, and thresholding the partial correlation matrix. This program requires R 3.4.4 or older version of R.

## Real Data Analyses

###### Brain Connectivity Analysis  
The corresponding datasets can be found under `RealDataAnalysis/Data/BrainConnectivity/AAL_YAN.csv`, `RealDataAnalysis/Data/BrainConnectivity/BrainRegions.csv` and `RealDataAnalysis/Data/BrainConnectivity/detrendedData`. The data is described and analyzed in [(Bu and Lederer, 2017)](https://arxiv.org/pdf/1704.02739.pdf). We refer to that paper for more details about the data. The program for Brain Connectivity Analysis can be found in `RealDataAnalysis/BrainConnectivity.R`. 


###### Human Microbiome Analysis
The program for Human Microbiome Analysis is included, and can be found in `RealDataAnalysis/AmericanGut.R`. The corresponding datasets can be found under `RealDataAnalysis/Data/AmericanGut/ag-cleaned_L2.txt`.

###### Atlantic Amphibians Abundance Analysis
The program for Atlantic Amphibians Abundance Analysis can be found under `RealDataAnalysis/AtlanticAmphibians.R`. The corresponding data sets `ATLANTIC_AMPHIBIANS_sites.csv` and `ATLANTIC_AMPHIBIANS_species.csv` can be found in [ATLANTIC AMPHIBIANS: a data set of amphibian communities from the Atlantic Forests of South America](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2392). The data is described in [(Vancine et al., 2018)](https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.2392&file=ecy2392-sup-0002-MetadataS1.pdf).


## Repository Authors
- Lu Yu — Ph.D. student in Statistics, University of Toronto 

- Tobias Kaufmann — Neuroscientist at NORMENT, University of Oslo and Oslo University Hospital

- Johannes Lederer — Professor in Mathematical Statistics, Ruhr-University Bochum


## Reference
[Yu, L., Kaufmann, T., and Lederer, J. (2019) False Discovery Rates in Biological Network.](https://arxiv.org/pdf/1907.03808.pdf)

[Bu, Y. and Lederer, J. (2017) Integrating additional knowledge into estimation of graphical models.](https://arxiv.org/pdf/1704.02739.pdf)

[Vancine, H., Duarte, K., de Souza, Y. et al. (2018) ATLANTIC AMPHIBIANS: a data set of amphibian communities from the Atlantic Forests of South America.](https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.2392&file=ecy2392-sup-0002-MetadataS1.pdf)





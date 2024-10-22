# False Discovery Rates in Biological Networks

This repository provides the data and implementations of the methods described in [(Yu et al., 2019)](https://arxiv.org/pdf/1907.03808.pdf). 

## Usage
The file `Knockoff.R` contains a function `GraphEstimation` to estimate a connectivity graph using the KO or KO+ method. The aforementioned paper contains detailed descriptions of these methods.

## Simulations
We provide an example code `Simulation.R` for a comparison in terms of power and FDR among KO, graphical lasso, neighborhood selection with the “and-rule” and the “or-rule”, thresholding the correlation matrix, and thresholding the partial correlation matrix. This program requires `R 3.4.4` or earlier.

## Real Data Analyses

###### Brain Connectivity Analysis  
The processed data is in `RealDataAnalysis/Data/BrainConnectivity/AAL_YAN.csv`, `RealDataAnalysis/Data/BrainConnectivity/BrainRegions.csv`, and `RealDataAnalysis/Data/BrainConnectivity/detrendedData`. 
The raw fMRI data was collected and provided by Dr. Dantao Peng, Dr. Yanlei Mu, and Dr. Xiao Zhang. 
The data preprocessing was conducted by Dr. Min Zhang.
The data is described and analyzed in [(Bu and Lederer, 2017)](https://arxiv.org/pdf/1704.02739.pdf). 
The programs for our statistical analysis are in `RealDataAnalysis/BrainConnectivity.R`. 

###### Human Microbiome Analysis
The processed data `RealDataAnalysis/Data/AmericanGut/ag-cleaned_L2.txt` is downloaded from the [American Gut Project](http://humanfoodproject.com/americangut/). 
We refer to [American Gut data repository](https://github.com/biocore/American-Gut) for more details about the data proprocessing. 
The programs for our statistical analysis are in `RealDataAnalysis/AmericanGut.R`.

###### Atlantic Amphibians Abundance Analysis
The data `ATLANTIC_AMPHIBIANS_sites.csv` and `ATLANTIC_AMPHIBIANS_species.csv` are taken from [ATLANTIC AMPHIBIANS: a data set of amphibian communities from the Atlantic Forests of South America](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2392). 
The data is described in [(Vancine et al., 2018)](https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.2392&file=ecy2392-sup-0002-MetadataS1.pdf)
and [ATLANTIC AMPHIBIANS data repository](https://github.com/mauriciovancine/ATLANTIC-AMPHIBIANS/tree/v1.0.0).
The program for our statistical analysis are in `RealDataAnalysis/AtlanticAmphibians.R`. 


## Repository Authors
- Lu Yu — Ph.D. student in Statistics, University of Toronto 

- Tobias Kaufmann — Neuroscientist at NORMENT, University of Oslo and Oslo University Hospital

- Johannes Lederer — Professor in Mathematical Statistics, Ruhr-University Bochum


## Reference
[Yu, L., Kaufmann, T., and Lederer, J. (2019) False Discovery Rates in Biological Network.](https://arxiv.org/pdf/1907.03808.pdf)

[Bu, Y. and Lederer, J. (2017) Integrating additional knowledge into estimation of graphical models.](https://arxiv.org/pdf/1704.02739.pdf)

[Vancine, H., Duarte, K., de Souza, Y. et al. (2018) ATLANTIC AMPHIBIANS: a data set of amphibian communities from the Atlantic Forests of South America.](https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.2392&file=ecy2392-sup-0002-MetadataS1.pdf)






## Paper ##

> Paper link on bioarxiv will go here


## Dependencies ##

All dependencies required to run the code can be found in Environment directory. Create the virtual environment and load it using [conda](https://docs.conda.io/en/latest/). Note this only applies for Python packages, packages in R will have to be downloaded separately. 

```
conda create -f /Environment/environment.yaml
source activate protein
```

## Analysis ##

All scripts used to run the analysis, which is written in Python and R, is under the `Analysis` directory. The rpy2 Python package is used to import R code into Python and can be found under `GetRCodeIntoPython.py`. Code in R is used for visualization and statistical analysis. All of the statistical analysis (gene set enrichment and regression modeling) is aggregated using `CalculateMetrics.py` and the R code imported can be found under `GetRegressionStats.R`. Input data and plotting of all figures can be used reproduce all of the analysis in `PlotAllFigures.py`, which uses visualization done in ggplot within `Plotting.R`.

All of the raw data imported and used for this analysis can be found under `GetData.py`. (Note that all of this data is already publicly available and links to get access to the data are provided in the manuscript.) All of the alternative splicing analysis used to examine gene silencing in high mutational load tumors can be found under `Splicing.py`. 

## Annotations ##

Proteostasis gene sets and CORUM complexes datatbases used are provided in the `GeneSets` directory. Other gene and drug category annotations used can be found under `GetAnnotations.py` and `DrugAnnotationMoa.R`, respectively. 

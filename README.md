# CyEMD: Differential analysis method developed for CyTOF data
R package providing a differential analysis method using the Earth Mover´s Distance to compare normalized distributions


## Installation
Currently the package is only available on GitHub. For installation use one of the following possibilities:

1. The package can be directly installed from GitHub running
   
   ```remotes::install_github("https://github.com/biomedbigdata/CyEMD_package.git")```

2. The package can be installed after cloning the repository running

   ```devtools::install("~/path/to/package")```


## Dependencies

R >= 4.2

Packages available on CRAN:
* data.table
* Rcpp
* RcppAlgos
* stats
* graphics

Packages available on Bioconductor:
* SummarizedExperiment
* CATALYST


## Cite   
Lis Arend, Judith Bernett, Quirin Manz, Melissa Klug, Olga Lazareva, Jan Baumbach, Dario Bongiovanni, Markus List, 
A systematic comparison of novel and existing differential analysis methods for CyTOF data, 
Briefings in Bioinformatics, Volume 23, Issue 1, January 2022, bbab471, https://doi.org/10.1093/bib/bbab471

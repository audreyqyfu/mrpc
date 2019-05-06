
MRPC builds on existing PC algorithms and learns a causal network with increased accuracy.  The inferred causal network contains directed and undirected edges, with the direction indicating causality.  For genomic data, MRPC determines edge direction under the principle of Mendelian randomization when genotype and molecular phenotype (e.g. gene expression) data are both available at the individual level. Nodes in the inferred network may be a genotype or a molecular phenotype.  
Development of the R package MRPC is at https://github.com/audreyqyfu/mrpc, and official releases are available on CRAN: https://cran.r-project.org/web/packages/MRPC/index.html.

References:

Md. Bahadur Badsha, Audrey Qiuyan Fu. Learning causal biological networks with the principle of Mendelian randomization. bioRxiv. 171348. doi:10.1101/171348.

Md. Bahadur Badsha, Evan A Martin, Audrey Qiuyan Fu. MRPC: An R package for accurate inference of causal graphs.  arXiv. 1806.01899.

## Installation

### 1. Installation of the most recent version from GitHub.

First install the R package devtools available on CRAN, if it is not already installed. This package provides function `install_github()` that enables installing packages directly from github with the following command.

Invoke R and then type with the following command:
```
R> install.packages ("devtools")
R> library (devtools)
#install R packages that MRPC depends on before running the next line 
#see details below
R>install_github ("audreyqyfu/mrpc")
```
MRPC depends on several R packages from CRAN and from Bioconductor.  It is likely that some of these packages are not installed on your computer.  If the R package is available on CRAN, you may use the following command line for installation (change _packagename_ to the name of the package to be installed, e.g, bnlearn, pcalg, etc.) before running function `install_github`:
```
R> install.packages("packagename")
```

The following Bioconductor packages also need to be installed before running function `install_github`:
```
R> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
R> BiocManager::install ('RBGL')
R> BiocManager::install ('Rgraphviz')
R> BiocManager::install ('GO.db')
R> BiocManager::install ('impute')
R> BiocManager::install ('preprocessCore')
```
### 2. Installation from the source of a released package.

Download the package source MRPC_xxx.tar.gz.  

In Terminal, navigate to the directory where the package is stored, and run the following command line:
```bash
$ R CMD INSTALL MRPC_xxx.tar.gz
```
Again, you may need to first install the Bioconductor packages that MRPC depends on using the instructions above.
Alternatively, you may also run the following command line in R, after changing the working directory to where MRPC_xxx.tar.gz is stored on your computer:
```
R> install.packages("MRPC_xxx.tar.gz", repos = NULL, type="source")
```
### 3. Installation from CRAN.

Official releases are available on CRAN.  To install,
```
R> install.packages("MRPC")
```
## Using MRPC
After installation, load the MRPC package into R:
```
R> library (MRPC)
```
Bring up the documentation of the package:
```
R> library (help=MRPC)
```


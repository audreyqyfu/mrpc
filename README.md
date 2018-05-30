The MRPC is the first machine learning algorithm that employs Principle of Mendelian randomization(PMR)and 
learns a causal biological network from integrating genotype and molecular phenotype (such as gene expression) data. 
MRPC, which incorporates Principle of MR (gMR) into PC algorithms and learns acausal graph where the nodes are genetic 
variants and molecular phenotypes (such as gene expression), and where the edges between nodes are undirected or directed, 
with the direction indicating causality. Crucially, by combining gMR with machine learning, our method is efficient and 
accurate.

Package implementation/installation: 

The R package MRPC, can be downloaded with git from github in the following link:
https://github.com/audreyqyfu/mrpc. To install the MRPC package, need to install the devtools package from CRAN first 
(if devtools already not install). Invoke R and then type with the following command:

R> install.packages(“devtools”)

R>library(devtools)

The devtools package provides install_github( ) that enables installing packages from github with the following command:

Based on our experience, some of the bioconductor packages needed for MRPC. These packages need to be installed before attempting to install MRPC.

R> source ('https://bioconductor.org/biocLite.R')

R> biocLite('Rgraphviz')

R> biocLite('GO.db')

R> biocLite('impute')

R> biocLite('preprocessCore')

When asked to update all/some/none [a/s/n], always select n.

R>install_github("audreyqyfu/mrpc")

Load the MRPC package.

R> library(MRPC)

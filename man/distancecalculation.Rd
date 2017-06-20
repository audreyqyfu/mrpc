\name{distancecalculation}
\alias{distancecalculation}

\title{
Structural Distance Calculation 
}
\description{
Structural distance calculation between two graph.If missing edges then penalty =1 and if missing direction then penalty =0.5. 
}
\usage{
distancecalculation(g1, g2, NQ)
}

\arguments{
  \item{g1}{
First graph object
}
  \item{g2}{
Second graph object
}
\item{NQ}{Need to specify for how many genetic variants (SNPs/indels/CNV/eQTL) in your input data. For example, if your dataset have one genetic variants (first column), then NQ=1, if 2 (1st and 2nd Column), then NQ=2..and so on.
}
}


\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

\examples{
## True model
n<-nrow (simu.data_M1)     #Number of row
V<-colnames(simu.data_M1)  #Column names

Rcor_R=RobustCor(simu.data_M1, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_1_GT <- MRPC(simu.data_M1,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest', labels= V, verbose = TRUE)

#Inferred graph 

n<-nrow (simu.data_M2)    #Number of row
V<-colnames(simu.data_M1) #Column names

Rcor_R=RobustCor(simu.data_M2, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_2_GT <- MRPC(simu.data_M2,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest', labels= V, verbose = TRUE)

#
Distance=distancecalculation(MRPC.fit_1_GT$graph,MRPC.fit_2_GT$graph,NQ=1)

}


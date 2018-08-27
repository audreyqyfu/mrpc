\name{SimulatedData}
\alias{SimulatedData}

\title{
Generating Simulated Data
}
\description{
Function for the generating simulated data
}
\usage{
SimulatedData(N, p,model,b0.1,b1.1,b1.2,b1.3,sd.1)
}

\arguments{

\item{N}{Number of observation}
\item{p}{Population frequency of the references allel. Real number between 0 to 1, which is the number of a particluar allel present}
\item{model}{To specify the model you want. For example, if you want to generate data for model 0 you would type 'model0' into the function.}
\item{b0.1}{Intercept of b0.1+b1.1*P1, where P1 is the 1st parent of the corresponding genes.}
\item{b1.1}{Slope of P1 for b0.1+b1.1*P1, where P1 is the 1st parent of the corresponding genes.}
\item{b1.2}{Slope of P2 for b0.1+b1.1*P1+b1.2*P2, where P2 is the 2nd parent of the corresponding genes.}
\item{b1.3}{Slope of P3 for b0.1+b1.1*P1+b1.2*P2+b1.3*P3, where P3 is the 3rd parent of the corresponding genes.}
\item{sd.1}{Standard deviation for corresponding data generated genes.}

}

\details{
In here, we used first column of the input matrix is the genotype of the expression quantitative trait loci (eQTL)/Copy number variation (CNVs) and the remaining columns are the gene expression data.
}

\value{
Matrix
}

\author{
Md Bahadur Badsha <mdbadsha@uidaho.edu>,Evan A Martin <mart9986@vandals.uidaho.edu> and Audrey Qiuyan Fu <audreyf@uidaho.edu>
}


\seealso{
\link{MRPC} for estimating a DAG using the modified pc (MRPC) algorithm; \link{ModiSkeleton} for estimating a undirected graph using modified \link{skeleton} function; \link{EdgeOrientation} for orientation rules for edges in MRPC algorithm; \link{Case_NP} the gene to have data generated has no parent; \link{Case_1P} the gene to have data generated has one parent; \link{Case_2P} the gene to have data generated has two parents.

}
\examples{

## Data pre-processing

#If you use only one genotype of the expression quantitative trait loci
#(eQTL)/Copy number variation (CNV), the 1st column of
#the input matrix will be #eQTL/CNV and the remaining
#columns are the gene expression data.


## Model 0
simu.data_M0<-SimulatedData(N = 10^3,p = 0.45
,'model0',b0.1 = 0,b1.1 = 1,b1.2 = 1,b1.3 = 1,sd.1 = 1)

## Model 1
simu.data_M1<-SimulatedData(N = 10^3,p = 0.45
,'model1',b0.1 = 0,b1.1 = 1,b1.2 = 1,b1.3 = 1,sd.1 = 1)

## Model 2
simu.data_M2<-SimulatedData(N = 10^3,p = 0.45
,'model2',b0.1 = 0,b1.1 = 1,b1.2 = 1,b1.3 = 1,sd.1 = 1)

## Model 3
simu.data_M3<-SimulatedData(N = 10^3,p = 0.45
,'model3',b0.1 = 0,b1.1 = 1,b1.2 = 1,b1.3 = 1,sd.1 = 1)

## Model 4
simu.data_M4<-SimulatedData(N = 10^3,p = 0.45
,'model4',b0.1 = 0,b1.1 = 1,b1.2 = 1,b1.3 = 1,sd.1 = 1)

## Multiple Parents Model
simu.data_multiparent<-SimulatedData(N = 10^3,p
=0.45, 'multiparent',b0.1 = 0,b1.1 = 1
,b1.2 = 1,b1.3=1,sd.1 = 1)

## Star Model
simu.data_starshaped<-SimulatedData(N = 10^3,p = 0.45, 
'starshaped',b0.1 = 0,b1.1 = 1,b1.2 = 1,b1.3 = 1,sd.1= 1)

## Layers Model
simu.data_layered<-SimulatedData(N = 10^3,p = 0.45,
'layered',b0.1 = 0,b1.1 = 1,b1.2 = 1,b1.3 =
1,sd.1 = 1)

}

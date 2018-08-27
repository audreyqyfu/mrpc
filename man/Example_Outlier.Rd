\name{Example_Outlier}
\alias{Example_Outlier}
\title{
Impact of outliers on graph inference
}
\description{
The data contain two genotype nodes, V1 and V2, and 3 phenotype nodes, T1,T2 and T3. The genotype nodes are discrete, whereas the phenotype nodes are continuous. We added 10 artificial outliers to our original data matrix randomly.First, we generated uniform distribution for outlying expression and add with the original data. Then we check the robustness by three methods (MRPC< mmhc and pc) in the presence of outliers. 
}

\author{
Md Bahadur Badsha <mdbadsha@uidaho.edu>,Evan A Martin <mart9986@vandals.uidaho.edu> and Audrey Qiuyan Fu <audreyf@uidaho.edu>
}

\examples{
library(pcalg) #for pc
library(bnlearn) #for mmhc

#Truth with and without outlier
tarmat=matrix(0,nrow=ncol(Data_withoutoutlier),ncol = ncol(Data_withoutoutlier))
colnames(tarmat)=colnames(Data_withoutoutlier)
rownames(tarmat)=colnames(Data_withoutoutlier)
tarmat[1,2]=1
tarmat[2,1]=1
tarmat[2,3]=1
tarmat[4,3]=1
tarmat[4,5]=1
Truth<-as(tarmat, "graphNEL")

#For withoutoutlier data
n<-nrow (Data_withoutoutlier)    #Number of row
V<-colnames(Data_withoutoutlier) #Column names

#Classical correlation
suffStat_C1= list(C = cor(Data_withoutoutlier), n = n)

#Robust correlation (Beta=0.005)
Rcor_R1=RobustCor(Data_withoutoutlier, 0.005) 
suffStat_R1= list(C = Rcor_R1$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_withoutoutlier<- MRPC(Data_withoutoutlier,suffStat_R1
,GV=2,FDR=0.05,indepTest ='gaussCItest',labels=V,verbose = TRUE
)

## Estimated graph by mmhc
data=data.frame(Data_withoutoutlier)
mmhc_withoutoutlier=mmhc(data)

## Estimated graph by existing pc with alpha=0.05
pc.fit_withoutoutlier<- pc(suffStat_C1,indepTest = gaussCItest
,alpha=0.05,labels =V,verbose = TRUE)

#For withoutlier data
n<-nrow (Data_withoutlier)    #Number of row
V<-colnames(Data_withoutlier) #Column names

#Classical correlation
suffStat_C1= list(C = cor(Data_withoutlier), n = n)

#Robust correlation (Beta=0.005)
Rcor_R1=RobustCor(Data_withoutlier, 0.005) 
suffStat_R1= list(C = Rcor_R1$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_withoutlier<- MRPC(Data_withoutlier,suffStat_R1,GV=2
,FDR=0.05,indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
data=data.frame(Data_withoutlier)
mmhc_withoutlier=mmhc(data)

## Estimated graph by existing pc with alpha=0.05
pc.fit_withoutlier<- pc(suffStat_C1,indepTest = gaussCItest
,alpha=0.05,labels =V,verbose = TRUE)

## Plot of the results
## Show estimated graph
#if (require(Rgraphviz)){
par(mfrow=c(2,4))
plot(Truth,main="Truth")
plot(MRPC.fit_withoutoutlier@graph,main="MRPC")
graphviz.plot(mmhc_withoutoutlier,main="mmhc")
plot(pc.fit_withoutoutlier,main="pc")
plot(Truth,main=" ")
plot(MRPC.fit_withoutlier@graph,main=" ")
graphviz.plot(mmhc_withoutlier,main=" ")
plot(pc.fit_withoutlier,main=" ")
#}
}

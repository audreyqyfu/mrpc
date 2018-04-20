\name{seqDiff}
\alias{seqDiff}

\title{
Deviation of the Two Sequnce 
}
\description{
Deviation of the inferred from the truth
}
\usage{
seqDiff(g1, g2)
}

\arguments{
  \item{g1}{
Adjacency matrix from the first graph object
}
  \item{g2}{
Adjacency matrix from the second graph object
}
}

\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

\examples{

#Truth graph (V1-->T1-->T2-->T3)
tarmat_s1=matrix(0,nrow=4,ncol = 4)
colnames(tarmat_s1)=c("V1","T1","T2","T3")
rownames(tarmat_s1)=colnames(tarmat_s1)
#Create an adjacency matrix for truth
tarmat_s1[1,2]=1
tarmat_s1[2,3]=1
tarmat_s1[3,4]=1

#Inferred graph (V1-->T1<--T2-->T3)
tarmat_s2=matrix(0,nrow=4,ncol = 4)
colnames(tarmat_s2)=c("V1","T1","T2","T3")
rownames(tarmat_s2)=colnames(tarmat_s2)
#Create an adjacency matrix for inferred
tarmat_s2[1,2]=1
tarmat_s2[3,2]=1
tarmat_s2[3,4]=1

#Deviation of the inferred from the truth
seqDiff(tarmat_s2,tarmat_s1)
         
}

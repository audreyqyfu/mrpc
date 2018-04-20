\name{Recall_Precision}
\alias{Recall_Precision}

\title{
Performance Evaluation by Recall and Precision in MRPC
}
\description{
Recall=(# edges correctly identified in inferred graph) / (# edges in true graph); #Precision = (# edges correctly identified in inferred graph) / (# edges in inferred graph). 
}
\usage{
Recall_Precision(g1, g2, GV, edge.presence, edge.direction)
}

\arguments{
  \item{g1}{
First graph object, we say from the truth graph
}
  \item{g2}{
Second graph object, we say from the inferred graph
}
  \item{GV}{
Need to specify for how many genetic variants (SNPs/indels/CNV/eQTL) in your input data for the consistence of MRPC with others inference, because only the MRPC implemented the gMR.For example, if your dataset have one genetic variants (first column), then GV=1, if 2 (1st and 2nd Column), then GV=2..and so on.
}
  \item{edge.presence}{
specify the weights for edge presence
}
  \item{edge.direction}{
specify the weights for edge direction
}
}
\details{
We consider it more important to be able to identify the presence of an edge than to also get the direct correct.  Therefore, we assign 1 to an edge with the correct direction and 0.5 to an edge with the wrong direction or no direction (Badsha et al., 2018).
}

\value{
\itemize{
\item \code{TP}: Total found edge in the inferred graph and edge exit in the true graph 
\item \code{FP}: Total found edge in the inferred graph but no edge exit in the true graph
\item \code{Recall}: Power, or sensitivity measures how many edges from the true graph a method can recover.
\item \code{Precision}: Measures how many correct edges are recovered in the inferred graph
}

}

\references{
1.Md. Bahadur Badsha, Audrey Qiuyan Fu: Learning causal biological networks with the principle of Mendelian randomization (Nat Comn; under reveiw 2018), 

}
\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}


\seealso{
\link{SimulatedData}: adjusted Structural Hamming Distance (aSHD)
}
\examples{

library(pcalg)   #for graph object
#Truth graph (V1-->T1-->T2-->T3)
tarmat_s1=matrix(0,nrow=4,ncol = 4)
colnames(tarmat_s1)=c("V1","T1","T2","T3")
rownames(tarmat_s1)=colnames(tarmat_s1)
#Create an adjacency matrix for truth
tarmat_s1[1,2]=1
tarmat_s1[2,3]=1
tarmat_s1[3,4]=1
#Graph object for truth
Truth<-as(tarmat_s1, "graphNEL")

#Inferred graph (V1-->T1<--T2-->T3)
tarmat_s2=matrix(0,nrow=4,ncol = 4)
colnames(tarmat_s2)=c("V1","T1","T2","T3")
rownames(tarmat_s2)=colnames(tarmat_s2)
#Create an adjacency matrix for truth
tarmat_s2[1,2]=1
tarmat_s2[3,2]=1
tarmat_s2[3,4]=1
#Graph objects for inferred
Inferred<-as(tarmat_s2, "graphNEL")
         
#Recall and Precision
Recall_Precision(Truth, Inferred, GV=1, edge.presence=1.0, edge.direction=0.5)
}      

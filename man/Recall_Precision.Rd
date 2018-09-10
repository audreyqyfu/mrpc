\name{Recall_Precision}
\alias{Recall_Precision}

\title{
Performance Evaluation by Recall and Precision in MRPC
}
\description{
Recall = (# edges correctly identified in inferred graph) / (# edges in true graph); Precision = (# edges correctly identified in inferred graph) / (# edges in inferred graph).For example, when the true graph is V-->T1-->T2, and the inferred graphs are i) V-->T1-->T2, and V-->T2; ii) V-->T1--T2; and iii) V-->T1<--T2, the number of correctly identified edges is then 2, 1.5, and 1.5, respectively. Recall is calculated to be 1, 0.75, and 0.75 respectively, whereas precision is 2/3 = 0.67, 1.5/2 = 0.75, and 1.5/2 = 0.75, respectively.  
}
\usage{
Recall_Precision(g1, g2, GV, edge.presence = 1.0, edge.direction = 0.5)
}

\arguments{
  \item{g1}{
First graph object, from the true graph
}
  \item{g2}{
Second graph object, from the inferred graph
}
  \item{GV}{
The number of genetic variants (SNPs/indels/CNV/eQTL) in the input data. For example, if the data has one genetic variant, first column, then GV = 1, if 2, 1st and 2nd Column, then GV = 2, and so on.
}
\item{edge.presence}{
The weight for an edge being present.
}
  \item{edge.direction}{
The weight for the edge direction.
}
}
\details{
We consider it more important to be able to identify the presence of an edge than to also get the direct correct. Therefore, we assign 1 to an edge with the correct direction and 0.5 to an edge with the wrong direction or no direction (Badsha et al., 2018).
}

\value{
A \link{list} of object that containing the following:
\itemize{
\item \code{Matrix}: Results store for TP and FP
\item \code{TP}: Total found edges in the inferred graph and edge exists in the true graph.
\item \code{FP}: Total found edges in the inferred graph but no edge exists in the true graph.
\item \code{Recall}: Power, or sensitivity measures how many edges from the true graph a method can recover.
\item \code{Precision}: Measures how many correct edges are recovered in the inferred graph.
}

}

\references{
1.Md. Bahadur Badsha, Audrey Qiuyan Fu: Learning causal biological networks with the principle of Mendelian randomization (Nat Comn; under reveiw 2018), 

}
\author{
Md Bahadur Badsha (\link{mdbadsha@uidaho.edu})
}


\seealso{
\link{aSHD}: adjusted Structural Hamming Distance (aSHD)
}
\examples{

# True model
# True graph (V1 --> T1 --> T2 --> T3)
tarmat_s1 <- matrix(0,
                    nrow = 4,
                    ncol = 4)
                    
colnames(tarmat_s1) <- c("V1", "T1", "T2", "T3")

rownames(tarmat_s1) <- colnames(tarmat_s1)

# Create an adjacency matrix for the true graph
tarmat_s1[1, 2] <- 1
tarmat_s1[2, 3] <- 1
tarmat_s1[3, 4] <- 1

# Graph object of the true graph
Truth <- as(tarmat_s1,
            "graphNEL")

# Inferred graph (V1 --> T1 <-- T2 --> T3)
tarmat_s2 <- matrix(0,
                    nrow = 4,
                    ncol = 4)
                    
colnames(tarmat_s2) <-c ("V1", "T1", "T2", "T3")

rownames(tarmat_s2) <- colnames(tarmat_s2)

# Create an adjacency matrix for the inferred graph
tarmat_s2[1, 2] <- 1
tarmat_s2[3, 2] <- 1
tarmat_s2[3, 4] <- 1

# Graph objects for the inferred graph
Inferred <- as(tarmat_s2,
               "graphNEL")
         
#Recall and Precision
Recall_Precision(Truth,
                 Inferred,
                 GV = 1,
                 edge.presence = 1.0,
                 edge.direction = 0.5)
}      

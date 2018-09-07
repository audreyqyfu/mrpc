\name{seqDiff}
\alias{seqDiff}

\title{
Deviation of the Two Sequnces
}
\description{
Deviation of the inferred graph from the true graph.
}
\usage{
seqDiff(g1, g2)
}

\arguments{
  \item{g1}{
Adjacency matrix from the first graph object.
}
  \item{g2}{
Adjacency matrix from the second graph object.
}
}

\author{
Md Bahadur Badsha <mdbadsha@uidaho.edu>, Evan A Martin <mart9986@vandals.uidaho.edu>, and Audrey Qiuyan Fu <audreyf@uidaho.edu>
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

# Deviation of the inferred graph from the true graph.
seqDiff(tarmat_s2,
        tarmat_s1)
         
}

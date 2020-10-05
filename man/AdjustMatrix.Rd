\name{AdjustMatrix}
\alias{AdjustMatrix}

\title{
Adjust the columns of the input matrix same as in the reference matrix
}

\description{
We adjusted the columns of the input matrix to have the same ordering as in the reference matrix.
}

\usage{
AdjustMatrix(reference, input)
}

\arguments{
  \item{reference}{The reference matrix.}
  \item{input}{The input matrix.}
}

\value{
Matrix
}

\author{
Md Bahadur Badsha (mbbadshar@gmail.com)
}

\examples{

# Reference matrix 
                  
reference <- matrix(0,nrow=4,ncol = 4)
colnames(reference) <- c("V1","T1","T2","T3")
rownames(reference) <- colnames(reference)
   
# Adjacency matrix for reference  
reference[1,2] <- 1
reference[2,3] <- 1
reference[3,4] <- 1

# Input matrix 
input <- matrix(0,
                nrow = 4,
                ncol = 4)
colnames(input) <- c("V1","T2","T3","T1")
rownames(input) <- colnames(input)
   
# Adjacency matrix for input
input[1,2] <- 1
input[2,3] <- 1
input[3,4] <- 1

# Adjust the columns of the input matrix same as in the reference matrix
AdjustMatrix <- AdjustMatrix (reference, input)

}


\name{mpinv}
\alias{mpinv}

\title{
Function of Calculate Inverse Matrix
}
\description{
Calculate inverse of the matrix, when the matrix is not square. This function is used to calculate the robust correlation matrix.
}
\usage{
mpinv(X)
}

\arguments{
  \item{X}{
Data Matrix
}
}

\value{
Matrix
}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}
\examples{
Inversematrix <- mpinv(simu.data_M0)
}

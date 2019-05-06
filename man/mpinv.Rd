\name{mpinv}
\alias{mpinv}

\title{
Calculate the inverse matrix
}
\description{
This function calculates the inverse of the non-square matrix as part of the calculation of the robust correlation matrix.
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
Inversematrix <- mpinv(simu_data_M0)
}

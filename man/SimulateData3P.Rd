\name{SimulateData3P}
\alias{SimulateData3P}

\title{
Simulate data for a node with three parents
}
\description{
Simulate data for a node with three parents
}
\usage{
  SimulateData3P(N, P1, P2, P3, b0.1, b1.1, b1.2, b1.3, sd.1)
}
\arguments{
  \item{N}{
    Number of observations.
  }
  \item{P1}{Data vector of the parent node, P1.}
  \item{P2}{Data vector of the parent node, P2.}
  \item{P3}{Data vector of the parent node, P3.}
  \item{b0.1}{Intercept of b0.1 + b1.1*P1 + b1.2*P2 + b1.3*P3, where P1, P2, and P3 are the parents of the corresponding node.}
  \item{b1.1}{Slope of P1 for b0.1 + b1.1*P1 + b1.2*P2 + b1.3*P3, where P1, P2, and P3 are the parents of the corresponding node.}
  \item{b1.2}{Slope of P2 for b0.1 + b1.1*P1 + b1.2*P2 + b1.3*P3, where P1, P2, and P3 are the parents of the corresponding node.}
  \item{b1.3}{Slope of P3 for b0.1 + b1.1*P1 + b1.2*P2 + b1.3*P3, where P1, P2, and P3 are the parents of the corresponding node.}
  \item{sd.1}{Standard deviation for corresponding data generated node.}
}

\value{
  Vector
}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}
\seealso{
  \link{SimulateData} for simulated data generating function.
}
\examples{
Data3P <- SimulateData3P(N = 10^3, 
                         P1 = 1, 
                         P2 = 1, 
                         P3 = 1, 
                         b0.1 = 0, 
                         b1.1 = 1, 
                         b1.2 = 1, 
                         b1.3 = 1, 
                         sd.1 = 1)
}

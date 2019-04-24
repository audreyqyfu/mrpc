\name{SimulateData2P}
\alias{SimulateData2P}

\title{
Simulate data for a node with two parents
}
\description{
Simulate data for a node with two parents
}
\usage{
SimulateData2P(N, P1, P2, b0.1, b1.1, b1.2, sd.1)
}
\arguments{
  \item{N}{
Number of observations
}
\item{P1}{Data vector of the parent node, P1.}
\item{P2}{Data vector of the parent node, P2.}
\item{b0.1}{Intercept of b0.1 + b1.1*P1 + b1.2*P2, where P1 and P2 are the parents of the corresponding node.}
\item{b1.1}{Slope of P1 for b0.1 + b1.1*P1+ b1.2*P2, where P1 and P2 are the parents of the corresponding node.}
\item{b1.2}{Slope of P2 for b0.1 + b1.1*P1 + b1.2*P2, where P1 and P2 are the parents of the corresponding node.}
\item{sd.1}{Standard deviation for corresponding data generated nodes.}
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
SimulateData2P(N = 10^3, P1 = 1, P2 = 1, 
        b0.1 = 0, b1.1 = 1,
        b1.2 = 1, sd.1 = 1)
}

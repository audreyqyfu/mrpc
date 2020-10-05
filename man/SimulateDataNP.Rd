\name{SimulateDataNP}
\alias{SimulateDataNP}

\title{
Simulate data for a node with no parent
}
\description{
Simulate data for a node with no parent
}
\usage{
SimulateDataNP(N, b0.1, sd.1)
}
\arguments{
  \item{N}{
Number of observations
}
  \item{b0.1}{Intercept of the corresponding simulated node.}

\item{sd.1}{Standard deviation for corresponding data generated node.}
}

\value{
Vector
}

\author{
Md Bahadur Badsha (mbbadshar@gmail.com)
}
\seealso{
\link{SimulateData} for simulated data generating function.
}
\examples{
DataNP <- SimulateDataNP(N = 10^3, 
                       b0.1 = 0, 
                       sd.1 = 1)
}

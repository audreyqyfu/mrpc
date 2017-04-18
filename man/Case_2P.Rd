\name{Case_2P}
\alias{Case_2P}

\title{
Case for Two Parents Generating Simulated Data
}
\description{
The gene to have data generated has two parents.
}
\usage{
Case_2P(N,P1,P2,b0.1,b1.1,b1.2,sd.1)
}
\arguments{
  \item{N}{
Number of observations
}
\item{P1}{P1 is the 1st parent of the corresponding genes.}
\item{P2}{P2 is the 2nd parent of the corresponding genes.}
\item{b0.1}{Intercept of b0.1+b1.1*P1, where P1 is the 1st parent of the corresponding genes.}
\item{b1.1}{Slope of P1 for b0.1+b1.1*P1, where P1 is the 1st parent of the corresponding genes.}
\item{b1.2}{Slope of P2 for b0.1+b1.1*P1+b1.2*P2, where P2 is the 2nd parent of the corresponding genes.}
\item{sd.1}{Standard deviation for corresponding data generated genes.}
}

\value{
Vector
}

\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}
\seealso{
\link{SimulatedData} for simulated data generating function.
}
\examples{
Case_2P(N=10^3,P1=1,P2=1,b0.1=0,b1.1=1,b1.2=1,sd.1=1)
}

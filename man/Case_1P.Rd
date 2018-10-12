\name{Case_1P}
\alias{Case_1P}

\title{
Case for One Parent Generating Simulated Data
}
\description{
The gene to have data generated has only one parent.
}
\usage{
Case_1P(N, P1, b0.1, b1.1, sd.1)
}
\arguments{
  \item{N}{
Number of observations
}
\item{P1}{Data vector of the parent gene P1.}

\item{b0.1}{Intercept of b0.1 + b1.1*P1, where P1 is the parent of the corresponding gene.}

\item{b1.1}{Slope of P1 for b0.1 + b1.1*P1, where P1 is the parent of the corresponding gene.}

\item{sd.1}{Standard deviation for corresponding data generated genes.}
}

\value{
Vector
}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}
\seealso{
\link{SimulatedData} for simulated data generating function.
}
\examples{
Case_1P(N = 10^3, P1 = 1, b0.1 = 0, b1.1 = 1, sd.1 = 1)
}

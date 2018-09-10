\name{Case_3P}
\alias{Case_3P}

\title{
  Case for Three Parent Generating Simulated Data
}
\description{
  The gene to have data generated has three parents.
}
\usage{
  Case_3P(N, P1, P2, P3, b0.1, b1.1, b1.2, b1.3, sd.1)
}
\arguments{
  \item{N}{
    Number of observations.
  }
  \item{P1}{Data vector of the parent gene, P1.}
  \item{P2}{Data vector of the parent gene, P2.}
  \item{P3}{Data vector of the parent gene, P3.}
  \item{b0.1}{Intercept of b0.1 + b1.1*P1 + b1.2*P2 + b1.3*P3, where P1, P2, and P3 are the parents of the corresponding gene.}
  \item{b1.1}{Slope of P1 for b0.1 + b1.1*P1 + b1.2*P2 + b1.3*P3, where P1, P2, and P3 are the parents of the corresponding gene.}
  \item{b1.2}{Slope of P2 for b0.1 + b1.1*P1 + b1.2*P2 + b1.3*P3, where P1, P2, and P3 are the parents of the corresponding gene.}
  \item{b1.3}{Slope of P3 for b0.1 + b1.1*P1 + b1.2*P2 + b1.3*P3, where P1, P2, and P3 are the parents of the corresponding gene.}
  \item{sd.1}{Standard deviation for corresponding data generated gene.}
}

\value{
  Vector
}

\author{
Md Bahadur Badsha (\link{mdbadsha@uidaho.edu})
}
\seealso{
  \link{SimulatedData} for simulated data generating function.
}
\examples{
  Case_3P(N = 10^3, P1 = 1, P2 = 1, P3 = 1, 
          b0.1 = 0, b1.1 = 1, b1.2 = 1, 
          b1.3 = 1, sd.1 = 1)
}

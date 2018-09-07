\name{Case_2P}
\alias{Case_2P}

\title{
Case for Two Parents Generating Simulated Data
}
\description{
The gene to have data generated has two parents.
}
\usage{
Case_2P(N, P1, P2, b0.1, b1.1, b1.2, sd.1)
}
\arguments{
  \item{N}{
Number of observations
}
\item{P1}{Data vector of the parent gene, P1.}
\item{P2}{Data vector of the parent gene, P2.}
\item{b0.1}{Intercept of b0.1 + b1.1*P1 + b1.2*P2, where P1 and P2 are the parents of the corresponding gene.}
\item{b1.1}{Slope of P1 for b0.1 + b1.1*P1+ b1.2*P2, where P1 and P2 are the parents of the corresponding gene.}
\item{b1.2}{Slope of P2 for b0.1 + b1.1*P1 + b1.2*P2, where P1 and P2 are the parents of the corresponding gene.}
\item{sd.1}{Standard deviation for corresponding data generated genes.}
}

\value{
Vector
}

\author{
Md Bahadur Badsha <mdbadsha@uidaho.edu>, Evan A Martin <mart9986@vandals.uidaho.edu>, and Audrey Qiuyan Fu <audreyf@uidaho.edu>
}
\seealso{
\link{SimulatedData} for simulated data generating function.
}
\examples{
Case_2P(N = 10^3, P1 = 1, P2 = 1, 
        b0.1 = 0, b1.1 = 1,
        b1.2 = 1, sd.1 = 1)
}

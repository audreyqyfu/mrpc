\name{Case_3P}
\alias{Case_3P}

\title{
  Case for Three Parent Generating Simulated Data
}
\description{
  The gene to have data generated has three parents.
}
\usage{
  Case_3P(N,P1,P2,P3,b0.1,b1.1,b1.2,b1.3,sd.1)
}
\arguments{
  \item{N}{
    Number of observations.
  }
  \item{P1}{P1 is the 1st parent of the corresponding genes.}
  \item{P2}{P2 is the 2nd parent of the corresponding genes.}
  \item{P3}{P3 is the 3rd parent of the corresponding genes.}
  \item{b0.1}{Intercept of b0.1+b1.1*P1, where P1 is the 1st parent of the corresponding genes.}
  \item{b1.1}{Slope of P1 for b0.1+b1.1*P1, where P1 is the 1st parent of the corresponding genes.}
  \item{b1.2}{Slope of P2 for b0.1+b1.1*P1+b1.2*P2, where P2 is the 2nd parent of the corresponding genes.}
  \item{b1.3}{Slope of P3 for b0.1+b1.1*P1+b1.2*P2+b1.3*P3, where P3 is the 3rd parent of the corresponding genes.}
  \item{sd.1}{Standard deviation for corresponding data generated genes.}
}

\value{
  Vector
}

\author{
Md Bahadur Badsha <mdbadsha@uidaho.edu>,Evan A Martin <mart9986@vandals.uidaho.edu> and Audrey Qiuyan Fu <audreyf@uidaho.edu>
}
\seealso{
  \link{SimulatedData} for simulated data generating function.
}
\examples{
  Case_3P(N=10^3,P1=1,P2=1,P3=1,b0.1=0,b1.1=1,b1.2=1,b1.3=1,sd.1=1)
}

\name{Case_NP}
\alias{Case_NP}

\title{
Case for No Parent Generating Simulated Data
}
\description{
The gene to have data generated has no parents.
}
\usage{
Case_NP(N, b0.1, sd.1)
}
\arguments{
  \item{N}{
Number of observations
}
  \item{b0.1}{Intercept of the corresponding simulated gene.}

\item{sd.1}{Standard deviation for corresponding data generated gene.}
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
Case_NP(N = 10^3, b0.1 = 0, sd.1 = 1)
}

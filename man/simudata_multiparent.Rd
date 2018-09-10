\name{simu.data_multiparent}
\alias{simu.data_multiparent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data for Multiple Parent Model
}
\description{
Simulate data for the multiple parent model, where the parents can be eQTL/CNVs along with several gene expression nodes.
}

\value{
Matrix
}

\author{
Md Bahadur Badsha (\link{mdbadsha@uidaho.edu})
}

\seealso{
\link{SimulatedData} for simulated data generating function.
}
\examples{
## Multiple Parent Model
simu.data_multiparent <- SimulatedData(N = 10^3, p = 0.45,
                                       'multiparent', b0.1 = 0,
                                       b1.1 = 1, b1.2 = 1,
                                       b1.3 = 1, sd.1 = 1)
}

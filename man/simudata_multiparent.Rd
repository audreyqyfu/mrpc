\name{simu.data_multiparent}
\alias{simu.data_multiparent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data for Multiple Parents Model
}
\description{
Generating simulated data for multiple parents model, where patrents can be the eQTL/CNVs with several gene expression.
}

\value{
Matrix
}

\author{
Md Bahadur Badsha <mdbadsha@uidaho.edu>,Evan A Martin <mart9986@vandals.uidaho.edu> and Audrey Qiuyan Fu <audreyf@uidaho.edu>
}

\seealso{
\link{SimulatedData} for simulated data generating function.
}
\examples{
## Multiple Parents Model
simu.data_multiparent<-SimulatedData(N = 10^3,p = 0.45,
'multiparent',b0.1 = 0,b1.1 = 1,b1.2 = 1,b1.3=1,sd.1 = 1)
  }

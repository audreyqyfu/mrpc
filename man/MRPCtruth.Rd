\name{MRPCtruth}
\alias{MRPCtruth}

\title{
Truth Models for MRPC Examples
}
\description{
Truth of the five basic topologies and three common topologies in biology, multi-parents, star and layered.
}

\author{
Md Bahadur Badsha (\link{mdbadsha@uidaho.edu})
}

\examples{
data("MRPCtruth") #load data

par(mfrow = c(2, 4))
plot(MRPCtruth$M0,
     main = "Model0")
plot(MRPCtruth$M1,
     main = "Model1")
plot(MRPCtruth$M2,
     main = "Model2")
plot(MRPCtruth$M3,
     main = "Model3")
plot(MRPCtruth$M4,
     main = "Model4")
plot(MRPCtruth$Multiparent,
     main = "Multiparent")
plot(MRPCtruth$Star,
     main = "Star")
plot(MRPCtruth$Layered,
     main = "Layered")

}

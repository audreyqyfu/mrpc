\name{MRPCtruth}
\alias{MRPCtruth}

\title{
Graphs used in simulation
}
\description{
Topologies of the five basic models and three common graphs in biology: namely the multi-parent graph, the star graph and the layered graph.  See details in Badsha and Fu (2018).
}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\references{
Md Bahadur Badsha and Audrey Qiuyan Fu (2018). Learning causal biological networks with the 
principle of Mendelian randomization. bioRxiv doi:10.1101/171348.     
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

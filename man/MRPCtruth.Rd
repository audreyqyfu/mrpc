\name{MRPCtruth}
\alias{MRPCtruth}

\title{
Truth Models for MRPC Examples
}
\description{
Truth of the five basic topology and three common (multi-parents, star and layered) topology in biology
}

\author{
Md Bahadur Badsha <mdbadsha@uidaho.edu>,Evan A Martin <mart9986@vandals.uidaho.edu> and Audrey Qiuyan Fu <audreyf@uidaho.edu>
}

\examples{
data("MRPCtruth") #load data
#if (require(Rgraphviz)){
par (mfrow=c(2,4))
plot(MRPCtruth$M0,main="Model0")
plot(MRPCtruth$M1,main="Model1")
plot(MRPCtruth$M2,main="Model2")
plot(MRPCtruth$M3,main="Model3")
plot(MRPCtruth$M4,main="Model4")
plot(MRPCtruth$Multiparent,main="Multiparent")
plot(MRPCtruth$Star,main="Star")
plot(MRPCtruth$Layered,main="Layered")
#}
}
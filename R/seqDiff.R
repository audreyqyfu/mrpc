
seqDiff<-function(g1, g2){
  #Convert binary to decimal for g1 (truth)
  G1=as.vector(t(g1))
  G11<- paste(G1,collapse="")
  S1=compositions::unbinary(G11)
  #Convert binary to decimal for g2 (inferred)
  G2=as.vector(t(g2))
  G22<- paste(G2,collapse="")
  S2=compositions::unbinary(G22)
  #Deviation of the inferred from the truth
  seqDiff=S2-S1
  return(seqDiff)
}
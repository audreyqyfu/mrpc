
#Data function for two papernts (P1,P2)
SimulateData2P<- function(N, P1,P2, b0.1, b1.1, b1.2,sd.1) {

  t <- rnorm(n = N,
             mean = b0.1 + b1.1*P1+b1.2*P2,
             sd = sd.1)

  return(t)

}

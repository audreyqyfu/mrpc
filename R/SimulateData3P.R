
#Data function for three papernts (P1,P2,P3)
SimulateData3P <- function(N,P1,P2,P3, b0.1, b1.1, b1.2,b1.3,sd.1) {

  t <- rnorm(n = N,
             mean = b0.1 + b1.1*P1+b1.2*P2+b1.3*P3,
             sd = sd.1)

  return(t)

}

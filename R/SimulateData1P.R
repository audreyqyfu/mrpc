#
#Data function one papernt (P1)
SimulateData1P <- function(N, P1,b0.1, b1.1, sd.1) {

  t <- rnorm(n = N,
             mean = b0.1 + b1.1*P1,
             sd = sd.1)

  return(t)

}

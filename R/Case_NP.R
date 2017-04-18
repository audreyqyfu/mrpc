#Data function no papernts
Case_NP<- function(N, b0.1,sd.1) {

  t <- rnorm(n = N,
             mean = b0.1,
             sd = sd.1)

  return(t)

}

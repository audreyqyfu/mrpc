
# N, b0.1 and sd.1 are the parameters to simulate data for no parent.
# Likewise N, b0.1, b1.1 and sd.1 are the parameters to simulate data for one parent.
# Likewise N, b0.1, b1.1,b1.2 and sd.1 are the parameters to simulate data for two parents.
# Likewise N, b0.1, b1.1,b1.2,b1.3 and sd.1 are the parameters to simulate data for three parents.
# To specify the model remember to use quotes. Ex. if you want to generate data for model 0 you would type 'model0' into the function.
#For example for Model0=simulateData(N = 10^3,p = 0.45,seed = 5,'model0',b0.1 = 0,b1.1 = 1,b1.2 = 1, sd.1 = 1)

SimulatedData<- function(N, p,model,b0.1, b1.1, b1.2,b1.3, sd.1) {

  #set.seed(seed)

  V <- c(sample(c(0, 1, 2),
                size = N,
                replace = TRUE,
                prob = c((1 - p)^2,
                         2*p*(1 - p),
                         p^2)))
  switch(model,

         model0 = {
           T1<- Case_1P(N=N,
                        P1=V,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)

           T2<- Case_NP(N=N,
                        b0.1 = b0.1,
                        sd.1 = sd.1)

           return(data.frame(V,
                             T1 = T1,
                             T2 = T2
           ))

         },

         model1 = {
           T1<- Case_1P(N=N,
                        P1=V,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)
           T2<- Case_1P(N=N,
                        P1=T1,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)
           return(data.frame(V,
                             T1 = T1,
                             T2 = T2
           ))
         },
         model2 = {
           T2<- Case_NP(N=N,
                        b0.1 = b0.1,
                        sd.1 = sd.1)
           T1<- Case_2P(N=N,
                        P1=V,
                        P2=T2,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        b1.2 = b1.2,
                        sd.1 = sd.1)
           return(data.frame(V,
                             T1 = T1,
                             T2 = T2
           ))

         },
         model3 = {
           T1<- Case_1P(N=N,
                        P1=V,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)
           T2<- Case_1P(N=N,
                        P1=V,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)
           return(data.frame(V,
                             T1 = T1,
                             T2 = T2
           ))

         },
         model4 = {
           T1.a<- Case_1P(N=N,
                          P1=V,
                          b0.1 = b0.1,
                          b1.1 = b1.1,
                          sd.1 = sd.1)
           T2.a<- Case_2P(N=N,
                          P1=V,
                          P2=T1.a,
                          b0.1 = b0.1,
                          b1.1 = b1.1,
                          b1.2 = b1.2,
                          sd.1 = sd.1)

           T2.b<- Case_1P(N=N,
                          P1=V,
                          b0.1 = b0.1,
                          b1.1 = b1.1,
                          sd.1 = sd.1)
           T1.b<- Case_2P(N=N,
                          P1=V,
                          P2=T2.b,
                          b0.1 = b0.1,
                          b1.1 = b1.1,
                          b1.2 = b1.2,
                          sd.1 = sd.1)
           # CoinToss creates a vector of 0s and 1s to determine which observations from T1.a, T1.b, T2.a, and T2.b are put into  T1.4 and T2.4.
           coinToss <- rbinom(n = N,
                              size = 1,
                              prob = 0.5)
           T1 <- rep(0, N)

           T1[which(coinToss == 0)] <- T1.a[which(coinToss == 0)]

           T1[which(coinToss == 1)] <- T1.b[which(coinToss == 1)]

           T2 <- rep(0, N)

           T2[which(coinToss == 0)] <- T2.a[which(coinToss == 0)]

           T2[which(coinToss == 1)] <- T2.b[which(coinToss == 1)]

           return(data.frame(V,
                             T1 = T1,
                             T2 = T2
           ))

         },
         
         multiparent= {
           
           T2<- Case_NP(N=N,
                        b0.1 = b0.1,
                        sd.1 = sd.1)
           T3<- Case_NP(N=N,
                        b0.1 = b0.1,
                        sd.1 = sd.1)
           T1<- Case_3P(N=N,
                        P1=V,
                        P2=T2,
                        P3=T3,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        b1.2 = b1.2,
                        b1.3 = b1.3,
                        sd.1 = sd.1)
           return(data.frame(V,
                             T1 = T1,
                             T2 = T2,
                             T3 = T3))
         },
         starshaped= {
           
           T1<- Case_1P(N=N,
                        P1=V,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)
           T2<- Case_1P(N=N,
                        P1=T1,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)
           T3<- Case_1P(N=N,
                        P1=T1,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)
           T4<- Case_1P(N=N,
                        P1=T1,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)
           T5<- Case_1P(N=N,
                        P1=T1,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)
           return(data.frame(V,
                             T1 = T1,
                             T2 = T2,
                             T3 = T3,
                             T4 = T4,
                             T5 = T5))
         },
         layered= {
           T1<- Case_1P(N=N,
                        P1=V,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)
           T2<- Case_1P(N=N,
                        P1=V,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)
           T3<- Case_1P(N=N,
                        P1=V,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)
           T4<- Case_1P(N=N,
                        P1=T1,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)
           T5<- Case_2P(N=N,
                        P1=T1,
                        P2=T2,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        b1.2 = b1.2,
                        sd.1 = sd.1)
           T6<- Case_1P(N=N,
                        P1=T2,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)
           T7<- Case_1P(N=N,
                        P1=T3,
                        b0.1 = b0.1,
                        b1.1 = b1.1,
                        sd.1 = sd.1)
           return(data.frame(V,
                             T1 = T1,
                             T2 = T2,
                             T3 = T3,
                             T4 = T4,
                             T5 = T5,
                             T6 = T6,
                             T7 = T7))
         },
         
         stop("Model not included or missing"))
}

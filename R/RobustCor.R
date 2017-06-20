RobustCor <- function(xx,Beta, plot=FALSE)
{
  #xx=as.matrix(xx,na.action(xx))
  #xx=na.omit(xx)
  #xx<- replace(xx, is.na(xx), 0)
  if (any(is.na(xx))) #if any NA
  {
    set.seed(103)
    xx=mice(xx,m=5,printFlag = F) #Multivariate Imputation by Chained Equations (MICE)
    xx=complete(xx) #Creates imputed data sets
  }
  xx=as.matrix(xx)
  nn <- dim(xx)[1] # sample size
  mm <- dim(xx)[2]
  #####===============Initialization=================================
  Mx<-NULL; Median<-NULL; Dist<-NULL; Data0<-NULL
  for (ii in 1:dim(xx)[2])
  {
      Median[ii] <- median(xx[,ii],na.rm = T)
    }

  #{Mode <- hist(xx[,ii], nclass = 20, plot = F)
  # Mx[ii] <- mean(Mode$mids[which(Mode$count==max(Mode$count))])
  #}
  for (jj in 1:dim(xx)[1])
  {
    Dist[jj]<-sqrt(sum((xx[jj,]-Median)^2,na.rm = T))
    }
  for (kk in 1:dim(xx)[1])
  {
    if (Dist[kk] <= as.numeric(quantile(Dist, p=.5,na.rm = T)))
    Data0 <- rbind(Data0, xx[kk,])
    }
  #plot(Dist)
  #readline('Pls Enter')
  #####==========End Initialization================================
  #source("mpinv.r")
  Mo <-as.numeric(colMeans(Data0,na.rm = T))
  Vo <- as.matrix(cov(Data0,use = "pairwise.complete.obs"))
  DiffTol = 0.005;
  DiffNorm = +10000;
  Iter = 0;
  #cat("Test=",dim(Vo),"\n")
  ##Wx <- NULL
  while (DiffNorm > DiffTol)
  {
    Wx <- NULL
    for (j1 in 1:nn)
    { 
    xx=replace(xx,is.na(xx),0)  
    zo <- as.numeric(xx[j1,]-Mo)
    #zz <- bta*(t(zo)%*%solve(Vo)%*%zo)
    #if(det(t(Vo)%*%zo)<0.000005)
    zz <- Beta*(t(zo)%*%mpinv(Vo)%*%zo)
    #if(det(t(Vo)%*%zo)>=0.000005)
    #zz <- Beta*(t(zo)%*%solve(Vo)%*%zo)
    Wx[j1] <- exp(-zz/2)
    }
    M.new <- matrix(0,nrow =mm)
    V.new <- matrix(0,nrow=mm,ncol =mm)
    for (j2 in 1:nn)
    { M.new <- M.new + (Wx[j2]*xx[j2,])/sum(Wx)
    V.new <- V.new + (1+Beta)*(Wx[j2]*(xx[j2,]-Mo)%*%t(xx[j2,]-Mo))/sum(Wx)
    }
    #norm1 <- sqrt(sum(M.new^2))+sqrt(sum(V.new^2))
    #norm2 <- sqrt(sum((M.new-Mo)^2))+sqrt(sum((V.new-Vo)^2))
    #DiffNorm <- norm2/norm1
    DiffNorm <- sqrt(sum((M.new-Mo)^2))/mm + sqrt(sum((V.new-Vo)^2))/mm
    Mo = M.new
    Vo = V.new
    Iter = Iter + 1
    #cat("Iter,Diff=",c(Iter,DiffNorm),"\n")
  }
  Cor.coef<-cov2cor(V.new)
  Wt.weight <-Wx
  if (plot) plot(Wt.weight)
  return(list(RR=Cor.coef, M=M.new, V=V.new, Wt=Wt.weight))
}



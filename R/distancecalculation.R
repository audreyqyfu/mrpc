distancecalculation=function (g1, g2,NQ) 
{
  if (is(g1, "pcAlgo")) 
    g1 <- g1@graph
  if (is(g2, "pcAlgo")) 
    g2 <- g2@graph
  if (is(g1, "graphNEL")) {
    m1 <- wgtMatrix(g1, transpose = FALSE)  #need library(pcalg) for wgtMatrix
    m1[m1 != 0] <- 1
  }
  if (is(g2, "graphNEL")) {
    m2 <- wgtMatrix(g2, transpose = FALSE)
    m2[m2 != 0] <- 1
  }
  
  Distannce_matrix=matrix(0, nrow(m1),ncol(m1)) #Output matrix
  ind1=which(m1==1,arr.ind = T)  #pullout the edges from 1st graph
  ind2=which(m2==1,arr.ind = T)  #pullout the edges from 2nd graph
  
  #For 1st graph
  for (i in seq_len(nrow(ind1))) {
    x <- ind1[i, 1]  #1st node
    y <- ind1[i, 2]  #2nd node
    #If missing edges then penalty =1
    if((m1[x,y]==1 || m1[y,x]!=1)  & (m2[x,y]!=1 & m2[y,x]!=1) & Distannce_matrix[x,y]==0 & Distannce_matrix[y,x]==0 )
    {
      Distannce_matrix[x,y]=1 
    }
    #If missing direction then penalty =0.5
    if(m1[x,y]==1 & m1[y,x]!=1 & m2[x,y]==1 & m2[y,x]==1)
      
    {
      Distannce_matrix[x,y]=0.5 
    }
    
    #if((x<=NQ & y<=NQ) & m1[x,y]==1 & m1[y,x]!=1 & (m2[x,y]==1 || m2[y,x]==1))
    
    #{
    # Distannce_matrix[x,y]=0.5 
    #}
    
  }
  #For 2nd graph
  for (i in seq_len(nrow(ind2))) {
    x <- ind2[i, 1] #1st node 
    y <- ind2[i, 2] #2nd node
    #If missing edge then penalty =1  
    if((m2[x,y]==1 || m2[y,x]!=1)  & (m1[x,y]!=1 & m1[y,x]!=1) & Distannce_matrix[x,y]==0 & Distannce_matrix[y,x]==0)
    {
      Distannce_matrix[x,y]=1 
    }
    
    #If missing direction then penalty =0.5
    if(m2[x,y]==1 & m2[y,x]!=1 & m1[y,x]==1 & Distannce_matrix[x,y]==0 & Distannce_matrix[y,x]==0)
      
    {
      Distannce_matrix[x,y]=0.5 
    }
  }
  #
  if(any(Distannce_matrix[1:NQ,]==.5) || any(Distannce_matrix[,1:NQ]==.5))
  {
    Distannce_matrix[1:NQ,]=0
    Distannce_matrix[,1:NQ]=0  
  }
  
  
  Distance=sum(Distannce_matrix)  
  return(Distance)
}

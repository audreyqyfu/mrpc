#The SHD, as is implemented in the R package pcalg (Kalisch et al., 2012) 
#and bnlearn(Scutari, 2010), counts how many differences exist between two
#directed graphs. This distance is 1 if an edge exists in one graph but missing
#in the other, or if the direction of an edge is different in the two graphs. 
#The larger this distance, the more different the two graphs are. We adjusted 
#the SHD to reduce the penalty on the wrong direction of an edge to 0.5.Details
#in help(aSHD)
aSHD=function (g1, g2,GV) 
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
  #If ordering nodes
   if(any(colnames(m1)!=colnames(m2)))
   {
     Order_node=match(colnames(m1),colnames(m2))
     m2<-m2[Order_node,Order_node] #new
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
    
    #if((x<=GV & y<=GV) & m1[x,y]==1 & m1[y,x]!=1 & (m2[x,y]==1 || m2[y,x]==1))
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
  #if(GV>1 & (any(Distannce_matrix[1:GV,1:GV]==.5) || any(Distannce_matrix[,1:GV]==.5)))
  if(GV>1)
  {
    W=which(Distannce_matrix[1:GV,1:GV]==0.5,arr.ind = T)
    Distannce_matrix[W]=0 
  }
  
  Distance=sum(Distannce_matrix)  
  return(Distance)
}

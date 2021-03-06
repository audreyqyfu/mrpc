#recall and precision calculation in MRPC
#Recall = (# edges correctly identified in inferred graph) / (# edges in true graph);
#Precision = (# edges correctly identified in inferred graph) / (# edges in inferred graph).
#we assign edge.presence=1 to an edge with the correct direction 
#and edge.direction=0.5 to an edge with the wrong direction or no direction
#Details please see help(RecallPrecision)
RecallPrecision <- function (g1, g2, GV, includeGV, edge.presence = 1.0, edge.direction = 0.5) 
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
    Order_node <- match(colnames(m1),colnames(m2))
    m2 <- m2[Order_node,Order_node] #new
  }
  
  #Output matrix
  Evaluation_matrix <- matrix(0, nrow = 1, ncol = 2) 
  colnames(Evaluation_matrix) <- c("TP", "FP")
  #index of evaluation matrix
  TP <- 1
  FP <- 2
  #library(combinat)
  #ind1=t(combinat::combn(ncol(m1), 2))
  
  if (includeGV) {     #Condition for recall and precision with GV
    #Calculate the edges from the true graph
    m11 <- m1
    for (i in 1:nrow(m11))
    {
      for (j in 1:ncol(m11))
      {
        if(m11[i,j]==m11[j,i])
        {
          m11[i,j] <- 0
        }
      }
    }
    #NTE=Number of edges in the truth graph 
    NTE <- length(which(m11==1)) 
    
    #Calculate the edges from the inferred graph
    m22 <- m2
    for (i in 1:nrow(m22))
    {
      for (j in 1:ncol(m22))
      {
        if(m22[i,j]==m22[j,i])
        {
          m22[i,j] <- 0
        }
      }
    }
    #NIE=Number of edges in the inferred graph
    NIE <- length(which(m22==1)) 
    ind1 <- t(combn(ncol(m1), 2))
    for (i in seq_len(nrow(ind1))) {
      x <- ind1[i, 1]  #1st node
      y <- ind1[i, 2]  #2nd node
      #H0:No edge
      #Ha: edge exit
      #TP: Reject H0,when Ha is true (H0 is false)
      #Found edge in the inferred graph and edge exit in the true graph 
      #FP: Reject H0,when H0 is true 
      #Found edge in the inferred graph but no edge exit in the true graph 
      
      #if x-->y  in true graph and x-->y in infrred graph then TP=1
      if ((m1[x, y]==1 & m1[y, x]!=1) & (m2[x, y]==1 & m2[y, x]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.presence
      }
      #if y-->x  in true graph and y-->x in infrred graph then TP=1
      if ((m1[y, x]==1 & m1[x, y]!=1) & (m2[y, x]==1 & m2[x, y]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.presence
      } 
      #if x<-->ytrue graph and x<-->y in infrred graph then TP=1 
      
      if((m1[x, y]==1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.presence
      }
      #if x-->y true graph and x<--y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]!=1) & (m2[x, y]!=1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }
      
      #if x<--y true graph and x-->y in infrred graph then TP=0.5
      if((m1[x, y]!=1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }
      
      #if x<-->ytrue graph and x-->y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }
      
      #if x<-->ytrue graph and x<--y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]==1) & (m2[x, y]!=1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }
      
      #if x-->ytrue graph and x<-->y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]!=1) & (m2[x, y]==1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }
      
      #if x<--ytrue graph and x<-->y in infrred graph then TP=0.5
      if((m1[x, y]!=1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }
      
      #if x y (no edge) true graph and x--y in infrred graph then FP=1 
      if((m1[x, y]==0 & m1[y, x]==0) & (m2[x, y]==1 || m2[y, x]==1))
      {
        Evaluation_matrix[FP] <- Evaluation_matrix[FP] + 1.0
      }
      
      #Recall
      if (NTE != 0) {
        Recall <- Evaluation_matrix[TP]/NTE
      }
      
      else {
        Recall <-NA
      }
      
      
      if (NIE != 0 ) {
        Precision <- Evaluation_matrix[TP]/NIE
      }
      
      
      else {
        Precision <- NA
      }
      
     
    }
  }
  ####Without GV###
  else
  { #Condition for recall and precision without GV
    
    
    if(GV==0)
    {
      m11 <- m1  #since GV=0 
      m22 <- m2  #since GV=0 
    }
    else
    {
      m11 <- m1[-c(1:GV), -c(1:GV)]  #exclude GV
      m22 <- m2[-c(1:GV), -c(1:GV)]  #exclude GV
    }
    #Calculate the edges from the true graph
    for (i in 1:nrow(m11))
    {
      for (j in 1:ncol(m11))
      {
        if(m11[i,j]==m11[j,i])
        {
          m11[i,j] <- 0
        }
      }
    }
    #NTE=Number of edges in the truth graph 
    NTE <- length(which(m11==1)) 
    
    #Calculate the edges from the inferred graph
    #m22 <- m2[-c(1:GV),-c(1:GV)]
    for (i in 1:nrow(m22))
    {
      for (j in 1:ncol(m22))
      {
        if(m22[i,j]==m22[j,i])
        {
          m22[i,j] <- 0
        }
      }
    }
    #NIE=Number of edges in the inferred graph
    NIE <- length(which(m22==1)) 
    if(GV==0)
    {
      m1 <- m1  
      m2 <- m2
    }
    else
    {
      m1 <- m1[-c(1:GV), -c(1:GV)] #exclude GV
      m2 <- m2[-c(1:GV), -c(1:GV)] #exclude GV
    }
    #m1 <- m1[-c(1:GV),-c(1:GV)] #exclude GV
    #m2 <- m2[-c(1:GV),-c(1:GV)] #exclude GV
    
    ind1 <- t(combn(ncol(m1), 2))
    for (i in seq_len(nrow(ind1))) {
      x <- ind1[i, 1]  #1st node
      y <- ind1[i, 2]  #2nd node
      #H0:No edge
      #Ha: edge exit
      #TP: Reject H0,when Ha is true (H0 is false)
      #Found edge in the inferred graph and edge exit in the true graph 
      #FP: Reject H0,when H0 is true 
      #Found edge in the inferred graph but no edge exit in the true graph 
      
      #if x-->y  in true graph and x-->y in infrred graph then TP=1
      if ((m1[x, y]==1 & m1[y, x]!=1) & (m2[x, y]==1 & m2[y, x]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.presence
      }
      #if y-->x  in true graph and y-->x in infrred graph then TP=1
      if ((m1[y, x]==1 & m1[x, y]!=1) & (m2[y, x]==1 & m2[x, y]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.presence
      } 
      #if x<-->ytrue graph and x<-->y in infrred graph then TP=1 
      
      if((m1[x, y]==1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.presence
      }
      #if x-->y true graph and x<--y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]!=1) & (m2[x, y]!=1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }
      
      #if x<--y true graph and x-->y in infrred graph then TP=0.5
      if((m1[x, y]!=1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }
      
      #if x<-->ytrue graph and x-->y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }
      
      #if x<-->ytrue graph and x<--y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]==1) & (m2[x, y]!=1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }
      
      #if x-->ytrue graph and x<-->y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]!=1) & (m2[x, y]==1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }
      
      #if x<--ytrue graph and x<-->y in infrred graph then TP=0.5
      if((m1[x, y]!=1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }
      
      #if x y (no edge) true graph and x--y in infrred graph then FP=1 
      if((m1[x, y]==0 & m1[y, x]==0) & (m2[x, y]==1 || m2[y, x]==1))
      {
        Evaluation_matrix[FP] <- Evaluation_matrix[FP] + 1.0
      }
      
      #Recall
      
      if (NTE != 0) {
        Recall <- Evaluation_matrix[TP]/NTE
      }
      
      else {
        Recall <-NA
      }
      
      
      if (NIE != 0 ) {
        Precision <- Evaluation_matrix[TP]/NIE
      }
      
      
      else {
        Precision <- NA
      }
      
      
    }
  }
  
  #return(Evaluation_matrix)
  return(list(Matrix = Evaluation_matrix,
              TP = Evaluation_matrix[TP],
              FP = Evaluation_matrix[FP],
              NTE = NTE,
              NIE = NIE,
              Recall = Recall,
              Precision = Precision))
}
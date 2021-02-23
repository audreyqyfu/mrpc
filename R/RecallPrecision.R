#recall and precision calculation in MRPC
#Recall = (# edges correctly identified in inferred graph) / (# edges in true graph);
#Precision = (# edges correctly identified in inferred graph) / (# edges in inferred graph).
#we assign edge.presence=1 to an edge with the correct direction 
#and edge.direction=0.5 to an edge with the wrong direction or no direction
#Details please see help(RecallPrecision)


RecallPrecision <- function (g1, g2, GV, includeGV, edge.presence = 1, edge.direction = 0.5) 
{
  if (is(g1, "pcAlgo")) 
    g1 <- g1@graph
  if (is(g2, "pcAlgo")) 
    g2 <- g2@graph
  if (is(g1, "graphNEL")) {
    m1 <- wgtMatrix(g1, transpose = FALSE)
    m1[m1 != 0] <- 1
  }
  if (is(g2, "graphNEL")) {
    m2 <- wgtMatrix(g2, transpose = FALSE)
    m2[m2 != 0] <- 1
  }
  if (any(colnames(m1) != colnames(m2))) {
    Order_node <- match(colnames(m1), colnames(m2))
    m2 <- m2[Order_node, Order_node]
  }
  Evaluation_matrix <- matrix(0, nrow = 1, ncol = 2)
  colnames(Evaluation_matrix) <- c("TP", "FP")
  TP <- 1
  FP <- 2
  if (includeGV) {
    m11 <- m1
    for (i in 1:nrow(m11)) {
      for (j in 1:ncol(m11)) {
        if (m11[i, j] == m11[j, i]) {
          m11[i, j] <- 0
        }
      }
    }
    NTE <- length(which(m11 == 1))
    m22 <- m2
    for (i in 1:nrow(m22)) {
      for (j in 1:ncol(m22)) {
        if (m22[i, j] == m22[j, i]) {
          m22[i, j] <- 0
        }
      }
    }
    NIE <- length(which(m22 == 1))
    
    
    
    ind1 <- t(combn(ncol(m1), 2))
    for (i in seq_len(nrow(ind1))) {
      x <- ind1[i, 1]
      y <- ind1[i, 2]
      if ((m1[x, y] == 1 & m1[y, x] != 1) & (m2[x, y] == 
                                             1 & m2[y, x] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.presence
      }
      if ((m1[y, x] == 1 & m1[x, y] != 1) & (m2[y, x] == 
                                             1 & m2[x, y] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.presence
      }
      if ((m1[x, y] == 1 & m1[y, x] == 1) & (m2[x, y] == 
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.presence
      }
      if ((m1[x, y] == 1 & m1[y, x] != 1) & (m2[x, y] != 
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.direction
      }
      if ((m1[x, y] != 1 & m1[y, x] == 1) & (m2[x, y] == 
                                             1 & m2[y, x] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.direction
      }
      if ((m1[x, y] == 1 & m1[y, x] == 1) & (m2[x, y] == 
                                             1 & m2[y, x] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.direction
      }
      if ((m1[x, y] == 1 & m1[y, x] == 1) & (m2[x, y] != 
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.direction
      }
      if ((m1[x, y] == 1 & m1[y, x] != 1) & (m2[x, y] == 
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.direction
      }
      if ((m1[x, y] != 1 & m1[y, x] == 1) & (m2[x, y] == 
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.direction
      }
      if ((m1[x, y] == 0 & m1[y, x] == 0) & (m2[x, y] == 
                                             1 || m2[y, x] == 1)) {
        Evaluation_matrix[FP] <- Evaluation_matrix[FP] + 
          1
      }
      
      
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
  else {
    if (GV == 0) {
      m11 <- m1
      m22 <- m2
    }
    else {
      m11 <- m1[-c(1:GV), -c(1:GV)]
      m22 <- m2[-c(1:GV), -c(1:GV)]
    }
    for (i in 1:nrow(m11)) {
      for (j in 1:ncol(m11)) {
        if (m11[i, j] == m11[j, i]) {
          m11[i, j] <- 0
        }
      }
    }
    NTE <- length(which(m11 == 1))
    for (i in 1:nrow(m22)) {
      for (j in 1:ncol(m22)) {
        if (m22[i, j] == m22[j, i]) {
          m22[i, j] <- 0
        }
      }
    }
    NIE <- length(which(m22 == 1))
    if (GV == 0) {
      m1 <- m1
      m2 <- m2
    }
    else {
      m1 <- m1[-c(1:GV), -c(1:GV)]
      m2 <- m2[-c(1:GV), -c(1:GV)]
    }
    ind1 <- t(combn(ncol(m1), 2))
    
    ########################
    for (i in seq_len(nrow(ind1))) {
      x <- ind1[i, 1]
      y <- ind1[i, 2]
      if ((m1[x, y] == 1 & m1[y, x] != 1) & (m2[x, y] == 
                                             1 & m2[y, x] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.presence
      }
      if ((m1[y, x] == 1 & m1[x, y] != 1) & (m2[y, x] == 
                                             1 & m2[x, y] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.presence
      }
      if ((m1[x, y] == 1 & m1[y, x] == 1) & (m2[x, y] == 
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.presence
      }
      if ((m1[x, y] == 1 & m1[y, x] != 1) & (m2[x, y] != 
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.direction
      }
      if ((m1[x, y] != 1 & m1[y, x] == 1) & (m2[x, y] == 
                                             1 & m2[y, x] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.direction
      }
      if ((m1[x, y] == 1 & m1[y, x] == 1) & (m2[x, y] == 
                                             1 & m2[y, x] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.direction
      }
      if ((m1[x, y] == 1 & m1[y, x] == 1) & (m2[x, y] != 
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.direction
      }
      if ((m1[x, y] == 1 & m1[y, x] != 1) & (m2[x, y] == 
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.direction
      }
      if ((m1[x, y] != 1 & m1[y, x] == 1) & (m2[x, y] == 
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + 
          edge.direction
      }
      if ((m1[x, y] == 0 & m1[y, x] == 0) & (m2[x, y] == 
                                             1 || m2[y, x] == 1)) {
        Evaluation_matrix[FP] <- Evaluation_matrix[FP] + 
          1
      }
      
      
      if (NTE != 0 ) {
        Recall <- Evaluation_matrix[TP]/NTE
      }
      
      
      else {
        Recall <-NA
      }
      
      
      if (NIE != 0) {
        Precision <- Evaluation_matrix[TP]/NIE
      }
      
      
      else {
        Precision <- NA
      }
      
 
    }
  }
  return(list(Matrix = Evaluation_matrix, TP = Evaluation_matrix[TP], 
              FP = Evaluation_matrix[FP], NTE = NTE, NIE = NIE, Recall = Recall, 
              Precision = Precision))
}


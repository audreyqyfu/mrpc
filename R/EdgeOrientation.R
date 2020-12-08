#This is the step 2 of MRPC to direction determination of the undirected edges.

EdgeOrientation <- function (gInput, GV, suffStat, FDR, alpha, indepTest,
                             FDRcontrol, tau = 0.5, lambda = 0.25,
                             verbose = FALSE) {
  
  g <- as(gInput@graph, "matrix") # g ia an adjacency from undirected graph (skleton)
  g1 <- g
  p <- nrow(g)
  
  # tarmat (adjacency matrix for directed graph) is updated in every step, 
  # and contains the output of final topology
  tarmat <- matrix(0, nrow(g),ncol(g)) #same row and column from g
  rownames(tarmat) <- rownames(g)     #same row names from g
  colnames(tarmat) <- colnames(g)     #same column names from g
  
  # extract all edges
  g[lower.tri(g)] <- 0  # Use only upper triangular because g is symmetric matrix
  edges <- which (g==1, arr.ind = TRUE)
  
  #Step-1 start
  if (GV>0) {
    # identify edges involving Vs
    edgesWithBothVs <- edges[which (edges[,1] <= GV & edges[,2] <= GV), ]
    edgesWithVFirst <- edges[which (edges[,1] <= GV & edges[,2] > GV), ]
    edgesWithVSecond <- edges[which (edges[,1] > GV & edges[,2] <= GV), ]
    
    # assign 1s to corresponding edges in tarmat
    # edges between two Vs are undirected (or bidirected)
    if (length (edgesWithBothVs) > 2) {
      tarmat[edgesWithBothVs] <- 1
      tarmat[edgesWithBothVs[,2:1]] <- 1
    } else {
      tarmat[edgesWithBothVs[1], edgesWithBothVs[2]] <- 1
      tarmat[edgesWithBothVs[2], edgesWithBothVs[1]] <- 1
    }
    # Edges involving one V go from V to the other node
    # Need to distinguish whether edgesWithVFirst or edgesWithVSecond is 
    # a matrix or a vector.
    if (length (edgesWithVFirst) > 2) {
      tarmat[edgesWithVFirst] <- 1      
    } else {
      tarmat[edgesWithVFirst[1], edgesWithVFirst[2]] <- 1
    }
    
    if (length(edgesWithVSecond) > 2) {
      tarmat[edgesWithVSecond[,2:1]] <- 1
    } else {
      tarmat[edgesWithVSecond[2], edgesWithVSecond[1]] <- 1
    }
  }
  #Step-1 end
  
  #Step-2 start
  #Start to orient v-structures
  if (verbose)
  cat("\n V-structures are as follows :\n")
  # extract edges involving at least one gene node
  #edgesWithGs <- edges[which (edges[,1]>GV | edges[,2]>GV), ]
  #ind <- edgesWithGs
  
  m <- gInput@test  #current test
  Alpha <- gInput@alpha  #alpha
  pval <- gInput@pval
  R <- gInput@R          #decision of test
  ind <- which(g1 == 1, arr.ind = TRUE)  #Pullout the all relation in adjacency matrix from undirected graph
  V <- colnames(g)
  
  # Update the addis vectors/values output from the ModiSkeleton function.
  w0 <- tau * lambda * FDR/2
  alphai <- gInput@alphai
  kappai <- gInput@kappai
  kappai_star <- gInput@kappai_star
  K <- gInput@K
  Ci <- gInput@Ci
  Ci_sum <- sum(Ci)
  Si <- gInput@Si
  Si_sum <- sum(Si)
  Ci_plus <- gInput@Ci_plus
  gammai <- gInput@gammai
  normalizer <- gInput@normalizer
  exponent <- gInput@exponent
  
  for (i in seq_len(nrow(ind))) {
    x <- ind[i, 1]
    y <- ind[i, 2]
    
    allZ <- setdiff(which(g1[y, ] == 1), x)
    for (z in allZ) {
      # Triplet x-y-z is directed x-->y<--z if x and z conditionally dependent given y
      if ((g1[x, z] == 0 & g1[x, y] == 1 & g1[y, z] == 1)  & !(tarmat[y, x] ==1) & !(tarmat[z, y] ==1) & !(tarmat[y, z] ==1) & 
          !(y %in% gInput@sepset[[x]][[z]] || y %in% gInput@sepset[[z]][[x]])) 
      {
        m <- m+1
        
        # Increase the length of the R, pval, and other addis vectors if
        # their length is greater than or equal to the current iteration.
        if (m >= length(R)) {
          
          R <- c(R, numeric(m / 2))
          pval <- c(pval, numeric(m / 2))
          kappai <- c(kappai, numeric(m / 2))
          kappai_star <- c(kappai_star, numeric(m / 2))
          Ci <- c(Ci, numeric(m / 2))
          Si <- c(Si, numeric(m / 2))
          Ci_plus <- c(Ci_plus, numeric(m / 2))
          gammai <- c(gammai, numeric(m / 2))
          alphai <- c(alphai, numeric(m / 2))
          
        }
        
        if(indepTest=="gaussCItest") #if indepTest=gaussCItest
        {
          pval[m] <- gaussCItest(x, z, y, suffStat)
        }
        if(indepTest=="disCItest") #if indepTest=disCItest
        {
          pval[m] <- disCItest(x, z, y, suffStat) #additional
        }
        #pval=disCItest(x, z, y, suffStat) #additional conditional test
        if (FDRcontrol == 'LOND') {
          
          # Calculate alpha using the LOND method.
          alphai[m] <- SeqFDR(m,FDR,a=2,R)
          
          Alpha <- alphai[m]
          
          # Calculate the alpha value using the ADDIS algorithm
        } else if (FDRcontrol == 'ADDIS') {
          
          # Run ADDIS on the current iteration and update all the vectors
          # and other values.
          run_addis <- addis(alpha = FDR,
                             tau = tau,
                             lambda = lambda,
                             iter = m,
                             w0 = w0,
                             pval = pval,
                             alphai = alphai,
                             gammai = gammai,
                             kappai = kappai,
                             kappai_star = kappai_star,
                             K = K,
                             Ci = Ci,
                             Si = Si,
                             Ri = R,
                             Ci_plus = Ci_plus,
                             Ci_sum = Ci_sum,
                             Si_sum = Si_sum,
                             normalizer = normalizer,
                             exponent = exponent)
          
          # Update all values and vectors output from the addis function.
          alphai[[m]] <- run_addis[[1]]
          gammai[[m]] <- run_addis[[2]]
          K <- run_addis[[3]]
          R[[m]] <- run_addis[[5]]
          Si[[m - 1]] <- run_addis[[6]]
          Ci[[m - 1]] <- run_addis[[7]]
          Ci_sum <- run_addis[[9]]
          Si_sum <- run_addis[[10]]
          
          # Only update the kappai and Ci_plus vectors if K is greater
          # than one. If K is zero then the first element in kappai will
          # remain zero until the first rejection.
          if (K != 0) {
            
            kappai[[K]] <- run_addis[[4]]
            kappai_star[[K]] <- run_addis[[11]]
            Ci_plus[1:K] <- run_addis[[8]]
            
          }
          
          # Update the Alpha value with the current alphai value.
          Alpha <- alphai[[m]]
          
        } else {
          
          Alpha <- alpha
          
        }
        if (verbose){
        cat("x=", x, " y=", z, " S=", y,"\n")
        cat("Test number =", m, "\n")
        cat("Additional pval value =", pval[m], "\n")
        cat("Alpha value =", Alpha, "\n")
        }
        if (pval[m] <= Alpha) {  #Reject H0 (H0:nodes are independent)
          R[m] <- 1
          if (verbose) {
        cat(V[x], "->", V[y], "<-", V[z], "\n")  #Printout the v-structures
        cat("Since pval<Alpha,additional test is rejected;", "Nodes", V[x] ,"and" ,V[z] ,"are dependent given", V[y], "\n")
          }
          tarmat[x, y] <- tarmat[z, y] <- 1 #directed x-->y<--z
        } 
        else {
          R[m] <- 0  #Accept H0
          #tarmat[x, y] <- tarmat[y, z] <- 1 #directed x-->y-->z
          if (verbose)
          cat("Since pval>Alpha,additional test is accepted;", "Nodes", V[x] ,"and" ,V[z] ,"are independent given", V[y], "\n")
        }
        
      }
    }
  }
  #Step-2 end to orient v-structures
  
  #Step-3
  #Orient remaining edges, weather gMR is applicable or not based on the following process:
  #Pullout all edges from g1 matrix and define 1st and 2nd node (goal to make a triplet for the remaining edges)
  #extract all edges with involving those 1st and 2nd nodes (from row and column of g)
  #ignore the nodes that already have direction
  #form a triplet
  #determine the direction using test results from part 1 based on gMR
  #Repeat until all undirected edges to directed
  # m <- m
  # Alpha <- Alpha
  # R <- R
  #start 
  #when data contain genetic variants
  if (any(tarmat == 1) & GV>0) #if at least one edge directed already made so far
  {
    WW1 <- unique(which(tarmat==1,arr.ind = T)[,2]) #pullout the all canidate genes for v-structure
    #if(WW1<GV){
    #WW1=WW1[-c(1:GV)] #ignor canidate genes for GV
    WW1 <- setdiff(WW1,1:GV)
    #}
    if(length(WW1)!=0)
    {
      #WW1=setdiff(WW1,1:GV)
      for (v1 in 1:length(WW1))
      {
        WW2 <- unique(which(tarmat[,WW1[v1]]==1,arr.ind = T)) #edge between v-structure  
        WW3 <- unique(which(g1[,WW1[v1]]==1,arr.ind = T)) #edges between canidate genes and others  
        WW3 <- setdiff(WW3,WW2)   #ignore if already direction
        if(length(WW3)!=0)
        {
          for (v2 in 1:length(WW3))
          {
            x <- WW2[1]
            y <- WW1[v1]
            z <- WW3[v2]
            
            
            if ((!is.na(x) & !is.na(y) & !is.na(z)) & (x!=0 & y!=0 & z!=0) & (x!=y & y!=z & z!=x))
            {
              #Case-2: If,y and z are adjacent, x and z conditionally independent given y,
              #then the edge direction will be y-->z
              if (g1[y, z] == 1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 & ((y %in% gInput@sepset[[x]][[z]]) || (y %in% gInput@sepset[[z]][[x]])))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 0
              }
              #Case-3: If,y and z are adjacent, x and z conditionally dependent given y,
              #then the edge direction will be z-->y.
              if (g1[y, z] == 1 & g1[x, z] != 1 & tarmat[x,y]==1 & tarmat[y, x]!=1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 &!(y %in% gInput@sepset[[x]][[z]]) & !(y %in% gInput@sepset[[z]][[x]]))
                
              {
                m <- m+1
                
                # Increase the length of the R, pval, and other addis vectors if
                # their length is greater than or equal to the current iteration.
                if (m >= length(R)) {
                  
                  R <- c(R, numeric(m / 2))
                  pval <- c(pval, numeric(m / 2))
                  kappai <- c(kappai, numeric(m / 2))
                  kappai_star <- c(kappai_star, numeric(m / 2))
                  Ci <- c(Ci, numeric(m / 2))
                  Si <- c(Si, numeric(m / 2))
                  Ci_plus <- c(Ci_plus, numeric(m / 2))
                  gammai <- c(gammai, numeric(m / 2))
                  alphai <- c(alphai, numeric(m / 2))
                  
                }
                
                if(indepTest=="gaussCItest") #if indepTest=gaussCItest for continuous data
                {
                  pval[m] <- gaussCItest(x, z, y, suffStat) #additional pval
                }
                if(indepTest=="disCItest") #if indepTest=disCItest for discrete data 
                {
                  pval[m] <- disCItest(x, z, y, suffStat) #additional pval
                }
                
                if (FDRcontrol == 'LOND') {
                  
                  # Calculate alpha using the LOND method.
                  alphai[m] <- SeqFDR(m,FDR,a=2,R)
                  
                  Alpha <- alphai[m]
                  
                  # Calculate the alpha value using the ADDIS algorithm
                } else if (FDRcontrol == 'ADDIS') {
                  
                  # Run ADDIS on the current iteration and update all the vectors
                  # and other values.
                  run_addis <- addis(alpha = FDR,
                                     tau = tau,
                                     lambda = lambda,
                                     iter = m,
                                     w0 = w0,
                                     pval = pval,
                                     alphai = alphai,
                                     gammai = gammai,
                                     kappai = kappai,
                                     kappai_star = kappai_star,
                                     K = K,
                                     Ci = Ci,
                                     Si = Si,
                                     Ri = R,
                                     Ci_plus = Ci_plus,
                                     Ci_sum = Ci_sum,
                                     Si_sum = Si_sum,
                                     normalizer = normalizer,
                                     exponent = exponent)
                  
                  # Update all values and vectors output from the addis function.
                  alphai[[m]] <- run_addis[[1]]
                  gammai[[m]] <- run_addis[[2]]
                  K <- run_addis[[3]]
                  R[[m]] <- run_addis[[5]]
                  Si[[m - 1]] <- run_addis[[6]]
                  Ci[[m - 1]] <- run_addis[[7]]
                  Ci_sum <- run_addis[[9]]
                  Si_sum <- run_addis[[10]]
                  
                  # Only update the kappai and Ci_plus vectors if K is greater
                  # than one. If K is zero then the first element in kappai will
                  # remain zero until the first rejection.
                  if (K != 0) {
                    
                    kappai[[K]] <- run_addis[[4]]
                    kappai_star[[K]] <- run_addis[[11]]
                    Ci_plus[1:K] <- run_addis[[8]]
                    
                  }
                  
                  # Update the Alpha value with the current alphai value.
                  Alpha <- alphai[[m]]
                  
                } else {
                  
                  Alpha <- alpha
                  
                }
                if (verbose){
                  cat("Additional pval value =", pval[m], "\n")
                  cat("Alpha value =", Alpha, "\n")    
                }

                if (pval[m] <= Alpha) {  #Reject H0 (H0:nodes are independent)
                  R[m] <- 1
                  if (verbose) {
                    V <- colnames(g)
                    cat("\n", V[x], "->", V[y], "<-", V[z], "\n")  #Printout the v-structures
                    cat("Since pval<Alpha,additional test is rejected", "Nodes", V[x] ,"and" ,V[z] ,"are dependent given", V[y], "\n")
                  }
                  tarmat[z, y] <- 1 #directed z-->y
                } 
                else {
                  tarmat[y, z] <- 1
                  R[m] <- 0  #Accept H0
                }
              }
              #Case-4:If, y & z have relation and x & y conditionally dependent given z and x & z conditionally independent given y,
              #then edge direction will be z<-->y
              if (g1[y, z] == 1 & g1[x, z] == 1 & g1[x, y] == 1 & !(z %in% gInput@sepset[[x]][[y]]) &!(y %in% gInput@sepset[[x]][[z]]))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 1
              }
            }
            
          }
          
        }
      }
      #if(length(WW3)!=0)
      #{
      #Orient remaining edges,
      T1 <- which(tarmat==1,arr.ind = T)
      G1 <- which(g==1,arr.ind = T)
      D1 <- setdiff(G1,T1)
      #
      if(length(D1)!=0)
      {
        for (d1 in 1:length(D1))
        {
          Rem1_row <- which(g[,D1[d1]]==1,arr.ind = T)
          Rem1_col <- which(g[D1[d1],]==1,arr.ind = T)
          Rem1 <- c(Rem1_row,Rem1_col)
          Rem2 <- which(tarmat[,Rem1]==1,arr.ind = T)
          if(length(Rem1)!=0 & length(Rem2)!=0)
          {
            for (d11 in 1:length(Rem1))
            {
              x <- Rem2[d11]
              y <- Rem1[d11]
              z <- D1[d1]
            }
            
            
            if ((!is.na(x) & !is.na(y) & !is.na(z)) & (x!=0 & y!=0 & z!=0) & (x!=y & y!=z & z!=x))
            {
              #Case-2: If,y and z are adjacent, x and z conditionally independent given y,
              #then the edge direction will be y-->z
              if (g1[y, z] == 1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 & ((y %in% gInput@sepset[[x]][[z]]) || (y %in% gInput@sepset[[z]][[x]])))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 0
              }
              
              #Case-3: If,y and z are adjacent, x and z conditionally dependent given y,
              #then the edge direction will be z-->y.
              if (g1[y, z] == 1 & g1[x, z] != 1 & tarmat[x,y]==1 & tarmat[y, x]!=1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 &!(y %in% gInput@sepset[[x]][[z]]) & !(y %in% gInput@sepset[[z]][[x]]))
                
              {
                m <- m+1
                
                # Increase the length of the R, pval, and other addis vectors if
                # their length is greater than or equal to the current iteration.
                if (m >= length(R)) {
                  
                  R <- c(R, numeric(m / 2))
                  pval <- c(pval, numeric(m / 2))
                  kappai <- c(kappai, numeric(m / 2))
                  kappai_star <- c(kappai_star, numeric(m / 2))
                  Ci <- c(Ci, numeric(m / 2))
                  Si <- c(Si, numeric(m / 2))
                  Ci_plus <- c(Ci_plus, numeric(m / 2))
                  gammai <- c(gammai, numeric(m / 2))
                  alphai <- c(alphai, numeric(m / 2))
                  
                }
                
                if(indepTest=="gaussCItest") #if indepTest=gaussCItest
                {
                  pval[m] <- gaussCItest(x, z, y, suffStat)
                }
                if(indepTest=="disCItest") #if indepTest=gaussCItest
                {
                  pval[m] <- disCItest(x, z, y, suffStat) #additional
                }                
                if (FDRcontrol == 'LOND') {
                  
                  # Calculate alpha using the LOND method.
                  alphai[m] <- SeqFDR(m,FDR,a=2,R)
                  
                  Alpha <- alphai[m]
                  
                  # Calculate the alpha value using the ADDIS algorithm
                } else if (FDRcontrol == 'ADDIS') {
                  
                  # Run ADDIS on the current iteration and update all the vectors
                  # and other values.
                  run_addis <- addis(alpha = FDR,
                                     tau = tau,
                                     lambda = lambda,
                                     iter = m,
                                     w0 = w0,
                                     pval = pval,
                                     alphai = alphai,
                                     gammai = gammai,
                                     kappai = kappai,
                                     kappai_star = kappai_star,
                                     K = K,
                                     Ci = Ci,
                                     Si = Si,
                                     Ri = R,
                                     Ci_plus = Ci_plus,
                                     Ci_sum = Ci_sum,
                                     Si_sum = Si_sum,
                                     normalizer = normalizer,
                                     exponent = exponent)
                  
                  # Update all values and vectors output from the addis function.
                  alphai[[m]] <- run_addis[[1]]
                  gammai[[m]] <- run_addis[[2]]
                  K <- run_addis[[3]]
                  R[[m]] <- run_addis[[5]]
                  Si[[m - 1]] <- run_addis[[6]]
                  Ci[[m - 1]] <- run_addis[[7]]
                  Ci_sum <- run_addis[[9]]
                  Si_sum <- run_addis[[10]]
                  
                  # Only update the kappai and Ci_plus vectors if K is greater
                  # than one. If K is zero then the first element in kappai will
                  # remain zero until the first rejection.
                  if (K != 0) {
                    
                    kappai[[K]] <- run_addis[[4]]
                    kappai_star[[K]] <- run_addis[[11]]
                    Ci_plus[1:K] <- run_addis[[8]]
                    
                  }
                  
                  # Update the Alpha value with the current alphai value.
                  Alpha <- alphai[[m]]
                  
                } else {
                  
                  Alpha <- alpha
                  
                }
                if (pval[m] <= Alpha) {  #Reject H0 (H0:nodes are independent)
                  R[m] <- 1
                  
                  if (verbose) {
                    V <- colnames(g)
                    cat("\n", V[x], "->", V[y], "<-", V[z], "\n")  #Printout the v-structures
                    cat("Since pval<Alpha,additional test is rejected", "Nodes", V[x] ,"and" ,V[z] ,"are dependent given", V[y], "\n")
                  }
                  tarmat[z, y] <- 1 #directed z-->y
                } 
                else {
                  tarmat[y, z] <- 1
                  R[m] <- 0  #Accept H0
                }
              }
              #Case-4:If, y & z have relation and x & y conditionally dependent given z and x & z conditionally independent given y,
              #then edge direction will be z<-->y
              if (g1[y, z] == 1 & g1[x, z] == 1 & g1[x, y] == 1 & !(z %in% gInput@sepset[[x]][[y]]) &!(y %in% gInput@sepset[[x]][[z]]))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 1
              }
            }
          }
        }
      }
      #Check the reamning edge orientation
      ind <- which(g1 == 1, arr.ind = TRUE)  #Pullout the all relation in adjacency matrix from undirected graph
      for (i in seq_len(nrow(ind))) {
        x <- ind[i, 1]
        y <- ind[i, 2]
        
        if(tarmat[x, y]==0 & tarmat[y, x]==0) #bidirected if still no edge in tarmat
        {
          tarmat[x, y] <- 1
          tarmat[y, x] <- 1
        }
      }
    }
    
  }#end when data contain genetic variants
  #Start to orient remaining edges, when data not necessary to contain genetic variants.
  #First identify all the canidate genes for v-structure
  #Then identify the edges (undirected) with canidate genes 
  #ignore the nodes that already have direction
  #Make a triplet
  #determine the direction using test results from part 1 based on gMR
  #Repeat until all undirected edges to directed
  
  #Start when data not necessary to contain genetic variants.
  if (any(tarmat == 1) & GV==0) #if at least one edge directed already made so far
  {
    WW1 <- unique(which(tarmat == 1,arr.ind = T)[,2]) #pullout the all canidate genes for v-structure
    #WW1=WW1[-c(1:GV)]
    if(length(WW1)!=0)
    {
      #WW1=setdiff(WW1,1:GV)
      for (v1 in 1:length(WW1))
      {
        WW2 <- unique(which(tarmat[,WW1[v1]]==1,arr.ind = T)) #edge between v-structure  
        WW3 <- unique(which(g1[,WW1[v1]]==1,arr.ind = T)) #edges between canidate genes and others  
        WW3 <- setdiff(WW3,WW2)   #ignore if already direction
        if(length(WW3)!=0)
        {
          for (v2 in 1:length(WW3))
          {
            x <- WW2[1]
            y <- WW1[v1]
            z <- WW3[v2]
            
            if ((!is.na(x) & !is.na(y) & !is.na(z)) & (x!=0 & y!=0 & z!=0) & (x!=y & y!=z & z!=x))
            {
              #Case-2: If,y and z are adjacent, x and z conditionally independent given y,
              #then the edge direction will be y-->z
              if (g1[y, z] == 1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 & ((y %in% gInput@sepset[[x]][[z]]) || (y %in% gInput@sepset[[z]][[x]])))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 0
              }
              #Case-3: If,y and z are adjacent, x and z conditionally dependent given y,
              #then the edge direction will be z-->y.
              if (g1[y, z] == 1 & g1[x, z] != 1 & tarmat[x,y]==1 & tarmat[y, x]!=1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 &!(y %in% gInput@sepset[[x]][[z]]) & !(y %in% gInput@sepset[[z]][[x]]))
                
              {
                m <- m+1
                
                # Increase the length of the R, pval, and other addis vectors if
                # their length is greater than or equal to the current iteration.
                if (m >= length(R)) {
                  
                  R <- c(R, numeric(m / 2))
                  pval <- c(pval, numeric(m / 2))
                  kappai <- c(kappai, numeric(m / 2))
                  kappai_star <- c(kappai_star, numeric(m / 2))
                  Ci <- c(Ci, numeric(m / 2))
                  Si <- c(Si, numeric(m / 2))
                  Ci_plus <- c(Ci_plus, numeric(m / 2))
                  gammai <- c(gammai, numeric(m / 2))
                  alphai <- c(alphai, numeric(m / 2))
                  
                }
                
                if(indepTest=="gaussCItest") #if indepTest=gaussCItest
                {
                  pval[m] <- gaussCItest(x, z, y, suffStat)
                }
                if(indepTest=="disCItest") #if indepTest=disCItest
                {
                  pval[m] <- disCItest(x, z, y, suffStat) #additional
                }
                
                if (FDRcontrol == 'LOND') {
                  
                  # Calculate alpha using the LOND method.
                  alphai[m] <- SeqFDR(m,FDR,a=2,R)
                  
                  Alpha <- alphai[m]
                  
                  # Calculate the alpha value using the ADDIS algorithm
                } else if (FDRcontrol == 'ADDIS') {
                  
                  # Run ADDIS on the current iteration and update all the vectors
                  # and other values.
                  run_addis <- addis(alpha = FDR,
                                     tau = tau,
                                     lambda = lambda,
                                     iter = m,
                                     w0 = w0,
                                     pval = pval,
                                     alphai = alphai,
                                     gammai = gammai,
                                     kappai = kappai,
                                     kappai_star = kappai_star,
                                     K = K,
                                     Ci = Ci,
                                     Si = Si,
                                     Ri = R,
                                     Ci_plus = Ci_plus,
                                     Ci_sum = Ci_sum,
                                     Si_sum = Si_sum,
                                     normalizer = normalizer,
                                     exponent = exponent)
                  
                  # Update all values and vectors output from the addis function.
                  alphai[[m]] <- run_addis[[1]]
                  gammai[[m]] <- run_addis[[2]]
                  K <- run_addis[[3]]
                  R[[m]] <- run_addis[[5]]
                  Si[[m - 1]] <- run_addis[[6]]
                  Ci[[m - 1]] <- run_addis[[7]]
                  Ci_sum <- run_addis[[9]]
                  Si_sum <- run_addis[[10]]
                  
                  # Only update the kappai and Ci_plus vectors if K is greater
                  # than one. If K is zero then the first element in kappai will
                  # remain zero until the first rejection.
                  if (K != 0) {
                    
                    kappai[[K]] <- run_addis[[4]]
                    kappai_star[[K]] <- run_addis[[11]]
                    Ci_plus[1:K] <- run_addis[[8]]
                    
                  }
                  
                  # Update the Alpha value with the current alphai value.
                  Alpha <- alphai[[m]]
                  
                } else {
                  
                  Alpha <- alpha
                  
                }
                if (verbose) {
                cat("Additional pval value =", pval[m], "\n")
                cat("Alpha value =", Alpha, "\n")
                }
                #Alpha=SeqFDR(m,FDR,a=2,R) #Alpha valued from sequential FDR test
                if (pval[m] <= Alpha) {  #Reject H0 (H0:nodes are independent)
                  R[m] <- 1
                  
                  if (verbose) {
                    V <- colnames(g)
                    cat("\n", V[x], "->", V[y], "<-", V[z], "\n")  #Printout the v-structures
                    cat("Since pval<Alpha,additional test is rejected", "Nodes", V[x] ,"and" ,V[z] ,"are dependent given", V[y], "\n")
                  }
                  tarmat[z, y] <- 1 #directed z-->y
                } 
                else {
                  tarmat[y, z] <- 1
                  R[m] <- 0  #Accept H0
                }
              }
              #Case-4:If, y & z have relation and x & y conditionally dependent given z and x & z conditionally independent given y,
              #then edge direction will be z<-->y
              if (g1[y, z] == 1 & g1[x, z] == 1 & g1[x, y] == 1 & !(z %in% gInput@sepset[[x]][[y]]) &!(y %in% gInput@sepset[[x]][[z]]))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 1
              }
            }
            
          }
          
        }
      }
      #if(length(WW3)!=0)
      #{
      #Orient remaining edges,
      T1 <- which(tarmat==1,arr.ind = T)
      G1 <- which(g==1,arr.ind = T)
      D1 <- setdiff(G1,T1)
      #
      if(length(D1)!=0)
      {
        for (d1 in 1:length(D1))
        {
          Rem1 <- which(g[,D1[d1]]==1,arr.ind = T)
          Rem2 <- which(tarmat[,Rem1]==1,arr.ind = T)
          if(length(Rem1)!=0 & length(Rem2)!=0)
          {
            for (d11 in 1:length(Rem1))
            {
            x <- Rem2[d11]
            y <- Rem1[d11]
            z <- D1[d1]
            }
            
            if ((!is.na(x) & !is.na(y) & !is.na(z)) & (x!=0 & y!=0 & z!=0) & (x!=y & y!=z & z!=x))
            {
              #Case-2: If,y and z are adjacent, x and z conditionally independent given y,
              #then the edge direction will be y-->z
              if (g1[y, z] == 1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 & ((y %in% gInput@sepset[[x]][[z]]) || (y %in% gInput@sepset[[z]][[x]])))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 0
              }
              
              #Case-3: If,y and z are adjacent, x and z conditionally dependent given y,
              #then the edge direction will be z-->y.
              if (g1[y, z] == 1 & g1[x, z] != 1 & tarmat[x,y]==1 & tarmat[y, x]!=1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 &!(y %in% gInput@sepset[[x]][[z]]) & !(y %in% gInput@sepset[[z]][[x]]))
                
              {
                m <- m+1
                
                # Increase the length of the R, pval, and other addis vectors if
                # their length is greater than or equal to the current iteration.
                if (m >= length(R)) {
                  
                  R <- c(R, numeric(m / 2))
                  pval <- c(pval, numeric(m / 2))
                  kappai <- c(kappai, numeric(m / 2))
                  kappai_star <- c(kappai_star, numeric(m / 2))
                  Ci <- c(Ci, numeric(m / 2))
                  Si <- c(Si, numeric(m / 2))
                  Ci_plus <- c(Ci_plus, numeric(m / 2))
                  gammai <- c(gammai, numeric(m / 2))
                  alphai <- c(alphai, numeric(m / 2))
                  
                }
                
                if(indepTest=="gaussCItest") #if indepTest=gaussCItest
                {
                  pval[m] <- gaussCItest(x, z, y, suffStat)
                }
                if(indepTest=="disCItest") #if indepTest=gaussCItest
                {
                  pval[m] <- disCItest(x, z, y, suffStat) #additional
                }                
                if (FDRcontrol == 'LOND') {
                  
                  # Calculate alpha using the LOND method.
                  alphai[m] <- SeqFDR(m,FDR,a=2,R)
                  
                  Alpha <- alphai[m]
                  
                  # Calculate the alpha value using the ADDIS algorithm
                } else if (FDRcontrol == 'ADDIS') {
                  
                  # Run ADDIS on the current iteration and update all the vectors
                  # and other values.
                  run_addis <- addis(alpha = FDR,
                                     tau = tau,
                                     lambda = lambda,
                                     iter = m,
                                     w0 = w0,
                                     pval = pval,
                                     alphai = alphai,
                                     gammai = gammai,
                                     kappai = kappai,
                                     kappai_star = kappai_star,
                                     K = K,
                                     Ci = Ci,
                                     Si = Si,
                                     Ri = R,
                                     Ci_plus = Ci_plus,
                                     Ci_sum = Ci_sum,
                                     Si_sum = Si_sum,
                                     normalizer = normalizer,
                                     exponent = exponent)
                  
                  # Update all values and vectors output from the addis function.
                  alphai[[m]] <- run_addis[[1]]
                  gammai[[m]] <- run_addis[[2]]
                  K <- run_addis[[3]]
                  R[[m]] <- run_addis[[5]]
                  Si[[m - 1]] <- run_addis[[6]]
                  Ci[[m - 1]] <- run_addis[[7]]
                  Ci_sum <- run_addis[[9]]
                  Si_sum <- run_addis[[10]]
                  
                  # Only update the kappai and Ci_plus vectors if K is greater
                  # than one. If K is zero then the first element in kappai will
                  # remain zero until the first rejection.
                  if (K != 0) {
                    
                    kappai[[K]] <- run_addis[[4]]
                    kappai_star[[K]] <- run_addis[[11]]
                    Ci_plus[1:K] <- run_addis[[8]]
                    
                  }
                  
                  # Update the Alpha value with the current alphai value.
                  Alpha <- alphai[[m]]
                  
                } else {
                  
                  Alpha <- alpha
                  
                }
                if (pval[m] <= Alpha) {  #Reject H0 (H0:nodes are independent)
                  R[m] <- 1
                  
                  if (verbose) {
                    V <- colnames(g)
                    cat("\n", V[x], "->", V[y], "<-", V[z], "\n")  #Printout the v-structures
                    cat("Since pval<Alpha,additional test is rejected", "Nodes", V[x] ,"and" ,V[z] ,"are dependent given", V[y], "\n")
                  }
                  tarmat[z, y] <- 1 #directed z-->y
                } 
                else {
                  tarmat[y, z] <- 1
                  R[m] <- 0  #Accept H0
                }
              }
              #Case-4:If, y & z have relation and x & y conditionally dependent given z and x & z conditionally independent given y,
              #then edge direction will be z<-->y
              if (g1[y, z] == 1 & g1[x, z] == 1 & g1[x, y] == 1 & !(z %in% gInput@sepset[[x]][[y]]) &!(y %in% gInput@sepset[[x]][[z]]))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 1
              }
            }
          }
        }
      }
      #Check the reaming edges
      ind <- which(g1 == 1, arr.ind = TRUE)  #Pullout the all relation in adjacency matrix from undirected graph
      for (i in seq_len(nrow(ind))) {
        x <- ind[i, 1]
        y <- ind[i, 2]
        if(tarmat[x, y]==0 & tarmat[y, x]==0) #bidirected if still no edge in tarmat
        {
          tarmat[x, y] <- 1
          tarmat[y, x] <- 1
        }
      }
    }
    
  }#end when data not necessary to contain genetic variants.
  
  #Step 4
  #Produce as a same skeleton if no edges involves with genetic variants and no v-structures
  if(all(tarmat==0))
  {
    tarmat <- g1
  }
  #Step 5
  #If found the direction with only genetic variants, no v-structures and 
  #all the other nodes still undirected, then remaining edges will be bidirected 
  if(GV>0 & (any(tarmat[1:GV,]==1) || any(tarmat[,1:GV]==1)) & all(tarmat[-c(1:GV),-c(1:GV)]==0))
  {
    tarmat1 <- g1
    tarmat1[1:GV,] <- tarmat[1:GV,]
    tarmat1[,1:GV] <- tarmat[,1:GV]
    tarmat <- tarmat1
  }
  
  gInput@graph <- as(tarmat, "graphNEL")
  gInput@R <- R[1:m]
  gInput@K <- K
  gInput@pval <- pval[1:m]
  gInput@kappai <- kappai[1:K]
  gInput@kappai_star <- kappai_star[1:K]
  gInput@Ci <- Ci[1:m]
  gInput@Si <- Si[1:m]
  gInput@Ci_plus <- Ci_plus[1:K]
  gInput@gammai <- gammai[1:m]
  
  gInput
  
}

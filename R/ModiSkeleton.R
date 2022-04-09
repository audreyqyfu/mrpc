#This is the step 1 of MRPC to draw the undirected graph

ModiSkeleton <- function (data, suffStat, FDR, alpha, indepTest = c("gaussCItest", "disCItest","citest"), labels, p,
                          method = c("stable", "original", "stable.fast"), m.max = Inf, fixedGaps = NULL,
                          fixedEdges = NULL, NAdelete = TRUE, FDRcontrol = c("LOND", "ADDIS", "NONE"),
                          tau = 0.5, lambda = 0.25, verbose = FALSE) {
  cl <- match.call()
  if (!missing(p))
    stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if (missing(labels)) {
    if (missing(p))
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  }
  else {
    stopifnot(is.character(labels))
    if (missing(p))
      p <- length(labels)
    else if (p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
  }
  seq_p <- seq_len(p)
  method <- match.arg(method)
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedGaps), c(p, p)))
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps)))
    stop("fixedGaps must be symmetric")
  else G <- !fixedGaps
  diag(G) <- FALSE
  if (any(is.null(fixedEdges))) {
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (!identical(fixedEdges, t(fixedEdges)))
    stop("fixedEdges must be symmetric")
  {
    sepset <- lapply(seq_p, function(.) vector("list", p))
    pMax <- matrix(-Inf, nrow = p, ncol = p)
    m <- 0L    #Current test number
    Alpha <- 0L
    K <- 0
    
    # Initialize vectors/scalars that are specific to the addis function.
    R <- pval <- kappai <- Ci <- Si <- Ci_plus <- numeric(dim(data)[2]^2)
    gammai <- kappai_star <- alphai <- numeric(dim(data)[2]^2)
    
    # Create objects for the numerator and exponent of the gamma series. This
    # is a p-series whose infinite sum is 1.
    normalizer <- 0.4374901658
    exponent <- 1.6
    
    # Initialize the gammai_sum to be zero. This will be updated after the first
    # rejection occurs.
    gammai_sum <- 0
    
    diag(pMax) <- 1
    done <- FALSE
    ord <- 0L
    
    n.edgetests <- numeric(1)
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord1 <- ord + 1L] <- 0
      done <- TRUE
      #ind <- which(G, arr.ind = TRUE)
      ind <- which(upper.tri(G), arr.ind = TRUE)
      ind <- ind[order(ind[, 1]), ]
      remEdges <- nrow(ind)
      # Order refers to the number of nodes being conditioned on.
      if (verbose)
        cat("Order=", ord, "; remaining edges:", remEdges,
            "\n", sep = "")
      if (method == "stable") {
        G.l <- split(G, gl(p, p))
      }
      
      # Loop through all possible edges.
      for (i in 1:remEdges) {
        
        # Print every 100th index of the edges considered and the number of all
        # possible edges. If a test is performed, the details of this test will
        # also be printed later (starting at line 161).
        if (verbose && (verbose >= 2 || i%%100 == 0))
          cat("|i=", i, "|iMax=", remEdges, "\n")
        x <- ind[i, 1]
        y <- ind[i, 2]
        
        if (G[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- if (method == "stable")
            
            nbrsBool1 <- G.l[[x]]
          nbrsBool2 <- G.l[[y]] 
          # G.l[[x]]
          #else G[, x]
          nbrsBool1[y] <- FALSE
          nbrs_x <- seq_p[nbrsBool1]
          #G.l[[y]]
          #else G[, y]
          nbrsBool2[x] <- FALSE
          nbrs_y <- seq_p[nbrsBool2]
          
          nbrs <- unique(union(nbrs_x,nbrs_y))
          
          #G.l[[x]]
          #else G[, x]
          #nbrsBool[y] <- FALSE
          #nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord)
              done <- FALSE
            S <- seq_len(ord)
            repeat {
              n.edgetests[ord1] <- n.edgetests[ord1] + 1
              
              m <- m+1  #Total number of the test
              
              # Increase the length of the R, pval, and other ADDIS vectors if
              # their length is greater than or equal to the current iteration.
              if (m >= length(R)) {
                
                R <- c(R, numeric(dim(data)[2]^2))
                pval <- c(pval, numeric(dim(data)[2]^2))
                kappai <- c(kappai, numeric(dim(data)[2]^2))
                kappai_star <- c(kappai_star, numeric(dim(data)[2]^2))
                Ci <- c(Ci, numeric(dim(data)[2]^2))
                Si <- c(Si, numeric(dim(data)[2]^2))
                Ci_plus <- c(Ci_plus, numeric(dim(data)[2]^2))
                gammai <- c(gammai, numeric(dim(data)[2]^2))
                alphai <- c(alphai, numeric(dim(data)[2]^2))
                
              }
              
              ##Start to calculate P-value using ci.test and gaussCItest
              
              if(indepTest == "citest") #if indepTest=ci.test
              {   
                x <- data[,ind[i,1]]
                y <- data[,ind[i,2]]
                z <- data[,nbrs[S]]
                if (length(S)==0) {
                  P <- ci.test(x, y)
                  pval[m] <- P$p.value
                } else {
                  P <- ci.test(x, y, z)
                  pval[m] <- P$p.value  #P-Value
                }
                x <- ind[i, 1]
                y <- ind[i, 2]
              }
              if(indepTest=="gaussCItest") #if indepTest=gaussCItest
              {
                pval[m] <- gaussCItest(x, y, nbrs[S], suffStat)
              }
              if(indepTest=="disCItest") #if indepTest=disCItest
              {
                pval[m] <- disCItest(x, y, nbrs[S], suffStat)
              }
              
              #End to calculate P-value using ci.test and gaussCItest
              
              if (verbose)
                cat("x=", x, " y=", y, " S=", nbrs[S],"\n")
              if (is.na(pval[m]))
                pval[m] <- as.numeric(NAdelete)
              if (pMax[x, y] < pval[m])
                pMax[x, y] <- pval[m]
              if (verbose)
                cat("Test number =", m, "\n")
              if (verbose)
                cat("pval =", pval[m], "\n")
              
              if (FDRcontrol == 'LOND') { #if want to control sequential FDR 
                
                # Calculate alpha using the LOND method.
                alphai[m] <- SeqFDR(m,FDR,a=2,R)
                
                Alpha <- alphai[m]
                
              } else if (FDRcontrol == 'ADDIS') {
                
                # Calculate alpha using the ADDIS algorithm.
                
                # Initialize all the vectors and values for the first iteration.
                if (m == 1) {
                  
                  # Calculate w0 from tau, lambda, and alpha. This value is used
                  # in calculating alpha_t at each iteration.
                  w0 <- tau * lambda * FDR/2
                  
                  # Calculate the sum of the candidate p-values.
                  Ci_sum <- 0
                  
                  # Calculate the sum of the selected tests for each p-value tested.
                  Si_sum <- 0
                  
                  # The total number of rejections so far.
                  K <- 0
                  
                  # Calculate the first element in the gamma sequence.
                  gammai[1] <- normalizer / 1^exponent
                  
                  # Calculate alphai for the first test.
                  alphai[1] <- w0 * gammai[1]
                  
                  # Update the Alpha value with alphai[1]
                  Alpha <- alphai[1]
                  
                  # Determine if the first test should be rejected
                  R[1] <- pval[1] <= alphai[1]
                  
                } else {
                  
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
                                     gammai_sum = gammai_sum,
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
                  gammai_sum <- run_addis[[12]]
                  
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
                  
                }
                
                
              } else if (FDRcontrol == 'NONE') {
                Alpha <- alpha #if want to use fixed significance level 
              }
              if (verbose)
                cat("Alpha value =", Alpha, "\n")
              
              if (pval[m] <= Alpha) {  #Reject H0 (H0:nodes are independent)
                R[m] <- 1
                if (verbose)
                  cat("Since pval<Alpha,test is rejected: Nodes are dependent", "\n")
              } else {
                R[m] <- 0  #Accept H0
                if (verbose)
                  cat("Since pval>Alpha,test is accepted:Nodes are independent", "\n")
              }
              if (pval[m] >= Alpha) {
                G[x, y] <- G[y, x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                
                break
              } else {
                nextSet <- getNextSet(length_nbrs, ord, S)
                if (nextSet$wasLast)
                  break
                S <- nextSet$nextSet
              }
              
            }
            
          }
        }
        
      }
      
      ord <- ord + 1L
    }
    for (i in 1:(p - 1)) {
      for (j in 2:p) {
        pMax[i, j] <- pMax[j, i] <- max(pMax[i,j], pMax[j, i])
      }
    }
    }
  Gobject <- if (sum(G) == 0) {
    new("graphNEL", nodes = labels)
  } else {
    colnames(G) <- rownames(G) <- labels
    as(G, "graphNEL")
  }
  
  #temp<-new("pcAlgo",graph = Gobject,call = cl, n = integer(0),
  # max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
  # sepset = sepset,pMax = pMax, zMin = matrix(NA, 1, 1))
  
  new("MRPCclass",
      graph = Gobject,
      call = cl,
      n = integer(0),
      max.ord = as.integer(ord - 1),
      n.edgetests = n.edgetests,
      sepset = sepset,
      pMax = pMax,
      zMin = matrix(NA, 1, 1),
      test = m,
      alpha = Alpha,
      R = R,
      K = K,
      pval = pval,
      normalizer = normalizer,
      exponent = exponent,
      alphai = alphai,
      kappai = kappai,
      kappai_star = kappai_star,
      Ci = Ci,
      Si = Si,
      Ci_plus = Ci_plus,
      gammai = gammai,
      gammai_sum = gammai_sum)
  
  #return(list(obj=temp,test=m,alpha=Alpha,R=R))
  #else
  #{
  #return(list(graph = Gobject,call = cl, n = integer(0),
  # max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
  # sepset = sepset,pMax = pMax, zMin = matrix(NA, 1, 1),test=m,alpha=Alpha,R=R))
  
  #}
}

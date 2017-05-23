#This is the part 1 of MRPC to draw the undirected graph
ModiSkeleton<-function (data,suffStat,FDR, indepTest = c("gaussCItest", "citest"), labels, p, method = c("stable",
                                                             "original", "stable.fast"), m.max = Inf, fixedGaps = NULL,
                  fixedEdges = NULL, NAdelete = TRUE, verbose = FALSE)
{
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
    R<-0L      #Rejection number
    Alpha <-0L
    pval<-0L

    diag(pMax) <- 1
    done <- FALSE
    ord <- 0L

    n.edgetests <- numeric(1)
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord1 <- ord + 1L] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      ind <- ind[order(ind[, 1]), ]
      remEdges <- nrow(ind)
      if (verbose)
        cat("Order=", ord, "; remaining edges:", remEdges,
            "\n", sep = "")
      if (method == "stable") {
        G.l <- split(G, gl(p, p))
      }

      for (i in 1:remEdges) {
        if (verbose && (verbose >= 2 || i%%100 == 0))
          cat("|i=", i, "|iMax=", remEdges, "\n")
        x <- ind[i, 1]
        y <- ind[i, 2]

        if (G[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- if (method == "stable")
            G.l[[x]]
          else G[, x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord)
              done <- FALSE
            S <- seq_len(ord)
            repeat {
              n.edgetests[ord1] <- n.edgetests[ord1] +
                1
              m=m+1  #Total number of the test

              ##Start to calculate P-value using ci.test and gaussCItest
              
              if(indepTest=="citest") #if indepTest=ci.test
                              {   
                                x=data[,ind[i,1]]
                                y=data[,ind[i,2]]
                                z=data[,nbrs[S]]
                                if (length(S)==0) {
                                    P<- ci.test(x, y)
                                    pval[m]=P$p.value
                                } else {
                                    P<- ci.test(x, y, z)
                                    pval[m]=P$p.value  #P-Value
                                }
                                x <- ind[i, 1]
                                y <- ind[i, 2]
              }
              if(indepTest=="gaussCItest") #if indepTest=gaussCItest
                {
                pval[m]<- gaussCItest(x, y, nbrs[S], suffStat)
              }
                                   
              #End to calculate P-value using ci.test and gaussCItest

              if (verbose)
                cat("x=", x, " y=", y, " S=", nbrs[S],"\n")
              if (is.na(pval[m]))
                pval <- as.numeric(NAdelete)
              if (pMax[x, y] < pval[m])
                pMax[x, y] <- pval[m]

            cat("Test number =", m, "\n")
            cat("pval =", pval[m], "\n")
            
            Alpha=SeqFDR(m,FDR,a=2,R) #Alpha valued from sequential FDR test
            cat("Alpha value =", Alpha, "\n")

            if (pval[m]<= Alpha) {  #Reject H0 (H0:nodes are independent)
              R[m]=1
              cat("Since pval<Alpha,test is rejected: Nodes are dependent", "\n")
            } else {
              R[m]=0  #Accept H0
              cat("Since pval>Alpha,test is accepted:Nodes are independent", "\n")
            }
            if (pval[m]>= Alpha) {
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
  
  temp<-new("pcAlgo",graph = Gobject,call = cl, n = integer(0),
            max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
            sepset = sepset,pMax = pMax, zMin = matrix(NA, 1, 1))

  return(list(obj=temp,test=m,alpha=Alpha,R=R))
}

#This is the main function of MRPC algorithm combine with ModiSkeleton and EdgeOrientation

MRPC <- function (data, suffStat, GV, FDR = 0.05, alpha = 0.05, indepTest = c("gaussCItest","disCItest", "citest"), labels, p, fixedGaps = NULL,
               fixedEdges = NULL, NAdelete = TRUE, m.max = Inf, u2pd = c("relaxed", "rand", "retry"),
               skel.method = c("stable", "original", "stable.fast"), conservative = FALSE, maj.rule = FALSE,
               solve.confl = FALSE, FDRcontrol = c("LOND", "ADDIS", "NONE"), tau = 0.5, lambda = 0.25, verbose = FALSE)
{
  cl <- match.call()
  
  if (FDRcontrol == "LOND") {
      cat ("Using the LOND method for online FDR control at FDR = ", FDR, "\n")
  } else if (FDRcontrol == "ADDIS") {
      cat ("Using the ADDIS method for online FDR control at FDR = ", FDR, "\n")
  } else if (FDRcontrol == "NONE") {
      cat ("Not applying error control. The type I error rate for each test is ", alpha, "\n")
  }
  
  if (!missing(p))
    stopifnot(is.numeric(p), length(p <- as.integer(p)) ==
                1, p >= 2)
  if (missing(labels)) {
    if (missing(p))
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  }
  else {
    stopifnot(is.character(labels))
    if (missing(p)) {
      p <- length(labels)
    }
    else if (p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else message("No need to specify 'p', when 'labels' is given")
  }
  u2pd <- match.arg(u2pd)
  skel.method <- match.arg(skel.method)
  if (u2pd != "relaxed") {
    if (conservative || maj.rule)
      stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")
    if (solve.confl)
      stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n can only be run with 'u2pd = relaxed'")
  }
  if (conservative && maj.rule)
    stop("Choose either conservative PC or majority rule PC!")
    
  if (verbose)
  cat ("Test for independence:", indepTest, "\n")
  skel <- ModiSkeleton(data, suffStat, FDR = FDR, alpha = alpha, indepTest = indepTest,
                       labels = labels, method = skel.method, fixedGaps = fixedGaps,
                       fixedEdges = fixedEdges, NAdelete = NAdelete, m.max = m.max, 
                       FDRcontrol = FDRcontrol, tau = tau, lambda = lambda,
                       verbose = verbose)
  skel@call <- cl
  
  #indepTest <- match.fun (indepTest)
  #if (GV) {
    #cat ("test for independence:", indepTest, "\n")
    #skel <- ModiSkeleton(data,suffStat,GV,FDR,indepTest,labels = labels,
                         #method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                         #NAdelete = NAdelete, m.max = m.max, verbose = verbose)
    #skel$call <- cl
    if (!conservative && !maj.rule) {
      switch(u2pd,relaxed = EdgeOrientation(skel, GV = GV, suffStat, FDR, alpha,
                                            FDRcontrol = FDRcontrol, indepTest,
                                            tau = tau, lambda = lambda,
                                            verbose = verbose))
    }
    else {
      pc. <- pc.cons.intern(skel, suffStat, match.fun(indepTest), alpha = alpha,
                            version.unf = c(2, 1), maj.rule = maj.rule,
                            verbose = verbose)
      
      EdgeOrientation(pc.$sk, FDRcontrol = FDRcontrol, tau = tau,
                      lambda = lambda, verbose = verbose)
    }
  #} #else {
    #skel <- ModiSkeleton(data,suffStat,GV,FDR,indepTest,labels = labels,
                        # method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                         #NAdelete = NAdelete, m.max = m.max, verbose = verbose)
    #skel$obj@call <- cl
    #if (!conservative && !maj.rule) {
     # switch(u2pd, rand = udag2pdag(skel), retry = udag2pdagSpecial(skel)$pcObj,
            # relaxed = udag2pdagRelaxed(skel$obj, verbose = verbose,
             #                          solve.confl = solve.confl))
    #}
    #else {
     # pc. <- pc.cons.intern(skel, suffStat, match.fun(indepTest), alpha=FDR,
                           # version.unf = c(2, 1), maj.rule = maj.rule, verbose = verbose)
      #udag2pdagRelaxed(pc.$sk, verbose = verbose, unfVect = pc.$unfTripl,
                      # solve.confl = solve.confl)
    #}
  #}
}

#This is the main function of MRPC algorithm combine with ModiSkeleton and EdgeOrientation

MRPC <- function (data, suffStat, GV, FDR = 0.05, alpha = 0.05, indepTest = c("gaussCItest","disCItest", "citest"), labels, p, fixedGaps = NULL,
               fixedEdges = NULL, NAdelete = TRUE, m.max = Inf, u2pd = c("relaxed", "rand", "retry"),
               skel.method = c("stable", "original", "stable.fast"), conservative = FALSE, maj.rule = FALSE,
               solve.confl = FALSE, FDRcontrol = c("LOND", "ADDIS", "NONE"), tau = 0.5, lambda = 0.25, verbose = FALSE)
{
  cl <- match.call()
  
  # check FDRcontrol method
  if (FDRcontrol == "LOND") {
      cat ("Using the LOND method for online FDR control at FDR = ", FDR, "\n")
  } else if (FDRcontrol == "ADDIS") {
      cat ("Using the ADDIS method for online FDR control at FDR = ", FDR, "\n")
  } else if (FDRcontrol == "NONE") {
      cat ("Not applying error control. The type I error rate for each test is ", alpha, "\n")
  } else {
      stop ("Choose one of the following three options for FDRcontrol: LOND, ADDIS, or NONE")
  }
  
  # provide either the node labels
  # or they are generated based on the number of nodes p
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
  
  # Step I: Infer the graph skeleton
  skel <- ModiSkeleton(data, suffStat, FDR = FDR, alpha = alpha, indepTest = indepTest,
                       labels = labels, method = skel.method, fixedGaps = fixedGaps,
                       fixedEdges = fixedEdges, NAdelete = NAdelete, m.max = m.max, 
                       FDRcontrol = FDRcontrol, tau = tau, lambda = lambda,
                       verbose = verbose)
  skel@call <- cl
  
  # Step II: Orient the edges
    EdgeOrientation(skel, GV = GV, suffStat, FDR, alpha,
                                          FDRcontrol = FDRcontrol, indepTest,
                                          tau = tau, lambda = lambda,
                                          verbose = verbose)
                                          
}

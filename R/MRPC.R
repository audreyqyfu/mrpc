
MRPC<-function (data,suffStat,NQ,FDR,indepTest = c("gaussCItest", "citest"),labels, p, fixedGaps = NULL,
               fixedEdges = NULL, NAdelete = TRUE, m.max = Inf,u2pd = c("relaxed", "rand", "retry"),
               skel.method = c("stable", "original","stable.fast"), conservative = FALSE, maj.rule = FALSE,
               solve.confl = FALSE, verbose = FALSE)
{
  cl <- match.call()
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
  
  #indepTestName <- as.character (quote(indepTest))

  cat ("test for independence:", indepTest, "\n")
  skel <- ModiSkeleton(data,suffStat,FDR,indepTest,labels = labels,
                       method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                       NAdelete = NAdelete, m.max = m.max, verbose = verbose)
  skel$call <- cl
  
  #indepTest <- match.fun (indepTest)
  if (NQ) {
    if (!conservative && !maj.rule) {
      switch(u2pd,relaxed = EdgesOrientation(skel,NQ=NQ,suffStat,FDR,verbose = verbose))
    }
    else {
      pc. <- pc.cons.intern(skel, suffStat, match.fun(indepTest), alpha = FDR,
                            version.unf = c(2, 1), maj.rule = maj.rule, verbose = verbose)
      EdgesOrientation(pc.$sk, verbose = verbose)
    }
  } else {
    if (!conservative && !maj.rule) {
      switch(u2pd, rand = udag2pdag(skel), retry = udag2pdagSpecial(skel)$pcObj,
             relaxed = udag2pdagRelaxed(skel$obj, verbose = verbose,
                                        solve.confl = solve.confl))
    }
    else {
      pc. <- pc.cons.intern(skel, suffStat, match.fun(indepTest), alpha=FDR,
                            version.unf = c(2, 1), maj.rule = maj.rule, verbose = verbose)
      udag2pdagRelaxed(pc.$sk, verbose = verbose, unfVect = pc.$unfTripl,
                       solve.confl = solve.confl)
    }
  }
}

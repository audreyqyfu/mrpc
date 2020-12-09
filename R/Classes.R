setClass("MRPCclass",
         slots = c(call = "call",
                   n = "integer",
                   max.ord = "integer",
                   n.edgetests = "numeric",
                   sepset = "list",
                   pMax = "matrix",
                   graph = "graph", 
                   zMin = "matrix",
                   test ="numeric",
                   alpha ="numeric",
                   R = "numeric",
                   K = "numeric",
                   pval = "numeric",
                   normalizer = "numeric",
                   exponent = "numeric",
                   alphai = "numeric",
                   kappai = "numeric",
                   kappai_star = "numeric",
                   Ci = "numeric",
                   Si = "numeric",
                   Ci_plus = "numeric",
                   gammai = "numeric",
                   gammai_sum = "numeric"))

##' auxiliary, hidden
show.MRPC.amat <- function(amat, zero.print, ...) {
  cat("\nAdjacency Matrix G:\n")
  print.table(amat, zero.print=zero.print, ...)
}

setMethod("summary", "MRPCclass",
          function(object, amat = TRUE, zero.print = ".", ...) {
            cat("Object of class 'MRPCclass', from Call:\n",
                paste(deparse(object@call), sep = "\n", collapse = "\n"),
                "\n\nNmb. edgetests during skeleton estimation:\n", sep = "")
            cat("===========================================\n")
            cat("Max. order of algorithm: ", object@max.ord,
                "\nNumber of edgetests from m = 0 up to m =", object@max.ord,
                ": ", object@n.edgetests)
            g <- object@graph
            nbrs <- vapply(g@edgeL, function(x) length(x$edges), 1L)
            cat("\n\nGraphical properties of skeleton:\n")
            cat("=================================\n")
            cat("Max. number of neighbours: ", max(nbrs),
                "at node(s)", which(nbrs==max(nbrs)),
                "\nAvg. number of neighbours: ", mean(nbrs),"\n")
            if(amat)
              show.MRPC.amat(as(g, "matrix"), zero.print=zero.print)
          })

print.MRPCclass <- function(x, amat = FALSE, zero.print = ".", ...) {
  cat("Object of class 'MRPCclass', from Call:\n",
      paste(deparse(x@call), sep = "\n", collapse = "\n"),
      "\n", sep="")
  A <- as(x@graph, "matrix")
  if(amat)
    show.MRPC.amat(A, zero.print=zero.print, ...)
  amat2 <- A + 2*t(A)
  ude <- sum(amat2 == 3)/2
  de <- sum(amat2 == 1)
  cat("Number of undirected edges: ", ude, "\n")
  cat("Number of directed edges:   ", de, "\n")
  cat("Total number of edges:      ", de + ude, "\n")
  invisible(x)
}
setMethod("show", "MRPCclass", function(object) print.MRPCclass(object))


setMethod("plot", signature(x = "MRPCclass"),
          function(x, y, main = NULL, zvalue.lwd = FALSE,
                   lwd.max = 7, labels = NULL, ...)
          {
            #check.Rgraphviz()
            
            if(is.null(main))
              main <- deparse(x@call)
            attrs <- nodeAttrs <- list()
            p <- numNodes(G <- x@graph)
            if (!is.null(labels)) {
              attrs$node <- list(shape = "ellipse", fixedsize = FALSE)
              names(labels) <- nodes(G)
              nodeAttrs$label <- labels
            }
            
            if (zvalue.lwd && numEdges(G) != 0) {
              lwd.mat <-
                if(is.matrix(Z <- x@zMin) && all(dim(Z) == p)) Z
              else ## from newer pc(): 'zMin' is deprecated there, but pMax corresponds:
                qnorm(x@pMax/2, lower.tail=FALSE)
              lwd.mat <- lwd.max * lwd.mat/max(lwd.mat)
              z <- Rgraphviz::agopen(G, name = "lwdGraph",
                                     nodeAttrs = nodeAttrs, attrs = attrs)
              for (i in seq_along(z@AgEdge)) {
                z@AgEdge[[i]]@lwd <- lwd.mat[as.integer(z@AgEdge[[i]]@head),
                                             as.integer(z@AgEdge[[i]]@tail)]
              }
              Rgraphviz::plot(z, main = main, ...)
              #plot(z, main = main, ...)
            } else {
              #plot(G, nodeAttrs = nodeAttrs, main = main,
                              #attrs = attrs, ...)
              Rgraphviz::plot(G, nodeAttrs = nodeAttrs, main = main,
                              attrs = attrs, ...)
              }
          })

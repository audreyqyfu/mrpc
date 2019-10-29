#Plot a graph with nodes in modules indicated by colors

PlotGraphWithModules <- function(Adj_directed,PlotDendrogramObj,GV=GV,node.size=8,arrow.size = 5,label.size = 3,alpha = 1,...) {
  
  
  AC <- PlotDendrogramObj$dynamicColors #All (genes and GV) colors
  CT <- table(AC)   #Color table
  Adj_symmetric_matrix <- PlotDendrogramObj$Adj_symmetric_matrix #symmetric matrix from PlotDendrogram
  G_Order <- vector("list")
  n <- vector("list")
  MO <- c()
  for (i in 1:length(CT)) {
    G_Order[[i]] <- (which(AC==rownames(CT)[i]))
    n[[i]] <- length(G_Order[[i]])
    MO <- c(MO, G_Order[[i]])
  }
  #New adjacency matrix by module order (MO)
  New_Mat_MO <- Adj_directed[MO,MO]
  
  MM <- 0:(ncol(Adj_directed)-1)
  NN1 <- as.matrix(cumsum(n))
  NN2 <- NN1#[-nrow(NN1),]
  Module <- CutModules(MM,c(0,NN2))
  
  col <- rownames(CT)
  names(col) <- levels(Module)
  #library(network) need for network function
  #Craete a graphical network used for the next step
  net <- network(New_Mat_MO, directed = TRUE)
  #Make two group (i) for phenotypes (e.g., gene expression)=circle and (ii) genotypes=Triangle
  #shape.palette = c("Genotype" = 17,"Phenotype" = 19)
  #char=colnames(Adj_symmetric_matrix)[GV+1:(ncol(Adj_symmetric_matrix)-1)] #all GV
  if(GV==0){
    char <- NULL
  }
  else{
    char <- colnames(Adj_symmetric_matrix)[1:GV] #all GV
  }
  #char=colnames(Adj_symmetric_matrix)[-c(1:GV)] #all GV
  net%v%"phono" <- ifelse((colnames(New_Mat_MO) %in% char),"Genotype","Phenotype")
  
  #library(GGally) #Need for ggnet2 function
  #Plot the final graph
  PlotGraphWithModulesObj <- ggnet2(net,color = Module,palette = col,node.size=node.size,arrow.size = arrow.size,label=TRUE,label.size = label.size,alpha = alpha,
                    shape.legend = "",edge.label.color = Module,shape = net%v%"phono",shape.palette = c("Genotype" = 17,"Phenotype" = 19),
                    color.legend = "Modules",legend.position = "bottom",arrow.gap = 0.0280)
  return(PlotGraphWithModulesObj)
}

DendroModuleGraph=function(Adj_directed,minModuleSize,GV) {

#Adj_directed is a matrix from directed graph
#Start to find Module based on library(WGCNA)
#convert to a symmetric matrix because Adj_directed is not symmetric
  Ad_Matrix=Adj_directed
  for (i in 1:nrow(Ad_Matrix))
    {
    for (j in 1:ncol(Ad_Matrix))
      {
      if(Ad_Matrix[i,j]==1)
        {
        Ad_Matrix[j,i]=1
      }
    }
  }
#Ad_Matrix is symmetric matrix after converting
  TOM= TOMsimilarity(Ad_Matrix); #Need Symmetric matrix
  dissTOM = 1-TOM
#Call the hierarchical clustering function from fastcluster
#This hclust that provides a much faster hierarchical clustering routine than the standard hclust function.
  geneTree = fastcluster::hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04);

  minModuleSize =minModuleSize; #Need to mention minimum module size
# Module identification using dynamic tree cut:
#The Dynamic Tree Cut may identify modules whose expression profiles are very similar.
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = TRUE,
                              minClusterSize = minModuleSize);

  Grouplist=table(dynamicMods)
  #print(Grouplist) #0=means unassigned genes

  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)

  Colorlist=table(dynamicColors)
  #print(Colorlist)

# Plot the dendrogram and colors underneath
#par(mfrow=c(1,2))
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = FALSE, guideHang = 0.05,
                      main = "Dendrogram and module colors of nodes")


  AC<- dynamicColors #All (genes and CNV) colors
  CT<- Colorlist   #Color table

  G_Order <- vector("list")
  n<-vector("list")
  MO=c()
  for (i in 1:length(CT)) {
    G_Order[[i]]=(which(AC==rownames(CT)[i]))
    n[[i]]=length(G_Order[[i]])
    MO=c(MO, G_Order[[i]])
  }
  #New adjacency matrix by module order (MO)
  New_Mat_MO=Adj_directed[MO,MO]

  MM=0:(ncol(Adj_directed)-1)
  NN1=as.matrix(cumsum(n))
  NN2=NN1#[-nrow(NN1),]
  Module<- Cut_Modules(MM,c(0,NN2))

  col= rownames(CT)
  names(col) = levels(Module)
  #library(network) need for network function
  #Craete a graphical network used for the next step
  net= network(New_Mat_MO, directed = TRUE)
  #Make two group (i) for phenotypes (e.g., gene expression)=circle and (ii) genotypes=Triangle
  #shape.palette = c("Genotype" = 17,"Phenotype" = 19)
  #char=colnames(Ad_Matrix)[GV+1:(ncol(Ad_Matrix)-1)] #all CNV
  char=colnames(Ad_Matrix)[1:GV] #all CNV
  #char=colnames(Ad_Matrix)[-c(1:GV)] #all CNV
  net%v%"phono"= ifelse((colnames(New_Mat_MO) %in% char),"Genotype","Phenotype")

  #library(GGally) #Need for ggnet2 function
  #Plot the final graph
  plotobj=ggnet2(net,color = Module,palette = col,node.size=5,arrow.size = 3,label=TRUE,label.size = 1,alpha = 1,
                 shape.legend = "Node type",edge.label.color = Module,shape = net%v%"phono",shape.palette = c("Genotype" = 17,"Phenotype" = 19),
                 color.legend = "Modules color",legend.position = "bottom",arrow.gap = 0.010)
  return(list(graph=plotobj,dynamicColors=dynamicColors,GroupMods=Grouplist,GroupModsColors=Colorlist,Adjmatrixdirected=Adj_directed))
}

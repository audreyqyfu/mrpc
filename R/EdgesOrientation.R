#This is the part 2 of MRPC to direction determination of the undirected edges.
#We cosider x (1st col) is a genetic variants (SNPs/indels/CNV/eQTL) and y and z (remaning col's) are the gene expression
#We consider two scenario (with 4 cases), because we're working triplets to orient the v-structure first and then orient remaing edges
#So,when more than three nodes, then the position of the genetic variant and genes are changed

EdgesOrientation<-function (gInput,NQ=NQ,verbose = FALSE)
{
  g <- as(gInput@graph, "matrix") # g ia an adjacency from undirected graph
  g1=g
  p <- nrow(g)

  # tarmat (adjacency matrix for directed graph) is updated in every step, 
  # and contains the output of final topology
  tarmat=matrix(0, nrow(g),ncol(g)) #same row and column from g
  rownames(tarmat)=rownames(g)     #same row names from g
  colnames(tarmat)=colnames(g)     #same column names from g

  # extract all edges
  g[lower.tri(g)] <- 0  # Use only upper triangular because g is symmetric matrix
  edges <- which (g==1, arr.ind=TRUE)
  
  #Step-1
  if (NQ>0) {
    # identify edges involving Vs
    edgesWithBothVs <- edges[which (edges[,1]<= NQ & edges[,2]<=NQ), ]
    edgesWithVFirst <- edges[which (edges[,1]<= NQ & edges[,2]>NQ), ]
    edgesWithVSecond <- edges[which (edges[,1]> NQ & edges[,2]<=NQ), ]

    # assign 1s to corresponding edges in tarmat
    # edges between two Vs are undirected (or bidirected)
    tarmat[edgesWithBothVs] <- 1
    if (length (edgesWithBothVs) > 2) {
      tarmat[edgesWithBothVs[,2:1]] <- 1
    } else {
      tarmat[edgesWithBothVs[2:1]] <- 1
    }
    # edges involving one V go from V to the other node
    tarmat[edgesWithVFirst] <- 1
    if (length(edgesWithVSecond) > 2) {
      tarmat[edgesWithVSecond[,2:1]] <- 1
    } else {
      tarmat[edgesWithVSecond[2:1]] <- 1
    }
  }
  
  #Step-2
  #Start to orient v-structures
  cat("\n V-structures are as follows :\n")
  edgesWithGs <- edges[which (edges[,1]>NQ | edges[,2]>NQ), ]
  ind <- edgesWithGs
  
  #ind <- which(g1 == 1, arr.ind = TRUE)  #Pullout the all relation in adjacency matrix from undirected graph
  for (i in seq_len(nrow(ind))) {
    x <- ind[i, 1]
    y <- ind[i, 2]
    allZ <- setdiff(which(g1[y, ] == 1), x)

    for (z in allZ) {
      # Triplet x-y-z is directed x-->y<--z if x and z conditionally dependent given y
      if ((g1[x, z] == 0 & g1[x, y] == 1) & !(tarmat[y, x] ==1) & !(tarmat[z, y] ==1) & !(tarmat[y, z] ==1) & !(y %in% gInput@sepset[[x]][[z]] ||
                              y %in% gInput@sepset[[z]][[x]])) {
        if (verbose) {
          V=colnames(g)
          cat("\n", V[x], "->", V[y], "<-", V[z], "\n")  #Printout the v-structures
        }
        tarmat[x, y] <- tarmat[z, y] <- 1
        #tarmat[y, x] <- tarmat[y, z] <- 0
      }
    }
  }
  #End to orient v-structures

  #Step-3 Orient remaining edges, weather MR is applicable or not based on the following process:
      #Pullout edges to extract parent and child/gene nodes (goal to make a triplet)
      #extract edges with involving parent and child/gene nodes
      #ignore the nodes that already have direction
      #form a triplet
      #determine the direction using test results from part 1 based on MR
      #Repeat step 3 until all undirected edges to directed
  #Start to orient remaining edges.
  repeat {

    S1=which(tarmat==1,arr.ind = T) #Pullout edges to extract parent and child/gene nodes from tarmat
    for (v in 1:nrow(S1)) {
      S_pn=S1[v,1] #extract the parent nodes
      S_gn=S1[v,2] #extract the gene/child nodes

      S11_cp=which(g[,S_pn]==1,arr.ind = T) #pullout edges with parent nodes in column
      S11_rp=which(g[S_pn,]==1,arr.ind = T) #pullout edges with parent nodes in row

      S11_cg=which(g[,S_gn]==1,arr.ind = T) #pullout edges with gene/child nodes in column
      S11_rg=which(g[S_gn,]==1,arr.ind = T) #pullout edges with gene/child nodes in row

      S11=c(S11_cg,S11_rg,S11_cp,S11_rp)  #cobine all
      S11<-S11[S11!= S_pn]                #ignore if parent nodes that already have direction
      S11<-S11[S11!= S_gn]                #ignore if gene/child nodes that already have direction

      if(length(S11)!=0) {
        for (d2 in 1:length(S11)) {
          #ignore if already have direction
          if(tarmat[S_gn,S11[d2]]==1||tarmat[S11[d2],S_gn]==1) {
            S111=integer(0)
          } else {
            S111=S11[d2]
          }

          if(length(S111)>0) {
            for (f1 in 1:length(S111)) {
              y=S111[f1]
              #for (f2 in 1:length(S_gn[])) {
              x=S_gn
              if(g[S_gn,S11[d2]]==0 & g[S11[d2],S_gn]==0) {
                x=S_pn
              }

              z<-0

            #Scenario-2 (gene is the first column of the triplets)
            #If Gene (other than genetic variant) is the first column
            if (x>NQ) {
              ny=x  #New y
              z=y
              y=ny
              d=tarmat[ ,c(y,z)]
              x1=which(d== 1, arr.ind = TRUE) # Which node have relation with present nodes y and z
              if (!empty(x1)) {
                x=x1[1,1]  #New x is the row of x1
              }

              if (all(d==0)) {
                d=tarmat[c(y,z), ]
                x1=which(d== 1, arr.ind = TRUE) # Which node have relation with present nodes y and z
                if (!empty(x1)) {
                  x=x1[1,2]   #New x is the column of x1
                }
              }
            }
            if ((x!=0 & y!=0 & z!=0) & (x!=y & y!=z & z!=x)) {
              #Case-2: If,y and z are adjacent, x and z conditionally independent given y,
              #then the edge direction will be y-->z
              if (g1[y, z] == 1 & (y %in% gInput@sepset[[x]][[z]]) || (y %in% gInput@sepset[[z]][[x]])) {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 0
              }
              
              #Case-3: If,y and z are adjacent, x and z conditionally dependent given y,
              #then the edge direction will be z-->y.
              if (x!=y & y!=z & z!=x & g1[y, z] == 1 & tarmat[y, z]!=1 & !(y %in% gInput@sepset[[x]][[z]]) & !(y %in% gInput@sepset[[z]][[x]])) {
                tarmat[z, y]  <- 1
                tarmat[y, z]  <- 0
              }
              
              #Case-4:If, y & z have relation and x & y conditionally dependent given z and x & z conditionally independent given y,
              #then edge direction will be z<-->y
              if (g1[y, z] == 1 & g1[x, z] == 1 & g1[x, y] == 1 & !(z %in% gInput@sepset[[x]][[y]])&!(y %in% gInput@sepset[[x]][[z]])) {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 1
              }
            }
          }
        }
        }
      }
    }
  #End to orient remaining edges.

#Repeat untill fulfill the given condition
    tar1 <- tarmat
    for (i in 1:nrow(tar1))
      {
      for (j in 1:ncol(tar1))
        {
        if(tar1[i,j]==tar1[j,i])
        {
          tar1[j,i]=0
        }
      }
    }

    if (length(which(tar1==1))==length(which(g==1))) #Condition for repeating
      break
  }
  repeat {
    old_tarmat <- tarmat
    if (all(tarmat == old_tarmat))
      break
  }
  
  gInput@graph<-as(tarmat, "graphNEL")
  gInput

}

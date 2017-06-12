#This is the part 2 of MRPC to direction determination of the undirected edges.
#We cosider x (1st col) is a genetic variants (SNPs/indels/CNV/eQTL) and y and z (remaning col's) are the gene expression
#We consider two scenario (with 4 cases), because we're working triplets to orient the v-structure first and then orient remaing edges
#So,when more than three nodes, then the position of the genetic variant and genes are changed

EdgesOrientation<-function (gInput,NQ=NQ,suffStat,FDR,verbose = FALSE)
{
  g <- as(gInput$graph, "matrix") # g ia an adjacency from undirected graph
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
  
  #Step-2
  #Start to orient v-structures
  cat("\n V-structures are as follows :\n")
  # extract edges involving at least one gene node
  #edgesWithGs <- edges[which (edges[,1]>NQ | edges[,2]>NQ), ]
  #ind <- edgesWithGs
  m=gInput$test  #current test
  Alpha=gInput$alpha  #alpha
  R=gInput$R          #decision of test
  ind <- which(g1 == 1, arr.ind = TRUE)  #Pullout the all relation in adjacency matrix from undirected graph
  for (i in seq_len(nrow(ind))) {
    x <- ind[i, 1]
    y <- ind[i, 2]
    
    allZ <- setdiff(which(g1[y, ] == 1), x)
    for (z in allZ) {
    # Triplet x-y-z is directed x-->y<--z if x and z conditionally dependent given y
      if ((g1[x, z] == 0 & g1[x, y] == 1 & g1[y, z] == 1)  & !(tarmat[y, x] ==1) & !(tarmat[z, y] ==1) & !(tarmat[y, z] ==1) & 
          !(y %in% gInput$sepset[[x]][[z]] || y %in% gInput$sepset[[z]][[x]])) 
        {
        m=m+1
        pval=gaussCItest(x, z, y, suffStat) #additional conditional test
        Alpha=SeqFDR(m,FDR,a=2,R) #Alpha valued from sequential FDR test
        
        cat("Additional pval value =", pval, "\n")
        cat("Alpha value =", Alpha, "\n")
        
        if (pval<= Alpha) {  #Reject H0 (H0:nodes are independent)
          R[m]=1
        if (verbose) {
            V=colnames(g)
            cat("\n", V[x], "->", V[y], "<-", V[z], "\n")  #Printout the v-structures
            cat("Since pval<Alpha,additional test is rejected", "Nodes", V[x] ,"and" ,V[z] ,"are dependent given", V[y], "\n")
             }
          tarmat[x, y] <- tarmat[z, y] <- 1 #directed x-->y<--z
          } 
        else {
          R[m]=0  #Accept H0
            }

    }
  }
  }
  #End to orient v-structures

  
  #Step-3
       #Orient remaining edges, weather MR is applicable or not based on the following process:
       #Pullout all edges from g1 matrix and define 1st and 2nd node (goal to make a triplet for the remaining edges)
       #extract all edges with involving those 1st and 2nd nodes (from row and column of g)
       #ignore the nodes that already have direction
       #form a triplet
       #determine the direction using test results from part 1 based on MR
       #Repeat until all undirected edges to directed
  m=m
  Alpha=Alpha
  R=R
  #Start to orient remaining edges.
  if (any(tarmat == 1)) #if at least one edge directed already made so far
  {
    edgeall=which(g1==1,arr.ind = T) #Pullout the all edges again from g1 matrix
    for (u in 1:nrow(edgeall))
    {
      FirstNode=edgeall[u,1]  #extract the 1st node
      SecondNode=edgeall[u,2] #extract the 2nd node
      
      edgesWithFirstNodeinCol=which(g[,FirstNode]==1,arr.ind = T) #pullout edges with 1st node in column from g
      edgesWithFirstNodeinRow=which(g[FirstNode,]==1,arr.ind = T) #pullout edges with 1stt node in row from g
      
      edgesWithSecondNodeinCol=which(g[,SecondNode]==1,arr.ind = T) #pullout edges with 2nd nodes in column from g
      edgesWithSecondNodeinRow=which(g[SecondNode,]==1,arr.ind = T) #pullout edges with gene/child nodes in row from g
      
      combineEdgesAll=c(edgesWithSecondNodeinCol,edgesWithSecondNodeinRow,edgesWithFirstNodeinCol,edgesWithFirstNodeinRow)  #cobine all
      combineEdgesAll<-combineEdgesAll[combineEdgesAll!= FirstNode]   #ignore if 1st nodes
      combineEdgesAll<-combineEdgesAll[combineEdgesAll!= SecondNode]  #ignore if 2nd nodes
      
      if(length(combineEdgesAll)!=0)
      {
        for (v in 1:length(combineEdgesAll))
        {
          if(tarmat[SecondNode,combineEdgesAll[v]]==1||tarmat[combineEdgesAll[v],SecondNode]==1)  #ignore if already have direction
          {
            edgesNew=integer(0)
          }
          else
          {
            edgesNew=combineEdgesAll[v]
          }
          
          if(length(edgesNew)>0)
          {
            for (w in 1:length(edgesNew))
            {
              y=edgesNew[w]
              x=SecondNode
              if(g[SecondNode,combineEdgesAll[v]]==0 & g[combineEdgesAll[v],SecondNode]==0)
              {
                x=FirstNode
              }
              
              z<-0
              
      #Scenario-2 (gene is the first column of the triplets)
              if (x>NQ) #If Gene (other than genetic variant) is the first column
              {
                ny=x  #New y
                z=y
                y=ny
                d=tarmat[ ,c(y,z)]  #pullout the relation in column with y and z in tarmat
                x1=which(d== 1, arr.ind = TRUE) # Which node have relation with present nodes y and z
                if (!empty(x1))
                {
                  x=x1[1,1]  #New x 
                }
                
                if (all(d==0))  #If no relation found in column with y and z in tarmat
                {
                  d=tarmat[c(y,z), ]  #pullout the relation in row with y and z in tarmat
                  x1=which(d== 1, arr.ind = TRUE) # Which node have relation with present nodes y and z
                  if (!empty(x1))
                  {
                    x=x1[1,2]   #New x 
                  }
                  
                }
                
              }
  if ((x!=0 & y!=0 & z!=0) & (x!=y & y!=z & z!=x))
              {
                #Case-2: If,y and z are adjacent, x and z conditionally independent given y,
                #then the edge direction will be y-->z
                if (g1[y, z] == 1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 & ((y %in% gInput$sepset[[x]][[z]]) || (y %in% gInput$sepset[[z]][[x]])))
                {
                  tarmat[y, z]  <- 1
                  tarmat[z, y]  <- 0
                }
                
                #Case-3: If,y and z are adjacent, x and z conditionally dependent given y,
                #then the edge direction will be z-->y.
                if (g1[y, z] == 1 & g1[x, z] != 1 & tarmat[x,y]==1 & tarmat[y, x]!=1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 &!(y %in% gInput$sepset[[x]][[z]]) & !(y %in% gInput$sepset[[z]][[x]]))
                
                {
                  m=m+1
                  pval=gaussCItest(x, z, y, suffStat) #additional conditional test
                  Alpha=SeqFDR(m,FDR,a=2,R) #Alpha valued from sequential FDR test
                  if (pval<= Alpha) {  #Reject H0 (H0:nodes are independent)
                    R[m]=1
                   
                if (verbose) {
                      V=colnames(g)
                      cat("\n", V[x], "->", V[y], "<-", V[z], "\n")  #Printout the v-structures
                      cat("Since pval<Alpha,additional test is rejected", "Nodes", V[x] ,"and" ,V[z] ,"are dependent given", V[y], "\n")
                      }
                    tarmat[z, y] <- 1 #directed z-->y
                  } 
                  else {
                    tarmat[y, z] <- 1
                    R[m]=0  #Accept H0
                  }
                }
                #Case-4:If, y & z have relation and x & y conditionally dependent given z and x & z conditionally independent given y,
                #then edge direction will be z<-->y
                if (g1[y, z] == 1 & g1[x, z] == 1 & g1[x, y] == 1 & !(z %in% gInput$sepset[[x]][[y]]) &!(y %in% gInput$sepset[[x]][[z]]))
                {
                  tarmat[y, z]  <- 1
                  tarmat[z, y]  <- 1
                }
              }
            }
          }
        }
      }
    }
  }
  #Step 4
  #Produce as a same skeleton if no edges involves with genetic variants and no v-structures
  if(all(tarmat==0))
  {
    tarmat=g1
  }
  #Step 5
  #If found the direction with only genetic variants, no v-structures and 
  #all the other nodes still undirected, then remaining edges will be bidirected 

  if((any(tarmat[1:NQ,]==1) || any(tarmat[,1:NQ]==1)) & all(tarmat[-c(1:NQ),-c(1:NQ)]==0))
  {
    tarmat1=g1
    tarmat1[1:NQ,]=tarmat[1:NQ,]
    tarmat1[,1:NQ]=tarmat[,1:NQ]
    tarmat=tarmat1
  }
  gInput$graph<-as(tarmat, "graphNEL")
  gInput

}

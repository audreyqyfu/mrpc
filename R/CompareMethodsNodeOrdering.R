# In this code we used the same data set but permute the T nodes 

# to generate multiple permuted data sets and apply the methods 

# to check the variation in inferred graphs. 


CompareMethodsNodeOrdering <- function (N, signal,model,n_data, n_nodeordering) {
  
  
  # Parameters
  
  p <- 0.45
  
  b0.1 <- 0
  
  b1.1 <- signal
  
  b1.2 <- signal
  
  b1.3 <- signal
  
  sd.1 <- 1
  
  switch(model,
         truth1 = {
           
           # Truth1 model (V1-->T1-->T2-->T3)
           
           tarmat_s1 <- matrix(0,nrow=4,ncol = 4)
           
           colnames(tarmat_s1) <- c("V1","T1","T2","T3")
           
           rownames(tarmat_s1) <- colnames(tarmat_s1)
           
           # Adjacency matrix from the true graph
           
           tarmat_s1[1,2] <- 1
           
           tarmat_s1[2,3] <- 1
           
           tarmat_s1[3,4] <- 1
           
           # The number of genetic variants in the graph.
           
           GV <- 1
           
           # The number of nodes
           
           n.nodes <- ncol (tarmat_s1)
           
           # Lists for the adjacency matrix output by each method.
           
           # MRPC
           
           Adjlist_MRPC <- vector (mode = 'list', length = n_data)
           
           # pc
           
           Adjlist_PC <- vector (mode = 'list', length = n_data)
           
           # pc.stable
           
           Adjlist_pc.stable <- vector (mode = 'list', length = n_data)
           
           # mmpc
           
           Adjlist_mmpc <- vector (mode = 'list', length = n_data)
           
           # mmhc
           
           Adjlist_mmhc <- vector (mode = 'list', length = n_data)
           
           # Matrix for the number of unique graphs for each data set (across n permutations).
           
           countMatrix <- matrix (nrow = n_data, ncol = 5)
           
           # name the columns according to the method.
           
           colnames (countMatrix) <- c('MRPC', 'pc', 'pc.stable', 'mmpc', 'mmhc')
           
           # Loop through the independent data sets.
           
           for (e in 1:n_data) {
             
             cat ("n_data=", e)
             
             # Simulate data for the true1.
             
             V1 <- c(sample(c(0, 1, 2),size = N,replace = TRUE,prob = c((1 - p)^2,2*p*(1 - p),p^2)))
             
             T1 <- SimulateData1P(N=N,P1=V1,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
             
             T2 <- SimulateData1P(N=N,P1=T1,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
             
             T3 <- SimulateData1P(N=N,P1=T2,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
             
             Data1 <- cbind(V1,T1,T2,T3)
             
             # MRPC
             
             Adjlist_MRPC[[e]] <- vector(mode = 'list', length = n_nodeordering)
             
             # pc
             
             Adjlist_PC[[e]] <- vector(mode = 'list', length = n_nodeordering)
             
             # pc.stable
             
             Adjlist_pc.stable[[e]] <- vector(mode = 'list', length = n_nodeordering)
             
             # mmpc
             
             Adjlist_mmpc[[e]] <- vector(mode = 'list', length = n_nodeordering)
             
             # mmhc
             
             Adjlist_mmhc[[e]] <- vector(mode = 'list', length = n_nodeordering)
             
             # Loop through all node permutations for one data set.
             
             for (j in 1:n_nodeordering) {
               
               cat ("n_nodeordering=", j)
               
               # Create a new ordering for the T nodes.
               
               temp.order <- c (GV, sample((GV + 1):n.nodes))
               
               # New data matrix with permuted nodes.
               
               Data2 <- Data1[, temp.order]
               
               n <- nrow (Data2)    # Sample size
               
               V <- colnames (Data2) # Node labels
               
               # Calculate Pearson correlation
               
               suffStat <- list (C = cor (Data2), n = n)
               
               # Infer the graph by MRPC
               
               MRPC_Inferred <- MRPC (Data2, suffStat, GV = GV, FDR = 0.05, alpha = 0.05, FDRcontrol = TRUE, indepTest = 'gaussCItest', labels = V, verbose = FALSE)
               
               # Adjacency matrix from the graph by MRPC
               
               G_MRPC <- as (MRPC_Inferred@graph, "matrix")
               
               # Save adjacency matrix
               
               Adjlist_MRPC[[e]][[j]] <- AdjustMatrix (tarmat_s1, G_MRPC)
               
               # Infer the graph by PC
               
               PC_Inferred <- pc (suffStat, alpha = 0.05, indepTest = gaussCItest, labels = V, verbose = F)
               
               # Adjacency matrix from the graph by pc
               
               G_PC <- as (PC_Inferred@graph, "matrix")
               
               # Save adjacency matrix
               
               Adjlist_PC[[e]][[j]] <- AdjustMatrix (tarmat_s1, G_PC)
               
               # arcs not to be included from gene expression to genotype
               
               to <- rep (colnames (Data2)[1:GV], each = (ncol (Data2) - GV))
               
               from <- rep (colnames (Data2)[(GV + 1):ncol (Data2)], GV)
               
               bl <- cbind (from, to)
               
               # Infer the graph by pc.stable
               
               pc.stable_Inferred <- pc.stable (data.frame (Data2), blacklist = bl, alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
               
               # Inferred graph object by pc.stable
               
               G_pc.stable <- amat (pc.stable_Inferred)
               
               # Save adjacency matrix
               
               Adjlist_pc.stable[[e]][[j]] <- AdjustMatrix (tarmat_s1, G_pc.stable)
               
               # Infer the graph by mmpc
               
               mmpc_Inferred <- mmpc (data.frame (Data2), blacklist = bl, alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
               
               # Inferred graph object by mmpc
               
               G_mmpc <- amat (mmpc_Inferred)
               
               # Save adjacency matrix
               
               Adjlist_mmpc[[e]][[j]] <- AdjustMatrix (tarmat_s1, G_mmpc)
               
               # Infer the graph by mmhc
               
               mmhc_Inferred <- mmhc (data.frame (Data2), blacklist = bl, debug = FALSE)
               
               # Inferred graph object by mmpc
               
               G_mmhc <- amat (mmhc_Inferred)
               
               # Save adjacency matrix
               
               Adjlist_mmhc[[e]][[j]] <- AdjustMatrix (tarmat_s1, G_mmhc)
               
             }
             
             # Calculate the number of unique graphs inferred by each method across all permutations.
             
             countMatrix[e, ] <- c(length (unique (Adjlist_MRPC[[e]])), length (unique (Adjlist_PC[[e]])), length (unique (Adjlist_pc.stable[[e]])), length (unique (Adjlist_mmpc[[e]])), length( unique (Adjlist_mmhc[[e]])))
             
           }
           
           return(countMatrix)
           
         },
         
         truth2 = {
           
           # Truth2 model (V1-->T1<--T2-->T3)
           
           tarmat_s2 <- matrix(0,nrow=4,ncol = 4)
           
           colnames(tarmat_s2) <- c("V1","T1","T2","T3")
           
           rownames(tarmat_s2) <- colnames(tarmat_s2)
           
           # Adjacency matrix from the true graph
           
           tarmat_s2[1,2] <- 1
           
           tarmat_s2[3,2] <- 1
           
           tarmat_s2[3,4] <- 1
           
           # The number of genetic variants in the graph.
           
           GV <- 1
           
           # The number of nodes
           
           n.nodes <- ncol (tarmat_s2)
           
           # Lists for the adjacency matrix output by each method.
           
           # MRPC
           
           Adjlist_MRPC <- vector (mode = 'list', length = n_data)
           
           # pc
           
           Adjlist_PC <- vector (mode = 'list', length = n_data)
           
           # pc.stable
           
           Adjlist_pc.stable <- vector (mode = 'list', length = n_data)
           
           # mmpc
           
           Adjlist_mmpc <- vector (mode = 'list', length = n_data)
           
           # mmhc
           
           Adjlist_mmhc <- vector (mode = 'list', length = n_data)
           
           # Matrix for the number of unique graphs for each data set (across n permutations).
           
           countMatrix <- matrix (nrow = n_data, ncol = 5)
           
           # name the columns according to the method.
           
           colnames (countMatrix) <- c('MRPC', 'pc', 'pc.stable', 'mmpc', 'mmhc')
           
           # Loop through the independent data sets.
           
           for (e in 1:n_data) {
             
             cat ("n_data=", e)
             
             # Simulate data for the true2.
             
             V1 <- c(sample(c(0, 1, 2),size = N,replace = TRUE,prob = c((1 - p)^2,2*p*(1 - p),p^2)))
             
             T2 <- SimulateDataNP(N=N,b0.1=b0.1,sd.1=sd.1)
             
             T1 <- SimulateData2P(N=N,P1=V1,P2=T2,b0.1=b0.1,b1.1=b1.1,b1.2=b1.2,sd.1=sd.1)
             
             T3 <- SimulateData1P(N=N,P1=T2,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
             
             # Combined
             Data1 <- cbind(V1,T1,T2,T3)
             
             # MRPC
             
             Adjlist_MRPC[[e]] <- vector(mode = 'list', length = n_nodeordering)
             
             # pc
             
             Adjlist_PC[[e]] <- vector(mode = 'list', length = n_nodeordering)
             
             # pc.stable
             
             Adjlist_pc.stable[[e]] <- vector(mode = 'list', length = n_nodeordering)
             
             # mmpc
             
             Adjlist_mmpc[[e]] <- vector(mode = 'list', length = n_nodeordering)
             
             # mmhc
             
             Adjlist_mmhc[[e]] <- vector(mode = 'list', length = n_nodeordering)
             
             # Loop through all node permutations for one data set.
             
             for (j in 1:n_nodeordering) {
               
               cat ("n_nodeordering=", j)
               
               # Create a new ordering for the T nodes.
               
               temp.order <- c (GV, sample((GV + 1):n.nodes))
               
               # New data matrix with permuted nodes.
               
               Data2 <- Data1[, temp.order]
               
               n <- nrow (Data2)    # Sample size
               
               V <- colnames (Data2) # Node labels
               
               # Calculate Pearson correlation
               
               suffStat <- list (C = cor (Data2), n = n)
               
               # Infer the graph by MRPC
               
               MRPC_Inferred <- MRPC (Data2, suffStat, GV = GV, FDR = 0.05, alpha = 0.05, FDRcontrol = TRUE, indepTest = 'gaussCItest', labels = V, verbose = FALSE)
               
               # Adjacency matrix from the graph by MRPC
               
               G_MRPC <- as (MRPC_Inferred@graph, "matrix")
               
               # Save adjacency matrix
               
               Adjlist_MRPC[[e]][[j]] <- AdjustMatrix (tarmat_s2, G_MRPC)
               
               # Infer the graph by PC
               
               PC_Inferred <- pc (suffStat, alpha = 0.05, indepTest = gaussCItest, labels = V, verbose = F)
               
               # Adjacency matrix from the graph by pc
               
               G_PC <- as (PC_Inferred@graph, "matrix")
               
               # Save adjacency matrix
               
               Adjlist_PC[[e]][[j]] <- AdjustMatrix (tarmat_s2, G_PC)
               
               # arcs not to be included from gene expression to genotype
               
               to <- rep (colnames (Data2)[1:GV], each = (ncol (Data2) - GV))
               
               from <- rep (colnames (Data2)[(GV + 1):ncol (Data2)], GV)
               
               bl <- cbind (from, to)
               
               # Infer the graph by pc.stable
               
               pc.stable_Inferred <- pc.stable (data.frame (Data2), blacklist = bl, alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
               
               # Inferred graph object by pc.stable
               
               G_pc.stable <- amat (pc.stable_Inferred)
               
               # Save adjacency matrix
               
               Adjlist_pc.stable[[e]][[j]] <- AdjustMatrix (tarmat_s2, G_pc.stable)
               
               # Infer the graph by mmpc
               
               mmpc_Inferred <- mmpc (data.frame (Data2), blacklist = bl, alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
               
               # Inferred graph object by mmpc
               
               G_mmpc <- amat (mmpc_Inferred)
               
               # Save adjacency matrix
               
               Adjlist_mmpc[[e]][[j]] <- AdjustMatrix (tarmat_s2, G_mmpc)
               
               # Infer the graph by mmhc
               
               mmhc_Inferred <- mmhc (data.frame (Data2), blacklist = bl, debug = FALSE)
               
               # Inferred graph object by mmpc
               
               G_mmhc <- amat (mmhc_Inferred)
               
               # Save adjacency matrix
               
               Adjlist_mmhc[[e]][[j]] <- AdjustMatrix (tarmat_s2, G_mmhc)
               
             }
             
             # Calculate the number of unique graphs inferred by each method across all permutations.
             
             countMatrix[e, ] <- c(length (unique (Adjlist_MRPC[[e]])), length (unique (Adjlist_PC[[e]])), length (unique (Adjlist_pc.stable[[e]])), length (unique (Adjlist_mmpc[[e]])), length( unique (Adjlist_mmhc[[e]])))
             
           }
           
           return(countMatrix)
           
         },
         
         stop("Model not included or missing"))
  
}

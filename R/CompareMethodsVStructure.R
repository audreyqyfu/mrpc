# Comparison of inference accuracy of different methods on data 
# with and without a v-structure

CompareMethodsVStructure <- function(N,signal,model,includeGV,ita) {

  #parameters setting
  p <- 0.45
  b0.1 <- 0
  b1.1 <- signal
  b1.2 <- signal
  b1.3 <- signal
  sd.1 <- 1
  
  # Initial Recall
  MRPC_Recall <- c()
  PC_Recall <- c()
  pc.stable_Recall <- c()
  mmpc_Recall <- c()
  mmhc_Recall <- c()
  
  # Initial Precision
  MRPC_Precision <- c()
  PC_Precision <- c()
  pc.stable_Precision <- c()
  mmpc_Precision <- c()
  mmhc_Precision <- c()
  
  switch(model,
         
         model1 = {
           # Truth for model 1 (V1-->T1-->T2) without v-structure 
           
           tarmat_1 <- matrix(0,nrow=3,ncol = 3)
           colnames(tarmat_1) <- c("V1","T1","T2")
           rownames(tarmat_1) <- c("V1","T1","T2")
           # Adjacency matrix
           tarmat_1[1,2] <- 1
           tarmat_1[2,3] <- 1
           
           # Adjacency matrix from the graph by MRPC
           Truth_1 <- as(tarmat_1, "graphNEL")
           
           # The number of genetic variants in the graph.
           GV <- 1
           # The number of nodes
           n.nodes <- ncol (tarmat_1)
           
           #Ietaration for model 1
           for (i in 1:ita) {
             #Data for model 1
             simu.data_1 <- SimulateData(N = N,p = p,'model1',b0.1 = b0.1,b1.1 = b1.1,b1.2 = b1.2,b1.3 = b1.3,sd.1 = sd.1)
             
             # Create a new ordering for the T nodes.
             temp.order <- c (GV, sample((GV + 1):n.nodes))
             
             # New data with permute
             simu.data_2 <- simu.data_1[,temp.order]
             n <- nrow (simu.data_2)    #Number of row
             V <- colnames(simu.data_2) #Column names
             
             # Calculate Pearson correlation
             suffStat<- list(C = cor(simu.data_2), n = n)
             
             # Infer the graph by MRPC
             MRPC_Inferred <- MRPC(simu.data_2,suffStat,GV=GV,FDR=0.05, FDRcontrol = TRUE,
                                   indepTest ='gaussCItest',labels=V,verbose = TRUE)
             # Recall and Precision by MRPC
             MRPC_Recall[i] <- RecallPrecision(Truth_1, MRPC_Inferred@graph, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Recall
             MRPC_Precision[i] <- RecallPrecision(Truth_1, MRPC_Inferred@graph, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Precision
             
             # Infer the graph by pc
             PC_Inferred <- pc(suffStat,alpha =0.05,
                               indepTest =gaussCItest,labels=V,verbose = TRUE)
             # Recall and Precision by pc
             PC_Recall[i] <- RecallPrecision(Truth_1, PC_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Recall
             PC_Precision[i] <- RecallPrecision(Truth_1,PC_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Precision
             
             # arcs not to be included from gene expression to genotype
             to <- rep (colnames (simu.data_2)[1:GV], each = (ncol (simu.data_2) - GV))
             from <- rep (colnames (simu.data_2)[(GV + 1):ncol (simu.data_2)], GV)
             bl <- cbind (from, to)
             
             # Infer the graph by pc.stable
             pc.stable.fit <- pc.stable(simu.data_2, blacklist=bl,alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
             
             # Inferred graph object by pc.stable
             G_pc.stable <- amat (pc.stable.fit)
             pc.stable_Inferred <- as(G_pc.stable,"graphNEL")

             
             # Recall and Precision by pc.stable
             pc.stable_Recall[i] <- RecallPrecision(Truth_1,pc.stable_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Recall
             pc.stable_Precision[i] <- RecallPrecision(Truth_1, pc.stable_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Precision
             
             
             # Infer the graph by mmpc
             mmpc.fit <- mmpc(simu.data_2, blacklist=bl,alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
             
             # Inferred graph object by mmpc
             G_mmpc <- amat (mmpc.fit)
             mmpc_Inferred <- as(G_mmpc,"graphNEL")
             
             # Recall and Precision by mmpc
             mmpc_Recall[i] <- RecallPrecision(Truth_1,mmpc_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Recall
             mmpc_Precision[i] <- RecallPrecision(Truth_1, mmpc_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Precision
                        
             # Infer the graph by mmhc
             mmhc.fit <- mmhc(simu.data_2,blacklist=bl) 
             
             # Inferred graph object by mmhc
             G_mmhc <- amat (mmhc.fit)
             mmhc_Inferred <- as(G_mmhc,"graphNEL")
             
             # Recall and Precision by mmhc
             mmhc_Recall[i] <- RecallPrecision(Truth_1, mmhc_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Recall
             mmhc_Precision[i] <- RecallPrecision(Truth_1, mmhc_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Precision
             
           }
         },
         
         model2 = {
           #Truth for model 2 (V1-->T1<--T2) with v-structure 
           
           tarmat_2 <- matrix(0,nrow=3,ncol = 3)
           colnames(tarmat_2) <- c("V1","T1","T2")
           rownames(tarmat_2) <- c("V1","T1","T2")
           # Adjacency matrix
           tarmat_2[1,2] <- 1
           tarmat_2[3,2] <- 1
           
           # Adjacency matrix from the graph by MRPC
           Truth_2 <- as(tarmat_2, "graphNEL")
           
           # The number of genetic variants in the graph.
           GV <- 1
           # The number of nodes
           n.nodes <- ncol (tarmat_2)
           
           #Ietaration for model 2
           for (i in 1:ita) {
             #Data for model 2
             simu.data_1 <- SimulateData(N = N,p = p,'model2',b0.1 = b0.1,b1.1 = b1.1,b1.2 = b1.2,b1.3 = b1.3,sd.1 = sd.1)
             # Create a new ordering for the T nodes.
             temp.order <- c (GV, sample((GV + 1):n.nodes))
             # New data with permute
             simu.data_2 <- simu.data_1[,temp.order]
             n <- nrow (simu.data_2)    # Number of row
             V <- colnames(simu.data_2) # Column names
             
             # Calculate Pearson correlation
             suffStat<- list(C = cor(simu.data_2), n = n)
             
             # Infer the graph by MRPC
             MRPC_Inferred <- MRPC(simu.data_2,suffStat,GV=GV,FDR=0.05, FDRcontrol = TRUE,
                                   indepTest ='gaussCItest',labels=V,verbose = TRUE)
             # Recall and Precision by MRPC
             MRPC_Recall[i] <- RecallPrecision(Truth_2, MRPC_Inferred@graph, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Recall
             MRPC_Precision[i] <- RecallPrecision(Truth_2, MRPC_Inferred@graph, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Precision
             
             # Infer the graph by pc
             PC_Inferred <- pc(suffStat,alpha =0.05,
                               indepTest =gaussCItest,labels=V,verbose = TRUE)
             # Recall and Precision by pc
             PC_Recall[i] <- RecallPrecision(Truth_2, PC_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Recall
             PC_Precision[i] <- RecallPrecision(Truth_2,PC_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Precision
             
             # arcs not to be included from gene expression to genotype
             to <- rep (colnames (simu.data_2)[1:GV], each = (ncol (simu.data_2) - GV))
             from <- rep (colnames (simu.data_2)[(GV + 1):ncol (simu.data_2)], GV)
             bl <- cbind (from, to)
             
             # Infer the graph by pc.stable
             pc.stable.fit <- pc.stable(simu.data_2, blacklist=bl,alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
             
             # Inferred graph object by pc.stable
             G_pc.stable <- amat (pc.stable.fit)
             pc.stable_Inferred <- as(G_pc.stable,"graphNEL")
             
             # Recall and Precision by pc.stable
             pc.stable_Recall[i] <- RecallPrecision(Truth_2,pc.stable_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Recall
             pc.stable_Precision[i] <- RecallPrecision(Truth_2, pc.stable_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Precision
             
             # Infer the graph by mmpc
             mmpc.fit <- mmpc(simu.data_2, blacklist=bl,alpha = 0.05, B = NULL, max.sx = NULL, debug = FALSE, undirected = FALSE)
             
             # Inferred graph object by mmpc
             G_mmpc <- amat (mmpc.fit)
             mmpc_Inferred <- as(G_mmpc,"graphNEL")
             
             
             # Recall and Precision by mmpc
             mmpc_Recall[i] <- RecallPrecision(Truth_2,mmpc_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Recall
             mmpc_Precision[i] <- RecallPrecision(Truth_2, mmpc_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Precision
             
             # Infer the graph by mmhc
             mmhc.fit <- mmhc(simu.data_2,blacklist=bl) 
             # Inferred graph object by mmhc
             G_mmhc <- amat (mmhc.fit)
             mmhc_Inferred <- as(G_mmhc,"graphNEL")
             
             # Recall and Precision by mmhc
             mmhc_Recall[i] <- RecallPrecision(Truth_2, mmhc_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Recall
             mmhc_Precision[i] <- RecallPrecision(Truth_2, mmhc_Inferred, GV=GV, includeGV=includeGV, edge.presence=1.0, edge.direction=0.5)$Precision
             
           }
         },
         stop("Model not included or missing")) 
  
  #Mean and sd Recall
  #MRPC
  Mean_Recall_MRPC <- mean(MRPC_Recall) #Mean
  SD_Recall_MRPC <- sd(MRPC_Recall)     #SD
  #pc
  Mean_Recall_PC <- mean(PC_Recall)    #mean
  SD_Recall_PC <- sd(PC_Recall)        #SD
  #pc.stable
  Mean_Recall_pc.stable <- mean(pc.stable_Recall) #mean
  SD_Recall_pc.stable <- sd(pc.stable_Recall)     #SD
  #mmpc
  Mean_Recall_mmpc <- mean(mmpc_Recall) #mean
  SD_Recall_mmpc <- sd(mmpc_Recall)     #SD
  #mmhc
  Mean_Recall_mmhc <- mean(mmhc_Recall) #mean
  SD_Recall_mmhc <- sd(mmhc_Recall)     #SD
  
  #Mean and sd Precision
  #MRPC
  Mean_Precision_MRPC <- mean(MRPC_Precision) #mean
  SD_Precision_MRPC <- sd(MRPC_Precision)     #SD
  #pc
  Mean_Precision_PC <- mean(PC_Precision)     #mean
  SD_Precision_PC <- sd(PC_Precision)         #SD
  #pc.stable
  Mean_Precision_pc.stable <- mean(pc.stable_Precision) #mean
  SD_Precision_pc.stable <- sd(pc.stable_Precision)     #SD
  #mmpc
  Mean_Precision_mmpc <- mean(mmpc_Precision) #mean
  SD_Precision_mmpc <- sd(mmpc_Precision)     #SD
  #mmhc
  Mean_Precision_mmhc <- mean(mmhc_Precision) #mean
  SD_Precision_mmhc <- sd(mmhc_Precision)     #SD
  
  #All outputs
  Outputs <- matrix(c(Mean_Recall_MRPC,SD_Recall_MRPC,Mean_Precision_MRPC,SD_Precision_MRPC,
                      Mean_Recall_PC,SD_Recall_PC,Mean_Precision_PC,SD_Precision_PC,
                      Mean_Recall_pc.stable,SD_Recall_pc.stable,Mean_Precision_pc.stable,SD_Precision_pc.stable,
                      Mean_Recall_mmpc,SD_Recall_mmpc,Mean_Precision_mmpc,SD_Precision_mmpc,
                      Mean_Recall_mmhc,SD_Recall_mmhc,Mean_Precision_mmhc,SD_Precision_mmhc
                      ),
                    nrow = 5,ncol = 4,byrow = T)
  
  colnames(Outputs) <- c("Mean_Recall","SD_Recall","Mean_Precision","SD_Pricision")
  rownames(Outputs) <- c("MRPC","pc","pc.stable","mmpc","mmhc")
  return(Outputs)
}
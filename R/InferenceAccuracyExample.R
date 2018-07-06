InferenceAccuracyExample<- function(N=1000,signal,model,ita=1000) {

#Truth for model 1 (V1-->T1-->T2) without v-structure 
tarmat_1=matrix(0,nrow=3,ncol = 3)
colnames(tarmat_1)=c("V1","T1","T2")
rownames(tarmat_1)=c("V1","T1","T2")

tarmat_1[1,2]=1
tarmat_1[2,3]=1
Truth_1<-as(tarmat_1, "graphNEL")

#Truth for model 2 (V1-->T1<--T2) with v-structure 
tarmat_2=matrix(0,nrow=3,ncol = 3)
colnames(tarmat_2)=c("V1","T1","T2")
rownames(tarmat_2)=c("V1","T1","T2")

tarmat_2[1,2]=1
tarmat_2[3,2]=1
Truth_2<-as(tarmat_2, "graphNEL")

#parameters setting
N=N     #sample size
p=0.45
signal=signal  # 0.2/0.5/1.0
b0.1 = 0
b1.1 =signal
b1.2 =signal
b1.3 =signal
sd.1 = 1

ita=ita # Iteration number
#Initial
MRPC_Recall=0
MRPC_Precision=0
mmhc_Recall=0
mmhc_Precision=0
PC_Recall=0
PC_Precision=0

switch(model,
       
model1 = {
#Ietaration for model 1
for (i in 1:ita) {
  #Data for model 1
  simu.data_1=SimulatedData(N = N,p = p,'model1',b0.1 = b0.1,b1.1 = b1.1,b1.2 = b1.2,b1.3 = b1.3,sd.1 = sd.1)
  #Permute
  GV=1
  temp.order<-c(GV,gtools::permute((GV+1):ncol(simu.data_1)))
  # New data with permute
  simu.data_2=simu.data_1[,temp.order]
  n<-nrow (simu.data_2)    #Number of row
  V<-colnames(simu.data_2) #Column names
  
  #Classical correlation (Beta=0)
  suffStat_C1= list(C = cor(simu.data_2), n = n)
  #Robust correlation (Beta=0.005)
  Rcor_R1=RobustCor(simu.data_2, 0.005) 
  suffStat_R1= list(C = Rcor_R1$RR, n = n)
  
  ## Estimated graph by MRPC using gaussCItest
  MRPC_Inferred<- MRPC(simu.data_2,suffStat_R1,GV=1,FDR=0.05,
                       indepTest ='gaussCItest',labels=V,verbose = TRUE)
  # Recall and Precision by MRPC
  MRPC_Recall[i]=Recall_Precision(Truth_1, MRPC_Inferred$graph, GV=1, edge.presence=1.0, edge.direction=0.5)$Recall
  MRPC_Precision[i]=Recall_Precision(Truth_1, MRPC_Inferred$graph, GV=1, edge.presence=1.0, edge.direction=0.5)$Precision
  
  ## Estimated graph by mmhc using gaussCItest
  M1=mmhc(simu.data_2) 
  mmhc_Inferred=graphviz.plot(M1)
  # Recall and Precision by mmhc
  mmhc_Recall[i]=Recall_Precision(Truth_1, mmhc_Inferred, GV=1, edge.presence=1.0, edge.direction=0.5)$Recall
  mmhc_Precision[i]=Recall_Precision(Truth_1, mmhc_Inferred, GV=1, edge.presence=1.0, edge.direction=0.5)$Precision
  ## Estimated graph by pc using gaussCItest
  PC_Inferred<- pc(suffStat_C1,alpha =0.05,
                   indepTest =gaussCItest,labels=V,verbose = TRUE)
  # Recall and Precision by mmhc
  PC_Recall[i]=Recall_Precision(Truth_1, PC_Inferred, GV=1, edge.presence=1.0, edge.direction=0.5)$Recall
  PC_Precision[i]=Recall_Precision(Truth_1,PC_Inferred, GV=1, edge.presence=1.0, edge.direction=0.5)$Precision
                }
                  },

model2 = {
  #Ietaration for model 2
  for (i in 1:ita) {
    #Data for model 2
    simu.data_1=SimulatedData(N = N,p = p,'model2',b0.1 = b0.1,b1.1 = b1.1,b1.2 = b1.2,b1.3 = b1.3,sd.1 = sd.1)
    #Permute
    GV=1
    temp.order<-c(GV,gtools::permute((GV+1):ncol(simu.data_1)))
    # New data with permute
    simu.data_2=simu.data_1[,temp.order]
    n<-nrow (simu.data_2)    #Number of row
    V<-colnames(simu.data_2) #Column names
    
    #Classical correlation (Beta=0)
    suffStat_C1= list(C = cor(simu.data_2), n = n)
    #Robust correlation (Beta=0.005)
    Rcor_R1=RobustCor(simu.data_2, 0.005) 
    suffStat_R1= list(C = Rcor_R1$RR, n = n)
    
    ## Estimated graph by MRPC using gaussCItest
    MRPC_Inferred<- MRPC(simu.data_2,suffStat_R1,GV=1,FDR=0.05,
                         indepTest ='gaussCItest',labels=V,verbose = TRUE)
    # Recall and Precision by MRPC
    MRPC_Recall[i]=Recall_Precision(Truth_2, MRPC_Inferred$graph, GV=1, edge.presence=1.0, edge.direction=0.5)$Recall
    MRPC_Precision[i]=Recall_Precision(Truth_2, MRPC_Inferred$graph, GV=1, edge.presence=1.0, edge.direction=0.5)$Precision
    
    ## Estimated graph by mmhc using gaussCItest
    M1=mmhc(simu.data_2) 
    mmhc_Inferred=graphviz.plot(M1)
    # Recall and Precision by mmhc
    mmhc_Recall[i]=Recall_Precision(Truth_2, mmhc_Inferred, GV=1, edge.presence=1.0, edge.direction=0.5)$Recall
    mmhc_Precision[i]=Recall_Precision(Truth_2, mmhc_Inferred, GV=1, edge.presence=1.0, edge.direction=0.5)$Precision
    ## Estimated graph by pc using gaussCItest
    PC_Inferred<- pc(suffStat_C1,alpha =0.05,
                     indepTest =gaussCItest,labels=V,verbose = TRUE)
    # Recall and Precision by mmhc
    PC_Recall[i]=Recall_Precision(Truth_2, PC_Inferred, GV=1, edge.presence=1.0, edge.direction=0.5)$Recall
    PC_Precision[i]=Recall_Precision(Truth_2,PC_Inferred, GV=1, edge.presence=1.0, edge.direction=0.5)$Precision
        }
          },
stop("Model not included or missing")) 
  
#Mean and sd Recall
Mean_Recall_MRPC=mean(MRPC_Recall) #Mean
SD_Recall_MRPC=sd(MRPC_Recall)     #SD

Mean_Recall_mmhc=mean(mmhc_Recall) #mean
SD_Recall_mmhc=sd(mmhc_Recall)     #SD

Mean_Recall_PC=mean(PC_Recall)    #mean
SD_Recall_PC=sd(PC_Recall)        #SD

#Mean and sd Precision
Mean_Precision_MRPC=mean(MRPC_Precision) #mean
SD_Precision_MRPC=sd(MRPC_Precision)     #SD

Mean_Precision_mmhc=mean(mmhc_Precision) #mean
SD_Precision_mmhc=sd(mmhc_Precision)     #SD

Mean_Precision_PC=mean(PC_Precision)     #mean
SD_Precision_PC=sd(PC_Precision)         #SD

#All outputs
Outputs=matrix(c(Mean_Recall_MRPC,SD_Recall_MRPC,Mean_Precision_MRPC,SD_Precision_MRPC,Mean_Recall_mmhc,SD_Recall_mmhc,Mean_Precision_mmhc,SD_Precision_mmhc,Mean_Recall_PC,SD_Recall_PC,Mean_Precision_PC,SD_Precision_PC),nrow = 3,ncol = 4,byrow = T)
colnames(Outputs)=c("Mean_Recal","SD_Recall","Mean_Precision","SD_Pricision")
rownames(Outputs)=c("MRPC","mmhc","PC")
return(Outputs)
}
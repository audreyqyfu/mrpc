#For software paper truth 1 
#Same data but different node ordering 
#Details please see help(SoftwarePaperSimData)
SoftwarePaperSimData<- function(N, model,signal,n_data,n_nodeordering) {
  
#parameters settings
N=N
p=0.45
signal=signal
b0.1 = 0
b1.1 =signal
b1.2 =signal
b1.3 =signal
sd.1 = 1

ita_data=n_data # Number of data sets
ita_node=n_nodeordering # Number of node ordering 

switch(model,
       
       truth1 = {
         #Truth 1 (V1-->T1-->T2-->T3)
         tarmat_s1=matrix(0,nrow=4,ncol = 4)
         colnames(tarmat_s1)=c("V1","T1","T2","T3")
         rownames(tarmat_s1)=colnames(tarmat_s1)
         #Create an adjacency matrix for truth
         tarmat_s1[1,2]=1
         tarmat_s1[2,3]=1
         tarmat_s1[3,4]=1
         #Graph
         Truth_s1<-as(tarmat_s1, "graphNEL")
         
         #Convert binary to decimal for truth
         library(compositions)  #for binary to decimal
         V1=as.vector(t(tarmat_s1))
         A1 <- paste(V1,collapse="")
         S_truth=unbinary(A1)
         
         #Output matrices
         Diff_MRPC=matrix(NA,ncol = ita_node,nrow = ita_data) #For MRPC
         Diff_mmhc=matrix(NA,ncol = ita_node,nrow = ita_data) #For mmhc
         Diff_PC=matrix(NA,ncol = ita_node,nrow = ita_data)   #For PC
         
         #Data
         for (i in 1:ita_data) {
           cat("ita_data=",i)
           V1 <- c(sample(c(0, 1, 2),size = N,replace = TRUE,prob = c((1 - p)^2,2*p*(1 - p),p^2)))
           T1=Case_1P(N=N,P1=V1,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
           T2=Case_1P(N=N,P1=T1,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
           T3=Case_1P(N=N,P1=T2,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
           
           #Data combine
           Data1=cbind(V1,T1,T2,T3)
           
           for (j in 1:ita_node) {
             cat("ita_node=",j)
             library(gtools) #Permute
             GV=1
             temp.order<-c(GV,permute((GV+1):ncol(Data1)))
             # New data with permute
             Data2=Data1[,temp.order]
             
             n<-nrow (Data2)    #Number of row
             V<-colnames(Data2) #Column names
             Rcor_R=RobustCor(Data2, 0.005) #Robust correlation (Beta=0.005)
             suffStat_R= list(C = Rcor_R$RR, n = n)
             
             Rcor_C=RobustCor(Data2, 0) #Classical correlation (Beta=0)
             suffStat_C= list(C = Rcor_C$RR, n = n)
             
             ## Estimated graph by MRPC using gaussCItest
             MRPC.fit_1 <- MRPC(Data2,suffStat_R,GV=1,FDR=0.05,
                                indepTest ='gaussCItest',labels=V,verbose = TRUE)
             #adjacency matrix from directed graph by MRPC
             Adj_MRPC=as(MRPC.fit_1$graph,"matrix")
             
             #If ordering nodes
             if(any(colnames(tarmat_s1)!=colnames(Adj_MRPC)))
             {
               Order_node=match(colnames(tarmat_s1),colnames(Adj_MRPC))
               Adj_MRPC<-Adj_MRPC[Order_node,Order_node] #new
             }
             
             #Convert binary to decimal for MRPC
             V2_MRPC=as.vector(t(Adj_MRPC))
             A2_MRPC<- paste(V2_MRPC,collapse="")
             S_MRPC=unbinary(A2_MRPC)
             #Difference between truth and inferred
             Diff_MRPC[i,j]=S_truth-S_MRPC
             
             ## Estimated graph by mmhc
             M_mmhc=mmhc(data.frame(Data2)) 
             G_mmhc=graphviz.plot(M_mmhc)
             #adjacency matrix from mmhc 
             Adj_mmhc=as(G_mmhc,"matrix")
             
             #If ordering nodes
             if(any(colnames(tarmat_s1)!=colnames(Adj_mmhc)))
             {
               Order_node=match(colnames(tarmat_s1),colnames(Adj_mmhc))
               Adj_mmhc<-Adj_mmhc[Order_node,Order_node] #new
             }
             #Convert binary to decimal for mmhc
             V2_mmhc=as.vector(t(Adj_mmhc))
             A2_mmhc<- paste(V2_mmhc,collapse="")
             S_mmhc=unbinary(A2_mmhc)
             #Difference between truth and inferred
             Diff_mmhc[i,j]=S_truth-S_mmhc
             
             ## Estimated graph by PC
             PC.fit_1 <- pc(suffStat_C,alpha =0.05,
                            indepTest =gaussCItest,labels=V,verbose = TRUE)
             #adjacency matrix from PC 
             Adj_PC=as(PC.fit_1@graph,"matrix")
             
             #If ordering nodes
             if(any(colnames(tarmat_s1)!=colnames(Adj_PC)))
             {
               Order_node=match(colnames(tarmat_s1),colnames(Adj_PC))
               Adj_PC<-Adj_PC[Order_node,Order_node] #new
             }
             #Convert binary to decimal for PC
             V2_PC=as.vector(t(Adj_PC))
             A2_PC<- paste(V2_PC,collapse="")
             S_PC=unbinary(A2_PC)
             #Difference between truth and inferred
             Diff_PC[i,j]=S_truth-S_PC
             
             # write results to csv file
             #row indicates the number of data and col indicates number of node ordering by MRPC,mmhc and PC respectively
             write.table (cbind (Diff_MRPC,Diff_mmhc,Diff_PC), file = "Results_1.csv",sep=",",row.names = F,col.names = F)
             
           }
         }  
       },
       truth2 = {
         #Truth 2 (V1-->T1<--T2-->T3)
         tarmat_s2=matrix(0,nrow=4,ncol = 4)
         colnames(tarmat_s2)=c("V1","T1","T2","T3")
         rownames(tarmat_s2)=colnames(tarmat_s2)
         #Create an adjacency matrix for truth
         tarmat_s2[1,2]=1
         tarmat_s2[3,2]=1
         tarmat_s2[3,4]=1
         #Graph
         Truth_s2<-as(tarmat_s2, "graphNEL")
         
         #Convert binary to decimal for truth
         library(compositions)  #for binary to decimal
         V2=as.vector(t(tarmat_s2))
         A2 <- paste(V2,collapse="")
         S_truth=unbinary(A2)
         
         #Output matrices
         Diff_MRPC=matrix(NA,ncol = ita_node,nrow = ita_data) #For MRPC
         Diff_mmhc=matrix(NA,ncol = ita_node,nrow = ita_data) #For mmhc
         Diff_PC=matrix(NA,ncol = ita_node,nrow = ita_data)   #For pc
         
         #Data
         for (i in 1:ita_data) {
           cat("ita_data=",i)
           #Data
           V1 <- c(sample(c(0, 1, 2),size = N,replace = TRUE,prob = c((1 - p)^2,2*p*(1 - p),p^2)))
           T2=Case_NP(N=N,b0.1=b0.1,sd.1=sd.1)
           T1=Case_2P(N=N,P1=V1,P2=T2,b0.1=b0.1,b1.1=b1.1,b1.2=b1.2,sd.1=sd.1)
           T3=Case_1P(N=N,P1=T2,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
           #Data combine
           Data1=cbind(V1,T1,T2,T3)
           
           for (j in 1:ita_node) {
             cat("ita_node=",j)
             library(gtools) #Permute
             GV=1
             temp.order<-c(GV,permute((GV+1):ncol(Data1)))
             # New data with permute
             Data2=Data1[,temp.order]
             
             n<-nrow (Data2)    #Number of row
             V<-colnames(Data2) #Column names
             Rcor_R=RobustCor(Data2, 0.005) #Robust correlation (Beta=0.005)
             suffStat_R= list(C = Rcor_R$RR, n = n)
             
             Rcor_C=RobustCor(Data2, 0) #Classical correlation (Beta=0)
             suffStat_C= list(C = Rcor_C$RR, n = n)
             
             ## Estimated graph by MRPC using gaussCItest
             MRPC.fit_1 <- MRPC(Data2,suffStat_R,GV=1,FDR=0.05,
                                indepTest ='gaussCItest',labels=V,verbose = TRUE)
             #adjacency matrix from directed graph by MRPC
             Adj_MRPC=as(MRPC.fit_1$graph,"matrix")
             
             #If ordering nodes
             if(any(colnames(tarmat_s2)!=colnames(Adj_MRPC)))
             {
               Order_node=match(colnames(tarmat_s2),colnames(Adj_MRPC))
               Adj_MRPC<-Adj_MRPC[Order_node,Order_node] #new
             }
             
             #Convert binary to decimal for MRPC
             V2_MRPC=as.vector(t(Adj_MRPC))
             A2_MRPC<- paste(V2_MRPC,collapse="")
             S_MRPC=unbinary(A2_MRPC)
             #Difference between truth and inferred
             Diff_MRPC[i,j]=S_truth-S_MRPC
             
             ## Estimated graph by mmhc
             M_mmhc=mmhc(data.frame(Data2)) 
             G_mmhc=graphviz.plot(M_mmhc)
             #adjacency matrix from mmhc 
             Adj_mmhc=as(G_mmhc,"matrix")
             
             #If ordering nodes
             if(any(colnames(tarmat_s2)!=colnames(Adj_mmhc)))
             {
               Order_node=match(colnames(tarmat_s2),colnames(Adj_mmhc))
               Adj_mmhc<-Adj_mmhc[Order_node,Order_node] #new
             }
             #Convert binary to decimal for mmhc
             V2_mmhc=as.vector(t(Adj_mmhc))
             A2_mmhc<- paste(V2_mmhc,collapse="")
             S_mmhc=unbinary(A2_mmhc)
             #Difference between truth and inferred
             Diff_mmhc[i,j]=S_truth-S_mmhc
             
             ## Estimated graph by PC
             PC.fit_1 <- pc(suffStat_C,alpha =0.05,
                            indepTest =gaussCItest,labels=V,verbose = TRUE)
             #adjacency matrix from PC 
             Adj_PC=as(PC.fit_1@graph,"matrix")
             
             #If ordering nodes
             if(any(colnames(tarmat_s2)!=colnames(Adj_PC)))
             {
               Order_node=match(colnames(tarmat_s2),colnames(Adj_PC))
               Adj_PC<-Adj_PC[Order_node,Order_node] #new
             }
             
             #Convert binary to decimal for PC
             V2_PC=as.vector(t(Adj_PC))
             #Convert 
             A2_PC<- paste(V2_PC,collapse="")
             S_PC=unbinary(A2_PC)
             #Difference between truth and inferred
             Diff_PC[i,j]=S_truth-S_PC
             # write results to csv file
             #row indicates the number of data and col indicates number of node ordering by MRPC,mmhc and PC respectively
             write.table (cbind (Diff_MRPC,Diff_mmhc,Diff_PC), file = "Results_2.csv",sep=",",row.names = F,col.names = F)
           }
         }
       },
       stop("Model not included or missing"))
}

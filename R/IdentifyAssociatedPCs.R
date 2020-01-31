# This is the code to identify PCs that are significantly associated with eQTL-gene sets

IdentifyAssociatedPCs <- function (PCs.matrix,no.PCs=10,data,fdr.level=0.05,corr.threshold=FALSE,corr.value=0.3) {
  
  # Compute the correlation and corresponding p values between the top PCs and the eQTLs and genes
  # library(psych) # to use corr.test
  # no.PCs <- no.PCs  
  corr.PCs <- corr.test(PCs.matrix[,1:no.PCs],data)
  
  # The correlation matrix
  corr.matrix <- corr.PCs$r
  
  # The p values
  Pvalues <- corr.PCs$p
  # Apply the q value method at FDR of 0.05
  # library(WGCNA) # qvalue
  # qobj <- qvalue(Pvalues, fdr.level=0.05,robust = TRUE) 
  qobj <- qvalue(Pvalues, fdr.level=fdr.level) 
  
  # Significant associations
  Significant.asso <- qobj$significant
  List.significant.asso1 <- which(Significant.asso, arr.ind = TRUE, useNames = TRUE)
  # 1st column contains the PCs
  # 2nd column contains the associated eQTLs or genes name
  if(length(List.significant.asso1)!=0)
  {
    
    if (corr.threshold)
    { 
      # cellnote <- ifelse(abs(corr.matrix)>0.30,'Y','')
      # High.cor.index <- which(cellnote=="Y",arr.ind = T) 
      
      High.cor.index <- which(abs(corr.matrix)>corr.value,arr.ind = T)
      
      # library(plyr) # join.keys
      
      Match.index <- with(join.keys(data.frame(List.significant.asso1),data.frame(High.cor.index)), which(x %in% y))
      List.significant.asso <- data.frame(List.significant.asso1)[Match.index,]
      #List.significant.asso <- List.significant.asso2
    }
    else
    {
      List.significant.asso <- List.significant.asso1
    }
    
    # eQTLs or genes that are significantly associated with selected PCs
    PC.list <- list()

    for (m in 1:no.PCs) {
      PC.list[[m]] <- colnames(data)[List.significant.asso[which(List.significant.asso[,1]==paste(m)),2]]
    }
    names(PC.list) <- paste("PC",1:no.PCs,sep="") 
  
    PCs.asso.list <- list()
    for (i in 1:length(colnames(data))) {
      toMatch <- colnames(data)[i]
      #lapply(l1, function(x) all(toMatch %in% x))
      pc <- names(PC.list)[sapply(PC.list, function(x) all(toMatch %in% x))]
      #pc <- unique(unlist(PC.list))[sapply(unique(unlist(PC.list)), function(x) all(toMatch %in% x))]
      #if(length(pc)!=0){
        PCs.asso.list[[i]] <- pc
      #}
      
    }
    AssociatedPCs <- sort(unique(na.omit(unlist(PCs.asso.list))))
    
    if(length(AssociatedPCs)!=0)
    {
      Match.index.pcmatrix <- as.vector(match(AssociatedPCs,colnames(PCs.matrix)))
      data.withPC <- cbind(data,PCs.matrix[,Match.index.pcmatrix])
      colnames(data.withPC)[(ncol(data)+1):(ncol(data)+length(Match.index.pcmatrix))] <-  AssociatedPCs
    }
    else
    {
      AssociatedPCs <- NULL
      PCs.asso.list <- NULL
      data.withPC <- data
    }
  }
  else
  {
    AssociatedPCs <- NULL
    PCs.asso.list <- NULL
    data.withPC <- data
  }

return(list(AssociatedPCs=AssociatedPCs,
            data.withPC=data.withPC,
            corr.PCs=corr.matrix,
            PCs.asso.list=PCs.asso.list,
            qobj=qobj))
}


getPurity <- function(tumor.data,normal.data = NULL,tumor.type = NULL)
{
  purity = c()
  
  # if only having a few tumor samples, you should specify tumor type and load pre-defined iDMCs. 
  if (is.vector(tumor.data)){
    # one tumor sample
    if (is.null(tumor.type)){
      stop("should specify tumor.type if having only one tumor sample")
    }
    load("/mnt/Storage/home/zhengxq/zhengxq/dna_methylation/gb_revision/call_purity/iDMC.RData")
    if ( ! any(names(iDMC) == tumor.type)) {
      ## tumor.type is not in iDMC
      stop(paste0(tumor.type," is not found! See abbr.txt for detail"))
    }
    
    tumor.CpGs = names(tumor.data)
    comm.CpGs = intersect(tumor.CpGs,rownames(iDMC[[tumor.type]]))
    
    idmc.dat = data.frame(tumor.data[comm.CpGs])
    idmc.dat$hyper = iDMC[[tumor.type]][comm.CpGs,]$hyper
    
    cat("Calculating tumor purity ...\n")
    beta.adj = c(idmc.dat[idmc.dat$hyper == TRUE,1],1-idmc.dat[idmc.dat$hyper == FALSE,1])
    purity = .get_peak(beta.adj)
  }
  
  else if (ncol(tumor.data) < 20) {
    # multiple tumor samples, but less than 20
    cat(paste0("The number of tumor samples is: ",ncol(tumor.data),"\n"))
    tumor.sample = colnames(tumor.data)
    tumor.CpGs = rownames(tumor.data)
    
    if (is.null(tumor.type)){
      stop("should specify tumor.type if less than 20 tumor samples")
    }
    load("/mnt/Storage/home/zhengxq/zhengxq/dna_methylation/gb_revision/call_purity/iDMC.RData")
    if ( ! any(names(iDMC) == tumor.type)) {
      ## tumor.type is not in iDMC
      stop(paste0(tumor.type," is not found! See abbr.txt for detail"))
    }
    
    comm.CpGs = intersect(tumor.CpGs,rownames(iDMC[[tumor.type]]))

    idmc.dat = tumor.data[comm.CpGs,]
    idmc.dat$hyper = iDMC[[tumor.type]][comm.CpGs,]$hyper
    
    cat("Calculating tumor purity ...\n")
    for(t in tumor.sample){
      beta.adj = c(idmc.dat[idmc.dat$hyper == TRUE,t],1-idmc.dat[idmc.dat$hyper == FALSE,t])
      pu = .get_peak(beta.adj)
      purity[t] = pu
    }
  }
  
  else{
    cat(paste0("The number of tumor samples is: ",ncol(tumor.data),"\n"))
    # if we have enough tumor samples, get iDMCs and call purity
    if(is.null(normal.data)){
      # If no normal sample, use universal normals from 21 TCGA cancer types 
      cat("Loading univeral normals ...\n")
      load("/mnt/Storage/home/zhengxq/zhengxq/dna_methylation/gb_revision/normal43.RData") # should change in R package
      normal.data = n.betamatrix
    }
    else if (ncol(normal.data) < 20){
      # If less than 20 normal samples, use universal normals from 21 TCGA cancer types 
      cat("Loading univeral normals ...\n")
      load("/mnt/Storage/home/zhengxq/zhengxq/dna_methylation/gb_revision/normal43.RData") # should change in R package
      normal.data = n.betamatrix
    }
    normal.sample = colnames(normal.data)
    normal.CpGs = rownames(normal.data)
    
    tumor.sample = colnames(tumor.data)
    tumor.CpGs = rownames(tumor.data)
    
    comm.CpGs = intersect(normal.CpGs,tumor.CpGs)
    
    cat("Getting iDMCs ...\n")
    iDMC = get_iDMC(tumor.data[comm.CpGs,],normal.data[comm.CpGs,])
    
    idmc.dat = cbind(tumor.data[iDMC,tumor.sample],normal.data[iDMC,normal.sample])
    idmc.dat$hyper = rowMeans(idmc.dat[,tumor.sample],na.rm = TRUE) > rowMeans(idmc.dat[,normal.sample],na.rm = TRUE)
    
    cat("Calculating tumor purity ...\n")
    for(t in tumor.sample){
      beta.adj = c(idmc.dat[idmc.dat$hyper == TRUE,t],1-idmc.dat[idmc.dat$hyper == FALSE,t])
      pu = .get_peak(beta.adj)
      purity[t] = pu
    }
  }

  return(purity)
}


.get_peak <- function(dat){
  d = density(dat,na.rm = TRUE,kernel = "gaussian")
  d$x[which.max(d$y)]
}


.myVar <- function(x){
  var(x,na.rm = TRUE)
}

get_iDMC <- function(tumor.data,normal.data){
  
  all.dat = cbind(tumor.data,normal.data)
  tumor.sample = colnames(tumor.data)
  normal.sample = colnames(normal.data)
  
  .myRanksum <- function(x){
    ## only works within function get_iDMC 
    if (length(na.omit(x[tumor.sample])) == 0 | length(na.omit(x[normal.sample])) == 0){
      return(NA)
    }
    else{
      pval = wilcox.test(as.numeric(x[tumor.sample]),as.numeric(x[normal.sample]))$p.value
      return(pval)
    }
  }

  ranksum.pval = apply(all.dat,1,.myRanksum)
  tumor.var = apply(all.dat[,tumor.sample],1,.myVar)
  
  out = data.frame(ranksum.pval)
  out = cbind(out,tumor.var)
  
  cDMC = out[out$tumor.var >= 0.005,]
  idx = order(cDMC$ranksum.pval,decreasing=FALSE)[1:1000]
  iDMC = rownames(cDMC)[idx]
  
  return(iDMC)
}



load_hormone_annotation <- function(hormSet){
  # read the ahd annotations for genetic evidence
  # ahd <- read.xlsx("../20141008 AHD List conversion/AHD_genetice_evidence.xlsx",sheetIndex=1,stringsAsFactors=F)
  # ahd$AGI <- toupper(trim(ahd$AGI))
  # save(ahd,file="ahd.RData")
  load("ahdwoCIPK14.RData")
  ahdColnames <- colnames(ahd)[1:9]
  ahdColnamesL <- colnames(ahd)[1:10] # Changed [1:9] to [1:10] 2016/12/09
  
  # read the ahd annotations for GO evidence
  # ahdGo <- read.xlsx("../20141008 AHD List conversion/AHD_gene_ontology_evidence.xlsx",sheetIndex=1,stringsAsFactors=F)
  # colnames(ahdGo) <- colnames(ahd)
  # ahdGo$AGI <- toupper(trim(ahdGo$AGI))
  # save(ahdGo,file="ahdGo.RData")
  load("ahdGo.RData")
  
  # read the tair go annotations
  # tairGo <- read.xlsx("../20160607 GO hormone annotation/TAIR_GO_evidence.xlsx",sheetIndex = 1,stringsAsFactors = F)
  # tairGo$AGI <- toupper(trim(tairGo$AGI))
  # save(tairGo,file="tairGo.RData")
  # load("tairGo.RData") # 2019/10/18 replaced by file ../20160607 GO hormone annotation/Results 20181212/TAIR_GO_HORMONE_ANNOTATION.RData
  
  # 
  # 2019/10/18 read the tair go annotations also used for combined interactions
  load("TAIR_GO_HORMONE_ANNOTATION.RData")
  tairGo <- ahd.df
  rm(ahd.df)
  
  # combine annotations
  if(hormSet == "AHD_Genetic"){
    # nothing to do
  }else if(hormSet == "AHD_GO"){
    ahd <- ahdGo
  }else if(hormSet == "AHD_ALL"){
    lociList <- unique(c(ahd$AGI,ahdGo$AGI))
    newMatrix <- c()
    for(myLocus in lociList){
      myIndexA <- which(ahd$AGI == myLocus)  
      myIndexB <- which(ahdGo$AGI == myLocus)  
      newVec <- c()
      if(length(myIndexA) == 1 & length(myIndexB) == 1){
        newVec <- ahd[myIndexA,2:9] | ahdGo[myIndexB,2:9]
      }else if(length(myIndexA) == 1 & length(myIndexB) == 0){
        newVec <- ahd[myIndexA,2:9]
      }else if(length(myIndexA) == 0 & length(myIndexB) == 1){
        newVec <- ahdGo[myIndexB,2:9]
      }else{
        stop("No index found")
      }
      
      hIndex <- which(newVec == 1)
      if(length(hIndex) == 1){
        summary <- ahdColnames[2:9][hIndex]
      }else{
        summary <- "multiple"
      }
      
      newMatrix <- rbind(newMatrix,data.frame(myLocus,newVec,summary))
      
    }
    
    ahd <- newMatrix
    colnames(ahd) <- ahdColnamesL
  }else if(hormSet == "TAIR_GO"){
    ahd <- tairGo
    colnames(ahd) <- ahdColnamesL
  }else if(hormSet == "AHDGen_TAIR"){
    lociList <- unique(c(ahd$AGI,tairGo$AGI))
    newMatrix <- c()
    for(myLocus in lociList){
      myIndexA <- which(ahd$AGI == myLocus)  
        
      myIndexC <- which(tairGo$AGI == myLocus)
      newVec <- c()
      if(length(myIndexA) == 1 & length(myIndexC) == 1){
        newVec <- ahd[myIndexA,2:9] | tairGo[myIndexC,2:9]
      }else if(length(myIndexA) == 1){
        newVec <- ahd[myIndexA,2:9]
      }else if(length(myIndexC) == 1){
        newVec <- tairGo[myIndexC,2:9]
      }else{
        stop("No index found")
      }
      
      hIndex <- which(newVec == 1)
      if(length(hIndex) == 1){
        summary <- ahdColnames[2:9][hIndex]
      }else{
        summary <- "multiple"
      }
      
      newMatrix <- rbind(newMatrix,data.frame(myLocus,newVec,summary))
      
    }
    
    ahd <- newMatrix
    colnames(ahd) <- ahdColnamesL
  }else if(hormSet == "ALL"){
    lociList <- unique(c(ahd$AGI,ahdGo$AGI,tairGo$AGI))
    newMatrix <- c()
    for(myLocus in lociList){
      myIndexA <- which(ahd$AGI == myLocus)  
      myIndexB <- which(ahdGo$AGI == myLocus)  
      myIndexC <- which(tairGo$AGI == myLocus)
      newVec <- c()
      if(length(myIndexA) == 1 & length(myIndexB) == 1  & length(myIndexC) == 1){
        newVec <- ahd[myIndexA,2:9] | ahdGo[myIndexB,2:9] | tairGo[myIndexC,2:9]
      }else if(length(myIndexA) == 1 & length(myIndexB) == 1){
        newVec <- ahd[myIndexA,2:9] | ahdGo[myIndexB,2:9]
      }else if(length(myIndexA) == 1 & length(myIndexC) == 1){
        newVec <- ahd[myIndexA,2:9] | tairGo[myIndexC,2:9]
      }else if(length(myIndexB) == 1 & length(myIndexC) == 1){
        newVec <- ahdGo[myIndexB,2:9] | tairGo[myIndexC,2:9]
      }else if(length(myIndexA) == 1){
        newVec <- ahd[myIndexA,2:9]
      }else if(length(myIndexB) == 1){
        newVec <- ahdGo[myIndexB,2:9]
      }else if(length(myIndexC) == 1){
        newVec <- tairGo[myIndexC,2:9]
      }else{
        stop("No index found")
      }
      
      hIndex <- which(newVec == 1)
      if(length(hIndex) == 1){
        summary <- ahdColnames[2:9][hIndex]
      }else{
        summary <- "multiple"
      }
      
      newMatrix <- rbind(newMatrix,data.frame(myLocus,newVec,summary))
      
    }
    
    ahd <- newMatrix
    colnames(ahd) <- ahdColnamesL
  }else{
    stop("Hormone set not found")
  }
  
  ### add a piechart column to the hormone list
  # piechart: valuelist=" 0,0,0,0,0,1,0,0 " colorlist="#00CCCC,#0066CC,#6600CC,#CC00CC,#CC0000,#CC6600,#CCCC00,#00CC00" labellist="abscisic_acid,auxin,brassinosteroid,cytokinin,ethylene,gibberellin,jasmonic_acid,salicylic_acid"
  # piechart: valuelist=" 0,0,0,0,0,1,0,0 " colorlist="#00CCCC,#0066CC,#6600CC,#CC00CC,#CC0000,#CC6600,#CCCC00,#00CC00" labellist="aba,aux,br,ck,et,ga,ja,sa"
  pieVec <- c()
  for(l in 1:length(ahd[,1])){
    # pieTemp <- paste("piechart: valuelist=\"", paste(ahd[l,2:9],collapse = ",") ,"\" colorlist=\"#00CCCC,#0066CC,#6600CC,#CC00CC,#CC0000,#CC6600,#CCCC00,#00CC00\" labellist=\"aba,aux,br,ck,et,ga,ja,sa\"",sep="")
    pieTemp <- paste("piechart: valuelist=\"", paste(ahd[l,2:9],collapse = ",") ,"\" colorlist=\"#00CCCC,#0066CC,#6600CC,#CC00CC,#CC0000,#CC6600,#CCCC00,#00CC00\" labellist=\"\"",sep="")
    pieVec <- c(pieVec,pieTemp)
  }
  
  ahd <- cbind(ahd,pieVec)
  ahd$summary <- sub("_"," ",ahd$summary)
  return(ahd)
}




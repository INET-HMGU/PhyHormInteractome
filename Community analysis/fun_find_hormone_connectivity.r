
suppressPackageStartupMessages(library(igraph))

find_hormone_connectivity <- function(my.graph,ahd,dbg){
  ahd <- ahd[,1:9]
  # find number of interactions between the hormones
  connectivityList <- c()
  connectivityListNorm <- c()
  hormoneNodes <- c()
  for(i in 2:9){
    if(dbg){print(paste("i",i))}
    myLoci <- ahd[which(ahd[,i] == 1),1]
    myLoci <- myLoci[which(myLoci %in% names(V(my.graph)))]
    hormoneNodes <- c(hormoneNodes,length(myLoci))
    myTempCon <- rep(0,8)
    if(length(myLoci) > 0){
      for(j in 1:length(myLoci)){
        if(dbg){print(paste("j",j))}
        # if(is.na(myLoci[j])){
          # do nothing
        # }else{
          myNei <- names(neighbors(my.graph,as.character(myLoci[j])))
          myIndex <- which(ahd$AGI %in% myNei)
          for(k in 1:8){
            if(dbg){print(paste("k",k))}
            # print(length(myIndex))
            myCons <- length(which(ahd[myIndex,(k+1)] == 1))
            myTempCon[k] <- myTempCon[k] + myCons
          }
        # }
        
      }
    
    }
    connectivityList <- rbind(connectivityList,myTempCon)
    connectivityListNorm <- rbind(connectivityListNorm,myTempCon/length(myLoci))
  }
  rownames(connectivityList) <- seq(1,8)
  connectivityList <- data.frame(connectivityList)
  colnames(connectivityList) <- colnames(ahd)[2:9]
  rownames(connectivityList) <- colnames(ahd)[2:9]
  
  if(dbg){print(connectivityList);print(connectivityListNorm)}
  
  # convert to matrix for cytoscape
  myNames <- colnames(ahd)[2:9]
  conCyto <- c()
  # conCytoNorm <- c()
  nodeSize <- c()
  nodeSize <- data.frame(colnames(ahd)[2:9],hormoneNodes)
  # nodeSizeNorm <- c()
  for(i in 1:8){
    for(j in i:8){
      # if(i == j){
        # nodeSize <- rbind(nodeSize,data.frame(myNames[i],connectivityList[i,j]))
        # nodeSizeNorm <- rbind(nodeSizeNorm,data.frame(myNames[i],connectivityListNorm[i,j]))
      # }else{
        conCyto <- rbind(conCyto,data.frame(myNames[i],myNames[j],connectivityList[i,j]))
        # conCytoNorm <- rbind(conCytoNorm,data.frame(myNames[i],myNames[j],connectivityListNorm[i,j]))
      # }
    }
  }
  
  colnames(conCyto) <- c("Horm1","Horm2","Ext_Connectivity")
  # colnames(conCytoNorm) <- c("Horm1","Horm2","Ext_Connectivity")
  
  colnames(nodeSize) <- c("Hormone","Number of Loci")
  # colnames(nodeSizeNorm) <- c("Hormone","Int_Connectivity")
  
  write.xlsx(conCyto,paste(result.dir,file.name," - communities.xlsx",sep=""),sheetName = "Hormone connectivity - Edges",append = T)
  write.xlsx(nodeSize,paste(result.dir,file.name," - communities.xlsx",sep=""),sheetName = "Hormone connectivity - Nodes",append = T)
  # write.xlsx(conCytoNorm,paste(result.dir,file.name," - communities.xlsx",sep=""),sheetName = "Hormone connectivity norm - ext",append = T)
  # write.xlsx(nodeSizeNorm,paste(result.dir,file.name," - communities.xlsx",sep=""),sheetName = "Hormone connectivity norm - int",append = T)
  
  
  
}

suppressPackageStartupMessages(library(gplots)) # for colorpanel()
suppressPackageStartupMessages(library(lattice)) # for levelplot()
suppressPackageStartupMessages(library(RColorBrewer)) # for colors in levelplot
suppressPackageStartupMessages(library(latticeExtra))

# function to find largest connected graph in network
giant.component <- function(graph) {
  cl <- components(graph)
  induced_subgraph(graph, which(cl$membership == which.max(cl$csize)))
}

find_shortest_paths_each <- function(my.graph,ahd,dbg){
  ahd <- ahd[,1:9]
  # remove all non-connected components from the biggest connected component
  # the non connected components would receive an artificial small mean shortest path
  my.graph.small <- giant.component(my.graph)
  
  # determine all nodes in the graph
  vNames <- vertex_attr(my.graph.small,"name")
  
  # find node annotated with 
  myShortestPaths <- list()
  myShortestPaths2 <- list()
  
  hormIndex <- c(2:9) # helper vector
  
  myCounter <- 1 # internal counter for outer loop
  for(myHorm in 2:9){
    # all AGIs annotated with the same hormone
    if(dbg){print(colnames(ahd)[myHorm])}
    myIndex <- which(ahd[,myHorm] == 1)
    myVertices <- ahd$AGI[myIndex] 
    myNames <- vNames[which(vNames %in% myVertices)] # which of the annotated AGIs are in the network
    
    # all AGIs annotated with another hormone
    ahdTemp <- ahd[-which(ahd$AGI %in% ahd$AGI[myIndex]),]
    # ahdTemp <- ahd[-myIndex,] # remove loci of proteins annotated with the current hormone
    myNamesList <- list() # store loci of differently annotated proteins
    myCounterHormIndex <- 1
    for(j in hormIndex[-myCounter]){
      if(dbg){print(j)}
      myIndex2 <- which(ahdTemp[,j] == 1)
      myVertices2 <- unique(ahdTemp$AGI[myIndex2])
      myNamesList[[myCounterHormIndex]] <- vNames[which(vNames %in% myVertices2)]
      myCounterHormIndex <- myCounterHormIndex + 1
    }
    
    # internal connections can only be calculated, if at least two nodes with the same annotation
    # are available
    if(length(myNames) < 2){
      myShortestPaths[[myHorm-1]] <- NA
      myShortestPaths2[[myHorm-1]] <- NA
    }else{
      myAllLengths <- c() # same annotated results
      myAllLengths2List <- list() # differently annotated results
      for(i in 1:(length(myNames))){
        # same annotation
        if(i < length(myNames)){
          for(j in (i+1):length(myNames)){
            # print(paste("from",myNames[i],"to",myNames[j], sep = ""))
            myShortestPath <- shortest_paths(graph = my.graph.small,from = myNames[i],to = myNames[j])
            myLength <- length(myShortestPath$vpath[[1]])
            # print(myLength)
            myAllLengths <- c(myAllLengths,myLength-1)
          }
        }
        # test differently annotated proteins
        myCounterNames <- 1
        # loop over the set of AGIs from the other hormones
        for(myNames2 in myNamesList){
          if(dbg){print(myNames2)}
          myAllLengths2 <- c()
          if(length(myNames2) > 0){
            for(j in 1:length(myNames2)){ # test length 
              # print(paste("from",myNames[i],"to",myNames2[j], sep = ""))
              myShortestPath <- shortest_paths(graph = my.graph.small,from = myNames[i],to = myNames2[j])
              myLength <- length(myShortestPath$vpath[[1]])
              # print(myLength)
              myAllLengths2 <- c(myAllLengths2,myLength-1)
            }
          }else{
            myAllLengths2 <- c(myAllLengths2,NA)
          }
          
          # test, if the list already contains a value
          if(length(myAllLengths2List) < 7){
            myAllLengths2List[[myCounterNames]] <- myAllLengths2
          }else{
            myAllLengths2List[[myCounterNames]] <- c(myAllLengths2List[[myCounterNames]],myAllLengths2)
          }
          myCounterNames <- myCounterNames + 1
        }
        
      }
      myShortestPaths[[myHorm-1]] <- myAllLengths
      myShortestPaths2[[myHorm-1]] <- myAllLengths2List
    }
    
    myCounter <- myCounter + 1
  }
  
  pdf(paste(result.dir,file.name," - Hormone Distance Boxplots.pdf",sep=""),width = 10, height = 10)
  par(mar = c(10,4,4,4),mfrow=c(2,1))
  boxplot(myShortestPaths, names = colnames(ahd)[2:9],las=2,main="Shortest Path length to proteins annotated with same hormone")
  stripchart(myShortestPaths,vertical=T,add=T,method = "jitter",pch = 20, col = "blue",cex = 0.5)
  # stripchart(myShortestPaths, vertical = TRUE,
  #            method = "jitter", jitter = 0.1, add = TRUE, pch = 20, col = 'blue')
  # 
  myShortestPathsComb <- lapply(myShortestPaths2,unlist)
  
  boxplot(myShortestPathsComb, names = colnames(ahd)[2:9],las=2,main="Shortest Path length to proteins annotated with different hormones")
  stripchart(myShortestPathsComb,vertical=T,add=T,method = "jitter",pch = 20, col = "blue",cex = 0.5)
  dev.off()
  # library("beanplot")
  # par(mar = c(10,4,4,4))
  # beanplot(myShortestPaths,col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6",
  #          names = colnames(ahd)[2:9],las=2)
  # 
  # library("vioplot")
  # vioplot(myShortestPaths[[1]],myShortestPaths[[2]],myShortestPaths[[3]],myShortestPaths[[4]],
  #         myShortestPaths[[5]],myShortestPaths[[6]],myShortestPaths[[7]],myShortestPaths[[8]],
  #         col = "grey",at = c(1:8),names = colnames(ahd)[2:9],las=2)
  
  
  # make a matrix with means and medians
  myHelpVec <- seq(1,8)
  myMeanMatrix <- matrix(data = NA, nrow = 8, ncol = 8)
  myMedianMatrix <- matrix(data = NA, nrow = 8, ncol = 8)
  for(i in 1:8){
    myMeanMatrix[i,i] <- mean(myShortestPaths[[i]])
    myMedianMatrix[i,i] <- median(myShortestPaths[[i]])
    
    myTmpVec <- myHelpVec[-i]
    for(j in 1:7){
      if(length(myShortestPaths2[[i]]) == 1){
        myMeanMatrix[i,myTmpVec[j]] <- NA
        myMedianMatrix[i,myTmpVec[j]] <- NA
      }else{
        myMeanMatrix[i,myTmpVec[j]] <- mean(myShortestPaths2[[i]][[j]])
        myMedianMatrix[i,myTmpVec[j]] <- median(myShortestPaths2[[i]][[j]])  
      }
      
      
    }
  }
  
  # myMeanMatrix <- data.frame(myMeanMatrix)
  # colnames(myMeanMatrix) <- colnames(ahd[2:9])
  # rownames(myMeanMatrix) <- colnames(ahd[2:9])
  # # 
  # # myMedianMatrix <- data.frame(myMedianMatrix)
  # colnames(myMedianMatrix) <- colnames(ahd[2:9])
  # rownames(myMedianMatrix) <- colnames(ahd[2:9])
  
  hormShortNames <- c("ABA","AUX","BR","CK","ET","GA","JA","SA")
  colnames(myMeanMatrix) <- hormShortNames
  rownames(myMeanMatrix) <- hormShortNames
  myMeanMatrix[lower.tri(myMeanMatrix)] <- NA
  
  ind <- which( ! is.na(myMeanMatrix) , arr.ind = TRUE ) 
  outMean <<- cbind(ind , myMeanMatrix[ ! is.na( myMeanMatrix ) ] )
  colnames(outMean) <- c("y","x","z")
  # 
  # myMedianMatrix <- data.frame(myMedianMatrix)
  colnames(myMedianMatrix) <- hormShortNames
  rownames(myMedianMatrix) <- hormShortNames
  myMedianMatrix[lower.tri(myMedianMatrix)] <- NA
  
  ind <- which( ! is.na(myMedianMatrix) , arr.ind = TRUE ) 
  outMedian <<- cbind(ind , myMedianMatrix[ ! is.na( myMedianMatrix ) ] )
  colnames(outMedian) <- c("y","x","z")
  
  # attach(outMean)
  # attach(outMedian)
  
  
  marker = list(color = colorRampPalette(brewer.pal(11,"Spectral"))(100))
  
  pdf(paste(result.dir,file.name," - Hormone Distance Matrix.pdf",sep=""),width = 4, height = 4)
  # print(levelplot(myMeanMatrix,scales=list(tck=0, x=list(rot=90)), xlab=NULL, ylab=NULL, 
  #           col.regions=colorpanel(16, "yellow", "blue")))
  # print(levelplot(myMedianMatrix,scales=list(tck=0, x=list(rot=90)), xlab=NULL, ylab=NULL, 
  #                 col.regions=colorpanel(16, "yellow", "blue")))
  
  colorRange <- 13
  
  # print(levelplot(myMeanMatrix,scales=list(tck=0, x=list(rot=90)), xlab=NULL, ylab=NULL,at = c(-1,1:colorRange),
  #                 col.regions=marker$color,colorkey = list(at = c(0:colorRange)), main=list('Mean Distance - all paths', side=1, line = 0.5))
  #       + layer(panel.text(x, y,  round(z,1), data=outMean)))
  
  l1 <- levelplot(myMeanMatrix,scales=list(tck=0, x=list(rot=90)), xlab=NULL, ylab=NULL,at = c(-1,1:colorRange),
            col.regions=marker$color,colorkey = list(at = c(0:colorRange)), main=list('Mean Distance - all paths', side=1, line = 0.5))
  plot(l1)
  l1 <- l1 + layer(panel.text(x, y,  round(z,1), data=outMean))
  plot(l1)
  
  # print(levelplot(myMedianMatrix,scales=list(tck=0, x=list(rot=90)), xlab=NULL, ylab=NULL,at = c(-1,1:colorRange),
  #                 col.regions=marker$color,colorkey = list(at = c(0:colorRange)), main=list('Median Distance - all paths', side=1, line = 0.5))
  #       + layer(panel.text(x, y,  z, data=outMedian)))
  
  l2 <- levelplot(myMedianMatrix,scales=list(tck=0, x=list(rot=90)), xlab=NULL, ylab=NULL,at = c(-1,1:colorRange),
                  col.regions=marker$color,colorkey = list(at = c(0:colorRange)), main=list('Median Distance - all paths', side=1, line = 0.5))
  plot(l2)
  l2 <- l2 + layer(panel.text(x, y,  z, data=outMedian))
  plot(l2)
  dev.off()
  
  write.xlsx(myMeanMatrix, paste(result.dir, file.name," - Shortest Paths - all paths.xlsx", sep = ""), sheetName = "Mean Values")
  write.xlsx(myMedianMatrix, paste(result.dir, file.name," - Shortest Paths - all paths.xlsx", sep = ""), sheetName = "Median Values", append=T)
  
}
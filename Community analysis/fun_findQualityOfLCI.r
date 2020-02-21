
library(sqldf)

findQuality <- function(interactionData, phoFilter, binaryMethods, pho){
  # interactionData <- intact.tair  
  # phoFilter <- FALSE
  print("----------------------")
  
  print(paste("Number of interactions for A. thaliana: ",length(interactionData[,1]),sep=""))
  
  intact.tair.unique <- interactionData[,1:2]
  intact.tair.unique <- data.frame(t(apply(intact.tair.unique[,],1,sort))) # sort systematic names by column
  # colnames(intact.pho.interactions) <- myColnames
  names(intact.tair.unique$X1) <- NULL # remove names of the first column
  names(intact.tair.unique$X2) <- NULL # remove names of the second column
  # sort systematic names by row
  if(length(intact.tair.unique[,1])>1){
    intact.tair.unique <- intact.tair.unique[order(intact.tair.unique$X1,
                                                   intact.tair.unique$X2),]
  }else{
    colnames(intact.tair.unique) <- c("X1","X2")
  }
  # select distinct interactions
  intact.tair.unique <- sqldf("select distinct `X1`,`X2` from `intact.tair.unique`")
  
  print(paste("Number of unique interactions: ",length(intact.tair.unique[,1]),sep=""))
  
  
  # select interactions, where both interaction partners are in pho
  if(phoFilter){
    intact.pho.interactions.index <- which(intact.tair.unique[,1] %in% pho$Locus_ID & intact.tair.unique[,2] %in% pho$Locus_ID)
    intact.pho.interactions <- intact.tair.unique[intact.pho.interactions.index,1:2]
    
    # if(length(intact.pho.interactions[,1]) == 0){
    #   next
    # }
    # myColnames <- colnames(intact.pho.interactions)
    
    intact.pho.interactions <- data.frame(t(apply(intact.pho.interactions,1,sort))) # sort systematic names by column
    # colnames(intact.pho.interactions) <- myColnames
    names(intact.pho.interactions$X1) <- NULL # remove names of the first column
    names(intact.pho.interactions$X2) <- NULL # remove names of the second column
    # sort systematic names by row
    if(length(intact.pho.interactions[,1])>1){
      intact.pho.interactions <- intact.pho.interactions[order(intact.pho.interactions$X1,
                                                               intact.pho.interactions$X2),]
    }else{
      colnames(intact.pho.interactions) <- c("X1","X2")
    }
    # select distinct interactions
    intactPhoInteractionsUnique <- sqldf("select distinct `X1`,`X2` from `intact.pho.interactions`")
    intactPhoInteractionsUnique <- intactPhoInteractionsUnique[order(intactPhoInteractionsUnique$X1,
                                                                     intactPhoInteractionsUnique$X2),]
    
    print(paste("Intact - number of unique interactions in search space: ",length(intactPhoInteractionsUnique[,1]),sep=""))
  }else{
    intactPhoInteractionsUnique <- intact.tair.unique
  }
  
  # if(length(intactPhoInteractionsUnique[,1]) == 0){
  #   next
  # }
  ### identify high quality, binary and other interactions
  interactionQuality <- c()
  pubmedCombined <- c()
  methodsCombined <- c()
  for(i in 1:length(intactPhoInteractionsUnique$X1)){
    intIndexes <- which((interactionData[,1] == intactPhoInteractionsUnique$X1[i] & 
                           interactionData[,2] == intactPhoInteractionsUnique$X2[i])|
                          (interactionData[,1] == intactPhoInteractionsUnique$X2[i] & 
                             interactionData[,2] == intactPhoInteractionsUnique$X1[i]))
    if(length(intIndexes)==0){
      stop("interaction not found")
    }
    # print(paste("loop",i))
    # print(intIndexes)
    
    # extract method and publication information for the current interactions
    intInfo <- data.frame(interactionData[intIndexes,3],interactionData[intIndexes,4])
    intInfoUnique <- sqldf("select distinct * from intInfo")
    colnames(intInfo) <- c("Method","PubmedID")
    # check if the interaction is hq or binary or nothing
    binaryIndex <- which(intInfo$Method %in% binaryMethods)
    noBinary <- length(unique(intInfo$Method)) # changed 20160201 from length(binaryIndex) to length(unique((intInfo$Method)))
    noPubmed <- length(unique(intInfo$PubmedID))
    noPubmedBinary <- length(unique(intInfo$PubmedID[binaryIndex]))
    noMethods <- length(unique(intInfo$Method))
    noMethodsBinary <- length(unique(intInfo$Method[binaryIndex]))
    
    if(noBinary > 0 & noMethods > 1){
      interactionQuality[i] <- "HQ"
    }else if(noBinary == 1 & noPubmedBinary > 1){
      interactionQuality[i] <- "HQ"
    }else if(noBinary == 1 & noPubmedBinary == 1){
      interactionQuality[i] <- "binary"
    }else{
      interactionQuality[i] <- "non-binary"
    }
    
    pubmedCombined <- c(pubmedCombined, paste(intInfo$PubmedID,collapse = ","))
    methodsCombined <- c(methodsCombined, paste(intInfo$Method, collapse = ","))
    
  }
  
  intactPhoInteractionsUnique <- cbind(intactPhoInteractionsUnique,interactionQuality,pubmedCombined,methodsCombined)
  
  ## find number of interactions in each category
  # intactHQ <- length(which(intactPhoInteractionsUnique$interactionQuality == "HQ"))
  # intactBinary <- length(which(intactPhoInteractionsUnique$interactionQuality == "binary"))
  # intactNonBinary <- length(which(intactPhoInteractionsUnique$interactionQuality == "non-binary"))
  ## print number of interactions in each category
  # print(paste("intact search space HQ:", intactHQ))
  # print(paste("intact search space Binary:", intactBinary))
  # print(paste("intact search space Non-Binary:", intactNonBinary))
  
  return(intactPhoInteractionsUnique)
}



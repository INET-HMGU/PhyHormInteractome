
calcEnrichement <- function(myInteractome,dbg){
  number.members <- c() # collect the number of the members of the communities
  number.members.annotated <- c()
  # hormone names
  ahdHormones <- colnames(ahd)[2:9]
  p.val.vec <- c() # collect the p-values
  p.val.matrix <- c() # collect the p-values
  p.val.matrix2 <- c() # collect the p-values
  p.val.matrix3 <- c() # collect the p-values
  myIndex <- which(ahd$AGI %in% unique(c(as.character(myInteractome$IntA),as.character(myInteractome$IntB))))
  allLoci <- unique(c(as.character(myInteractome$IntA),as.character(myInteractome$IntB)))
  notAnnotatedLoci <- allLoci[which(!allLoci %in% ahd$AGI)]
  ahdInteractome <- ahd[myIndex,]
  # loop over the communities
  for(i in 1:length(unique(eb.com.list$membership))){
    if(dbg){
      print(paste("Community",i,"##################################################"))
    }
    # determine AGI of the members of the current community
    my.members <- eb.com.list$names[which(eb.com.list$membership == i)]
    number.members <- c(number.members,length(my.members))
    my.members.ahd.index <- which(ahdInteractome$AGI %in% my.members)
    # loop over the hormones and test their enrichment
    p.val.vec <- c() # collect the p-values
    p.val.vec2 <- c() # collect the p-values
    p.val.vec3 <- c() # collect the p-values
    for(j in 1:8){
      # annotated with current hormone in the community
      myIndex <- which(ahdInteractome[my.members.ahd.index,j+1] == 1)
      hormXcomIndex <- my.members.ahd.index[myIndex]
      HormXCom <- length(myIndex)
      
      ### not annotated with the current hormone in the community
      newMatrix <- c()
      if(HormXCom > 0){
        newMatrix <- ahdInteractome[my.members.ahd.index[-myIndex],]
      }else{
        newMatrix <- ahdInteractome[my.members.ahd.index,]
      }
      # test, if the matrix containing hormone annotation information is empty
      if(dim(newMatrix)[1] == 0){
        NotHormXCom <- 0
      }else{
        NotHormXCom <- length(unique(which(newMatrix == 1,arr.ind = T)[,1]))  
      }
      # in the community, but with different or no annotation
      notAnnotatedInCom <- length(which(notAnnotatedLoci %in% my.members))
      NotHormNotAnnoXCom <- NotHormXCom + notAnnotatedInCom
      
      # annotated with the current hormone in the complete interactome
      myIndex <- which(ahdInteractome[,j+1] == 1)
      HormXPi <- length(myIndex)
      # annotated with the currect hormone not in the current community
      myIndexNc <- which(ahdInteractome[,j+1] == 1)
      if(length(hormXcomIndex) > 0){
        removeIndex <- which(myIndexNc %in% hormXcomIndex)
        if(length(removeIndex) > 0){
          myIndexNc <- myIndexNc[-removeIndex]  
        }
      }
      HormXNotCom <- length(myIndexNc)
      
      
      # not annotated with the current hormone in the complete interactome
      if(length(myIndex) == 0){
        newMatrix <- ahdInteractome
      }else{
        newMatrix <- ahdInteractome[-myIndex,]  
      }
      
      if(dim(newMatrix)[1] == 0){
        NotHormXPi <- 0
      }else{
        NotHormXPi <- length(unique(which(newMatrix == 1,arr.ind = T)[,1]))
      }
      # not annotated with the current hormone and not in the current community
      NotHormXNotCom <- length(which(ahdInteractome[-my.members.ahd.index,j+1] == 0))
      
      # not annotated with the current hormone and not annotated at all 
      # and not in the current community
      removeIndex <- which(notAnnotatedLoci %in% my.members)
      notAnnotatedLoci2 <- c()
      if(length(removeIndex) > 0){
        notAnnotatedLoci2 <- notAnnotatedLoci[-removeIndex]
      }
      NotHormNotAnnoXNotCom <- length(notAnnotatedLoci2) + NotHormXNotCom
      
      ### build the matrix
      # using values of the complete interactome -> wrong
      x <- matrix(c(HormXCom,NotHormXCom,HormXPi,NotHormXPi),2,2) # erstellt eine 2x2 Felder Tafel
      dimnames(x) <-  list(c("HormX", "NotHormX"), c("Community", "PI")) # Namen zuweisen
      # using values of the communities
      y <- matrix(c(HormXCom,NotHormXCom,HormXNotCom,NotHormXNotCom),2,2)
      dimnames(y) <-  list(c("HormX", "NotHormX"), c("Community", "NotCommunity")) # Namen zuweisen
      # using values of the communities plus non-annotated proteins
      z <- matrix(c(HormXCom,NotHormNotAnnoXCom,HormXNotCom,NotHormNotAnnoXNotCom),2,2)
      dimnames(z) <-  list(c("HormX", "NotHormX"), c("Community", "NotCommunity")) # Namen zuweisen
      
      
      if(dbg){
        print(paste("Community", i, ", Hormone ", ahdHormones[j]))
        # print(x)
        # print("V2 ----------")
        # print(y)
        print("Version 3 - using values of the communities plus non-annotated proteins")
        print(z)
      }
      
      my.p.val <- NA
      my.p.val2 <- NA
      my.p.val3 <- NA
      if(HormXCom == 0){
        my.p.val <- 1
        my.p.val2 <- 1
        my.p.val3 <- 1
      }else{
        my.p.val <- fisher.test(x)$p.value
        my.p.val2 <- fisher.test(y)$p.value
        my.p.val3 <- fisher.test(z)$p.value
      }
      p.val.vec <- c(p.val.vec, my.p.val)
      p.val.vec2 <- c(p.val.vec2, my.p.val2)
      p.val.vec3 <- c(p.val.vec3, my.p.val3)
      
      if(dbg){
        # print(paste("P-Value:",my.p.val))
        # print(paste("P-Value2:",my.p.val2))
        print(paste("P-Value3:",my.p.val3))
      }
      
    }
    p.val.matrix <- rbind(p.val.matrix,p.val.vec)
    p.val.matrix2 <- rbind(p.val.matrix2,p.val.vec2)
    p.val.matrix3 <- rbind(p.val.matrix3,p.val.vec3)
  }
  
  my.min.values <- apply(p.val.matrix,1,min) # minimum p-value per community
  p.val.matrix <- cbind(p.val.matrix,number.members,my.min.values) # add vector containing minimum p-value to p-value matrix
  
  my.min.values2 <- apply(p.val.matrix2,1,min) # minimum p-value per community
  p.val.matrix2 <- cbind(p.val.matrix2,number.members,my.min.values2) # add vector containing minimum p-value to p-value matrix
  
  my.min.values3 <- apply(p.val.matrix3,1,min) # minimum p-value per community
  p.val.matrix3 <- cbind(p.val.matrix3,number.members,my.min.values3) # add vector containing minimum p-value to p-value matrix
  
  
  # correct p-values for multiple testing
  # holm correction
  correctedPvaluesHolm <- sapply(p.val.matrix[,10],p.adjust,method = "holm", n=8)
  correctedPvaluesBH <- sapply(p.val.matrix[,10],p.adjust,method = "BH", n=8)
  correctedPvaluesFDR <- sapply(p.val.matrix[,10],p.adjust,method = "fdr", n=8)
  correctedPvaluesBY <- sapply(p.val.matrix[,10],p.adjust,method = "BY", n=8)
  correctedPvaluesBon <- sapply(p.val.matrix[,10],p.adjust,method = "bonferroni", n=8)
  p.val.matrix <- cbind(p.val.matrix,correctedPvaluesHolm,correctedPvaluesBH,correctedPvaluesFDR,correctedPvaluesBY,correctedPvaluesBon)
  
  rownames(p.val.matrix) <- seq(1,length(p.val.matrix[,1])) # set rownames
  colnames(p.val.matrix) <- c(colnames(ahd)[2:9],"Community.Size","min.value","Corrected.P-Value.Holm",
                              "Corrected.P-Value.BH","Corrected.P-Value.fdr","Corrected.P-Value.BY",
                              "Corrected.P-Value.Bonferroni") # set colnames
  # version restricted on communities
  correctedPvaluesHolm <- sapply(p.val.matrix2[,10],p.adjust,method = "holm", n=8)
  correctedPvaluesBH <- sapply(p.val.matrix2[,10],p.adjust,method = "BH", n=8)
  correctedPvaluesFDR <- sapply(p.val.matrix2[,10],p.adjust,method = "fdr", n=8)
  correctedPvaluesBY <- sapply(p.val.matrix2[,10],p.adjust,method = "BY", n=8)
  correctedPvaluesBon <- sapply(p.val.matrix2[,10],p.adjust,method = "bonferroni", n=8)
  p.val.matrix2 <- cbind(p.val.matrix2,correctedPvaluesHolm,correctedPvaluesBH,correctedPvaluesFDR,correctedPvaluesBY,correctedPvaluesBon)
  
  rownames(p.val.matrix2) <- seq(1,length(p.val.matrix2[,1])) # set rownames
  colnames(p.val.matrix2) <- c(colnames(ahd)[2:9],"Community.Size","min.value","Corrected.P-Value.Holm",
                              "Corrected.P-Value.BH","Corrected.P-Value.fdr","Corrected.P-Value.BY",
                              "Corrected.P-Value.Bonferroni") # set colnames
  
  # version restricted on communities including non-annotated proteins
  correctedPvaluesHolm <- sapply(p.val.matrix3[,10],p.adjust,method = "holm", n=8)
  correctedPvaluesBH <- sapply(p.val.matrix3[,10],p.adjust,method = "BH", n=8)
  correctedPvaluesFDR <- sapply(p.val.matrix3[,10],p.adjust,method = "fdr", n=8)
  correctedPvaluesBY <- sapply(p.val.matrix3[,10],p.adjust,method = "BY", n=8)
  correctedPvaluesBon <- sapply(p.val.matrix3[,10],p.adjust,method = "bonferroni", n=8)
  p.val.matrix3 <- cbind(p.val.matrix3,correctedPvaluesHolm,correctedPvaluesBH,correctedPvaluesFDR,correctedPvaluesBY,correctedPvaluesBon)
  
  rownames(p.val.matrix3) <- seq(1,length(p.val.matrix3[,1])) # set rownames
  colnames(p.val.matrix3) <- c(colnames(ahd)[2:9],"Community.Size","min.value","Corrected.P-Value.Holm",
                               "Corrected.P-Value.BH","Corrected.P-Value.fdr","Corrected.P-Value.BY",
                               "Corrected.P-Value.Bonferroni") # set colnames
  
  return(list(p.val.matrix,p.val.matrix2,p.val.matrix3))
}



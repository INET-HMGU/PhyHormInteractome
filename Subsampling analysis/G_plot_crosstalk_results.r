# Stefan Altmann
# this script was executed using R3.6.1 and Windows 10
# count and plot number of solid (Type I) and possible (Type II) crosstalks
# afterwards to monte carlo subsampling

rm(list = ls())

setwd("~/../Desktop/Github/Subsampling")

library(ggplot2)
library(ggthemes)

library(openxlsx)

library(igraph)


# select hormone annotation set for analysis
myHormoneAnnotationSetList <-
  c("AHDGen_TAIR", "AHD_Genetic","TAIR_GO10","AHDGen_TAIR10")

# result directory
result.dir <- "Results"
if(!dir.exists(result.dir)){
  dir.create(result.dir)
}

# store data
annoData <- list()
intData <- c()

# read interaction data
intData <- read.xlsx("Extended Table 2 - InteractionData.xlsx", sheet = "Interaction datasets", startRow = 2)

# read intact interactions
intDataIntAct <- read.csv2("IntAct Interactions.csv")

# read biogrid interactions
intDataBioGrid <- read.csv2("BioGRID Interactions.csv")

# uncomment set of networks to be analyzed for crosstalks
# selectedNetwork <- c("PvP")
selectedNetwork <- c("PhIMAIN","Intact_SS_BMBS","Biogrid_SS_BMBS")

networkResults <- list()


################
detCrosstalks <- function(myInteractions){
  # empty annotation matrix for "solid" crosstalks
  annoMatrixS <- matrix(data = 0, nrow = hc, ncol = hc)
  # empty annotation matrix for "possible" crosstalks
  annoMatrixP <- matrix(data = 0, nrow = hc, ncol = hc)
  
  # iterate over interactions
  for(i in 1:length(myInteractions[,1])){
    indexA <- which(annoData[[1]]$AGI == as.character(myInteractions[i,1]))
    indexB <- which(annoData[[1]]$AGI == as.character(myInteractions[i,2]))
    
    if(length(indexA) == 1 & length(indexB) == 1){
      annoA <- annoData[[1]][indexA,2:(1+hc)]
      annoB <- annoData[[1]][indexB,2:(1+hc)]
      
      allA <- which(annoA == 1)
      allB <- which(annoB == 1)
      allExp <- expand.grid(allA, allB)
      
      if (any(annoA) & any(annoB)) {
        if (any(xor(annoA, annoB))) {
          if (any(annoA & annoB)) {
            # crosstalkResult <- "possible"
            for(j in 1:length(allExp[,1])){
              annoMatrixP[allExp[j, 1], allExp[j, 2]] <-
                annoMatrixP[allExp[j, 1], allExp[j, 2]] + 1
            }
          } else{
            # crosstalkResult <- "solid"
            for(j in 1:length(allExp[,1])){
              ### for debug only
              # if((allExp[j,1] == 1 & allExp[j, 2] == 2)|
              #    (allExp[j,1] == 2 & allExp[j, 2] == 1)){
              #   print(paste(intData[[1]]$IntA[index1[i]], ":", intData[[1]]$IntB[index1[i]]))
              # }
              annoMatrixS[allExp[j, 1], allExp[j, 2]] <-
                annoMatrixS[allExp[j, 1], allExp[j, 2]] + 1
            }
          }
        } else{
          # strictly same hormone annotation
          if(length(allExp[,1]) == 1){
            annoMatrixS[allExp[1, 1], allExp[1, 2]] <-
              annoMatrixS[allExp[1, 1], allExp[1, 2]] + 1
          }
          # for(j in 1:length(allExp[,1])){
          #   annoMatrixP[allExp[j, 1], allExp[j, 2]] <-
          #     annoMatrixP[allExp[j, 1], allExp[j, 2]] + 1
          # }
        }
      } else{
        # crosstalkResult <- "no"
      }
    }else{
      # do nothing
    }
    
  }
  
  colnames(annoMatrixS) <- colnames(annoData[[1]])[2:(1+hc)]
  rownames(annoMatrixS) <- colnames(annoData[[1]])[2:(1+hc)]
  
  colnames(annoMatrixP) <- colnames(annoData[[1]])[2:(1+hc)]
  rownames(annoMatrixP) <- colnames(annoData[[1]])[2:(1+hc)]
  
  sumMatrixS <- matrix(data = 0, nrow = hc, ncol = hc)
  
  for(i in 1:hc){
    sumMatrixS[i,i] <- annoMatrixS[i,i]
  }
  
  for(i in 1:(hc-1)){
    for(j in (i+1):hc){
      sumMatrixS[i,j] <- annoMatrixS[i,j] + annoMatrixS[j,i]
    }
  }
  
  
  
  sumMatrixP <- matrix(data = 0, nrow = hc, ncol = hc)
  
  for(i in 1:hc){
    sumMatrixP[i,i] <- annoMatrixP[i,i]
  }
  
  for(i in 1:(hc-1)){
    for(j in (i+1):hc){
      sumMatrixP[i,j] <- annoMatrixP[i,j] + annoMatrixP[j,i]
    }
  }
  
  # set hormone names as column and row names
  hormNames <- c("ABA", "AUX", "BR", "CK", "ET", "GA", "JA", "SA", "KAR", "SL")
  
  colnames(sumMatrixS) <- hormNames[1:hc]
  rownames(sumMatrixS) <- hormNames[1:hc]
  colnames(sumMatrixP) <- hormNames[1:hc]
  rownames(sumMatrixP) <- hormNames[1:hc]
  
  return(list(annoMatrixS, annoMatrixP, sumMatrixS, sumMatrixP))
  
}

##################




for(s in 1:length(selectedNetwork)){
  
  print(selectedNetwork[s])
  
  if(selectedNetwork[s] == "Intact_SS_BMBS"){
    myNetworkIndex <- which(colnames(intDataIntAct) == "Intact_SS")
  }else if(selectedNetwork[s] == "Biogrid_SS_BMBS"){
    myNetworkIndex <- which(colnames(intDataBioGrid) == "Biogrid_SS")  
  }else{
    myNetworkIndex <- which(colnames(intData) == selectedNetwork[s])  
  }
  
  annoSet <- c("HormAnno-AHDGen,tairGO",
               "HormAnno-AHDGen",
               "HormAnno-tairGO10",
               "HormAnno-AHDGen,tairGO10")
  
  hormoneResults <- list()
  
  for(k in 1:4){
    print(annoSet[k])
    annoData[[1]] <- read.csv2(file = paste(annoSet[k],".csv",sep = ""))
    
    # depending on the hormone annotations, 8 or 10 columns must be read
    hc <- 0
    if(k < 3){
      hc <- 8
    }else{
      hc <- 10
    }
    
    # select interactions for selected network
    if(selectedNetwork[s] == "Intact_SS_BMBS"){
      index1 <- which(intDataIntAct[,myNetworkIndex] == 1 & 
                        (intDataIntAct$IntAct_Quality == "binary" | intDataIntAct$IntAct_Quality == "HQ"))
    }else if(selectedNetwork[s] == "Biogrid_SS_BMBS"){
      index1 <- which(intDataBioGrid[,myNetworkIndex] == 1 & 
                        (intDataBioGrid$BioGRID_Quality == "binary" | intDataBioGrid$BioGRID_Quality == "HQ"))
    }else{
      index1 <- which(intData[,myNetworkIndex] == 1)
    }
    
    # determine number of crosstalks
  
    if(selectedNetwork[s] == "Intact_SS_BMBS"){
      mylist <- detCrosstalks(intDataIntAct[index1,1:2])
    }else if(selectedNetwork[s] == "Biogrid_SS_BMBS"){
      mylist <- detCrosstalks(intDataBioGrid[index1,1:2])
    }else{
      mylist <- detCrosstalks(intData[index1,1:2])
    }
    
    
    
    annoMatrixS <- mylist[[1]] 
    annoMatrixP <- mylist[[2]] 
    sumMatrixS <- mylist[[3]] 
    sumMatrixP <- mylist[[4]] 
    
    # store values for later use
    sumMatSReal <- sumMatrixS
    sumMatPReal <- sumMatrixP

    ##### in grid
    
    hormNames <- c("ABA", "AUX", "BR", "CK", "ET", "GA", "JA", "SA", "KAR", "SL")
    
    # convert matrix to data frame for plotting
    zzz <- c()
    for(i in 1:hc){
      for(j in 1:hc){
        xxx <- data.frame(sumMatrixS[i,j], "solid", hormNames[i],hormNames[j])
        colnames(xxx) <- c("value","type","horm1", "horm2")
        zzz <- rbind(zzz, xxx)
        
        
        xxx <- data.frame(sumMatrixP[i,j], "possible", hormNames[i],hormNames[j])
        colnames(xxx) <- c("value","type","horm1", "horm2")
        zzz <- rbind(zzz, xxx)
        
      }
    }
    # for later use
    zzzReal <- zzz
    
    # write results to file
    wb <- createWorkbook()
    addWorksheet(wb, sheetName = "solid")
    writeData(wb, sheet = "solid", sumMatrixS,colNames = T, rowNames = T)
    addWorksheet(wb, sheetName = "possible")
    writeData(wb, sheet = "possible", sumMatrixP, colNames = T, rowNames = T)
    
    saveWorkbook(wb, file = paste(result.dir,"/Crosstalks-",selectedNetwork[s],"-",annoSet[k],".xlsx", sep = ""), overwrite = T)
    
    # plot results
    pdf(paste(result.dir,"/Crosstalks-",selectedNetwork[s],"-",annoSet[k],".pdf", sep = ""), width = 6, height = 8)
    
    bp4 <- ggplot(zzz, aes(x="", y=value, fill = type)) +
      facet_grid(horm1~horm2) +
      geom_bar(width = .8, stat = "identity", position = position_dodge(width = .7)) +
      scale_fill_manual(values = c("#CCA303", "#6495ED")) +
      scale_x_discrete(breaks = NULL) +
      theme_gray()
    print(bp4)
    
    dev.off()
    
    ################ randomization test
    doRandomization <- T
    if(doRandomization){
      nr <- 1000 # number of random networks to analyze
      # collecting results of annoMatrixS of random network
      rams <- array(rep(0, nr*hc*hc), dim=c(hc, hc, nr))
      # collecting results of annoMatrixP of random network
      ramp <- array(rep(0, nr*hc*hc), dim=c(hc, hc, nr))
      # collecting results of sumMatrixS of random network
      rsms <- array(rep(0, nr*hc*hc), dim=c(hc, hc, nr))
      # collecting results of sumMatrixP of random network
      rsmp <- array(rep(0, nr*hc*hc), dim=c(hc, hc, nr))
      
      for(l in 1:nr){
        print(l)
        mySample <- sample.int(length(index1), 100)
        # select 100 random interactions
        my.graph.sample <- graph_from_data_frame(intData[index1[mySample],1:2], directed = FALSE)
        
        # test for crosstalks
        mylist <- detCrosstalks(as_data_frame(my.graph.sample))
        annoMatrixS <- mylist[[1]]
        annoMatrixP <- mylist[[2]]
        sumMatrixS <- mylist[[3]]
        sumMatrixP <- mylist[[4]]
        rams[1:hc,1:hc,l] <- annoMatrixS
        ramp[1:hc,1:hc,l] <- annoMatrixP
        rsms[1:hc,1:hc,l] <- sumMatrixS
        rsmp[1:hc,1:hc,l] <- sumMatrixP
      }
      
      ### calculate distribution for each hormone combination and for all combinations
      crosstalk.table.s <- c()
      crosstalk.table.p <- c()
      crosstalk.table.colnames <- c()
      for(i in 1:(hc-1)){
        for(j in (i+1):hc){
          crosstalk.table.s <- cbind(crosstalk.table.s, rsms[i,j,])
          crosstalk.table.p <- cbind(crosstalk.table.p, rsmp[i,j,])
          crosstalk.table.colnames <- c(crosstalk.table.colnames, paste(hormNames[i], hormNames[j], sep = "-"))
        }
      }
      
      colnames(crosstalk.table.s) <- crosstalk.table.colnames
      colnames(crosstalk.table.p) <- crosstalk.table.colnames
      
      
      
      # write results to file
      wb <- createWorkbook()
      
      addWorksheet(wb, sheetName = "PCP Type 1")
      writeData(wb, sheet = "PCP Type 1", crosstalk.table.s,colNames = T, rowNames = T)
      
      addWorksheet(wb, sheetName = "PCP Type 2")
      writeData(wb, sheet = "PCP Type 2", crosstalk.table.p, colNames = T, rowNames = T)
      
      saveWorkbook(wb, file = paste(result.dir,"/Crosstalks-",selectedNetwork[s],"-",annoSet[k],".xlsx", sep = ""), overwrite = T)
      
      
      # # prepare randomization results for plotting
      # # convert matrix to data frame for plotting
      # zzz <- c()
      # for(i in 1:hc){
      #   for(j in 1:hc){
      #     xxx <- data.frame(mean(rsms[i,j,]),sd(rsms[i,j,]), "solid", hormNames[i],hormNames[j])
      #     colnames(xxx) <- c("value","sd","type","horm1", "horm2")
      #     zzz <- rbind(zzz, xxx)
      #     
      #     
      #     xxx <- data.frame(mean(rsmp[i,j,]), sd(rsmp[i,j,]), "possible", hormNames[i],hormNames[j])
      #     colnames(xxx) <- c("value","sd","type","horm1", "horm2")
      #     zzz <- rbind(zzz, xxx)
      #     
      #   }
      # }
      # 
      # pdf(paste(result.dir,"/Crosstalks-",selectedNetwork[s],"-",annoSet[k],"RandomNetwork.pdf", sep = ""), width = 6, height = 8)
      # 
      # # print results
      # bp4 <- ggplot(zzz, aes(x="", y=value, fill = type)) +
      #   facet_grid(horm1~horm2) +
      #   geom_bar(width = .8, stat = "identity", position = position_dodge(width = .7)) +
      #   geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0, position=position_dodge(.6)) +
      #   scale_fill_manual(values = c("#CCA303", "#6495ED")) +
      #   scale_x_discrete(breaks = NULL) +
      #   theme_gray()
      # print(bp4)
      # 
      # zzz <- cbind(zzz, rep("Random", length(zzz[,1])))
      # colnames(zzz)[6] <- "source"
      # zzzReal <- cbind(zzzReal[,1],rep(0, length(zzzReal[,1])), zzzReal[,2:4], rep("Observed", length(zzzReal[,1])))
      # colnames(zzzReal)[6] <- "source"
      # colnames(zzzReal)[2] <- "sd"
      # colnames(zzzReal)[1] <- "value"
      # zzzComp <- rbind(zzz, zzzReal)
      # bp5 <- ggplot(zzzComp, aes(x=source, y=value, fill = type)) +
      #   facet_grid(horm1~horm2) +
      #   geom_bar(width = .8, stat = "identity", position = position_dodge()) +
      #   geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0, position=position_dodge(.7)) +
      #   scale_fill_manual(values = c("#CCA303", "#6495ED")) +
      #   theme_gray()
      # print(bp5)
      # 
      # dev.off()
      
      
      
    }
    
    
  
  }
  
  
}

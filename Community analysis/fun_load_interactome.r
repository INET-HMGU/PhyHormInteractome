

suppressPackageStartupMessages(library(gdata))
suppressPackageStartupMessages(library(sqldf))

load_interactome <- function(interactomesToUse, restrictAllToSs, dbg, writeToFile,SSExtended){
  
  
  myInteractome <- c()
  # read PvP
  if("PvP" %in% interactomesToUse){
    # path to pvp interactome
    myInteractomePath <- "Extended Table 2 - InteractionData.xlsx"
    # read the interactome
    myInteractomeTemp <- openxlsx::read.xlsx(myInteractomePath, sheet = "Interaction datasets", startRow = 2)
    # extract interactions of PhIMAIN
    myInteractome <- rbind(myInteractome,myInteractomeTemp[which(myInteractomeTemp$PhIMAIN == 1),1:2])
  }
  # read PvA
  if("PvA" %in% interactomesToUse){
    # path to pva interactome
    myInteractomePath <- "Extended Table 2 - InteractionData.xlsx"
    # read the interactome
    myInteractomeTemp <- openxlsx::read.xlsx(myInteractomePath, sheet = "Interaction datasets", startRow = 2)
    # extract interactions of PvA
    myInteractome <- rbind(myInteractome,myInteractomeTemp[which(myInteractomeTemp$PvA == 1),1:2])
  }
  # read PvTF
  if("PvTF" %in% interactomesToUse){
    # path to pva interactome
    myInteractomePath <- "Extended Table 2 - InteractionData.xlsx"
    # read the interactome
    myInteractomeTemp <- openxlsx::read.xlsx(myInteractomePath, sheet = "Interaction datasets", startRow = 2)
    # extract interactions of PvTF
    myInteractome <- rbind(myInteractome,myInteractomeTemp[which(myInteractomeTemp$PvTF == 1),1:2])
  }
  # read AI1Main
  if("AI1Main" %in% interactomesToUse){
    # ai1 <- xlsx::read.xlsx("Evidence for Network Evolution in an Arabidopsis Interactome Map/Arabidopsisinteractome_SOM_Tables/Arabidopsisinteractome_SOM_TableS4.xlsx",
    #                  sheetIndex = 1, stringsAsFactors = F)
    # save(ai1,file="ai1.RData")
    load("ai1.RData")
    mainIndex <- which(ai1$AI.1MAIN == 1)
    ai1Main <- ai1[mainIndex,c(1:2)]
    ai1Main <- data.frame(t(apply(ai1Main,1,sort)),row.names = NULL)
    colnames(ai1Main) <- c("IntA","IntB")
    names(ai1Main$X1) <- NULL
    names(ai1Main$X2) <- NULL
    myInteractomeTemp <- ai1Main
    myInteractome <- rbind(myInteractome,myInteractomeTemp)
  }
  
  # read AI1Repeat
  if("AI1Repeat" %in% interactomesToUse){
    # ai1 <- xlsx::read.xlsx("Evidence for Network Evolution in an Arabidopsis Interactome Map/Arabidopsisinteractome_SOM_Tables/Arabidopsisinteractome_SOM_TableS4.xlsx",
    #                  sheetIndex = 1, stringsAsFactors = F)
    load("ai1.RData")
    repeatIndex <- which(ai1$AI.1REPEAT== 1)
    ai1Repeat <- ai1[repeatIndex,c(1:2)]
    ai1Repeat <- data.frame(t(apply(ai1Repeat,1,sort)),row.names = NULL)
    colnames(ai1Repeat) <- c("IntA","IntB")
    names(ai1Repeat$X1) <- NULL
    names(ai1Repeat$X2) <- NULL
    myInteractomeTemp <- ai1Repeat
    myInteractome <- rbind(myInteractome,myInteractomeTemp)
  }
  
  if(any(grepl("intact",interactomesToUse,ignore.case = T)) | any(grepl("biogrid",interactomesToUse,ignore.case = T))){
    ### check literature interactions
    # load literature interactions
    ### IntAct ############################################################################
    # load intact data
    source("fun_load_intact.r")
    
    ### BioGrid ##############################################################
    # load biogrid
    source("fun_load_biogrid.r")
    
    # find quality function
    source("fun_findQualityOfLCI.r")
    
    
    ### load the experimental systems categorization: binary vs. not binary detection method
    
    # binary methods biogrid
    biogrid.exp.sys <- xlsx::read.xlsx("experimental systems/biogrid.experimental.systems.final.xlsx",sheetIndex=1)
    biogrid.binary <- biogrid.exp.sys$Experimental.System[which(biogrid.exp.sys$Binary.Final == "Yes")]
    
    # binary methods intact
    intact.exp.sys <- xlsx::read.xlsx("experimental systems/intact.experimental.systems.final.xlsx",sheetIndex=1)
    intact.binary <- intact.exp.sys$Experimental.System[which(intact.exp.sys$Binary.Final == "Yes")]
    
    # remove psi-mi prefix of methods
    intact.binary.tmp <- c()
    for(i in 1:length(intact.binary)){
      method.split <- unlist(strsplit(as.character(intact.binary[i]),"\\("))
      intact.binary.tmp <- c(intact.binary.tmp,substr(method.split[2],1,nchar(method.split[2])-1))
    }
    intact.binary <- intact.binary.tmp
  }
  
  if(any(grepl("intact",interactomesToUse,ignore.case = T)) | 
     any(grepl("biogrid",interactomesToUse,ignore.case = T)) | 
     restrictAllToSs){


    # load search space
    pho <- xlsx::read.xlsx("Extended Table 1 - Search space.xlsx", sheetIndex = 2)
    colnames(pho) <- c("Locus_ID", "Gene_family", "AHD2")

    
    pho$Locus_ID <- toupper(trim(pho$Locus_ID))
    print("###")
    print("LENGTH OF SS")
    print(length(pho$Locus_ID))
    print("###")
  }
  
  # read intact
  if("IntactSS_BM" %in% interactomesToUse){
    # load intact data
    # intact <- loadIntact()
    # save(intact,file="intact.RData")
    load("intact.RData")

    # keep only columns interactor A, interactor B, detection method, pubmed id, interaction type
    intact.tair <- intact[,c(5,6,7,9,12)]
    intact.tair$Publication.Identifier.s. <- gsub("pubmed:","",intact.tair$Publication.Identifier.s.)
    # determine quality to find binary multiple (BM) interactions
    # parameter: interaction set, only interactions in pho search space, set of binary methods
    # print("finding quality")
    myInteractionsIntact <- findQuality(intact.tair, FALSE, intact.binary, pho) 
    # print("finished finding quality")
    
    # myInteractionsIntact <- findQuality(intact.tair, FALSE, intact.binary)
    myInteractomeTemp <- myInteractionsIntact[which(myInteractionsIntact$interactionQuality == "HQ"),1:2]
    colnames(myInteractomeTemp) <- c("IntA","IntB")
    
    indexA <- which(myInteractomeTemp$IntA %in% pho$Locus_ID)
    indexB <- which(myInteractomeTemp$IntB %in% pho$Locus_ID)
    indexC <- intersect(indexA,indexB)
    
    myInteractomeTemp <- myInteractomeTemp[indexC,]
    myInteractomeTemp <- data.frame(t(apply(myInteractomeTemp,1,sort)))
    colnames(myInteractomeTemp) <- c("IntA","IntB")
    
    myInteractionsTemp <- sqldf("select distinct `IntA`,`IntB` from `myInteractomeTemp`")
    
    
    myInteractome <- rbind(myInteractome,myInteractomeTemp)
  }
  
  # read intact
  if("IntactSS_BS" %in% interactomesToUse){
    # load intact data
    # intact <- loadIntact()
    # save(intact,file="intact.RData")
    load("intact.RData")
    
    # keep only columns interactor A, interactor B, detection method, pubmed id, interaction type
    intact.tair <- intact[,c(5,6,7,9,12)]
    intact.tair$Publication.Identifier.s. <- gsub("pubmed:","",intact.tair$Publication.Identifier.s.)
    # determine quality to find binary multiple (BM) interactions
    # parameter: interaction set, only interactions in pho search space, set of binary methods
    myInteractionsIntact <- findQuality(intact.tair, FALSE, intact.binary,pho) 
    
    # myInteractionsIntact <- findQuality(intact.tair, FALSE, intact.binary)
    myInteractomeTemp <- myInteractionsIntact[which(myInteractionsIntact$interactionQuality == "binary"),1:2]
    colnames(myInteractomeTemp) <- c("IntA","IntB")
    
    indexA <- which(myInteractomeTemp$IntA %in% pho$Locus_ID)
    indexB <- which(myInteractomeTemp$IntB %in% pho$Locus_ID)
    indexC <- intersect(indexA,indexB)
    
    myInteractomeTemp <- myInteractomeTemp[indexC,]
    myInteractomeTemp <- data.frame(t(apply(myInteractomeTemp,1,sort)))
    colnames(myInteractomeTemp) <- c("IntA","IntB")
    
    myInteractionsTemp <- sqldf("select distinct `IntA`,`IntB` from `myInteractomeTemp`")
    
    
    myInteractome <- rbind(myInteractome,myInteractomeTemp)
  }
  
  if("IntAct_BINARY_Complete" %in% interactomesToUse){
    # load intact data
    # intact <- loadIntact()
    # save(intact,file="intact.RData")
    load("intact.RData")
    
    # keep only columns interactor A, interactor B, detection method, pubmed id, interaction type
    intact.tair <- intact[,c(5,6,7,9,12)]
    intact.tair$Publication.Identifier.s. <- gsub("pubmed:","",intact.tair$Publication.Identifier.s.)
    # determine quality to find binary multiple (BM) interactions
    # parameter: interaction set, only interactions in pho search space, set of binary methods
    myInteractionsIntact <- findQuality(intact.tair, FALSE, intact.binary, pho) 
    
    # myInteractionsIntact <- findQuality(intact.tair, FALSE, intact.binary)
    myInteractomeTemp <-
      myInteractionsIntact[which(
        myInteractionsIntact$interactionQuality == "binary" |
          myInteractionsIntact$interactionQuality == "HQ"
      ), 1:2]
    colnames(myInteractomeTemp) <- c("IntA","IntB")
    
    myInteractomeTemp <- data.frame(t(apply(myInteractomeTemp,1,sort)))
    colnames(myInteractomeTemp) <- c("IntA","IntB")
    
    myInteractionsTemp <- sqldf("select distinct `IntA`,`IntB` from `myInteractomeTemp`")
    
    myInteractome <- rbind(myInteractome,myInteractomeTemp)
  }
  
  if("BiogridSS_BM" %in% interactomesToUse){
    # biogrid <- loadBiogrid()
    # save(biogrid,file="biogrid.RData")
    load("biogrid.RData")
    biogrid.tair <- biogrid[, c(6,7,12,15,13)]
    myInteractionsBiogrid <- findQuality(biogrid.tair,FALSE, biogrid.binary, pho)
    
    myInteractomeTemp <- biogrid.tair[which(myInteractionsBiogrid$interactionQuality == "HQ"),1:2]
    
    colnames(myInteractomeTemp) <- c("IntA","IntB")
    
    indexA <- which(myInteractomeTemp$IntA %in% pho$Locus_ID)
    indexB <- which(myInteractomeTemp$IntB %in% pho$Locus_ID)
    indexC <- intersect(indexA,indexB)
    
    myInteractomeTemp <- myInteractomeTemp[indexC,]
    myInteractomeTemp <- data.frame(t(apply(myInteractomeTemp,1,sort)))
    colnames(myInteractomeTemp) <- c("IntA","IntB")
    
    myInteractionsTemp <- sqldf("select distinct `IntA`,`IntB` from `myInteractomeTemp`")
    
    myInteractome <- rbind(myInteractome,myInteractomeTemp)
  }
  
  if("BiogridSS_BS" %in% interactomesToUse){
    # biogrid <- loadBiogrid()
    # save(biogrid,file="biogrid.RData")
    load("biogrid.RData")
    biogrid.tair <- biogrid[, c(6,7,12,15,13)]
    myInteractionsBiogrid <- findQuality(biogrid.tair,FALSE, biogrid.binary, pho)
    
    myInteractomeTemp <- biogrid.tair[which(myInteractionsBiogrid$interactionQuality == "binary"),1:2]
    
    colnames(myInteractomeTemp) <- c("IntA","IntB")
    
    indexA <- which(myInteractomeTemp$IntA %in% pho$Locus_ID)
    indexB <- which(myInteractomeTemp$IntB %in% pho$Locus_ID)
    indexC <- intersect(indexA,indexB)
    
    myInteractomeTemp <- myInteractomeTemp[indexC,]
    myInteractomeTemp <- data.frame(t(apply(myInteractomeTemp,1,sort)))
    colnames(myInteractomeTemp) <- c("IntA","IntB")
    
    myInteractionsTemp <- sqldf("select distinct `IntA`,`IntB` from `myInteractomeTemp`")
    
    myInteractome <- rbind(myInteractome,myInteractomeTemp)
  }
  
  if("BioGRID_BINARY_Complete" %in% interactomesToUse){
    # biogrid <- loadBiogrid()
    # save(biogrid,file="biogrid.RData")
    load("biogrid.RData")
    biogrid.tair <- biogrid[, c(6,7,12,15,13)]
    myInteractionsBiogrid <- findQuality(biogrid.tair,FALSE, biogrid.binary, pho)
    
    myInteractomeTemp <-
      biogrid.tair[which(
        myInteractionsBiogrid$interactionQuality == "binary" |
          myInteractionsBiogrid$interactionQuality == "HQ"
      ), 1:2]
    
    colnames(myInteractomeTemp) <- c("IntA","IntB")
    
    myInteractomeTemp <- data.frame(t(apply(myInteractomeTemp,1,sort)))
    colnames(myInteractomeTemp) <- c("IntA","IntB")
    
    myInteractionsTemp <- sqldf("select distinct `IntA`,`IntB` from `myInteractomeTemp`")
    
    myInteractome <- rbind(myInteractome,myInteractomeTemp)
  }
  
  # restrict the selected interactome of all selected networks to search space
  if(restrictAllToSs & !("PvP" %in% interactomesToUse)){
    indexA <- which(myInteractome$IntA %in% pho$Locus_ID)
    indexB <- which(myInteractome$IntB %in% pho$Locus_ID)
    indexC <- intersect(indexA,indexB)
    
    ### find loci not in search space
    # indexAMissing <- which(!myInteractome$IntA %in% pho$Locus_ID)
    # indexBMissing <- which(!myInteractome$IntB %in% pho$Locus_ID)
    # print(myInteractome$IntA[indexAMissing])
    # print(myInteractome$IntB[indexBMissing])
    
    myInteractome <- myInteractome[indexC,]
    myInteractome <- data.frame(t(apply(myInteractome,1,sort)))
    colnames(myInteractome) <- c("IntA","IntB")
    
    myInteractome <- sqldf("select distinct `IntA`,`IntB` from `myInteractome`")
  }else{
    myInteractome <- data.frame(t(apply(myInteractome,1,sort)))
    colnames(myInteractome) <- c("IntA","IntB")
    
    myInteractome <- sqldf("select distinct `IntA`,`IntB` from `myInteractome`")
  }
  
  # remove self interactions = loops from the network
  myIndex <- which(as.character(myInteractome$IntA) == as.character(myInteractome$IntB))
  myInteractomeSelf <- c()
  myInteractomeFull <- myInteractome
  if(length(myIndex) > 0){
    myInteractomeSelf <- myInteractome[myIndex,]
    myInteractome <- myInteractome[-myIndex,]
  }
  
  myInteractomeList <- list(myInteractome,myInteractomeFull)
  
  
  # write interactions to the excel file
  if(writeToFile){
    write.xlsx(myInteractomeFull,paste(result.dir,file.name," - communities.xlsx",sep=""),sheetName = "Combined interactions",append = T)
  }
  
  return(myInteractomeList)
}

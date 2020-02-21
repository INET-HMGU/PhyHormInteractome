# load the ClueGO results for import in cytoscape as network attributes

suppressPackageStartupMessages(library(stringr))

load_cluego <- function(){
  
  # read the mapping file Symbol <-> AGI
  genes <- read.table("./ClueGOResults7/ClueGOResults7 Genes With Corresponding Functions.txt",
             header = T, sep = "\t", comment.char = "")
  
  # read the file containing the go groups
  groups <- read.table("./ClueGOResults7/ClueGOResults7 Functional Groups With Genes.txt",
                       header = T, sep = "\t", comment.char = "")
  
  # extract the AGIs from the alias column
  AGIs <- c()
  for(i in 1:length(genes$NAME)){
    tmp <- grepl("AT[1-5,C,M]G[0-9]{5}",genes$Aliases[i])
    if(!tmp){
      tmp2 <- grepl("AT[1-5,C,M]G[0-9]{5}",genes$NAME[i])
    }
    
    tmpAGI <- ""
    if(tmp){
      tmpAGI <- str_extract(genes$Aliases[i],"AT[1-5,C,M]G[0-9]{5}")
    }else if(tmp2){
      tmpAGI <- str_extract(genes$NAME[i],"AT[1-5,C,M]G[0-9]{5}")
    }else{
      stop("AGI not found")
    }
    
    AGIs <- c(AGIs,tmpAGI)
    
    # print(paste(i,tmpAGI))
  }
  
  
  # find corresponding groups for each AGI
  goGroupNumber <- rep("",length(AGIs))
  goGroupFunction <- rep("",length(AGIs))
  for(i in 1:length(groups$Function)){
    myGenes <- unlist(strsplit(as.character(groups$Group.Genes[i]),"\\|"))
    myIndex <- which(genes$NAME %in% myGenes)
    for(j in 1:length(myIndex)){
      if(goGroupNumber[myIndex[j]] == ""){
        goGroupNumber[myIndex[j]] <- as.character(groups$Groups[i])
      }else{
        goGroupNumber[myIndex[j]] <- paste(goGroupNumber[myIndex[j]],groups$Groups[i],sep=",")
      }
      
      if(goGroupFunction[myIndex[j]] == ""){
        goGroupFunction[myIndex[j]] <- as.character(groups$Function[i])
      }else{
        goGroupFunction[myIndex[j]] <- paste(goGroupFunction[myIndex[j]],groups$Function[i],sep=",")
      }
    }
  }
  
  cluegoData <- data.frame(AGIs,goGroupNumber,goGroupFunction)
  
  write.xlsx(cluegoData,paste(result.dir,file.name," - communities.xlsx",sep=""), sheetName = "ClueGo", append = T)
}



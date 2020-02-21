# Stefan Altmann
# Analyze network for communities and test, if communities are enriched in a hormone annotation
# this script was executed using R3.6.1 and Windows 10

rm(list=ls())

setwd("~/../Desktop/Github/Community analysis")

# load libraries
library(xlsx)
library(igraph)
library(sqldf)
library(gdata)
library(doParallel)
library(foreach)
library(utils)
library(iterators)
library(doParallel)
library(snow)
library(tcltk)


interactomesToUse <- c()

##################################### BEGIN PARAMETER

# parameters
dbg <- T # show debug information in main script
dbg2 <- F # show debug information in loaded functions
doRandomize <- T # do randomization test
randomNetworks <- 1000 # number of randomized networks 
SSExtended <- T # if T, then search space comprising 1252 loci is used, otherwise ss comprising 1227 loci is used
restrictAllToSs <- F # restricted to search space, if analyzing LCI_IntA or LCI_BioG
noCores <- 5 # number of cores to use for randomization

##################################### Network
# select the required networks!
myTestList <- c(1) # PhI_MAIN
# myTestList <- c(5) # LCI_IntA - Intact binary interactions
# myTestList <- c(15) # LCI_BioG - BioGRID binary interactions
# myTestList <- c(5,15) # run both LCI networks

sink.output <- F # sink output to file


outMean <- c()
outMedian <- c()

##################################### Community detection algorithm to use
# select the community finding algorithm
method <- "Edge Betweenness"
# method <- "Walktrap"
# method <- "Infomap"
# method <- "Fast Greedy" # 3
# method <- "Label Prop"
# method <- "Leading Eigen" # not run, R is crashing
# method <- "Louvain"
# method <- "Optimal" # not run, to slow
# method <- "Spinglass" # not run, does not work with unconnected graphs

# result directory
dateDir <- "Results/"
if(!dir.exists(dateDir)){
  dir.create(dateDir)  
}

# set number of steps, required only for a subset of community detection algorithm
walkSteps <- 5

# result directory for selected algorithm and parameter
methodDir <- paste(dateDir, method, " Step length ", walkSteps, "/", sep = "")
if(!dir.exists(methodDir)){
  dir.create(methodDir)  
}

##################################### Hormone annotations
# all interesting combinations of hormone annotations
# hormoneAnnotations <- c("AHD_Genetic","AHD_GO","AHD_ALL","TAIR_GO","ALL")
# hormoneAnnotations <- c("AHD_Genetic","TAIR_GO","AHDGen_TAIR")
# hormoneAnnotations <- c("ALL")
hormoneAnnotations <- c("AHDGen_TAIR")
# hormoneAnnotations <- c("AHD_Genetic") # 4
# hormoneAnnotations <- c("AHD_ALL")
# hormoneAnnotations <- c("AHD_Genetic","AHDGen_TAIR")
# hormoneAnnotations <- c("ALL")
# hormoneAnnotations <- c("AHD_GO","AHD_ALL","AHDGen_TAIR","TAIR_GO","ALL")
# hormoneAnnotations <- c("TAIR_GO")

# write results in functions to file
writeToFile <- T

##################################### END PARAMETER

# function to calculate the enrichment of communities
source("fun_calc_enrichment.r")

# function to load interactomes
source("fun_load_interactome.r")

# function to load hormone annotations
source("fun_load_hormone_annotation.r")

# function to load ClueGo Results
source("fun_load_cluego.r")

# function to find communities and calculate other network measures
source("fun_find_communities.r")

# function to find the proportion of interactions within / between communities
source("fun_find_interaction_proportion.r")

# calculate some network measure for network comparison
source("fun_calc_network_measures.r")

# calculate the length of the shortest paths between all combinations of hormone annotations
source("fun_find_shortest_paths_each_hormone.r")

# calc shortest paths, where only unique shortest paths are kept
# i.e. the shortest path is not allowed to contain a protein annotated with the same hormone
source("fun_find_unique_shortest_paths.r")

# calc shortest paths, where only unique shortest paths are kept
# i.e. the shortest path is not allowed to contain a protein annotated with the same hormone
# extension 2017/08/23: exlclude nodes, which have multiple annotations
source("fun_find_unique_shortest_paths_wo_multiple.r")

# calculate the hormone connectivity
source("fun_find_hormone_connectivity.r")

# iterate over different sets of hormone annotations
for(hormSet in hormoneAnnotations){
  # hormSet <- "AHDGen_TAIR"
  # output to file
  result.dir <- paste(methodDir, "Results - ", hormSet,"/",sep="")
  if(!dir.exists(result.dir)){
    dir.create(result.dir)  
  }
  
  if(restrictAllToSs){
    result.dir <- paste(result.dir,"/SS/",sep="")
    if(!dir.exists(result.dir)){
      dir.create(result.dir)  
    }
  }else{
    result.dir <- paste(result.dir,"/all/",sep="")
    if(!dir.exists(result.dir)){
      dir.create(result.dir)  
    }
  }
  
  if(dbg){
    print(paste(hormSet, collapse = ","))
  }
  
  # iterate of selected interactomes
  for(k in myTestList){
    # print loop iteration
    if(dbg){print(paste("Network collection",k))}
    
    ### select list of interactomes to use
    if(k == 1){
      interactomesToUse <- c("PvP")
    }else if(k == 5){
      interactomesToUse <- c("IntactSS_BM","IntactSS_BS")
    }else if(k == 15){
      interactomesToUse <- c("BiogridSS_BM","BiogridSS_BS")
    }else{
      stop("Set not found")
    }
    
    # file name and excel file set up
    file.name <- ""
    if(restrictAllToSs){
      file.name <- paste("eb - ",paste(interactomesToUse,collapse = ",")," (SS)",sep="")
    }else{
      file.name <- paste("eb - ",paste(interactomesToUse,collapse = ",")," (all)",sep="")
    }
    write.xlsx("readme",paste(result.dir,file.name," - communities.xlsx",sep=""),sheetName = "Read me",append = F)
    
    ### output to file ###########################
    if(sink.output){
      out.file.name <- paste(result.dir,"output-",paste(interactomesToUse,collapse=","),".txt",sep="")
      out.file <- file(description=out.file.name, open="w+")
      sink(file = out.file, append = TRUE, type = "output", split = FALSE)
    }
    
    ### load the hormone annotations
    if(dbg){print("Loading hormone annotations")}
    ahd <- load_hormone_annotation(hormSet)
    
    ### load the interactomes and make the graph
    if(dbg){print("Loading Interactomes")}
    myInteractome <- load_interactome(interactomesToUse, restrictAllToSs, dbg2, writeToFile, SSExtended)
    
    my.graph <- graph.data.frame(myInteractome[[1]][,1:2], directed=FALSE, vertices=NULL)
    my.graph.full <- graph.data.frame(myInteractome[[2]][,1:2], directed=FALSE, vertices=NULL)
    
    # set colors per node
    colorlist=c("#00CCCC","#0066CC","#6600CC","#CC00CC","#CC0000","#CC6600","#CCCC00","#00CC00","#CCCCCC","#FFFFFF")
    hormonelist=c("abscisic acid","auxin","brassinosteroid","cytokinin","ethylene","gibberellin","jasmonic acid","salicylic acid","multiple","none")
    
    for(myV in 1:length(V(my.graph.full))){
      # print(myV)
      myIndex <- which(ahd$AGI == names(V(my.graph.full)[myV]))
      if(length(myIndex) == 1){
        V(my.graph.full)$hormone[myV] <- ahd$summary[myIndex]
        myIndex2 <- which(hormonelist == ahd$summary[myIndex])
        V(my.graph.full)$color[myV] <- colorlist[myIndex2]
      }else{
        V(my.graph.full)$hormone[myV] <- "none"
        V(my.graph.full)$color[myV] <- colorlist[10]
      }
    }
    
    
    pdf(paste(result.dir,file.name," - network.pdf",sep=""),width=16,height = 10)
    
    l <- layout.fruchterman.reingold(my.graph.full, niter = 10)
    
    plot(my.graph.full, layout=l,vertex.size=5,vertex.label.cex = 0.5)
    
    dev.off()
    
    
    
    # find connectivity between and within proteins with same / different hormone annotation
    if(dbg){print("Finding hormone connectivity")}
    find_hormone_connectivity(my.graph,ahd,dbg2)
    
    # find communities
    if(dbg){print("Finding communities")}
    eb.com.list <- find_communities(my.graph.full,my.graph.full,method,walkSteps,TRUE, writeToFile)
    
    # calculate the shortest path between each combination of hormone annotations
    if(dbg){print("Finding shortest paths for each hormone combination")}
    find_shortest_paths_each(my.graph,ahd,dbg2)
    # source("fun_find_shortest_paths_each_hormone.r")
    # calculate the unique shortest path between each combination of hormone annotations
    if(dbg){print("Finding unique shortest paths for each hormone combination")}
    find_unique_shortest_paths(my.graph,ahd,dbg2)
    
    # calculate the unique shortest path between each combination of hormone annotations
    # exlcuding multiple annotated nodes
    # if(dbg){print("Finding unique shortest paths for each hormone combination wo multiple annotated nodes")}
    # find_unique_shortest_paths_wo_multiple(my.graph,ahd,dbg2,1)
    # find_unique_shortest_paths_wo_multiple(my.graph,ahd,dbg2,2)
    # find_unique_shortest_paths_wo_multiple(my.graph,ahd,dbg2,3)
    # find_unique_shortest_paths_wo_multiple(my.graph,ahd,dbg2,4)
    # find_unique_shortest_paths_wo_multiple(my.graph,ahd,dbg2,5)
    # find_unique_shortest_paths_wo_multiple(my.graph,ahd,dbg2,6)
    
    # calc network measures
    if(dbg){print("Calculating network measures")}
    calc_network_measures(my.graph,"linear",dbg2)
    calc_network_measures(my.graph,"power law",dbg2)
    
    
    ### ENRICHMENT ##########################################################################
    if(dbg){print("Determining enrichment")}
    
    # p.val.matrix.list <- calcEnrichement(myInteractome[[1]],dbg2) # 2016/12/16
    p.val.matrix.list <- calcEnrichement(myInteractome[[2]],dbg2) # 2016/12/16
    # p.val.matrix <- p.val.matrix.list[[1]]
    # p.val.matrix2 <- p.val.matrix.list[[2]]
    p.val.matrix3 <- p.val.matrix.list[[3]]
    
    # find index of enriched communities including non-annotated proteins
    # use bonferroni corrected p-values
    no.sig.enriched.real <- which(p.val.matrix3[,15] < 0.05)
    
    # append a column, which contains the hormone, in which the community is enriched
    hormoneEnrichment <- rep("none",length(p.val.matrix3[,1]))
    hormNames <- colnames(p.val.matrix3)[1:8]
    for(lineCount in no.sig.enriched.real){
      hormIndex <- which(p.val.matrix3[lineCount,1:8] == p.val.matrix3[lineCount,10])
      hormoneEnrichment[lineCount] <- hormNames[hormIndex]
    }
    
    p.val.matrix3 <- data.frame(p.val.matrix3,hormoneEnrichment)
    
    # write p-values to file
    write.xlsx(p.val.matrix3,paste(result.dir,file.name," - communities.xlsx",sep=""), sheetName = "P-Values3", append = T)
    
    ### check, which community is enriched and check in which hormone it is enriched
    enriched.community <- c()
    enriched.hormone <- c()
    for (mem in eb.com.list$membership) {
      if (p.val.matrix3[mem, 15] < 0.05) {
        enriched.community <- c(enriched.community, T)
        enriched.hormone <-
          c(enriched.hormone, as.character(p.val.matrix3[mem, 16]))
      } else{
        enriched.community <- c(enriched.community, F)
        enriched.hormone <- c(enriched.hormone, "none")
      }
    }
    
    # combine values
    community.data <-
      cbind(eb.com.list$names,
            eb.com.list$membership,
            enriched.community,
            enriched.hormone)
    
    colnames(community.data) <-
      c("Locus", "Community", "Enriched", "Enriched Hormone")
    
    # write values
    write.xlsx(
      community.data,
      paste(result.dir, file.name, " - communities.xlsx", sep = ""),
      sheetName = "Community membership - enriched",
      append = T,
      row.names = F
    )
    
    ### calculate the proportion of interactions between / within communities
    # if(dbg){print("Finding interaction ratio of communities")}
    # find_interaction_proportion(my.graph, eb.com.list, p.val.matrix3, dbg2)
    
    
    ######################### RANDOM Network analysis for enrichment ##########################
    if(doRandomize){
      if(dbg){print("Randomizing network")}
      # for randomized network
      no.sig.enriched.random <- c()
      
      # registerDoParallel(cores=6)  
      # Start a cluster
      cl <- makeCluster(noCores)
      registerDoParallel(cl)
      # run network randomization in parallel
      # no.sig.enriched.random <- foreach(loop=1:randomNetworks, .combine = f(randomNetworks), .packages = 'igraph') %dopar% {
      no.sig.enriched.random <- foreach(i = icount(randomNetworks), .combine = c, .packages = c('igraph','tcltk')) %dopar% {
        # if(dbg){print(loop)}
        if(!exists("pb")) pb <- tkProgressBar(as.character(i), min=1, max=randomNetworks)
        setTkProgressBar(pb, i)
        # my.graph.random <- rewire(my.graph,niter=50000)
        my.graph.random <- rewire(my.graph, with = keeping_degseq(loops = FALSE, niter = vcount(my.graph) * 10))
        # communities - edge betweenness community for random graph
        # eb.com.list <- edge.betweenness.community(my.graph.random,weights=NULL)
        # eb.com.list <- cluster_edge_betweenness(my.graph.random,weights=NULL)
        # eb.com.list <- walktrap.community(my.graph.random,steps = walkSteps)
        # eb.com.list <- cluster_infomap(my.graph.random,nb.trials = walkSteps)
        if(method == "Walktrap"){
          eb.com.list <- cluster_walktrap(my.graph.random,steps = walkSteps)
        }else if(method == "Infomap"){
          eb.com.list <- cluster_infomap(my.graph.random,nb.trials = walkSteps)
        }else if(method == "Edge Betweenness"){
          eb.com.list <- cluster_edge_betweenness(my.graph.random,weights=NULL,edge.betweenness = T,merges = T, bridges = T)
        }else if(method == "Fast Greedy"){
          eb.com.list <- cluster_fast_greedy(my.graph.random)
        }else if(method == "Label Prop"){
          eb.com.list <- cluster_label_prop(my.graph.random)
        }else if(method == "Leading Eigen"){
          eb.com.list <- cluster_leading_eigen(my.graph.random)
        }else if(method == "Louvain"){
          eb.com.list <- cluster_louvain(my.graph.random)
        }else if(method == "Optimal"){
          eb.com.list <- cluster_optimal(my.graph.random)
        }else if(method == "Spinglass"){
          eb.com.list <- cluster_spinglass(my.graph.random)
        }else{
          stop("Method for community detection not found")
        }
        
        # calculate the enrichtment
        randomInteractions <- as_data_frame(my.graph.random, what = "edges")
        colnames(randomInteractions) <- c("IntA","IntB")
        p.val.matrix.list <- calcEnrichement(randomInteractions,T)
        p.val.matrix <- p.val.matrix.list[[3]]
        
        no.sig.enriched <- which(p.val.matrix[,11] < 0.05)
        length(no.sig.enriched)
      }
      
      #Stop the cluster
      stopCluster(cl)
      
      # plot curve and boxplot in pdf      
      pdf(paste(result.dir,file.name," - enriched.communities.pdf",sep=""),width = 4,height = 4)
      
      # calculate the p-value of the enrichment
      realEnrichCom <- length(no.sig.enriched.real)
      
      my.table <- data.frame(table(no.sig.enriched.random))
      randomCounter <- length(which(no.sig.enriched.random >= realEnrichCom))
      myPvalue <- NA
      myPvaluePrint <- ""
      if(randomCounter > 0){
        myPvalue <-  randomCounter / randomNetworks
        myPvalue <- signif(myPvalue,2)
        myPvaluePrint <- as.character(myPvalue)
      }else{
        myPvalue <-  1 / randomNetworks
        myPvaluePrint <- paste("<",as.character(myPvalue),sep="")
      }
      
      # curve
      plot(seq(0,(length(my.table$Freq)-1)),my.table$Freq,xlim=c(0,max(no.sig.enriched.random,realEnrichCom)+2),type="b",ylim=c(0,max(my.table$Freq)),
           ylab="Frequency",xlab="No. enriched communities",main="")
      arrows(realEnrichCom,randomNetworks/5,realEnrichCom,0,col = "red")
      mtext(text = paste("Repeats ",randomNetworks, ", P-Value: ",myPvaluePrint,sep=""), side = 3)
      
      # boxplot
      boxplot(no.sig.enriched.random,ylim=c(0,max(no.sig.enriched.random,realEnrichCom)+2))
      par(new=T)
      boxplot(data.frame(realEnrichCom),ylim=c(0,max(no.sig.enriched.random,realEnrichCom)+2),border="red")
      # close pdf
      dev.off()
      
    }
    
    ### end redirection to file ###############################################################
    # sink(type = "message")
    # sink()
    # try to close all connections
    closeAllConnections()
    # check number of open connections
    sink.number()
    # close existing output streams
    if(sink.number() > 0){
      close.files <- sink.number()
      for(i in 1:close.files){
        sink(type = "message")
        sink()
      }
    }
    
  } # end loop over interactomes
  
  # moveDir <- paste("Results - ", hormSet, sep = "")
  # move directory
  # file.copy(moveDir,paste(method,"Step length",walkSteps), recursive = T)
  # unlink(moveDir, recursive = T)
  
  
} # end loop over annotation sets

print("Finished")



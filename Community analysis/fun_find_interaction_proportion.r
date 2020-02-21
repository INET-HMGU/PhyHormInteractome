
# piechart: valuelist="1,0,0,0,0,0,0,0" 
# colorlist=c("#00CCCC","#0066CC","#6600CC","#CC00CC","#CC0000","#CC6600","#CCCC00","#00CC00") 
# labellist="abscisic



suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages((library(ggrepel)))

# find and return the biggest connected component in the graph
giant.component <- function(graph) {
  cl <- components(graph)
  induced_subgraph(graph, which(cl$membership == which.max(cl$csize)))
}

# multiple ggplots on one page
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  suppressPackageStartupMessages(library(grid))
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# find the proportion of internal/external interactions per community
find_interaction_proportion <- function(my.graph, eb.com.list, p.val.matrix, dbg){
  # remove all non-connected components from the biggest connected component
  # the non connected components would receive an artificial small mean shortest path
  my.graph.small <- giant.component(my.graph)
  
  # ratio absolute
  ratioInternExtern <- c()
  # Community size
  communitySize <- c()
  # community number
  communityNumber <- c()
  for(i in 1:length(eb.com.list)){
    if(dbg){print(paste(i,"#################"))}
    if(all(eb.com.list[[i]] %in% names(V(my.graph.small)))){
      
      comSize <- length(eb.com.list[[i]])
      my.subgraph <- induced.subgraph(my.graph.small,eb.com.list[[i]])
      internEdges <- length(E(my.subgraph))
      
      externEdges <- 0
      for(j in 1:length(eb.com.list[[i]])){
        neighborhood(my.graph,1,nodes = eb.com.list[[i]][j])
        my.neighbors <- neighbors(my.graph.small,eb.com.list[[i]][j])
        externEdges <- externEdges + length(which(!names(my.neighbors) %in% eb.com.list[[i]]))
      }
      
      if(dbg){
        print(paste("Internal edges:", internEdges))
        print(paste("External edges:", externEdges))
      }
      
      ratioInternExtern <- c(ratioInternExtern,internEdges/externEdges)
      communitySize <- c(communitySize,comSize)
      communityNumber <- c(communityNumber,i)
    }
  }
  
  
  # determine colors of points
  colorlist=c("#00CCCC","#0066CC","#6600CC","#CC00CC","#CC0000","#CC6600","#CCCC00","#00CC00") 
  
  pSmall <- data.frame(p.val.matrix[communityNumber,])
  
  if(ncol(pSmall) == 1){
    sigHormoneIndex <- data.frame(which(pSmall[1:8,1] < 0.05,arr.ind = T))
  }else{
    temp <- which(pSmall[,1:8] < 0.05,arr.ind = T)
    rownames(temp) <- NULL
    sigHormoneIndex <- data.frame(temp)  
  }
  
  
  colorVec <- rep("#CCCCCC",length(communitySize))
  
  # this loop does not work correct, if one community is found and this community is 
  # enriched for two hormones
  for(i in 1:length(sigHormoneIndex[,1])){
    colorVec[sigHormoneIndex$row[i]] <- colorlist[sigHormoneIndex$col[i]]
  }
  
  logRatio <- data.frame(log2(ratioInternExtern))
  # p1 <- ggplot(logRatio,aes(logRatio)) + 
  #   geom_histogram(binwidth = 1) + scale_x_continuous(expand = c(0.05,0.05)) +
  #   labs(title = "Ratio of internal to external edges", x = "log intern/extern edges ratio", y = "Number of communities")
  p1 <- ggplot(logRatio,aes(log2.ratioInternExtern.)) + 
    geom_histogram(binwidth = 1) + scale_x_continuous(expand = c(0.05,0.05)) +
    labs(title = "Ratio of internal to external edges", 
         x = "log intern/extern edges ratio", 
         y = "Number of communities")
  
  dat <- data.frame(communitySize,log2(ratioInternExtern))
  colnames(dat) <- c("x","y")
  p2 <- ggplot(dat, aes(x,y)) + geom_point(colour = colorVec, size = communitySize) + 
    geom_text_repel(aes(x,y),label = communityNumber, box.padding = unit(2, "lines")) + 
    labs(title = "Community size vs internal/external edges ratio", x = "Community size", y = "log intern/extern edges ratio")
  # plot(p2)
  pdf(paste(result.dir,file.name," - Internal to External Interactions Ratio.pdf",sep=""), width = 10, height = 5)
  multiplot(p1, p2, cols=2)
  dev.off()
  
  pdf(
    paste(
      result.dir,
      file.name,
      " - Internal to External Interactions Ratio LARGE 1.pdf",
      sep = ""
    ),
    width = 25,
    height = 25
  )
  print(p2)
  dev.off()
  
  pdf(
    paste(
      result.dir,
      file.name,
      " - Internal to External Interactions Ratio LARGE 2.pdf",
      sep = ""
    ),
    width = 10,
    height = 10
  )
  print(p2)
  dev.off()
}
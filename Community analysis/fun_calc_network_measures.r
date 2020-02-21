#

suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(QuACN))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
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

calc_network_measures <- function(my.graph,curve_model,dbg){
  # convert igraph to graphnel format
  my.graphNel <- igraph.to.graphNEL(my.graph)
  
  # degree
  dd <- data.frame(QuACN::degreeDistribution(my.graphNel))
  dd <- cbind(dd,dd$Freq / sum(dd$Freq))
  colnames(dd) <- c("Degree","Frequency","Fraction")
  
  # clustering coefficient
  deg <- graph::degree(my.graphNel)
  lcc <- QuACN::localClusteringCoeff(my.graphNel,deg)
  
  uniqueDegrees <- sort(unique(deg))
  ccMeans <- c()
  for(i in 1:length(uniqueDegrees)){
    myMean <- mean(lcc[which(deg == uniqueDegrees[i])])
    ccMeans <- c(ccMeans,myMean)
    if(dbg){print(myMean)}
  }
  
  pdf(paste(result.dir,file.name," - Network measures (",curve_model,").pdf",sep=""), width = 6, height = 6)
  
  
  ### degree distribution
  
  
  
  dat <- data.frame(as.numeric(as.character(dd$Degree)),dd$Fraction)
  
  colnames(dat) <- c("x","y")
  toRemove <- which(dat$y == 0)
  if(length(toRemove) > 0){
    dat <- dat[-toRemove,]  
  }
  
  # myModel <- lm(y~x, data = dat,na.action = na.exclude)
  # myModel <- lm(log(y)~log(x), data = dat,na.action = na.exclude)
  if(curve_model == "linear"){
    myModel <- lm(y~x, data = dat,na.action = na.exclude)
    
    q1a <- ggplot(dat, aes(x,y)) + geom_point(colour = 'blue', size = 2) + 
      labs(title = "Degree distribution", x = "k", y = "P(k)") + 
      stat_smooth(method = "lm", formula = as.formula(myModel), size = 1, se = FALSE, colour = "grey", na.rm = T)
    
  }else if(curve_model == "power law"){
    myModel <- nls(y~b*x^z,data=dat)  
    
    q1a <- ggplot(dat, aes(x,y)) + geom_point(colour = 'blue', size = 2) + 
      labs(title = "Degree distribution", x = "k", y = "P(k)") + 
      stat_smooth(method = "nls", formula = as.formula(myModel), size = 1, se = FALSE, colour = "grey", na.rm = T)
    
  }else{
    print("curve model not valid. Select power law or linear")
    stop()
  }
  
  
  
  
  # q1b <- q1a + coord_trans(x="log10",y="log10") + 
  #   labs(title = "Degree distribution (log 10)", x = "log k", y = "log P(k)")
  
  # q1c <- q1a + scale_y_continuous(trans=log10_trans(), na.value = NA_real_) +
  #   scale_x_continuous(trans=log10_trans(), na.value = NA_real_, ) +
  #   labs(title = "Degree distribution (log 10)", x = "log k", y = "log P(k)")
  
  # doLog <- T
  # if(all(dat$y) == 0){
  #   doLog <- F
  # }
  
  # if(doLog)
  q1c <- q1a + scale_y_log10() + 
    scale_x_log10() + 
    labs(title = "Degree distribution (log 10)", x = "log k", y = "log P(k)")
  
  ### clustering coefficient
  

  dat <- data.frame(uniqueDegrees,ccMeans)
  colnames(dat) <- c("x","y")
  toRemove <- which(dat$y == 0)
  if(length(toRemove) > 0){
    dat <- dat[-toRemove,]  
  }
  dat$y[which(dat$y == -Inf)] <- NA
  # myModel <- lm(y~x, data = dat,na.action = na.exclude)
  
  if(curve_model == "linear"){
  q2a <- ggplot(dat, aes(x,y)) + geom_point(colour = 'red', size = 2) + 
    labs(title = "Clustering coefficient", x = "k", y = "C(k)") + 
    stat_smooth(method = "lm", formula = as.formula(myModel), size = 1, se = FALSE, colour = "grey", na.rm = T)
  }else if(curve_model == "power law"){
    q2a <- ggplot(dat, aes(x,y)) + geom_point(colour = 'red', size = 2) + 
      labs(title = "Clustering coefficient", x = "k", y = "C(k)") + 
      stat_smooth(method = "nls", formula = as.formula(myModel), size = 1, se = FALSE, colour = "grey", na.rm = T)
  }else{
    # nothing to do  
  }
  
  if(all(is.infinite(log10(dat$y)))){
    q2c <- q2a
    }else{
      q2c <- q2a + scale_y_continuous(trans=log10_trans(), na.value = NA_real_) + 
        scale_x_continuous(trans=log10_trans(), na.value = NA_real_) + 
        labs(title = "Clustering coefficient (log 10)", x = "log k", y = "log C(k)")
  }
  
  multiplot(q1a,q2a,q1c,q2c,cols=2)

  dev.off()
}
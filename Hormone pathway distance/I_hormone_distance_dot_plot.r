# 20180427
# create nice plots for hormone distance
# this script was last executed using R3.6.2 and Windows 10

rm(list = ls())

setwd("~/../Desktop/Github/Hormone distance")

library(xlsx)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

library(igraph)

# adjust directories, from where to read data
dirs <- c("../Community analysis/Results/Edge Betweenness Step length 5/Results - AHDGen_TAIR/")

# adjust networks, which should be plotted
networks <-
  c(
    "all/eb - PvP (all)",
    "SS/eb - IntactSS_BM,IntactSS_BS (SS)",
    "SS/eb - BiogridSS_BM,BiogridSS_BS (SS)"
  )

doPlot <- function(dir, network){
  # debug only
  # dir <- dirs[1]
  # network <- networks[1]
  
  # read the community interactions
  edges <- read.xlsx(paste(dir, network, " - communities.xlsx", sep = ""), sheetName = "Hormone connectivity - Edges")
  nodes <- read.xlsx(paste(dir, network, " - communities.xlsx", sep = ""), sheetName = "Hormone connectivity - Nodes")
  interactions <- read.xlsx(paste(dir, network, " - communities.xlsx", sep = ""), sheetName = "Combined interactions")
  
  myGraph <- graph.data.frame(interactions[,2:3], directed = F)
  myDiameter <- igraph::diameter(myGraph, directed = F, unconnected = T)
  myDiameter <- as.data.frame(myDiameter)
  colnames(myDiameter) <- c("Size")
  
  # ggplot(myDiameter, aes(x = 1, y = 1, size = Size)) + geom_point(shape = 1) + scale_size_continuous(range = c(1,50))
  # + geom_point(shape = "9")
  
  # read the path lengths
  # pathlengthMedian <- read.xlsx(paste(dir, network, " - Shortest Paths - unique path.xlsx", sep = ""), sheetName = "Median Values")
  pathlengthMean <- read.xlsx(paste(dir, network, " - Shortest Paths - unique path.xlsx", sep = ""), sheetName = "Mean Values")
  
  # melt matrix to data frame
  pl <- reshape2::melt(pathlengthMean)
  pl <- pl[-which(is.na(pl$value)),]
  
  # convert long hormone names to abbreviations
  longname <- c("abscisic_acid", "auxin", "brassinosteroid", "cytokinin", 
                "ethylene", "gibberellin", "jasmonic_acid", "salicylic_acid")
  shortname <- c("ABA", "AUX", "BR", "CK", "ET", "GA", "JA", "SA")
  
  for(i in 1:8){
    nodes$Hormone <- sub(longname[i], shortname[i], nodes$Hormone)  
    edges$Horm1 <- sub(longname[i], shortname[i], edges$Horm1)
    edges$Horm2 <- sub(longname[i], shortname[i], edges$Horm2)  
  }
  
  # extract path lenght in correct order and add to vector
  pathlength <- c()
  for(i in 1:length(edges$Horm1)){
    myIndex <- which((pl$NA. == edges$Horm1[i] & pl$variable == edges$Horm2[i]) |
                       (pl$NA. == edges$Horm2[i] & pl$variable == edges$Horm1[i]))
    pathlength <- c(pathlength, pl$value[myIndex])
  }
  
  # add path length values to edges data frame
  edges <- cbind(edges, pathlength)
  colnames(edges) <- c("No", "Hormone1", "Hormone2", "Connections", "Pathlength")
  
  # make plot
  # p <- ggplot(edges, aes(x = Hormone1, y = Hormone2, size = Connections)) + 
  #   geom_point(aes(color = pathlength)) + 
  #   scale_size_area(max_size = 15) + 
  #   scale_colour_gradientn(colours = brewer.pal(12, "Paired"), limits = c(0,12))
  # 
  # p <- ggplot(edges, aes(x = Hormone1, y = Hormone2, size = Connections)) + 
  #   geom_point(aes(color = pathlength)) + 
  #   scale_size_area(max_size = 15) + 
  #   scale_colour_gradientn(colours = topo.colors(13), limits = c(0,12)) + 
  #   scale_size()
  
  uplimitCon <- max(40,max(edges$Connections))
  p <- ggplot(edges, aes(x = Hormone1, y = Hormone2, size = Connections)) + 
    geom_point(aes(color = pathlength)) + 
    scale_size_continuous(range = c(2,15), limits = c(0,uplimitCon)) + 
    scale_colour_gradientn(colours = topo.colors(13), limits = c(0,12))
  
  
  
  # p
  
  # pEmpty <- ggplot(edges, aes(x = Hormone1, y = Hormone2, size = Connections)) +
  #   geom_blank() +
  #   theme(axis.text = element_blank(),
  #         axis.title = element_blank(),
  #         line = element_blank(),
  #         panel.background = element_blank())
  
  colnames(nodes)[3] <- "Number.Loci"
  # pTop <- ggplot(nodes, aes(x = Hormone, y = rep(1,8), size = Number.Loci)) +
  #   geom_point() +
  #   scale_size_area(max_size = 15)
  
  uplimitNodes <- max(40, max(nodes$Number.Loci))
  pTop <- ggplot(nodes, aes(x = Hormone, y = rep(1,8), size = Number.Loci)) +
    geom_point() +
    scale_size_continuous(range = c(2,15), limits = c(0,uplimitNodes))
  # pTop
  
  # https://deanattali.com/2015/03/29/ggExtra-r-package/
  
  # add histogram of connectivity
  # pins <- ggplot(edges, aes(Connections)) + 
  #   geom_histogram() + labs(x = NULL, y = NULL) 
  
  pins <- ggplot(edges, aes(pathlength, fill = "")) +
    geom_density(col=2, alpha=0.2)+ labs(x = NULL, y = NULL) +
    scale_fill_discrete(guide=FALSE) +
    xlim(range(c(0,12))) +
    theme_classic()
  
  # pins
  g <- ggplotGrob(pins)
  
  p <- p + annotation_custom(
    grob = g,
    xmin = 5,
    xmax = 8,
    ymin = 1,
    ymax = 4
  )
  
  
  # p
  # grid.arrange(pTop, pEmpty, p,  
  #              ncol = 2, nrow = 2, widths = c(3, 2), heights = c(1, 3))

  pdf(paste(dir, network, " - connections and pathlength 4.pdf", sep = ""), width = 6, height = 7)
  grid.arrange(pTop, p,  
               ncol = 1, nrow = 2, heights = c(1,3))
  dev.off()
}

for(dir in dirs){
  for(network in networks){
    doPlot(dir, network)
  }
}









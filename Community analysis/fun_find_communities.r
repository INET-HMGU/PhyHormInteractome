
find_communities <- function(my.graph,my.graph.full,method,walkSteps,dbg, writeToFile){
  
  if(method == "Walktrap"){
    eb.com.list <- cluster_walktrap(my.graph,steps = walkSteps)
  }else if(method == "Infomap"){
    eb.com.list <- cluster_infomap(my.graph,nb.trials = walkSteps)
  }else if(method == "Edge Betweenness"){
    eb.com.list <- cluster_edge_betweenness(my.graph,weights=NULL,edge.betweenness = T,merges = T, bridges = T)
  }else if(method == "Fast Greedy"){
    eb.com.list <- cluster_fast_greedy(my.graph)
  }else if(method == "Label Prop"){
    eb.com.list <- cluster_label_prop(my.graph)
  }else if(method == "Leading Eigen"){
    eb.com.list <- cluster_leading_eigen(my.graph)
  }else if(method == "Louvain"){
    eb.com.list <- cluster_louvain(my.graph)
  }else if(method == "Optimal"){
    eb.com.list <- cluster_optimal(my.graph)
  }else if(method == "Spinglass"){
    eb.com.list <- cluster_spinglass(my.graph)
  }else{
    stop("method not found")
  }
  # write data to xlsx
  myData <- data.frame(eb.com.list$names,eb.com.list$membership)
  colnames(myData) <- c("Locus","Community")
  if(writeToFile){
  write.xlsx(myData,paste(result.dir,file.name," - communities.xlsx",sep=""),sheetName = "Community membership",append = T)
  # write hormone annotation to excel file
  write.xlsx(ahd,paste(result.dir,file.name," - communities.xlsx",sep=""),sheetName = "Hormone annotation",append = T)
  
  }
  
  # write edge betweenness of removed edges to excel file
  myEdges <- get.edgelist(my.graph)
  
  if(length(grep("edge.betweenness",names(eb.com.list))) > 0){
    removedEdges <- eb.com.list$removed.edges
    myEdges <- myEdges[removedEdges,]
    myEdgeData <- data.frame(myEdges,eb.com.list$edge.betweenness)
    if(writeToFile){
      write.xlsx(myEdgeData,paste(result.dir,file.name," - communities.xlsx",sep=""),sheetName = "Edge betweenness (Comm)",append = T)
      }
  }
  
  # find the interactions between the communities
  # and determine the size of the communities
  if(dbg){print("Determining interactions between communities")}
  community.interactions <- c()
  for(i in 1:length(eb.com.list$membership)){
    my.neighbors.index <- unlist(neighborhood(my.graph,1,V(my.graph)[which(V(my.graph)$name == eb.com.list$names[i])]))
    my.neighbors <- V(my.graph)$name[my.neighbors.index]
    my.connections <- eb.com.list$membership[which(eb.com.list$name %in% my.neighbors)]
    my.source <- rep(eb.com.list$membership[i],length(my.connections))
    community.interactions <- rbind(community.interactions,cbind(my.source,my.connections))
  }
  community.interactions <- data.frame(community.interactions)
  remove <- which(community.interactions$my.source == community.interactions$my.connections)
  community.interactions <- community.interactions[-remove,]
  community.interactions.sorted <- data.frame(t(apply(community.interactions,1,sort)))
  community.interactions.sorted.all <- community.interactions.sorted
  suppressPackageStartupMessages(library(plyr))
  community.interactions.sorted <- ddply(community.interactions.sorted,.(X1,X2),unique)
  colnames(community.interactions.sorted) <- c("CommunityA","CommunityB")
  
  # count occurrences of each community interaction
  commIntCounter <- c()
  for(i in 1:length(community.interactions.sorted$CommunityA)){
    cic <- length(which(community.interactions.sorted.all$X1 == community.interactions.sorted$CommunityA[i] &
            community.interactions.sorted.all$X2 == community.interactions.sorted$CommunityB[i] ))
    commIntCounter <- c(commIntCounter, cic)
  }
  community.interactions.sorted <- cbind(community.interactions.sorted,commIntCounter/2)
  
  if(writeToFile){
    write.xlsx(community.interactions.sorted,paste(result.dir,file.name," - communities.xlsx",sep=""), 
             sheetName = "Community interactions",append = T)
  }
  
  community.size <- data.frame(table(eb.com.list$membership))
  colnames(community.size) <- c("Community", "Size")
  
  if(writeToFile){
    write.xlsx(data.frame(table(eb.com.list$membership)),paste(result.dir,file.name," - communities.xlsx",sep=""), 
             sheetName = "Community Size", append = T)
  }
  
  # determine betweenness of nodes independent of communities
  myNodeBetweenness <- betweenness(my.graph.full)
  if(writeToFile){
    write.xlsx(data.frame(myNodeBetweenness),paste(result.dir,file.name," - communities.xlsx",sep=""), 
             sheetName = "Node Betweenness", append = T)
  }
  # determine betweenness of edges independent of communities
  myEdgeBetweenness <- edge_betweenness(my.graph.full)
  myEdgeBetweenness <- data.frame(get.edgelist(my.graph.full),myEdgeBetweenness)
  # check, if the interaction is contained in lci
  load("intact.RData")
  intact.small <- intact[,c(5,6)]
  colnames(intact.small) <- c("IntA","IntB")
  confirmed <- rep(0,length(myEdgeBetweenness$X1))
  
  myCommonVec <- c()
  for(i in 1:length(myEdgeBetweenness$X1)){
    myIndex <- which((myEdgeBetweenness$X1[i] == as.character(intact.small$IntA)) & (myEdgeBetweenness$X2[i] == as.character(intact.small$IntB)) |
                       (myEdgeBetweenness$X2[i] == as.character(intact.small$IntA)) & (myEdgeBetweenness$X1[i] == as.character(intact.small$IntB)) )
    if(length(myIndex) > 0){
      myCommonVec <- c(myCommonVec,i)
    }
  }
  
  confirmed[unique(myCommonVec)] <- 1
  
  myEdgeBetweenness <- cbind(myEdgeBetweenness,confirmed)
  colnames(myEdgeBetweenness) <- c("IntA","IntB","EdgeBetweenness","Confirmed")
  
  if(writeToFile){
    write.xlsx(myEdgeBetweenness,paste(result.dir,file.name," - communities.xlsx",sep=""), 
             sheetName = "Edge Betweenness", append = T,row.names = F)
  }
  # determine the topology of the network
  # stop("determine the topology of the network")
  
  # calculate the shortest paths between proteins annotated with same hormone and proteins annotated with different hormones 
  # stop("calculate the shortest paths between proteins with annotation")
  
  # calculate the ratio of connections within communities vs. connections between communities for later use
  # stop("calculate the ratio of connections within communities vs. connections between communities for later use")
  
  ######################### new 2018-11-18
  ### determine various measures of the community network
  comInt <- graph_from_data_frame(community.interactions.sorted)
  comDegrees <- igraph::degree(comInt, v = V(comInt))
  meanDegree <- mean(comDegrees)
  clusteringCoeff <- igraph::transitivity(comInt, type = "global")
  assort <- assortativity_degree(comInt)
  
  CommProp <- c()
  CommProp <- rbind(CommProp,c("Mean Degree",meanDegree))
  CommProp <- rbind(CommProp,c("Clustering Coefficient",clusteringCoeff))
  CommProp <- rbind(CommProp,c("Assortativity",assort))
  print("############### Community network properities")
  print(paste("mean degree:", meanDegree))
  print(paste("clustering coefficient:", clusteringCoeff))
  print(paste("assortativity:", assort))
  print("#############################################")
  write.xlsx(CommProp,paste(result.dir,file.name," - communities.xlsx",sep=""),sheetName = "Community Net Properties",append = T)
  
  write.xlsx(comDegrees,paste(result.dir,file.name," - communities.xlsx",sep=""),sheetName = "Community Degrees",append = T)
  
  
  return(eb.com.list)
}
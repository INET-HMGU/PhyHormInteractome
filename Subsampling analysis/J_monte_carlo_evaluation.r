# Stefan Altmann
# this script was executed using R3.6.1 and Windows 10
# count and plot number of solid (Type I) and possible (Type II) crosstalks 1000 times
# in 100 randomly selected interactions from PhIMAIN, BioGRID_SS_BMBS, and IntAct_SS_BMBS
# draw distributions for all crosstalks and specific hormone combinations and compare them

rm(list = ls())

setwd("~/../Desktop/Github/Subsampling")

library(openxlsx)

library(ggplot2)

annoSet <- c("HormAnno-AHDGen,tairGO",
             "HormAnno-AHDGen",
             "HormAnno-tairGO10",
             "HormAnno-AHDGen,tairGO10")

# annoSet <- c("HormAnno-AHDGen,tairGO10")

selectedNetwork <- c("PhIMAIN","Intact_SS_BMBS","Biogrid_SS_BMBS")

# result directory
result.dir <- "Results"

pcp.types <- c(1,2)

for(pcp.type in pcp.types){
# read the results of the monte carlo simulation
for(k in 1:length(annoSet)){
  # store simulation results for each network
  myResults <- list()
  sumPerRun <- list()
  
  for(s in 1:length(selectedNetwork)){
    # read results
    myResults[[s]] <- read.xlsx(paste(result.dir,"/Crosstalks-",selectedNetwork[s],"-",annoSet[k],".xlsx", sep = ""),
              sheet = pcp.type, rowNames = T, colNames = T)
    
    # calc sum per row, which is sum per random network
    sumPerRun[[s]] <- apply((myResults[[s]]), 1, sum)
  }
  
  # set names
  names(sumPerRun) <- selectedNetwork
  # convert to data frame
  df <- data.frame(matrix(unlist(sumPerRun), nrow=length(sumPerRun[[1]]), byrow=F))
  colnames(df) <- selectedNetwork
  # stack data frame
  dfs <- stack(df)
  
  pdf(paste(result.dir, "/PCPs Type ",pcp.type," Distributions - ", annoSet[k], ".pdf", sep = ""))
  p <- ggplot(dfs, aes(x=values)) + geom_density(aes(group=ind, colour=ind, fill=ind), alpha=0.3)
  print(p)
  
  p <- ggplot(dfs, aes(y=values,x=ind)) + geom_boxplot(aes(group=ind, colour=ind, fill=ind), alpha=0.3)
  print(p)
  
  dev.off()
  
  tpi <- t.test(df$PhIMAIN, df$Intact_SS_BMBS)
  print("#########")
  print(tpi)
  print(tpi$p.value)
  
  tpb <- t.test(df$PhIMAIN, df$Biogrid_SS_BMBS)
  print("#########")
  print(tpb)
  print(tpb$p.value)
  
  # calculate p-value per hormone combination
  pdf(paste(result.dir, "/PCPs Type ",pcp.type," Distributions per Hormone Combination - ", annoSet[k], ".pdf", sep = ""))
  
  hormCombs <- names(myResults[[1]])
  pValues.pi <- c()
  pValues.bp <- c()
  for(i in 1:length(hormCombs)){
    p.horm <- myResults[[1]][,i] # pvp
    i.horm <- myResults[[2]][,i] # intact
    b.horm <- myResults[[3]][,i] # biogrid
    
    wpi <- wilcox.test(p.horm, i.horm)
    # print("#########")
    # print(wpi)
    # print(wpi$p.value)
    pValues.pi <- c(pValues.pi, wpi$p.value)
    
    wpb <- wilcox.test(p.horm, b.horm)
    # print("#########")
    # print(wpb)
    # print(wpb$p.value)
    pValues.bp <- c(pValues.bp, wpb$p.value)
    
    df <- data.frame(p.horm, i.horm, b.horm)
    colnames(df) <- selectedNetwork
    # stack data frame
    dfs <- stack(df)
    
    p <- ggplot(dfs, aes(x=values)) + geom_density(aes(group=ind, colour=ind, fill=ind), alpha=0.3)
    p <- p + ggtitle(hormCombs[i])
    print(p)
    
    p <- ggplot(dfs, aes(y=values,x=ind)) + geom_boxplot(aes(group=ind, colour=ind, fill=ind), alpha=0.3)
    p <- p + ggtitle(hormCombs[i])
    print(p)
    
  }
  
  dev.off()
  
  
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "PVal PhIMAIN-Intact")
  openxlsx::writeData(wb, sheet = "PVal PhIMAIN-Intact", data.frame(hormCombs, pValues.pi), colNames = T, rowNames = T)
  openxlsx::addWorksheet(wb, sheetName = "PVal PhIMAIN-BioGrid")
  openxlsx::writeData(wb, sheet = "PVal PhIMAIN-BioGrid", data.frame(hormCombs, pValues.bp), colNames = T, rowNames = T)
  
  openxlsx::saveWorkbook(wb, file = paste(result.dir,"/PCPs Type ",pcp.type," P-values - ",annoSet[k],".xlsx", sep = ""), overwrite = T)
  
  
}

}

# Stefan Altmann
# calculate GO enrichment for each community using GOstats package and hypergeometric test
# this script was executed using R3.6.1 and Windows 10

rm(list = ls())

setwd("~/../Desktop/Github/Community GO analysis")  

result.dir <- "Results"
if(!dir.exists(result.dir)){
  dir.create(result.dir)
}

if(!dir.exists(paste(result.dir, "/GOstats", sep = ""))){
  dir.create(paste(result.dir, "/GOstats", sep = ""))
}

library(openxlsx)

# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("GOstats")
# 
# source("https://bioconductor.org/biocLite.R")
# biocLite("org.At.tair.db")


### R code from vignette source 'GOstatsForUnsupportedOrganisms.Rnw'

###################################################
### code chunk number 1: available Schemas
###################################################
# source("https://bioconductor.org/biocLite.R")
# biocLite("AnnotationForge")
library("AnnotationForge")
# available.dbschemas()


###################################################
### code chunk number 2: Acquire annotation data
###################################################
# source("https://bioconductor.org/biocLite.R")
# biocLite("GO.db")
library(GO.db)
# library("org.Hs.eg.db")
# source("https://bioconductor.org/biocLite.R")
# biocLite("org.At.tair.db")
library("org.At.tair.db")
# frame = toTable(org.Hs.egGO)
frame = toTable(org.At.tairGO)
# remove IEP
toRemove <- which(frame$Evidence == "IEP")
frame <- frame[-toRemove,]
# make data frame
goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)
head(goframeData)


###################################################
### code chunk number 3: transformGOFrame
###################################################
# goFrame=GOFrame(goframeData,organism="Homo sapiens")
goFrame=GOFrame(goframeData,organism="Arabidopsis thaliana")
goAllFrame=GOAllFrame(goFrame)


###################################################
### code chunk number 4: Make GSC
###################################################
# source("https://bioconductor.org/biocLite.R")
# biocLite("GSEABase")
library("GSEABase")
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())


###################################################
### code chunk number 5: <make parameter
###################################################
# source("https://bioconductor.org/biocLite.R")
# biocLite("GOstats")
library("GOstats")
# universe = Lkeys(org.Hs.egGO)
# universe = Lkeys(org.At.tairGO)

# load the community data
# community.locus <- read.xlsx("../20170516 Communities using PvP complete/Edge Betweenness Step length 5/Results - AHDGen_TAIR/restricted to SS/eb - PvP (restricted to SS) - communities.xlsx", 
#                              sheetName = "Community membership - enriched", stringsAsFactors = F)
# read sheet "Community membership - enriched"
# community.locus <- read.xlsx("../20170516 Communities using PvP complete/Results 2018-04-27/Edge_Betweenness Step length 5/Results - AHD_Genetic/SS/eb - PvP (SS) - communities.xlsx",
#                              sheet = 13)

community.locus <- xlsx::read.xlsx(paste("../Community analysis/Results/",
                                   "Edge Betweenness Step length 5/Results - AHDGen_TAIR/all/eb - PvP (all) - communities.xlsx", sep=""),
                             sheetName = "Community membership - enriched", stringsAsFactors = F) 

communities <- unique(community.locus$Community)

wb <- createWorkbook()
# write.xlsx("README",file = paste(result.dir,"/GOstats/Over.xlsx", sep = ""))
addWorksheet(wb, sheetName = "README")
writeData(wb, sheet = "README", "README")


for(i in communities){
  print(i)
  # all genes
  allGenes <- factor(as.integer(community.locus$Community %in% i))
  names(allGenes) <- community.locus$Locus
  
  # exclude communities with less than three nodes
  noNodes <- length(which(allGenes == 1))
  if(noNodes < 3){
    next
  }
  
  universe <- names(allGenes) # universe = pvp
  # universe <- Lkeys(org.At.tairGO) # universe = complete arabidopsis
  genes = names(allGenes)[which(allGenes == 1)]
  params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params", 
                               geneSetCollection=gsc, 
                               geneIds = genes, 
                               universeGeneIds = universe, 
                               ontology = "BP", 
                               pvalueCutoff = 0.05, 
                               conditional = TRUE, 
                               testDirection = "over")
  
  params2 <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params", 
                               geneSetCollection=gsc, 
                               geneIds = genes, 
                               universeGeneIds = universe, 
                               ontology = "BP", 
                               pvalueCutoff = 1, 
                               conditional = TRUE, 
                               testDirection = "over")
  
  ###################################################
  ### code chunk number 6: call HyperGTest
  ###################################################
  Over <- hyperGTest(params)
  print(head(summary(Over)))
  mySummary <- summary(Over)
  
  if(nrow(mySummary) > 0){
    Over2 <- hyperGTest(params2)
    mySummary2 <- summary(Over2)
  
    PValue.adj <- p.adjust(mySummary$Pvalue, method = "BH", n = length(mySummary2$Pvalue))
  
    mySummary <- cbind(mySummary, PValue.adj)
  
  
    # write.xlsx(summary(Over),"GOstats/Over.xlsx",sheetName = paste("Com",i),append = T)
    addWorksheet(wb, sheetName = paste("Com",i))
    writeData(wb, sheet = paste("Com",i), mySummary)
    
  }

}

saveWorkbook(wb, file = paste(result.dir,"/GOstats/Over.xlsx", sep = ""), overwrite = T)

# calc similarity between GO terms
# source("https://bioconductor.org/biocLite.R")
# biocLite("GOSemSim")
# library(GOSemSim)
# library(GO.db)
# atGO <- godata('org.At.tair.db', ont = "BP")
# 
# oSum <- summary(Over)
# for(i in 1:(length(oSum$GOBPID)-1)){
#   for(j in (i+1):length(oSum$GOBPID)){
#     mySim <- goSim(oSum$GOBPID[i],oSum$GOBPID[j], atGO, measure = "Wang")
#     print(mySim)
#   }
# }



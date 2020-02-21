
library(gdata)

loadBiogrid <- function(){
  # load all biogrid interactions and extract a. thaliana interactions
  biogrid.all <- read.csv("C:/Users/Stefan Altmann/Uni/Data/BioGRID 3.4.129/BIOGRID-ALL-3.4.129.tab2/BIOGRID-ALL-3.4.129.tab2.txt",sep="\t",header=T)
  biogrid.all.ath <- biogrid.all[which(biogrid.all$Organism.Interactor.A == "3702" & biogrid.all$Organism.Interactor.B == "3702"),]
  ### unique list of interactions in biogrid
  # convert systematic names toupper characters
  biogrid.all.ath$Systematic.Name.Interactor.A <- trim(toupper(as.character(biogrid.all.ath$Systematic.Name.Interactor.A)))
  biogrid.all.ath$Systematic.Name.Interactor.B <- trim(toupper(as.character(biogrid.all.ath$Systematic.Name.Interactor.B)))
  
  # remove interactions with genetic experimental system type
  genetic.type.index <- which(biogrid.all.ath$Experimental.System.Type == "genetic")
  biogrid.all.ath <- biogrid.all.ath[-genetic.type.index,]
  # remove Protein-RNA interactions
  protein.rna.index <- which(biogrid.all.ath$Experimental.System == "Protein-RNA")
  biogrid.all.ath <- biogrid.all.ath[-protein.rna.index,]
  
  # return the interactions
  return(biogrid.all.ath)
}
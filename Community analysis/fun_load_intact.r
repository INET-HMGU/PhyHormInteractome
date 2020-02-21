
removePsi <- function(x){
  intact.binary.tmp <- c()
  for(i in 1:length(x)){
    method.split <- unlist(strsplit(as.character(x[i]),"\\("))
    intact.binary.tmp <- c(intact.binary.tmp,substr(method.split[2],1,nchar(method.split[2])-1))
  }
  return(intact.binary.tmp)  
}


loadIntact <- function(){
  ### INTACT ###########################################################################
  intact <- read.csv2("C:/Users/Stefan Altmann/Uni/Data/Intact 2015-10/intact/intact.txt",
                      sep="\t",header=T, quote="", stringsAsFactors = F)
  # copy variable for further use
  # intact.copy <- intact
  # identify interaction partners in A. thaliana
  intact.int.a.index <- grep("taxid:3702",intact$Taxid.interactor.A)
  intact.int.b.index <- grep("taxid:3702",intact$Taxid.interactor.B)
  intact.int.common <- which(intact.int.a.index %in% intact.int.b.index)
  intact.int.a.index <- intact.int.a.index[intact.int.common]
  # extract interactions between arabidopsis proteins
  intact <- intact[intact.int.a.index,]
  # remove protein complexes
  intact.complex.index <- grep("spoke",intact$Expansion.method.s.)
  intact <- intact[-intact.complex.index,]
  # check for interacting proteins - no protein-RNA/DNA interactions, keep only protein/peptide - protein/peptide interactions
  # intA.protein.index <- which(intact$Type.s..interactor.A == "psi-mi:\"MI:0326\"(protein)" | intact$Type.s..interactor.A == "psi-mi:\"MI:0327\"(peptide)")
  # intB.protein.index <- which(intact$Type.s..interactor.B == "psi-mi:\"MI:0326\"(protein)" | intact$Type.s..interactor.B == "psi-mi:\"MI:0327\"(peptide)")
  intA.protein.index <- which(intact$Type.s..interactor.A == "psi-mi:\"MI:0326\"(protein)")
  intB.protein.index <- which(intact$Type.s..interactor.B == "psi-mi:\"MI:0326\"(protein)")
  intAB.protein.index <- intA.protein.index[which(intA.protein.index %in% intB.protein.index)]
  intact <- intact[intAB.protein.index,]

  intact$Interaction.detection.method.s. <- removePsi(intact$Interaction.detection.method.s.)
  
  # keep only TAIR AGI
  require(gsubfn)
  removeIndexes <- c()
  # duplicates <- c()
#   aliasIntA <- c()
#   aliasIntB <- c()
  for(i in 1:length(intact$Alias.es..interactor.A)){
    
    # intact$Alias.es..interactor.A[i] <- unique(unlist(strapply(toupper(intact$Alias.es..interactor.A[i]), "AT[1-5,C,M]G[0-9]{5}")))
    xxx <- unique(unlist(strapply(toupper(intact$Alias.es..interactor.A[i]), "AT[1-5,C,M]G[0-9]{5}")))
    yyy <- unique(unlist(strapply(toupper(intact$Alias.es..interactor.B[i]), "AT[1-5,C,M]G[0-9]{5}")))
    zzz <- unique(unlist(strapply(intact$Publication.Identifier.s.[i], "pubmed:[0-9]*")))
    if(is.null(xxx) | is.null(yyy)){
      removeIndexes <- c(removeIndexes,i)
    }else{
      intact$Alias.es..interactor.A[i] <- paste(xxx,collapse=" ")
      intact$Alias.es..interactor.B[i] <- paste(yyy,collapse=" ")
      if(!is.null(zzz)){
        gsub("pubmed:","",zzz)
        intact$Publication.Identifier.s.[i] <- paste(zzz,collapse = " ")
      }
      
    }
  }
  
  intact <- intact[-removeIndexes,]
  
  # remove interactions annotated with two AGIs
  ncharIntA <- sapply(as.character(intact$Alias.es..interactor.A),nchar)
  ncharIntB <- sapply(as.character(intact$Alias.es..interactor.B),nchar)
  removeAmbiguous <- c(which(ncharIntA > 9),which(ncharIntB > 9))
  intact <- intact[-removeAmbiguous,]
  
  
  return(intact)
}
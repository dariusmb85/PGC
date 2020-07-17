library(data.table)
library(glmnet)
library(xtable)

GetName <- function(Fname, ToP){
  analysisName <- str_split(Fname,"-FDR.Anno")
  analysisName <- analysisName[[1]][1]
  switch (ToP,
    "tissue" ={analysisName <- str_split(analysisName,"Pheno_")
    },
    "pheno" ={analysisName <- str_split(analysisName,"Tissue_")
    })
  analysisName <- analysisName[[1]][2]
  return(analysisName)
}

AddName <- function(metaFile, ToP){
  Tissue <-GetName(metaFile,ToP)
  dat <- fread(metaFile, header = TRUE)
  dat <- cbind(Tissue, dat)
  return(dat)
}

makeTable <- function(files,ToP){
  bigDat <- AddName(files,ToP)
  myDat <- bigDat[,.(Tissue, phom, FDRphom,
                   phet, FDRphet, chr, 
                   gene_name)]
  return(myDat)
}


###########Console Script#####################
list_files <- list.files(pattern = "*Anno.txt")
myDatList <- lapply(list_files, makeTable,ToP = 'pheno')
myDat <- rbindlist(myDatList)
myDat <- myDat[order(phom)]
myDat <- myDat[phet < 0.00005]
print(xtable(myDat,digits = c(0,0,-2,-2,-2,-2,0,0)),
      floating=FALSE,latex.environments=NULL,
      include.rownames = FALSE, booktabs = TRUE)

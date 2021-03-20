library(data.table)
library(qqman)
#library(dplry)

make_manhattan <- function(anno_pval){
  lFiles <- list()
  lFiles <- list.files(pattern = "_DISCARD-FDR.Anno.txt$", 
                            full.names = TRUE, recursive = TRUE)
  datlist <- list()
  datlist <- lapply(lFiles, function(x) fread(x, header = TRUE, 
                                              data.table=TRUE))
  plot_manhattan(datlist,anno_pval)
}

plot_manhattan <- function(dlist,apval){
  for(i in 1:length(dlist)){
    manDat <- dlist[[i]]
    manDat <- manDat[complete.cases(manDat),]
    tp <- manDat$Bonfphet[[1]]
    otp <- manDat$Bonfphom[[1]]
    if(tp < otp){
      pval <- 'phet'
    }else{
      pval <- 'phom'
    }
    manhattan(manDat,chr="chr",bp='bp',p = pval,
              snp='gene_name', annotatePval = apval)
  }
}
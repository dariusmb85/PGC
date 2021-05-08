library(data.table)
library(glmnet)
library(xtable)


###########Console Script#####################
list_files <- list.files(pattern = "*Anno.txt")
myDat <- fread(list_files)
myDat <- myDat[order(Bonfphom)]
myDat2 <- myDat[Bonfphom < 0.05]
myDat2$Tissue <- gsub("^_","",myDat2$Tissue)
myDat2$rs <- NULL
myDat2$bp <- NULL
print(xtable(myDat2,digits = c(0,0,-2,-2,-2,-2,-2,-2,0,0)),
      floating=FALSE,latex.environments=NULL,
      include.rownames = FALSE, booktabs = TRUE)

library(data.table)
library(glmnet)
library(xtable)
library(tidyverse)


###########Console Script#####################
list_files <- list.files(pattern = "*Anno.txt")
myDat <- fread(list_files[2])


myDat <- myDat[order(Bonfphet)]
myDat2 <- myDat[Bonfphet < 0.05]
if(dim(myDat2)[1] < 10){
  myDat2 <- myDat[1:10]
}
myDat2$Tissue <- gsub("^_","",myDat2$Tissue)
# myDat2$rs <- NULL
myDat2$bp <- NULL
print(xtable(myDat2,digits = c(0,0,0,-2,-2,-2,-2,-2,-2,0,0)),
      floating=FALSE,latex.environments=NULL,
      include.rownames = FALSE, booktabs = TRUE)

myDat <- myDat[order(Bonfphom)]
myDat3 <- myDat[Bonfphom < 0.05]
if(dim(myDat3)[1] < 10){
  myDat3 <- myDat[1:10]
}
myDat3$Tissue <- gsub("^_","",myDat3$Tissue)
# myDat3$rs <- NULL
myDat3$bp <- NULL
print(xtable(myDat3,digits = c(0,0,0,-2,-2,-2,-2,-2,-2,0,0)),
      floating=FALSE,latex.environments=NULL,
      include.rownames = FALSE, booktabs = TRUE)

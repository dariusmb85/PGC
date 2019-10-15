#' @name predXUtils
#' @title Download the PrediXcan database
#' @param downloadDest character of length 1, the download directory
#' @param rmTar logical of length 1, Should orginal tar file be removed?
#' @param predDB character of length 1, file path to a predXcan .db file

NULL

#' @describeIn predXUtils Download predXcan to the given directory
#' @export

pgcDownloadPrediXcan <- function(downloadDest = "PrediXcan", rmTar = TRUE) {
  
  if (!dir.exists(downloadDest)) {
    stop("The given download destination does not exist")
  }
  predUrl <- "https://s3.amazonaws.com/predictdb2/deprecated"
  f <- "GTEx-V6p-HapMap-2016-09-08.tar.gz"
  tf <- file.path(downloadDest, f)
  download.file(paste(predUrl, f, sep = "/"), destfile = tf)
  untar(tf)
  if (rmTar){
    file.remove(tf)
  }
  sub(".tar.gz$", '', tf)
}

#' @describeIn predXUtils Combine tissues weights into single data.table object
#' @param predDir character of length 1, directory where db files have been stored
#' @import data.table
#' @export

pgcCombinePredxFiles <- function(predDir, SubsetByTissue = FALSE) {
  
  #creats list of all database filenames
  fnames <- list.files(path = predDir, pattern = ".db$", full.names = TRUE)
  #creates empty list that we will add to
  datList <- lapply(fnames, pgcGetPredxWeights)
  dat <- rbindlist(datList)
  saveRDS(dat,"combinedWeightFiles.rds")
  dat[]
  
  if (SubsetByTissue == TRUE){
    tl <- pgcListTissues(dat)
    print(tl)
    print("Type Tissue name to subset new weights file\n")
    tiss <- readline("ex.'Brain_Cortex' for Brain Cortex or 'Brain' for every Brain tissue:")
    dat <- pgcSubsetTissues(tiss)
    saveRDS(dat,"combinedWeightFiles.rds")
  }
  dat[]
}

#' @describeIn predXUtils Returns list of tissues involved in all prediXcan .db
#' @param combWts character of length 1, directory where db files have been stored
#' @import data.table
#' @export


pgcListTissues <- function(combWts){
  
  readline("Here is s alist of tissues: Press any key")
  tissList <- unique(as.data.table(combWts)[[6]])
  tissList
  
}

#' @describeIn predXUtils Reads weight table from the given predXcan .db file
#' @importFrom DBI dbConnect dbSendQuery dbClearResult dbDisconnect
#' @importFrom RSQLite SQLite
#' @import data.table
#' @export

pgcGetPredxWeights <- function(predDB) {
  
  con <- dbConnect(SQLite(), dbname = predDB) # connects to database
  # grabs data from weights table
  res <- dbSendQuery(con, "SELECT * FROM weights;") 
  dat <- as.data.table(dbFetch(res, n = -1))
  dat[ , tissue := gsub("TW_|_0.5.db", "", basename(predDB))]
  dbClearResult(res) # empties our return variable
  dbDisconnect(con) # closes our db connection
  dat[]
  
}

#' @describeIn predXUtils Subsets combined weight table for a specific tissue or tissues
#' @param tis character of length 1, Tissue you want to subset combined weights file by
#' @import plyr data.table
#' @export

pgcSubsetTissues <- function(tis) {

  combWts <- readRDS("combinedWeightFiles.rds")
  dat <- combWts[grep(tis, combWts$tissue), ]
  dat[]
  
}

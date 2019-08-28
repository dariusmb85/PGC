#' @name predXUtils
#' @title Download the PrediXcan database
#' @param predDir character of length 1, the download directory
#' @param predDB character of length 1, file path to a predXcan .db file

NULL

#' @describeIn predXUtils Download predXcan to the given directory
#' @export

pgcDownloadPrediXcan <- function(predDir = "PrediXcan", overwrite = FALSE) {
  if (!dir.exists(predDir)) {
    dir.create(predDir)
  } else if (!overwrite) {
    unlink(predDir, recursive = TRUE, force = TRUE)
    dir.create(predDir)
  } else {
    stop("")
  }
  predUrl <- "https://s3.amazonaws.com/predictdb2"
  u <- "deprecated/GTEx-V6p-HapMap-2016-09-08.tar.gz"
  download.file(paste(predUrl, u, sep = "/"), destfile = predDir)
}

#' @describeIn predXUtils Combine tissues weights into single data.table object
#' @import data.table
#' @export

pgcCombinePredxFiles <- function(predDir) {
  
  #creats list of all database filenames
  fnames <- list.files(path = predDir, pattern = ".db$", full.names = TRUE)
  #creates empty list that we will add to
  datList <- lapply(fnames, pgcGetPredxWeights)
  dat <- rbindlist(datList)
  dat[]
  
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

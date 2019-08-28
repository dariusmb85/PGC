#' @name pgcCombinePredxFiles
#' @title Combine Predixcan database files
#' @param destDir character of length 1, the download directory
#' @details
#' \code{pgcCombinePredxFiles} relies on the directory containing... 
#' @importFrom DBI dbConnect dbSendQuery dbClearResult dbDisconnect
#' @importFrom RSQLite RSQLite
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
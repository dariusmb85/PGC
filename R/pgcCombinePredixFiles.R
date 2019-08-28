#' @name pgcCombinePredxFiles
#' @title Combine Predixcan database files
#' @param destDir character of length 1, the download directory
#' @import RSQLite
#' @export

pgcCombinePredxFiles <- function{
  #creats list of all database filenames
  fnames <- list.files(pattern = ".db$", full.names = TRUE, recursive = TRUE)
  #creates empty list that we will add to
  datlist <- list()
  #loops through filenames
  for (f in fnames) {
    
    tmp <- dbConnect(SQLite(), dbname = f)#connects to database
    res <- dbSendQuery(tmp, "SELECT * FROM weights;")#grabs data from weights table
    dat <- dbFetch(res,n = -1)#returns data and places into array dat
    dbClearResult(res)#empties our return variable
    dbDisconnect(tmp)#closes our db connection
    dat$tissue <- basename(dirname(f))#creates a column for tissue by directory name
    datlist[[basename(f)]] <- dat#adds iteration of data to our list
    
  }
  dat <- do.call(rbind, datlist)#combines all list in datlist into one data frame
  return(dat)
}

#' @name pgcDownloadPrediXcan
#' @title Download the PrediXcan database
#' @param destDir character of length 1, the download directory
#' @export

pgcDownloadPrediXcan <- function(destDir = "PrediXcan", overwrite = FALSE) {
  if (!dir.exists(destDir)) {
    dir.create(destDir)
  } else if (!overwrite) {
    unlink(destDir, recursive = TRUE, force = TRUE)
    dir.create(destDir)
  } else {
    stop("")
  }
  predUrl <- "https://s3.amazonaws.com/predictdb2"
  u <- "deprecated/GTEx-V6p-HapMap-2016-09-08.tar.gz"
  download.file(paste(predUrl, u, sep = "/"), destfile = destDir)
}
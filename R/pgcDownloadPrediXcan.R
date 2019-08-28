#' @name pgcDownloadPrediXcan
#' @title Download the PrediXcan database
#' @param destDir character of length 1, the download directory
#' @export

pgcDownloadPrediXcan <- function(destDir = "PrediXcan") {
  dir.exists(destDir) 
  dir.create(destDir)
  download.file("https://s3.amazonaws.com/predictdb2/deprecated/GTEx-V6p-HapMap-2016-09-08.tar.gz", destfile = destDir)
}
#' @name pgcDownloadPrediXcan
#' @title Download the PrediXcan database
#' @param destDir character of length 1, the download directory
#' @export

pgcDownloadPrediXcan <- function(destDir = "PrediXcan") {
  dir.exists(destDir)  # Changes example
  dir.create(destDir)
  download.file(, destfile = destDir)
}
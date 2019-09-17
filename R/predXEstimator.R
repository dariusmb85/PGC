#' @name predXEstimateExpr
#' @title Expression Estimator
#' @param downloadDest character of length 1, the download directory
#' @param rmTar logical of length 1, Should orginal tar file be removed?
#' @param predDB character of length 1, file path to a predXcan .db file


# library(data.table)
# getUcscQuery <- function(qstr, db = "hg38") {
#   if(!curl::has_internet()) stop("Must be connected to the internet.")
#   ucsc <- dbConnect(MySQL(), 
#                     username = "genome", 
#                     dbname = db,
#                     host = "genome-mysql.soe.ucsc.edu",
#                     port = 3306,
#                     password = "")
#   dat <- suppressWarnings(dbGetQuery(ucsc, qstr))
#   dat <- as.data.table(dat)
#   dbDisconnect(ucsc)
#   warning("UCSC uses 0-based start and 1-based end positions.", call. = FALSE)
#   dat[]
# }
# flds <- c("name", "chrom", "chromStart", "chromEnd", 
#           "refUCSC","observed", "class")
# qs <- sprintf("SELECT %s FROM snp151Common", paste(flds, collapse = ", "))
# snp <- getUcscQuery(qs)



#########################################

# combDat <- readRDS("combinedWeightFiles.rds")
# combDat[grep("^Brain", combDat$tissue), ]
# hdr <- names(comDat)
# 
# avsnp150 <- unique(as.data.frame(combDat)[[1]])
# avsnp150 <- as.data.table(avsnp150)
# 
# 
# m <- list()
# 
# names(combDat)
# rs <- combDat[[1]][1]
# trans <- combDat[[2]][1]
# coef <- combDat[[3]][1]
# efallele <- combDat[[5]][1]
# tiss <- combDat[[6]][1]

#make forest plots by running plottingCode.R

source('~/GitHub/PGC/inst/plottingCode.R')
fls <- list.files(pattern = "*datObj",recursive = FALSE)
l <- lapply(fls, readRDS)
names(l) <- sub(".datObj", "", basename(fls))
pdf("AcrossTissue_FP_LOAD_Supplemental.pdf")
l <- lapply(l, fixInf)
lapply(l,pltCPA)
dev.off()



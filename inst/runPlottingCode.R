#make forest plots by running plottingCode.R

source('~/GitHub/PGC/inst/plottingCode.R')
fls <- list.files(pattern = "*datObj",recursive = TRUE)
l <- lapply(fls, readRDS)
names(l) <- sub(".datObj", "", basename(fls))
pdf("ForestPlot_Supplemental.pdf")
lapply(l,pltCPA)
dev.off()



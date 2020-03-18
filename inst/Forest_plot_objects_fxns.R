################################################################################
# Outputs: Build forest plotname nomenclature including the input metadata and

# Finetune how tissue names are presented to the Forest plotter
pruneNameForPlot <- function(itissue){
  if (itissue %in% "Anterior_cingulate_cortex_BA24") return ("Anterior_cingulate")
  if (itissue %in% "Nucleus_accumbens_basal_ganglia") return ("Nucleus_accum")
  if (itissue %in% "Putamen_basal_ganglia") return ("Putamen")
  if (itissue %in% "Frontal_Cortex_BA9") return ("Frontal_Cortex")
  if (itissue %in% "Cerebellar_Hemisphere") return ("Cerebellar_Hemi")
  if (itissue %in% "Caudate_basal_ganglia") return ("Caudate")
  return (itissue)
}

###################################
#Needed because of the epacts runs
buildSNPfiles <- function(phe){
  switch(phe,
         "Marijuana_Strong_Desire"={
           snpF <- paste0("/projects/sequence_analysis/vol3/UCSFplink/epacts/", phe,
                          "_filtered/MarajuanaStrDsr.all.epacts")
         },
         "MJ"={
           snpF <- paste0("/projects/sequence_analysis/vol3/UCSFplink/epacts/", phe,
                          "_filtered/MJ_DX_DSM4.all.epacts")
         },
         "NIC"={
           snpF <- paste0("/projects/sequence_analysis/vol3/UCSFplink/epacts/", phe,
                          "_filtered/NIC_DX_DSM4.all.epacts")
         },
         {
           snpF <- paste0("/projects/sequence_analysis/vol3/UCSFplink/epacts/",phe,
                          "_filtered/",phe,".all.epacts")
         })
}

###############################
# annoSNProws <- function(snpR,avsnp){
#   setkey(avsnp, Chr, BEG, END)
#   setkey(snpR, Chr, BEG, END)
#   snpRanno<- snpR[avsnp]
#   return(snpRanno)
# }

################################################################################
# Build a MAP from GTEX variants to rsIDS. Have variants as INDEX. Use None for 
# missing rsIDs.
# Independant of the tissues
buildVariantToRSIDmap <- function() {
  filepath <- file.path("","projects","sequence_analysis","vol3","predix_Scan",
                        "GTEx-V6p_flowOver","gtex_data","")
  filename <- paste(filepath,
                    "GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt",
                    sep="")
  data <- fread(filename, stringsAsFactors = FALSE, h = TRUE)
  setkey(data, rs_id_dbSNP147_GRCh37p13)
  mapdata <- data[!("."), .(variant_id, rs_id_dbSNP147_GRCh37p13)] 
  #mapdata[mapdata == "."] <- NA
  #mapdata <- na.omit(mapdata)
  # Takes a long time to update rownames(densemapdata) = densemapdata[,1]
  mapdata[]
}

convertVariantToRSID <- function(isnp, densemapdata) {
  #variant <- mapdata[mapdata$variant_id==isnp, 1][1]
  setkey(densemapdata, variant_id)
  rsid <- as.data.table(densemapdata[J(isnp), rs_id_dbSNP147_GRCh37p13])
  # if (rsid == "character(0)"){
  #   rsid <- ""
  # }
  return (rsid)
}

################################################################################
convertRSIDToVariant <- function(ivar, densemapdata) {
  #variant <- mapdata[mapdata$variant_id==isnp, 1][1]
  setkey(densemapdata, rs_id_dbSNP147_GRCh37p13)
  var <- toString(densemapdata[J(ivar), variant_id])
  if (var == "character(0)"){
    var <- ""
  }
  return (var)
}

################################################################################ 
# Find tissue-specific GTEX betas and se.
buildGTEXfilenames <- function(tissue) {
  tiss <- tissue
  if(tiss == "Frontal_Cortex"){
    tiss <- "Frontal_Cortex_BA9"
  }
  
  basedir <- '/projects/sequence_analysis/vol3/predix_Scan/GTEx-V6p_flowOver/gtex_data'
  filename <- paste('Brain', paste(tiss, 'allpairs.txt', sep= '.'), sep= "_")
  return (paste(basedir, filename, sep= '/'))
}
################################################################################

pred_df <- function(tis, coh, basedir){
  filename <- file.path(basedir, coh, paste0(tis,
                                             "_Univariate_Covariates_Output_", coh,
                                             ".assoc.txt"))
  trans <- fread(filename)
  setkey(trans, rs)
  row <- trans[J(transcript), .(beta, se, p_score)]
  df_row <- data.table(tis,row)
  setnames(df_row,c('Tissue', 'weight', 'se', 'pval')) 
  df_row[]
}

################################################################################

pred_gtex <- function(tis, trans, dmd, snps){
  filename <- buildGTEXfilenames(tis)
  gtexdata <- fread(filename, h=TRUE, stringsAsFactors=FALSE)
  setkey(gtexdata,gene_id)
  gtex <- gtexdata[J(trans), .(gene_id, variant_id, slope,
                               slope_se, pval_nominal)]
  gtex <- data.table(tis, gtex)
  rm(gtexdata)
  setkey(gtex, variant_id)
  setkey(dmd, variant_id)
  gtexNew <- dmd[gtex, .(tis, rs_id_dbSNP147_GRCh37p13,
                         slope, slope_se, pval_nominal)]
  setkey(gtexNew, rs_id_dbSNP147_GRCh37p13)
  gtexNew <- gtexNew[snps]
  gtexNew[]
}

#################################################################################

AlzCohortBuilder <- function(coh,phe){
  newCoh <- c()
  for(i in 1:length(coh)){
    for(j in 1:length(phe)){
      cohorts <- paste0(coh[i],"_",phe[j])
      newCoh <- c(newCoh,cohorts)
    }#Commit
  }
  return (newCoh)
}
################################################################################

################################################################################
library('data.table')
source('/home/dariusmb/GitHub/PGC/inst/Forest_plot_objects_fxns.R')
source('/home/dariusmb/GitHub/PGC/R/allClasses.R')
basedir <- file.path("","projects","sequence_analysis","vol5",
                     "dariusmb","PGC_v8","output")
# basedir <- file.path("","projects","sequence_analysis","vol3",
#                      "predix_Scan","GTEx-V6p_flowOver","GEMMApipe_3"
#                      ,"output")

cohorts <- c("Alc_Dep","Family_Alcoholism","Marijuana_Strong_Desire","MJ","NIC")
#cohorts <- c("ADGC","ADNI","GENADA","MERGE","NIALOAD","ROSMAP","TARCC","TGEN2","WASHU")

#pheno <- c("AAO","LOAD")

# if(cohorts[1] == "ADGC"){
#   cohorts <- AlzCohortBuilder(cohorts,pheno)
# }
############v6p
# tissues <- c('Anterior_cingulate_cortex_BA24','Cerebellar_Hemisphere','Cortex',
#           'Hippocampus','Nucleus_accumbens_basal_ganglia','Caudate_basal_ganglia'
#           ,'Cerebellum','Frontal_Cortex_BA9','Hypothalamus','Putamen_basal_ganglia')
############v8
tissues <- c('Anterior_cingulate_cortex_BA24','Caudate_basal_ganglia',
             'Cerebellar_Hemisphere','Cerebellum','Cortex','Frontal_Cortex_BA9',
             'Hippocampus','Hypothalamus','Nucleus_accumbens_basal_ganglia',
             'Putambasal_ganglia')

# tissues <- c('Amygdala','Anterior_cingulate_cortex_BA24','Cerebellar_Hemisphere','Cortex',
#           'Hippocampus','Nucleus_accumbens_basal_ganglia','Caudate_basal_ganglia'
#           ,'Cerebellum','Frontal_Cortex_BA9','Hypothalamus','Putamen_basal_ganglia',
#           'Substantia_nigra')


combBrain <- file.path("","projects","sequence_analysis","vol5","dariusmb",
                       "PGC_v8","CombinedWeightsBrain.csv")

# combBrain <- file.path("","projects","sequence_analysis","vol3","predix_Scan",
#                       "GTEx-V6p_flowOver","GEMMApipe_3","NewCombinedDBEdit.csv")

################################################################################

options(echo = FALSE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
transcriptIn <- args[1]
transNumID   <- args[2]
resp         <- args[3]
analysis     <- args[4]
CorT         <- args[5]
transcriptIn <- eval(parse(text = args[1]))
transNumID   <- eval(parse(text = args[2]))
resp         <- eval(parse(text = args[3]))
analysis     <- eval(parse(text = args[4]))
CorT         <- eval(parse(text = args[5]))

# print(args)
#print(transcriptIn)
#print(transcriptIn[transNumID])

#print(transcript[transNumID])
transcriptIn <- transcriptIn[transNumID]
respIn <- resp[transNumID]

##Testing
transcriptIn <- c('ENSG00000141127')
translength <- length(transcriptIn)

###Testing X Tissues
tissues <- tissues[1:2]
respIn <- 2
#respIn <- 9
cohorts <- cohorts[respIn]

##Testing X Cohorts
# transcriptIn <- c('ENSG00000127824.9')
# cohorts <- cohorts[1:3]
# respIn  <- 9
# tissues <- tissues[respIn]
# CorT    <- "C"

# Now we need to fetch data from disparate sources and build a final dataframe of the form
#    cpassoc Shom-Shet phom-phet
# expression beta se
#        snp beta se
#       gwas beta se
# Add separate data table for subplot construction of the forest plots

#####FIX buildVariantToRSIDmap() too slow ########
#densemapdata <- buildVariantToRSIDmap()
filepath <- file.path("","projects","sequence_analysis","vol3",
                      "predix_Scan","GTEx-V6p_flowOver","gtex_data",
                      "")
gtexFile <- paste(filepath,
                  "GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt",
                  sep="")

df_cpassoc <- data.table()
df_gtex <- data.table()
df_gwas <- data.table()
df_predixcan <- data.table()
df_PrdxWts <- data.table()

print(paste('Processing next cohort', cohorts, sep="="))

transcript <- transcriptIn
  
# First
# grab the CPASSOC Across cohort result for the current transcript 
# SHow contribution of COHORT to the overall CPASSOC

switch (CorT,
        "T"={
          #Across Tissue
          print('Process xTissues CPASSOC')
          cpFile <- file.path("MetaAnalysis_AcrossTissue",
                              "FDR_Anno",paste0("MetaxTissue-FDR.Anno.txt"))
        },
        "C"={
          #Across Cohort
          print('Process xPheno CPASSOC')
          cpFile <- file.path("MetaAnalysis_AcrossPheno",
                              "FDR_Anno",paste0("MetaxPheno-FDR.Anno.txt"))
          print(cpFile)
        }
)

# cpFile <- file.path("Across_Tissue_AAO","FDR",
#                     "FDR_Anno",paste0("MetaxTissue_",
#                                       cohorts, "-FDR.Anno.txt"))

ccData <- fread(cpFile, stringsAsFactors = FALSE, h=TRUE)
setkey(ccData, rs)
transcriptName <- ccData[transcript, gene_name]
transcriptName <- transcriptName[1]
if ( is.na(transcriptName) ){
  transcriptName <- transcript
}
df_cpassoc <- ccData[transcript, .(gene_name, phom,
                                   phet)]
df_cpassoc <- df_cpassoc[phom == min(df_cpassoc$phom)]
setnames(df_cpassoc, c('Gene Name', 'sHomPval',
                       'sHetPval'))
rm(ccData)
print(df_cpassoc)

# Second
# grab the individual cohort gene-expression beta/se 
# info for all tissues:
# GEMMA associations

print("Process Tissues")
print(tissues)
print(cohorts)
print(transcript)

switch(CorT,
       "T"={
         df_predixcan <- lapply(tissues, pred_df,coh = cohorts,
                                basedir = basedir, CorT = CorT)
       },
       "C"={
         df_predixcan <- lapply(cohorts, pred_df, tis = tissues, 
                                basedir = basedir, CorT = CorT)
       }
)

# #Across Tissues
# df_predixcan <- lapply(tissues, pred_df,coh = cohorts,
#                        basedir = basedir, CorT = CorT)
# #Across Cohorts
# df_predixcan <- lapply(cohorts, pred_df, tis = tissues, 
#                        basedir = basedir, CorT = CorT)

df_predixcan <- rbindlist(df_predixcan)
setnames(df_predixcan, c('Tissue','weight','se','pval'))
print(df_predixcan)

# For each SNP of this transcript. Grab all the known 
# beta/se from the SNP-phenotype association
# Not sure if all these SNPs were used in all the 
# tissue calculations - probably not

# Third phenotype ~ SNPs
# Find all the unique SNPs accross all tissues. Then find
# the Phenotype_SNPs GWAS data and build a data frame
# Keep the object SNPlist for the GTEX data object as well

setwd(file.path("","projects","sequence_analysis","vol3",
      "UCSFplink","avsnp150"))
file_list <-list.files(pattern= "*hg19_multianno.txt$",
                       full.names = TRUE)
files <- lapply(file_list, function(x) fread(x, header = TRUE, 
                                             data.table=TRUE))
annoVarsnp <- rbindlist(files)
rm(files)
setnames(annoVarsnp,old = c("Start", "End"), 
         new = c("BEG", "END"))
setwd(basedir)

switch (CorT,
  "T"={
    SNPfiles <- buildSNPfiles(phe = cohorts)
    SNProws <- fread(SNPfiles)
  },
  "C"={
    #Across Cohorts
    SNPfiles <- lapply(cohorts, buildSNPfiles)
    SNProws <- lapply(SNPfiles,fread)
    SNProws <- rbindlist(SNProws)
  }
)


rm(SNPfiles)
setnames(SNProws, old = "#CHROM", new = "Chr")
setkey(annoVarsnp, Chr, BEG, END)
setkey(SNProws, Chr, BEG, END)
SNProws <- SNProws[annoVarsnp]
rm(annoVarsnp)

#dat[ , tissue := gsub("TW_|_0.5.db", "", basename(predDB))]

fname <- combBrain
weights <- fread(fname)
weights$gene <- gsub("[.].*", "", weights$gene)

SNPlist <- unique(weights[weights[, gene] == transcript,
                          rsid])
setkey(SNProws, avsnp150)
rowbetas <- SNProws[J(SNPlist), 
                    .(avsnp150, BETA, SEBETA , PVALUE)]
df_PrdxWts <- weights[weights[,gene] == transcript,
                      .(tissue,rsid, weight)]
rm(SNProws)
rm(weights)
print('Process SNP data')
df_gwas <- rowbetas[, .(avsnp150, BETA,
                        SEBETA, PVALUE)]
rm(rowbetas)
pred <- cohorts
df_gwas <- cbind(df_gwas, pred)
setnames(df_gwas, c('SNP', 'weight', 'se',
                    'pval', 'predictor'))
rownames(df_gwas) <- NULL
print(head(df_gwas))

# Fourth
# For the items in SNPlist, for each tissue, grab the 
# GTEX raw beta,se terms
print("Process GTEX Tissues")
print(cohorts)
print(tissues)
print(transcript)

gtexDat <- fread(gtexFile, stringsAsFactors = FALSE, h = TRUE)
setkey(gtexDat, rs_id_dbSNP147_GRCh37p13)
densemapdata <- gtexDat[!("."), .(variant_id,
                                  rs_id_dbSNP147_GRCh37p13)] 
rm(gtexDat)

gtex <- lapply(tissues, pred_gtex, trans = transcript,
               dmd = densemapdata, snps = SNPlist)
df_gtex <- rbindlist(gtex)

setwd("forest_plot_dataFiles")
setnames(df_gtex, c("Tissue", 'SNP', 'weight',
                    'se', 'pval'))

setnames(df_PrdxWts, c("Tissue","SNP","weight"))
setkey(df_gtex, Tissue, SNP)
setkey(df_PrdxWts, Tissue, SNP)

df_PrdxWts <- df_PrdxWts[df_gtex, .(Tissue, SNP, weight)]

# df_cpassoc
# df_predixcan
# head(df_gwas)
# head(df_gtex)
# 2496260
switch (CorT,
  "C" ={
    datObj <- cpRes(transcript, pred = tissues, "Across_Cohort", transcriptName, 
                    df_cpassoc, df_predixcan, df_gwas, df_gtex, df_PrdxWts)
    saveRDS(datObj, file = paste(transcript, tissues, transcriptName, 
                                 analysis,".datObj", sep = "_"))
  },
  "T"={
    datObj <- cpRes(transcript, cohorts, "Across_Tissue", transcriptName, 
                    df_cpassoc, df_predixcan, df_gwas, df_gtex, df_PrdxWts)
    saveRDS(datObj, file = paste(transcript, cohorts, transcriptName, 
                                 analysis,".datObj", sep = "_"))
  }
)
# datObj <- cpRes(transcript, cohorts, "Across_Tissue", transcriptName, 
#                 df_cpassoc, df_predixcan, df_gwas, df_gtex, df_PrdxWts)
# saveRDS(datObj, file = paste(transcript, cohorts, transcriptName, 
#                              analysis,".datObj", sep = "_"))
setwd(basedir)

#print('finished')
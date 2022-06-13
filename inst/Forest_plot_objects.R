################################################################################
library('data.table')
source('~/GitHub/PGC/inst/Forest_plot_objects_fxns.R')
source('~/GitHub/PGC/R/allClasses.R')
library('tibble')
library('stringr')
library('dplyr')

################################################################################

options(echo = FALSE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
transcriptIn <- args[1]
transNumID   <- args[2]
resp         <- args[3]
analysis     <- args[4]
CorT         <- args[5]
trait        <- args[6]
AlzPheno     <- args[7]
transcriptIn <- eval(parse(text = args[1]))
transNumID   <- eval(parse(text = args[2]))
resp         <- eval(parse(text = args[3]))
analysis     <- eval(parse(text = args[4]))
CorT         <- eval(parse(text = args[5]))
trait        <- eval(parse(text = args[6]))
AlzPheno     <- eval(parse(text = args[7]))
print(args)

#print(transcriptIn)
#print(transcriptIn[transNumID])

#Testing

###Testing X Tissues
# transcriptIn <- c("ENSG00000127364")
# translength <- length(transcriptIn)
# tissues <- tissues[2:6]
# respIn <- 5
# cohorts <- cohorts[respIn]
# AlzPheno <- "AAO"
# trait <- "Alz"
# analysis <- "Across_Tissue"
# CorT <- "T"

##Testing X Cohorts
# transcriptIn <- c('ENSG00000104936','ENSG00000159905','ENSG00000100243','ENSG00000142252')
# translength <- length(transcriptIn)
# 
# transNumID<- 3
# transcript <- transcriptIn[transNumID]
# # cohorts <- cohorts[1:4]
# resp    <- c(8,6,8,10)
# CorT    <- "C"
# trait   <- "Alz"
# AlzPheno<- "AAO"
# 


################ Start of Script ######################

basedir <- file.path("","projects","sequence_analysis","vol5",
                     "dariusmb","PGC_AD","output")

# cohorts <- c("Alc_Dep","Family_Alcoholism","Marijuana_Strong_Desire","MJ","NIC")
cohorts <- c("ADGC","ADNI","GENADA","NIALOAD",
             "ROSMAP","TARCC","WASHU")

pheno <- c("AAO","LOAD")

if(cohorts[1] == "ADGC"){
  cohorts <- AlzCohortBuilder(cohorts,pheno)
}

############v8_All_tissues
# tissues <- c('Amygdala','Anterior_cingulate_cortex_BA24','Cerebellar_Hemisphere','Cortex',
#           'Hippocampus','Nucleus_accumbens_basal_ganglia','Caudate_basal_ganglia'
#           ,'Cerebellum','Frontal_Cortex_BA9','Hypothalamus','Putamen_basal_ganglia',
#           'Substantia_nigra')

##############v8_All_Tissues_AD
tissues <- c('Amygdala','Anterior_cingulate_cortex_BA24','Caudate_basal_ganglia',
             'Cerebellar_Hemisphere','Cerebellum','Cortex','Frontal_Cortex_BA9',
             'Hippocampus','Hypothalamus','Nucleus_accumbens_basal_ganglia',
             'Putambasal_ganglia','Spinal_cord_cervical_c-1','Substantia_nigra')


# Addiction Analysis
# combBrain <- file.path("","projects","sequence_analysis","vol5","dariusmb",
#                        "PGC_v8","CombinedWeightsBrain.csv")

# Alz Analysis
combBrain <- file.path("","projects","sequence_analysis","vol5","dariusmb",
                       "PGC_AD","CombinedWeightsBrain.csv")


# Now we need to fetch data from disparate sources and build a final dataframe of the form
#    cpassoc Shom-Shet phom-phet
# expression beta se
#        snp beta se
#       gwas beta se
# Add separate data table for subplot construction of the forest plots

filepath <- file.path("","projects","sequence_analysis","vol5",
                      "dariusmb","gtex_data",
                      "")
gtexFile <- paste(filepath,
                  "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt",
                  sep="")

df_cpassoc <- data.table()
df_gtex <- data.table()
df_gwas <- data.table()
df_predixcan <- data.table()
df_PrdxWts <- data.table()

print(paste('Processing next cohort', cohorts, sep="="))

transcriptIn <- transcriptIn[transNumID]
respIn <- resp[transNumID]

transcript <- transcriptIn

#TEST COMMIT

# First
# grab the CPASSOC Across cohort result for the current transcript 
# SHow contribution of COHORT to the overall CPASSOC
cpFile <- get_Cpassoc_Data(trait = trait, cort = CorT, pheno = AlzPheno)

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
print(basedir)

if(AlzPheno == "AAO"){
  switch(CorT,
         "T"={
           cohorts <- cohorts[respIn]
           df_predixcan <- lapply(tissues, pred_df,coh = cohorts,
                                  basedir = basedir, CorT = CorT)
         },
         "C"={
           tissues <- tissues[respIn]
           df_predixcan <- lapply(cohorts, pred_df, tis = tissues, 
                                  basedir = basedir, CorT = CorT)
         }
  )
  df_predixcan <- rbindlist(df_predixcan)
}else if(AlzPheno == "LOAD"){
  switch(CorT,
         "T"={
           cohorts <- cohorts[respIn]
           df_predixcan <- lapply(tissues, pred_df,coh = cohorts,
                                  basedir = basedir, CorT = CorT)
         },
         "C"={
           tissues <- tissues[respIn]
           df_predixcan <- lapply(cohorts, pred_df, tis = tissues, 
                                  basedir = basedir, CorT = CorT)
         }
  )
  df_predixcan <- rbindlist(df_predixcan)
  df_predixcan <- df_predixcan[Tissue %like% "_LOAD"]
}


print(df_predixcan)

# For each SNP of this transcript. Grab all the known 
# beta/se from the SNP-phenotype association
# Not sure if all these SNPs were used in all the 
# tissue calculations - probably not

# Third phenotype ~ SNPs
# Find all the unique SNPs accross all tissues. Then find
# the Phenotype_SNPs GWAS data and build a data frame
# Keep the object SNPlist for the GTEX data object as well

# Chr	Start	End	Ref	Alt	avsnp150  ##avnsp150 original header

annoVarsnp <- annoVarFile(trait,cohorts) %>% 
  unique(by="avsnp150")
setwd(basedir)

# switch (CorT,         ###########Addiction Analysis
#   "T"={
#     SNPfiles <- buildSNPfiles(phe = cohorts)
#     SNProws <- fread(SNPfiles)
#   },
#   "C"={
#     #Across Cohorts
#     SNPfiles <- lapply(cohorts, buildSNPfiles)
#     SNProws <- lapply(SNPfiles,fread)
#     SNProws <- rbindlist(SNProws)
#   }
# )

##################Alz Analysis
switch (CorT,
        "T"={
          SNPfiles <- buildSNPfiles_AD(phe = cohorts)
          SNProws <- fread(SNPfiles)
        },
        "C"={
          #Across Cohorts
          SNProws_aao <- list()
          SNProws_load <- list()
          SNPfiles <- lapply(cohorts, buildSNPfiles_AD)
          SNProws <- lapply(SNPfiles,fread)
          # SNProws <- rbindlist(SNProws)
          for(i in 1:length(SNProws)){
            if(i %% 2 == 1){
              SNProws_aao <- cbind(SNProws_aao, SNProws[i])
            }else if(i %% 2 == 0){
              SNProws_load <- cbind(SNProws_load, SNProws[i])
            }
          }
          SNProws_aao <- rbindlist(SNProws_aao)
          SNProws_load <- rbindlist(SNProws_load)
        }
        
)
if(AlzPheno == "LOAD"){
  SNProws <- SNProws_load[SNP %like% "rs"]
  setnames(SNProws, old = c("CHR","BP"), new = c("Chr","BEG"))
  SNProws <- mutate(SNProws,"SE" = abs(OR/ qnorm(P/2)))
}else if(AlzPheno == "AAO"){
  SNProws <- SNProws[SNP %like% "rs"]
  setnames(SNProws, old = c("CHR","BP"), new = c("Chr","BEG"))
}


rm(SNPfiles)
# setnames(SNProws, old = "#CHROM", new = "Chr")
# setkey(annoVarsnp, Chr, BEG, END)
# setkey(SNProws, Chr, BEG, END)
# SNProws <- SNProws[annoVarsnp]
# rm(annoVarsnp)


setkey(annoVarsnp, Chr, BEG)
setkey(SNProws, Chr, BEG)
SNProws <- SNProws[unique(annoVarsnp)]
rm(annoVarsnp)


fname <- combBrain
weights <- fread(fname)
weights$gene <- gsub("[.].*", "", weights$gene)

SNPlist <- unique(weights[weights[, gene] == transcript,
                          rsid])
setkey(SNProws, avsnp150)
if(AlzPheno == "LOAD"){
  rowbetas <- SNProws[J(SNPlist),
                      .(avsnp150, OR, SE , P)]
  df_gwas <- rowbetas[, .(avsnp150, OR, SE, P)]
}else if(AlzPheno == "AAO"){
  rowbetas <- SNProws[J(SNPlist),
                      .(avsnp150, BETA, SE, P)]
  df_gwas <- rowbetas[, .(avsnp150, BETA, SE, P)]
  
}
df_PrdxWts <- weights[weights[,gene] == transcript,
                      .(tissue,rsid, weight)]
rm(SNProws)
rm(weights)
print('Process SNP data')
# df_gwas <- rowbetas[, .(avsnp150, BETA,
#                         SEBETA, PVALUE)]
# rm(rowbetas)
pred <- cohorts
df_gwas <- cbind(df_gwas, pred)
setnames(df_gwas, c('SNP', 'weight', 'se',
                    'pval', 'predictor'))
# rownames(df_gwas) <- NULL
print(head(df_gwas))
print(head(df_PrdxWts))
# Fourth
# For the items in SNPlist, for each tissue, grab the 
# GTEX raw beta,se terms
print("Process GTEX Tissues")
print(cohorts)
print(tissues)
print(transcript)

gtexDat <- fread(gtexFile, stringsAsFactors = FALSE, h = TRUE)
setkey(gtexDat, rs_id_dbSNP151_GRCh38p7)
densemapdata <- gtexDat[!("."), .(variant_id,
                                  rs_id_dbSNP151_GRCh38p7)] 
rm(gtexDat)

gtex <- lapply(tissues, pred_gtex, trans = transcript,
               dmd = densemapdata, snps = SNPlist)
df_gtex <- rbindlist(gtex)

setwd("forest_plot_dataFiles")
setnames(df_gtex, c("Tissue", 'SNP', 'weight',
                    'se', 'pval'))

setnames(df_PrdxWts, c("Tissue","SNP","weight"))
df_PrdxWts$Tissue <- gsub("^Brain_*", "", df_PrdxWts$Tissue)
setkey(df_gtex, Tissue, SNP)
setkey(df_PrdxWts, Tissue, SNP)


df_PrdxWts <- df_PrdxWts[df_gtex, .(Tissue, SNP, weight)]

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

setwd(basedir)

#print('finished')
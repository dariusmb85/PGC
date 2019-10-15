################################################################################
library('data.table')
source('Forest_plot_objects_fxns.R')

basedir <- '/projects/sequence_analysis/vol3/UCSFplink/GemmaGTEx_V6p/output'
cohorts <- c("Alc_Dep","Family_Alcoholism","MarijuanaStrDsr","MJ","NIC")
tissues <- c('Anterior_cingulate_cortex_BA24','Cerebellar_Hemisphere','Cortex',
          'Hippocampus','Nucleus_accumbens_basal_ganglia','Caudate_basal_ganglia'
          ,'Cerebellum','Frontal_Cortex','Hypothalamus','Putamen_basal_ganglia')
################################################################################

options(echo = FALSE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
transcriptIn <- args[1]
transNumID   <- args[2]
resp         <- args[3]
forestMetaPlotname <- "AcrossTissue"
transcriptIn <- eval(parse(text = args[1]))
transNumID   <- eval(parse(text = args[2]))
resp         <- eval(parse(text = args[3]))


# print(args)
#print(transcriptIn)
#print(transcriptIn[transNumID])

#print(transcript[transNumID])
transcriptIn <- transcriptIn[transNumID]
respIn <- resp[transNumID]

##Testing
# transcriptIn <- c('ENSG00000127399.10')#,'ENSG00000084072.12','ENSG00000165502.6',
                  # 'ENSG00000213402.2','ENSG00000197958.8','ENSG00000184983.5',
                  # 'ENSG00000196071.3','ENSG00000237765.2','ENSG00000186230.6',
                  # 'ENSG00000164989.11','ENSG00000183605.12','ENSG00000176390.10',
                  # 'ENSG00000157895.7','ENSG00000186470.9','ENSG00000186468.8',
                  # 'ENSG00000103021.5','ENSG00000084207.11')

#translength <- length(transcriptIn)

###Testing
# tissues <- tissues[1:2]
# respIn <- 5
cohorts <- cohorts[respIn]


# Now we need to fetch data from disparate sources and build a final dataframe of the form
#    cpassoc Shom-Shet phom-phet
# expression beta se
#        snp beta se
#       gwas beta se
# Add separate data table for subplot construction of the forest plots

setwd("/projects/sequence_analysis/vol3/UCSFplink/avsnp150/")
file_list <-list.files(pattern= "*hg19_multianno.txt$",
                       full.names = TRUE)
files <- lapply(file_list, function(x) fread(x, header = TRUE, 
                                             data.table=TRUE))
annoVarsnp <- rbindlist(files)
rm(files)
setnames(annoVarsnp,old = c("Start", "End"), new = c("BEG", "END"))
setwd(basedir)

#####FIX buildVariantToRSIDmap() too slow ########
#densemapdata <- buildVariantToRSIDmap()
filepath <- file.path("","projects","sequence_analysis","vol3","predix_Scan",
                      "GTEx-V6p_flowOver","gtex_data","")
filename <- paste(filepath,
                  "GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt",
                  sep="")
gtexDat <- fread(filename, stringsAsFactors = FALSE, h = TRUE)
setkey(gtexDat, rs_id_dbSNP147_GRCh37p13)
densemapdata <- gtexDat[!("."), .(variant_id,
                                  rs_id_dbSNP147_GRCh37p13)] 
rm(gtexDat)

df_cpassoc <- data.table()
df_gtex <- data.table()
df_gwas <- data.table()
df_predixcan <- data.table()


print(paste('Processing next cohort', cohorts, sep="="))

transcript <- transcriptIn
  
# First
# grab the CPASSOC Across cohort result for the current transcript 
# SHow contribution of COHORT to the overall CPASSOC
print('Process CPASSOC')
cpFile <- file.path("MetaAnalysis_AcrossTissue_Restricted","FDR","FDR_Anno",paste0("MetaxTissue_",
                cohorts, "-FDR.Anno.txt"))
ccData <- fread(cpFile, stringsAsFactors = FALSE, h=TRUE)
setkey(ccData, rs)
transcriptName <- ccData[transcript, gene_name]  
if ( is.na(transcriptName) ){
  transcriptName <- transcript
}
df_cpassoc <- ccData[transcript, .(gene_name, Shom, phom,
                                   Shet, phet)] 
setnames(df_cpassoc, c('Gene Name', 'sHom', 'sHomPval',
                       'sHet', 'sHetPval'))
rm(ccData)
print(df_cpassoc)
# Second
# grab the individual cohort gene-expression beta/se info for all tissues:
# GEMMA associations

print("Process Tissues")
print(tissues)
print(cohorts)
print(transcript)
df_predixcan <- lapply(tissues, pred_df, coh = cohorts, basedir = basedir)
df_predixcan <- rbindlist(df_predixcan)
print(df_predixcan) 
# For each SNP of this transcript. Grab all the known beta/se from the SNP-phenotype association
# Not sure if all these SNPs were used in all the tissue calculations - probably not

# Third phenotype ~ SNPs
# Find all the unique SNPs accross all tissues. Then find the Phenotype_SNPs GWAS data and build a data frame
# Keep the object SNPlist for the GTEX data object as well

SNPfiles <- buildSNPfiles(phe = cohorts)
SNProws <- fread(SNPfiles)
rm(SNPfiles)
setnames(SNProws, old = "#CHROM", new = "Chr")
setkey(annoVarsnp, Chr, BEG, END)
setkey(SNProws, Chr, BEG, END)
SNProws <- SNProws[annoVarsnp]
rm(annoVarsnp)

fname <- '/projects/sequence_analysis/vol3/prediXcan/GTEx-V6p-HapMap-2016-09-08/CombinedDB_Brain.txt'
weights <- fread(fname)#rewrite for data.table syntax

#Maybe want to filter to the top 10 SNPs
SNPlist <- unique(weights[weights[, gene] == transcript, rsid])
setkey(SNProws, avsnp150)
rowbetas <- SNProws[J(SNPlist), 
                    .(avsnp150, BETA, SEBETA , PVALUE)]
rm(SNProws)
rm(weights)
print('Process SNP data')
df_gwas <- rowbetas[, .(avsnp150, BETA,
                        SEBETA, PVALUE),
                    nomatch = 0]
rm(rowbetas)
pred <- cohorts
df_gwas <- cbind(df_gwas, pred)
setnames(df_gwas, c('SNP', 'beta', 'se',
                    'pval', 'predictor'))
rownames(df_gwas) <- NULL
print(head(df_gwas))
# Fourth
# For the items in SNPlist, for each tissue, grab the GTEX raw beta,se terms
# Second
# grab the individual cohort gene-expression beta/se info for all tissues: 
# GEMMA associations

print("Process GTEX Tissues")
print(cohorts)
print(tissues)
print(transcript)
gtex <- lapply(tissues, pred_gtex, trans = transcript,
               dmd = densemapdata, snps = SNPlist)
df_gtex <- rbindlist(gtex)

setwd("forest_plot_dataFiles")
setnames(df_gtex, c("Tissue", 'SNP', 'beta',
                    'se', 'pval'))
print(head(df_gtex))
# df_cpassoc
# df_predixcan
# head(df_gwas)
# head(df_gtex)
# 2298884
datObj <- cpRes(transcript, cohorts, "Across_Tissue",transcriptName, 
                df_cpassoc, df_predixcan,df_gwas,df_gtex)
saveRDS(datObj, file = paste(transcript,cohorts,transcriptName,"Across_Tissue.datObj",sep = "_"))
setwd(basedir)
# for (itissue in tissues){
# filename <- buildGTEXfilenames(itissue)
# print(filename)
# gtexdata <- fread(filename, h=TRUE, stringsAsFactors=FALSE)
# setkey(gtexdata,gene_id)
# gtex <- gtexdata[J(transcript), .(gene_id, variant_id, slope, slope_se, pval_nominal)]
# rm(gtexdata)
#     # Build a SNPlist known to contain hits in the current gtex$variant_id object
#     foundVars <- gtex[, 2]
#     for (i in seq(1,dim(foundVars)[1])){
#       invar <- toString(foundVars[i,1])
#       isnp <- convertVariantToRSID(invar, densemapdata)
#       if (isnp %in% SNPlist){
#         print(isnp)
#         #var <- convertRSIDToVariant(isnp, densemapdata)
#         if ( invar !=""){
#           setkey(gtex, variant_id)
#           beta         <- gtex[J(invar), .(slope)]
#           se           <- gtex[J(invar), .(slope_se)]
#           pval         <- gtex[J(invar), .(pval_nominal)]
#           df_row       <- data.table(itissue, isnp, beta, se, pval)
#           names(df_row)<-c('Tissue', 'SNP', 'beta', 'se', 'pval')
#           #df_gtex = rbind(df_gtex,df_ro, stringsAsFactors=FALSE)
#           print(paste('Got a name ', isnp, sep="="))
#           print(paste(itissue, isnp, sep=" "))
#           print(beta)
#           print(se)
#           df_gtex <- rbind(df_gtex, df_row, stringsAsFactors=FALSE)
#         }
#       }
#     }
#   }
# # }
# 
# 
# names(df_gtex) <- c('Tissue', 'SNP', 'beta', 'se', 'pval')
# rownames(df_gtex) <- NULL
# 
# mylist <- list(transcript = transcript,cpassoc = list(df_cpassoc),
#                predix = list(df_predixcan), gwas = list(df_gwas),
#                gtex = list(df_gtex))
# #lapply(mylist, function(x) write.table( data.table(x), 'test.csv',
#  #                                       append= TRUE, sep='\t' ))
# 
# setwd("forest_plot_dataFiles")
# 
# saveRDS(df_cpassoc, paste(transcript, 'cpassoc_data.rds', sep= "-"))
# saveRDS(df_predixcan, paste(transcript, 'predixcan_data.rds', sep = "-"))
# saveRDS(df_gwas, paste(transcript, 'gwas_data.rds', sep = "-"))
# saveRDS(df_gtex, paste(transcript, 'gtex_data.rds', sep = "-"))
# saveRDS(mylist, paste(transcript, 'all-data.rds', sep = "-"))
# 
# # write.table(df_cpassoc, paste(transcript, 'cpassoc_data.tsv', sep= "-"),
# #             sep='\t', quote=FALSE, row.names = FALSE)
# # write.table(df_predixcan, paste(transcript, 'predixcan_data.tsv', sep = "-"),
# #             sep='\t', quote=FALSE, row.names = FALSE)
# # write.table(df_gwas, paste(transcript, 'gwas_data.tsv', sep = "-"),
# #             sep='\t', quote=FALSE, row.names = FALSE)
# # write.table(df_gtex, paste(transcript, 'gtex_data.tsv', sep = "-"),
# #             sep='\t', quote=FALSE, row.names = FALSE)
# 
# lapply(mylist, function(x) write.table( data.table(x), paste(transcript,
#                                                              'all-data.tsv',
#                                                              sep = "-"),
#                                         append= TRUE, sep=',', row.names = FALSE ))
# setwd(basedir)
print('finished')
#}



setClass("cpRes",
         slots = c(trans = "character",
                   resp = "character",
                   pred = "character",
                   gene = "character",
                   cpAssoc = "data.table",
                   stAssoc = "data.table",
                   gwas = "data.table",
                   gtex = "data.table",
                   predWts = "data.table"))

cpRes <- function(trans, resp, pred, gene, cpAssoc,
                  stAssoc, gwas, gtex, predWts){
  cpOb<- new("cpRes", trans = trans, resp = resp, pred = pred,
             gene = gene, cpAssoc = cpAssoc, stAssoc = stAssoc,
             gwas = gwas, gtex = gtex, predWts = predWts)
  cpOb
}

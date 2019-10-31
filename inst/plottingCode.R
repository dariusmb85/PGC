
## requies 'dlfUtils'
require(dlfUtils)

pltForest <- function(grp, snp, beta, se) {
  lVal <- beta - 2*se
  hVal <- beta + 2*se
  high <- lVal > 0
  low  <- hVal < 0
  col  <- rep("gray", length(grp))
  col[high] <- "#D3560E"
  col[low]  <- "#521168"
  
  nSnp <- length(unique(snp))
  nGrp <- length(unique(grp))
  snpInd <- as.integer(as.factor(snp))
  grpInd <- as.integer(as.factor(grp))
  x <- (grpInd + (snpInd - 1)*(nGrp + 1))/(nGrp + 1)
  par(mar = c(c(1, 4, 1, 1)))
  plot.new()
  # plot.window(xlim = c(0, (nGrp + 1)*nSnp), ylim = range(c(0, lVal, hVal)))
  plot.window(xlim = c(0, nSnp), ylim = range(c(0, lVal, hVal), na.rm = TRUE))
  axis(side = 2)
  segments(x0 = x, x1 = x, y0 = lVal, y1 = hVal, col = col)
  abline(h = 0, lty = "dashed")
}

pltStd <- function(lbl, pval, beta, se) {
  
  nLbl <- length(lbl)
  lVal <- beta - 2*se
  hVal <- beta + 2*se
  high <- lVal > 0
  low  <- hVal < 0
  col  <- rep("gray", nLbl)
  col[high] <- "#D3560E"
  col[low]  <- "#521168"
  
  lblInt <- as.integer(as.factor(lbl))
  
  par(mar = c(4, 25, 1, 1))
  plot.new()
  plot.window(ylim = c(0, nLbl), xlim = range(c(0, lVal, hVal), na.rm = TRUE))
  axis(side = 1)
  title(xlab = "Standard Association Model")
  segments(y0 = lblInt, y1 = lblInt, x0 = lVal, x1 = hVal, col = col)
  abline(v = 0, lty = "dashed")
  text(x = line2user(1, 2), 
       y = lblInt, 
       labels = formatC(pval, format = "e", digits = 2),
       adj = c(1, 0.5),
       xpd = NA)
  text(x = line2user(8, 2), 
       y = lblInt, 
       labels = lbl,
       adj = c(1, 0.5),
       xpd = NA)
  
}

pltCPA <- function(x) {
  
  layout(matrix(c(1, 2, 3, 4), ncol = 1), heights = c(2, 6, 4, 4))
  par(mar = rep(0, 4))
  plot.new()
  plot.window(xlim = 0:1, ylim = c(0, 3))
  tg <- sprintf("%s (%s)", x@trans, x@gene)
  text(x = 0, y = 3, adj = c(0, 1), tg, cex = 2)
  rp <- sprintf("%s ~ %s", x@resp, x@pred)
  text(x = 0, y = 2, adj = c(0, 0.5), rp, cex = 1.5)
  text(x = 0.5, y = 2, adj = c(1, 0.5), "sHom", cex = 1.5)
  text(x = 0.5, y = 1, adj = c(1, 0.5), "sHet", cex = 1.5)
  text(x = 0.65, y = 3, adj = c(0, 1), "Estimate", cex = 1.5)
  text(x = 0.65, y = 2, adj = c(0, 0.5), signif(x@cpAssoc$sHom, 4), cex = 1.5)
  text(x = 0.65, y = 1, adj = c(0, 0.5), signif(x@cpAssoc$sHet, 4), cex = 1.5)
  text(x = 0.8, y = 3, adj = c(0, 1), "p-value", cex = 1.5)
  shomP <- formatC(x@cpAssoc$sHomPval, format = "e", digits = 2)
  text(x = 0.8, y = 2, adj = c(0, 0.5), shomP, cex = 1.5)
  shetP <- formatC(x@cpAssoc$sHetPval, format = "e", digits = 2)
  text(x = 0.8, y = 1, adj = c(0, 0.5), shetP, cex = 1.5)
  
  with(x@stAssoc[order(Tissue)], pltStd(Tissue, pval, beta, se))
  with(x@gwas[order(SNP, predictor)], pltForest(predictor, SNP, beta, se))
  title(ylab = "GWAS Model")
  with(x@gtex[order(SNP, Tissue)], pltForest(Tissue, SNP, beta, se))
  title(ylab = "GTEx Model")
  
}



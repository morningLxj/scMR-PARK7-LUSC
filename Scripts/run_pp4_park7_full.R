suppressWarnings({ library(data.table); library(coloc); library(dplyr) })
eqtl_full_file <- "My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz"
ilcco_file <- "My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv"
out_txt <- "My_MR_Project/PARK7_PP4_FULL.txt"
eq <- fread(eqtl_full_file, sep="\t")
eq <- eq[, .(GENE, RSID, CHR, POS, SPEARMANS_RHO, P_VALUE)]
eg <- eq[GENE == "PARK7"]
if (nrow(eg) > 0) {
  eg[, z := abs(qnorm(pmin(pmax(1 - P_VALUE/2, .Machine$double.eps), 1 - .Machine$double.eps)))]
  eg[, se_eqtl := abs(SPEARMANS_RHO / z)]
  eg <- eg[is.finite(se_eqtl)]
  eg <- eg[, .(rsid = as.character(RSID), CHR = as.character(CHR), POS = as.integer(POS), beta_eqtl = SPEARMANS_RHO, se_eqtl = se_eqtl, p_eqtl = P_VALUE)]
  ilc <- fread(ilcco_file, sep="\t")
  setnames(ilc, c("hm_rsid","hm_beta","standard_error","p_value","hm_chrom","hm_pos"), c("rsid","beta_gwas","se_gwas","p_gwas","CHR","POS"))
  ilc[, rsid := as.character(rsid)]
  ilc[, CHR := as.character(CHR)]
  ilc[, POS := as.integer(POS)]
  m1 <- suppressWarnings(inner_join(eg, ilc[, .(rsid, beta_gwas, se_gwas, p_gwas, CHR, POS)], by = "rsid"))
  if (nrow(m1) < 20) {
    m2 <- suppressWarnings(inner_join(eg[, .(CHR, POS, beta_eqtl, se_eqtl, p_eqtl)], ilc[, .(CHR, POS, beta_gwas, se_gwas, p_gwas)], by = c("CHR","POS")))
    merged <- as.data.frame(m2)
  } else {
    merged <- as.data.frame(m1)
  }
  merged <- merged[is.finite(merged$beta_eqtl) & is.finite(merged$se_eqtl) & is.finite(merged$beta_gwas) & is.finite(merged$se_gwas), ]
  cat(paste0("Matched SNPs:", nrow(merged), "\n"))
  if (nrow(merged) > 20) {
    res <- suppressWarnings(coloc.abf(dataset1 = list(beta = merged$beta_eqtl, varbeta = (merged$se_eqtl)^2, type = "quant", N = 982, sdY = 1), dataset2 = list(beta = merged$beta_gwas, varbeta = (merged$se_gwas)^2, type = "cc", s = 0.3, N = 50000)))
    sm <- as.data.frame(t(res$summary))
    pp4 <- as.numeric(sm[["PP.H4.abf"]])
    dt <- data.table(metric=c("PARK7_PP4_FULL","matched"), value=c(pp4, nrow(merged)))
    write.table(dt, file = out_txt, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  } else {
    dt <- data.table(metric=c("PARK7_PP4_FULL","matched"), value=c(NA, nrow(merged)))
    write.table(dt, file = out_txt, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
} else {
  exp_file <- file.path("My_MR_Project/Exposure", "EXPOSURE_cd8et_PARK7_mr.tsv.gz")
  if (file.exists(exp_file)) {
    d <- fread(exp_file, sep = "\t")
    d[, SNP := as.character(SNP)]
    o <- fread(ilcco_file, sep = "\t")
    setnames(o, c("hm_rsid","hm_beta","standard_error"), c("SNP","by","se_by"))
    m <- merge(d[, .(SNP, beta, se)], o[, .(SNP, by, se_by)], by = "SNP")
    if (nrow(m) >= 10) {
      res <- suppressWarnings(coloc.abf(dataset1 = list(beta = m$beta, varbeta = (m$se)^2, type = "quant", N = 982, sdY = 1), dataset2 = list(beta = m$by, varbeta = (m$se_by)^2, type = "cc", s = 0.3, N = 50000)))
      sm <- as.data.frame(t(res$summary))
      pp4 <- as.numeric(sm[["PP.H4.abf"]])
      dt <- data.table(metric=c("PARK7_PP4_FULL","matched"), value=c(pp4, nrow(m)))
      write.table(dt, file = out_txt, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    } else {
      dt <- data.table(metric=c("PARK7_PP4_FULL","matched"), value=c(NA, nrow(m)))
      write.table(dt, file = out_txt, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  } else {
    dt <- data.table(metric=c("PARK7_PP4_FULL","matched"), value=c(NA, 0))
    write.table(dt, file = out_txt, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}

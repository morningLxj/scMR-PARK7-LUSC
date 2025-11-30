suppressWarnings({
  library(data.table)
  library(coloc)
})
plink_bin <- "My_MR_Project/plink.exe"
if (!file.exists(plink_bin)) plink_bin <- "plink.exe"
bfile <- "My_MR_Project/Reference/g1000_eur"
mr_csv <- "My_MR_Project/CD8ET_TopN_MR_Results.csv"
eqtl_file <- "My_MR_Project/Exposure/cd8et_esnp_table.tsv.gz"
out_file <- "My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv"
exp_dir <- "My_MR_Project/Exposure"
out_csv <- "D:/2026YJ/My_MR_Project/CD8ET_ILCCO_coloc_sensitivity.csv"
mr <- fread(mr_csv)
mr[, p_ivw := as.numeric(p_ivw)]
mr <- mr[p_ivw < 5e-8]
targets <- gsub("cd8et_", "", mr$gene)
out <- fread(out_file, sep = "\t")
out <- out[, .(hm_rsid, hm_chrom, hm_pos, hm_beta, standard_error, hm_effect_allele, hm_other_allele, hm_effect_allele_frequency)]
out[, hm_chrom := as.character(hm_chrom)]
out[, hm_pos := as.integer(hm_pos)]
res_all <- NULL
e_all <- fread(eqtl_file, sep = "\t", select = c("GENE","RSID","CHR","POS","SPEARMANS_RHO","P_VALUE","A2_FREQ_ONEK1K"))
for (g in targets) {
  eg <- e_all[GENE == g]
  if (nrow(eg) == 0) next
  uniq <- unique(eg$CHR)
  chr <- as.character(uniq[which.max(tabulate(match(eg$CHR, uniq)))])
  center <- as.integer(median(eg$POS))
  start <- center - 1000000L
  end <- center + 1000000L
  e_reg <- EG[CHR == chr & POS >= start & POS <= end]
  o_reg <- out[hm_chrom == chr & hm_pos >= start & hm_pos <= end]
  if (nrow(e_reg) > 0) {
    tmp_assoc <- tempfile(pattern = paste0("assoc_", g), fileext = ".txt")
    fwrite(e_reg[, .(SNP = RSID, P = P_VALUE)], tmp_assoc, sep = "\t")
    tmp_out <- tempfile(pattern = paste0("clump_", g))
    args <- c("--bfile", bfile, "--clump", tmp_assoc, "--clump-p1", "1", "--clump-p2", "1", "--clump-kb", "10000", "--clump-r2", "0.001", "--out", tmp_out)
    suppressWarnings(system2(plink_bin, args = args, stdout = TRUE, stderr = TRUE))
    clumped_file <- paste0(tmp_out, ".clumped")
    if (file.exists(clumped_file)) {
      cl <- fread(clumped_file, fill = TRUE)
      if ("SNP" %in% names(cl)) e_reg <- e_reg[RSID %in% unique(cl$SNP)]
    }
  }
  z_corr <- NA_real_
  peak_shared <- 0
  peak_dist <- -1L
  if (nrow(e_reg) > 0 & nrow(o_reg) > 0) {
    e_reg[, z_eqtl := qnorm(pmin(pmax(1 - P_VALUE/2, .Machine$double.eps), 1 - .Machine$double.eps)) * sign(SPEARMANS_RHO)]
    e_reg[, se_eqtl := abs(SPEARMANS_RHO / qnorm(pmin(pmax(1 - P_VALUE/2, .Machine$double.eps), 1 - .Machine$double.eps)))]
    o_reg[, z_ilcco := hm_beta / standard_error]
    setnames(e_reg, "RSID", "SNP")
    setnames(o_reg, "hm_rsid", "SNP")
    j <- merge(e_reg[, .(SNP, CHR, POS, SPEARMANS_RHO, se_eqtl, A2_FREQ_ONEK1K)], o_reg[, .(SNP, hm_beta, standard_error, hm_effect_allele_frequency, hm_pos)], by = "SNP")
    z_eq <- j$SPEARMANS_RHO / j$se_eqtl
    z_gw <- j$hm_beta / j$standard_error
    keep_z <- which(is.finite(z_eq) & is.finite(z_gw))
    if (length(keep_z) >= 3) {
      med_eq <- median(z_eq[keep_z]); mad_eq <- mad(z_eq[keep_z]); thr_eq <- med_eq + 8 * mad_eq
      med_gw <- median(z_gw[keep_z]); mad_gw <- mad(z_gw[keep_z]); thr_gw <- med_gw + 8 * mad_gw
      mask <- (abs(z_eq) <= abs(thr_eq)) & (abs(z_gw) <= abs(thr_gw))
      z_eq <- z_eq[mask]
      z_gw <- z_gw[mask]
    }
    z_corr <- if (nrow(j) >= 3) suppressWarnings(cor(z_eq, z_gw)) else NA_real_
    if (nrow(j) >= 5) {
      dfj <- as.data.frame(j)
      dfj$z_eq <- z_eq
      dfj$z_gw <- z_gw
      fit <- suppressWarnings(lm(z_gw ~ z_eq, data = dfj))
      cd <- suppressWarnings(cooks.distance(fit))
      thr <- 4 / nrow(dfj)
      keep <- which(is.finite(cd) & cd <= thr)
      if (length(keep) >= 3) {
        z_corr <- suppressWarnings(cor(dfj$z_eq[keep], dfj$z_gw[keep]))
      }
    }
    peak_e <- e_reg[which.max(abs(z_eqtl))]
    peak_o <- o_reg[which.max(abs(z_ilcco))]
    peak_shared <- as.integer(peak_e$SNP == peak_o$SNP)
    peak_dist <- abs(as.integer(peak_e$POS) - as.integer(peak_o$hm_pos))
    pp4 <- NA_real_
    if (nrow(j) >= 10) {
      d1 <- list(beta = j$SPEARMANS_RHO, varbeta = (j$se_eqtl)^2, type = "quant", MAF = j$A2_FREQ_ONEK1K)
      d2 <- list(beta = j$hm_beta, varbeta = (j$standard_error)^2, type = "quant", MAF = j$hm_effect_allele_frequency)
      csum <- suppressWarnings(coloc.abf(d1, d2)$summary)
      if (!is.null(csum) && "PP.H4.abf" %in% names(csum)) pp4 <- as.numeric(csum["PP.H4.abf"])
    }
  } else {
    exp_file <- file.path(exp_dir, paste0("EXPOSURE_cd8et_", g, "_mr.tsv.gz"))
    if (file.exists(exp_file)) {
      d <- fread(exp_file, sep = "\t")
      d[, SNP := as.character(SNP)]
      o <- out
      setnames(o, c("hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error","hm_chrom","hm_pos"), c("SNP","ea_outcome","nea_outcome","by","se_by","CHR","POS"))
      m <- merge(d, o[, .(SNP, ea_outcome, nea_outcome, by, se_by, CHR, POS)], by = "SNP")
      if (nrow(m) >= 2) {
        idx <- which(m$effect_allele != m$ea_outcome)
        if (length(idx) > 0) m$by[idx] <- -m$by[idx]
        z_eqtl_ins <- m$beta / m$se
        z_ilcco_ins <- m$by / m$se_by
        z_corr <- suppressWarnings(cor(z_eqtl_ins, z_ilcco_ins))
        pe <- m[which.max(abs(z_eqtl_ins))]
        po <- m[which.max(abs(z_ilcco_ins))]
        peak_shared <- as.integer(pe$SNP == po$SNP)
        peak_dist <- as.integer(abs(as.integer(pe$POS) - as.integer(po$POS)))
        chr <- as.character(modes <- names(sort(table(m$CHR), decreasing = TRUE))[1])
        center <- as.integer(median(m$POS))
      }
    }
  }
  exp_file <- file.path(exp_dir, paste0("EXPOSURE_cd8et_", g, "_mr.tsv.gz"))
  sens_Q <- NA_real_
  sens_Q_p <- NA_real_
  eg_int <- NA_real_
  eg_int_se <- NA_real_
  eg_int_p <- NA_real_
  conc <- NA_real_
  if (file.exists(exp_file)) {
    d <- fread(exp_file, sep = "\t")
    d[, SNP := as.character(SNP)]
    o <- out
    setnames(o, c("hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error"), c("SNP","ea_outcome","nea_outcome","by","se_by"))
    m <- merge(d, o[, .(SNP, ea_outcome, nea_outcome, by, se_by)], by = "SNP")
    if (nrow(m) >= 2) {
      idx <- which(m$effect_allele != m$ea_outcome)
      if (length(idx) > 0) {
        m$by[idx] <- -m$by[idx]
      }
      bx <- m$beta
      by <- m$by
      se_bx <- m$se
      se_by <- m$se_by
      ratio <- by / bx
      var_ratio <- (se_by^2 / (bx^2)) + ((by^2) * (se_bx^2) / (bx^4))
      w <- 1 / var_ratio
      ivw <- sum(w * ratio) / sum(w)
      q <- sum(w * (ratio - ivw)^2)
      sens_Q <- q
      sens_Q_p <- 1 - pchisq(q, df = nrow(m) - 1)
      fit <- suppressWarnings(lm(by ~ bx, weights = 1/(se_by^2)))
      co <- coef(fit)
      vc <- vcov(fit)
      eg_int <- co[1]
      eg_int_se <- sqrt(vc[1,1])
      eg_int_p <- 2 * pnorm(abs(eg_int/eg_int_se), lower.tail = FALSE)
      conc <- mean((bx * by) > 0)
    }
  }
  shared <- if (exists("j")) nrow(j) else if (exists("m")) nrow(m) else 0L
  res <- data.table(gene = g, chr = chr, center_pos = center, n_eqtl_region = if (exists("m")) nrow(m) else nrow(e_reg), n_ilcco_region = if (exists("m")) nrow(o_reg) else nrow(o_reg), n_shared_rsids = shared, z_corr = z_corr, peak_shared = peak_shared, peak_dist_bp = peak_dist, Q = sens_Q, Q_p = sens_Q_p, egger_intercept = eg_int, egger_intercept_se = eg_int_se, egger_intercept_p = eg_int_p, concordance = conc, pp4 = if (exists("pp4")) pp4 else NA_real_)
  res_all <- rbind(res_all, res, fill = TRUE)
}
if (is.null(res_all)) {
  schema <- data.table(gene=character(), chr=character(), center_pos=integer(), n_eqtl_region=integer(), n_ilcco_region=integer(), n_shared_rsids=integer(), z_corr=double(), peak_shared=integer(), peak_dist_bp=integer(), Q=double(), Q_p=double(), egger_intercept=double(), egger_intercept_se=double(), egger_intercept_p=double(), concordance=double(), pp4=double())
  fwrite(schema, out_csv)
} else {
  fwrite(res_all, out_csv)
}
writeLines(paste0("saved ", out_csv))

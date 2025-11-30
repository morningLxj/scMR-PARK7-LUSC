suppressWarnings({
  library(data.table)
})
plink_bin <- "My_MR_Project/plink.exe"
if (!file.exists(plink_bin)) plink_bin <- "plink.exe"
bfile <- "My_MR_Project/Reference/g1000_eur"
exp_dir <- "My_MR_Project/Exposure"
out_file <- "My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv"
mr_out_csv <- "My_MR_Project/CD8ET_TopN_MR_Results.csv"
exps <- list.files(exp_dir, pattern = "EXPOSURE_.*_mr.tsv.gz", full.names = TRUE)
log_file <- "My_MR_Project/mr_run_log.txt"
cat(paste0("exposure_files_count=", length(exps), "\n"), file = log_file, append = FALSE)
out <- tryCatch(
  fread(out_file, sep = "\t"),
  error = function(e) {
    writeLines(paste0("out_read_error=", e$message), con = log_file)
    NULL
  }
)
if (is.null(out)) {
  writeLines("no MR results produced")
  quit(save = "no")
}
setnames(out, c("hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error","hm_effect_allele_frequency","p_value"), c("SNP","ea_outcome","nea_outcome","by","se_by","eaf_outcome","p_outcome"))
mr_all <- NULL
clump_rsids <- function(rsids, pvals, tag) {
  tmp_assoc <- tempfile(pattern = paste0("assoc_",tag), fileext = ".txt")
  fwrite(data.table(SNP = rsids, P = pvals), tmp_assoc, sep = "\t")
  tmp_out <- tempfile(pattern = paste0("clump_",tag))
  args <- c("--bfile", bfile, "--clump", tmp_assoc, "--clump-p1", "1", "--clump-p2", "1", "--clump-kb", "10000", "--clump-r2", "0.001", "--out", tmp_out)
  system2(plink_bin, args = args, stdout = TRUE, stderr = TRUE)
  clumped_file <- paste0(tmp_out, ".clumped")
  if (!file.exists(clumped_file)) return(rsids)
  cl <- fread(clumped_file, fill = TRUE)
  if (!"SNP" %in% names(cl)) return(rsids)
  unique(cl$SNP)
}
for (f in exps) {
  d <- fread(f, sep = "\t")
  gene <- gsub("EXPOSURE_", "", basename(f))
  gene <- gsub("_mr.tsv.gz", "", gene)
  rs <- d$SNP
  pv <- d$pval
  rs_clump <- tryCatch({clump_rsids(rs, pv, gene)}, error = function(e) rs)
  d <- d[d$SNP %in% rs_clump]
  if (nrow(d) < 1) next
  o <- out[out$SNP %in% d$SNP]
  if (nrow(o) == 0) next
  m <- merge(d, o, by = "SNP")
  cat(paste0("merged_rows_", gene, "=", nrow(m), "\n"), file = log_file, append = TRUE)
  idx <- which(m$effect_allele != m$ea_outcome)
  if (length(idx) > 0) {
    m$by[idx] <- -m$by[idx]
    tmp <- m$ea_outcome[idx]
    m$ea_outcome[idx] <- m$nea_outcome[idx]
    m$nea_outcome[idx] <- tmp
  }
  bx <- m$beta
  by <- m$by
  se_bx <- m$se
  se_by <- m$se_by
  ratio <- by / bx
  var_ratio <- (se_by^2 / (bx^2)) + ((by^2) * (se_bx^2) / (bx^4))
  w <- 1 / var_ratio
  ivw <- sum(w * ratio) / sum(w)
  se_ivw <- sqrt(1 / sum(w))
  p_ivw <- 2 * pnorm(abs(ivw / se_ivw), lower.tail = FALSE)
  q <- sum(w * (ratio - ivw)^2)
  q_p <- 1 - pchisq(q, df = nrow(m) - 1)
  fit <- suppressWarnings(lm(by ~ bx, weights = 1/(se_by^2)))
  co <- coef(fit)
  vc <- vcov(fit)
  egger_intercept <- co[1]
  egger_intercept_se <- sqrt(vc[1,1])
  egger_intercept_p <- 2 * pnorm(abs(egger_intercept/egger_intercept_se), lower.tail = FALSE)
  egger_slope <- co[2]
  conc <- mean((bx * by) > 0)
  res <- data.table(gene = gene, n_iv = nrow(m), beta_ivw = ivw, se_ivw = se_ivw, p_ivw = p_ivw, Q = q, Q_p = q_p, egger_intercept = egger_intercept, egger_intercept_se = egger_intercept_se, egger_intercept_p = egger_intercept_p, egger_slope = egger_slope, concordance = conc)
  mr_all <- rbind(mr_all, res, fill = TRUE)
}
if (length(exps) == 0) {
  eqtl_file <- file.path(exp_dir, "cd8et_eqtl_table.tsv.gz")
  d_all <- tryCatch(
    fread(eqtl_file, sep = "\t", nrows = 1000000, select = c("GENE","RSID","A2","A1","A2_FREQ_ONEK1K","SPEARMANS_RHO","P_VALUE")),
    error = function(e) {
      cat(paste0("fallback_read_error=", e$message, "\n"), file = log_file, append = TRUE)
      NULL
    }
  )
  if (is.null(d_all)) {
    writeLines("no MR results produced")
    quit(save = "no")
  }
  d_sub <- d_all[P_VALUE < 1e-6]
  cat(paste0("fallback_rows_p_lt_1e6=", nrow(d_sub), "\n"), file = log_file, append = TRUE)
  cnt <- d_sub[, .N, by = GENE][order(-N)]
  targets <- head(cnt$GENE, 10)
  cat(paste0("fallback_targets_count=", length(targets), "\n"), file = log_file, append = TRUE)
  process_gene <- function(g) {
    d <- d_sub[GENE == g, .(SNP = RSID,
                             effect_allele = A2,
                             other_allele = A1,
                             eaf = A2_FREQ_ONEK1K,
                             beta = SPEARMANS_RHO,
                             pval = P_VALUE)]
    if (nrow(d) == 0) return(NULL)
    d[, se := abs(beta / qnorm(1 - pval/2))]
    d <- d[!is.na(se)]
    if (nrow(d) < 1) return(NULL)
    rs <- d$SNP
    pv <- d$pval
    rs_clump <- tryCatch({clump_rsids(rs, pv, g)}, error = function(e) rs)
    d <- d[d$SNP %in% rs_clump]
    if (nrow(d) < 1) return(NULL)
    o <- out[out$SNP %in% d$SNP]
    if (nrow(o) == 0) return(NULL)
    m <- merge(d, o, by = "SNP")
    cat(paste0("merged_rows_fallback_", g, "=", nrow(m), "\n"), file = log_file, append = TRUE)
    idx <- which(m$effect_allele != m$ea_outcome)
    if (length(idx) > 0) {
      m$by[idx] <- -m$by[idx]
      tmp <- m$ea_outcome[idx]
      m$ea_outcome[idx] <- m$nea_outcome[idx]
      m$nea_outcome[idx] <- tmp
    }
    bx <- m$beta
    by <- m$by
    se_bx <- m$se
    se_by <- m$se_by
    ratio <- by / bx
    var_ratio <- (se_by^2 / (bx^2)) + ((by^2) * (se_bx^2) / (bx^4))
    w <- 1 / var_ratio
    ivw <- sum(w * ratio) / sum(w)
    se_ivw <- sqrt(1 / sum(w))
    p_ivw <- 2 * pnorm(abs(ivw / se_ivw), lower.tail = FALSE)
    q <- sum(w * (ratio - ivw)^2)
    q_p <- 1 - pchisq(q, df = nrow(m) - 1)
    fit <- suppressWarnings(lm(by ~ bx, weights = 1/(se_by^2)))
    co <- coef(fit)
    vc <- vcov(fit)
    egger_intercept <- co[1]
    egger_intercept_se <- sqrt(vc[1,1])
    egger_intercept_p <- 2 * pnorm(abs(egger_intercept/egger_intercept_se), lower.tail = FALSE)
    egger_slope <- co[2]
    conc <- mean((bx * by) > 0)
    data.table(gene = paste0("cd8et_", g), n_iv = nrow(m), beta_ivw = ivw, se_ivw = se_ivw, p_ivw = p_ivw, Q = q, Q_p = q_p, egger_intercept = egger_intercept, egger_intercept_se = egger_intercept_se, egger_intercept_p = egger_intercept_p, egger_slope = egger_slope, concordance = conc)
  }
  for (g in targets) {
    res <- process_gene(g)
    if (!is.null(res)) mr_all <- rbind(mr_all, res, fill = TRUE)
  }
}
if (is.null(mr_all)) {
  writeLines("no MR results produced")
  quit(save = "no")
}
fwrite(mr_all, mr_out_csv)
writeLines(paste0("saved MR results to ", mr_out_csv))

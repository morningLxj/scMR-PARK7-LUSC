suppressWarnings({ library(data.table); library(coloc); library(dplyr) })
eqtl_file <- "My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz"
ilcco_csv <- "My_MR_Project/ILCCO_processed.csv"
ilcco_tsv <- "My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv"
finngen_gz <- "My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz"
out_csv <- "My_MR_Project/PARK7_PP4_FULL.csv"
dbg_csv <- "My_MR_Project/PARK7_merged_debug.csv"
log_txt <- "My_MR_Project/robust_pp4_log.txt"
try({ if (file.exists(out_csv)) unlink(out_csv) }, silent = TRUE)
try({ if (file.exists(log_txt)) unlink(log_txt); file.create(log_txt) }, silent = TRUE)
eq <- tryCatch(fread(eqtl_file, sep = "\t"), error = function(e) { cat(paste0("eqtl_read_error=", e$message, "\n"), file = log_txt, append = TRUE); data.table() })
setnames(eq, names(eq), toupper(names(eq)))
if (!("A2_FREQ_ONEK1K" %in% names(eq))) {
  eq[, A2_FREQ_ONEK1K := NA_real_]
}
cat("step=read_eq\n")
setnames(eq, names(eq), toupper(names(eq)))
eg <- eq[GENE == "PARK7", .(RSID, CHR, POS, SPEARMANS_RHO, P_VALUE, A2_FREQ_ONEK1K)]
if (nrow(eg) == 0) {
  eg <- eq[CHR == "1" & as.integer(POS) >= 6000000L & as.integer(POS) <= 10000000L, .(RSID, CHR, POS, SPEARMANS_RHO, P_VALUE)]
}
cat(paste0("eQTL rows=", nrow(eg), "\n"), file = log_txt, append = TRUE)
eg[, z := abs(qnorm(pmin(pmax(1 - P_VALUE/2, .Machine$double.eps), 1 - .Machine$double.eps)))]
eg[, se_eqtl := abs(SPEARMANS_RHO / z)]
eg <- eg[is.finite(se_eqtl)]
eg <- eg[, .(rsid = as.character(RSID), chr = as.character(CHR), pos = as.integer(POS), beta_eqtl = SPEARMANS_RHO, se_eqtl = se_eqtl, p_eqtl = P_VALUE, maf = as.numeric(A2_FREQ_ONEK1K))]
ilc <- NULL
if (file.exists(ilcco_csv)) {
  ilc <- tryCatch(fread(ilcco_csv), error = function(e) { cat(paste0("ilcco_csv_read_error=", e$message, "\n"), file = log_txt, append = TRUE); NULL })
  if (!is.null(ilc)) {
    setnames(ilc, names(ilc), tolower(names(ilc)))
    if (!("rsid" %in% names(ilc))) { cat("ilcco_csv_missing_rsid\n", file = log_txt, append = TRUE) ; ilc <- NULL }
    if (!is.null(ilc)) ilc <- ilc[, .(rsid = as.character(rsid), chr = as.character(ifelse("chromosome" %in% names(ilc), chromosome, NA)), pos = as.integer(ifelse("position" %in% names(ilc), position, NA)), beta_gwas = as.numeric(beta), se_gwas = as.numeric(se), p_gwas = as.numeric(ifelse("p" %in% names(ilc), p, NA)))]
  }
}
if (is.null(ilc)) {
  ilc <- tryCatch(fread(ilcco_tsv, sep = "\t"), error = function(e) { cat(paste0("ilcco_tsv_read_error=", e$message, "\n"), file = log_txt, append = TRUE); NULL })
  if (!is.null(ilc)) {
    setnames(ilc, c("hm_rsid","hm_chrom","hm_pos","hm_beta","standard_error","p_value"), c("rsid","chr","pos","beta_gwas","se_gwas","p_gwas"))
    ilc[, rsid := as.character(rsid)]
    ilc[, chr := as.character(chr)]
    ilc[, pos := as.integer(pos)]
  } else {
    cat("ilcco_unavailable\n", file = log_txt, append = TRUE)
    dt <- data.table(status = "failed", reason = "no_ilcco", n_snps = 0L)
    write.csv(dt, out_csv, row.names = FALSE)
  }
}
m1 <- tryCatch(suppressWarnings(inner_join(eg, ilc[, .(rsid, chr, pos, beta_gwas, se_gwas, p_gwas)], by = "rsid")), error = function(e) { cat(paste0("merge_rsid_error=", e$message, "\n"), file = log_txt, append = TRUE); NULL })
merged <- if (is.null(m1) || nrow(m1) < 10) tryCatch(suppressWarnings(inner_join(eg[, .(chr, pos, beta_eqtl, se_eqtl, p_eqtl, maf)], ilc[, .(chr, pos, beta_gwas, se_gwas, p_gwas)], by = c("chr","pos"))), error = function(e) { cat(paste0("merge_chrpos_error=", e$message, "\n"), file = log_txt, append = TRUE); NULL }) else m1
if (is.null(merged)) { dt <- data.table(status = "failed", reason = "merge_error", n_snps = 0L); tryCatch({ fwrite(dt, out_csv) }, error = function(e) { write.csv(dt, out_csv, row.names = FALSE) }); merged <- data.table() }
n_before <- nrow(merged)
merged <- merged %>% filter(is.finite(beta_eqtl), is.finite(se_eqtl), is.finite(beta_gwas), is.finite(se_gwas))
merged <- merged %>% filter(se_eqtl > 0, se_gwas > 0)
n_after_var <- nrow(merged)
var_removed <- n_before - n_after_var
cat(paste0("merged_rows=", nrow(merged), "\n"), file = log_txt, append = TRUE)
cat(paste0("removed_zero_var=", var_removed, "\n"), file = log_txt, append = TRUE)
cat(paste0("step=merge_done n=", nrow(merged), "\n"))
if (!("snp" %in% names(merged))) {
  if ("rsid" %in% names(merged) && any(!is.na(merged$rsid))) {
    merged$snp <- as.character(merged$rsid)
  } else if (all(c("chr","pos") %in% names(merged))) {
    merged$snp <- paste0(as.character(merged$chr), ":", as.integer(merged$pos))
  }
}
if ("p_eqtl" %in% names(merged) || "p_gwas" %in% names(merged)) {
  n_before_dedup <- nrow(merged)
  merged <- merged %>% mutate(p_min = pmin(p_eqtl, p_gwas, na.rm = TRUE)) %>% arrange(snp, p_min) %>% distinct(snp, .keep_all = TRUE)
  dedup_removed <- n_before_dedup - nrow(merged)
  cat(paste0("dedup_removed=", dedup_removed, "\n"), file = log_txt, append = TRUE)
}
if ("maf" %in% names(merged)) {
  n_before_maf <- nrow(merged)
  merged <- merged %>% filter(is.finite(maf), maf > 0, maf < 1)
  maf_removed <- n_before_maf - nrow(merged)
  cat(paste0("removed_maf_out_of_range=", maf_removed, "\n"), file = log_txt, append = TRUE)
} else {
  maf_removed <- NA_integer_
}
if ("p_eqtl" %in% names(merged)) {
  merged <- merged %>% arrange(p_eqtl)
  if (nrow(merged) > 3000) { merged <- merged[1:3000,] }
}
if (nrow(merged) > 10) {
  d1 <- list(beta = merged$beta_eqtl, varbeta = (merged$se_eqtl)^2, type = "quant", N = 982, sdY = 1, snp = merged$snp)
  d2 <- list(beta = merged$beta_gwas, varbeta = (merged$se_gwas)^2, type = "cc", s = 0.3, N = 50000, snp = merged$snp)
  cat("step=coloc_start\n")
  err_message <- NULL
  res <- tryCatch(suppressWarnings(coloc.abf(dataset1 = d1, dataset2 = d2)), error = function(e) { err_message <<- e$message; cat(paste0("coloc_error=", e$message, "\n"), file = log_txt, append = TRUE); NULL })
  if (is.null(res)) {
    dt <- data.table(status = "failed", reason = if (is.null(err_message)) "coloc_error" else err_message, n_snps = nrow(merged), n_distinct_snp = length(unique(merged$snp)), removed_zero_var = var_removed, removed_maf = maf_removed)
    write.csv(dt, out_csv, row.names = FALSE)
    cat("status=failed\n", file = log_txt, append = TRUE)
    exp_file <- file.path("My_MR_Project/Exposure", "EXPOSURE_cd8et_PARK7_mr.tsv.gz")
    if (file.exists(exp_file)) {
      d <- tryCatch(fread(exp_file, sep = "\t"), error = function(e) { NULL })
      if (!is.null(d)) {
        d[, SNP := as.character(SNP)]
        o <- copy(ilc)
        setnames(o, "rsid", "SNP")
        m <- merge(d[, .(SNP, beta, se)], o[, .(SNP, beta_gwas, se_gwas)], by = "SNP")
        m <- m[is.finite(beta) & is.finite(se) & is.finite(beta_gwas) & is.finite(se_gwas)]
        m <- m[se > 0 & se_gwas > 0]
        m <- m[!is.na(SNP)]
        if (nrow(m) >= 2) {
          setorder(m, SNP, se)
          m <- unique(m, by = "SNP")
        }
        if (nrow(m) >= 10) {
          d1 <- list(beta = m$beta, varbeta = (m$se)^2, type = "quant", N = 982, sdY = 1, snp = m$SNP)
          d2 <- list(beta = m$beta_gwas, varbeta = (m$se_gwas)^2, type = "cc", s = 0.3, N = 50000, snp = m$SNP)
          r2_err <- NULL
          r2 <- tryCatch(suppressWarnings(coloc.abf(dataset1 = d1, dataset2 = d2)), error = function(e) { r2_err <<- e$message; NULL })
          if (!is.null(r2)) {
            sm2 <- as.data.frame(t(r2$summary))
            sm2$n_snps <- nrow(m)
            write.csv(sm2, out_csv, row.names = FALSE)
            cat("status=fallback_ok\n", file = log_txt, append = TRUE)
          } else {
            dt2 <- data.table(status = "failed", reason = if (is.null(r2_err)) "fallback_coloc_error" else r2_err, n_snps = nrow(m))
            write.csv(dt2, out_csv, row.names = FALSE)
            cat("status=fallback_failed\n", file = log_txt, append = TRUE)
          }
        }
      }
    }
  } else {
    sm <- as.data.frame(t(res$summary))
    sm$n_snps <- nrow(merged)
    write.csv(sm, out_csv, row.names = FALSE)
    tryCatch({ fwrite(head(merged, 50), dbg_csv) }, error = function(e) { write.csv(head(merged, 50), dbg_csv, row.names = FALSE) })
    cat("step=write_ok\n", file = log_txt, append = TRUE)
    cat("status=ok\n", file = log_txt, append = TRUE)
  }
} else {
  dt <- data.table(status = "failed", n_snps = nrow(merged))
  write.csv(dt, out_csv, row.names = FALSE)
  cat("status=failed\n", file = log_txt, append = TRUE)
  exp_file <- file.path("My_MR_Project/Exposure", "EXPOSURE_cd8et_PARK7_mr.tsv.gz")
  if (file.exists(exp_file)) {
    d <- tryCatch(fread(exp_file, sep = "\t"), error = function(e) { NULL })
    if (!is.null(d)) {
      d[, SNP := as.character(SNP)]
      o <- ilc
      setnames(o, "rsid", "SNP")
      m <- merge(d[, .(SNP, beta, se)], o[, .(SNP, beta_gwas, se_gwas)], by = "SNP")
      m <- m[is.finite(beta) & is.finite(se) & is.finite(beta_gwas) & is.finite(se_gwas)]
      m <- m[se > 0 & se_gwas > 0]
      m <- m[!is.na(SNP)]
      if (nrow(m) >= 2) {
        setorder(m, SNP, se)
        m <- unique(m, by = "SNP")
      }
      if (nrow(m) >= 10) {
        d1 <- list(beta = m$beta, varbeta = (m$se)^2, type = "quant", N = 982, sdY = 1, snp = m$SNP)
        d2 <- list(beta = m$beta_gwas, varbeta = (m$se_gwas)^2, type = "cc", s = 0.3, N = 50000, snp = m$SNP)
        r2_err <- NULL
        r2 <- tryCatch(suppressWarnings(coloc.abf(dataset1 = d1, dataset2 = d2)), error = function(e) { r2_err <<- e$message; NULL })
        if (!is.null(r2)) {
          sm2 <- as.data.frame(t(r2$summary))
          sm2$n_snps <- nrow(m)
          tryCatch({ fwrite(sm2, out_csv) }, error = function(e) { write.csv(sm2, out_csv, row.names = FALSE) })
          cat("status=fallback_ok\n", file = log_txt, append = TRUE)
        } else {
          dt2 <- data.table(status = "failed", reason = if (is.null(r2_err)) "fallback_coloc_error" else r2_err, n_snps = nrow(m))
          tryCatch({ fwrite(dt2, out_csv) }, error = function(e) { write.csv(dt2, out_csv, row.names = FALSE) })
          cat("status=fallback_failed\n", file = log_txt, append = TRUE)
        }
      }
    }
  }
}
if (!file.exists(out_csv)) {
  dt <- data.table(status = "failed", reason = "no_output_written", n_snps = if (exists("merged")) nrow(merged) else 0L)
  write.csv(dt, out_csv, row.names = FALSE)
}

fg_out_csv <- "My_MR_Project/PARK7_PP4_FINNGEN.csv"
fg_num_txt <- "My_MR_Project/PARK7_PP4_FINNGEN_NUMBER.txt"
if (file.exists(finngen_gz)) {
  fg <- tryCatch(fread(finngen_gz, sep = "\t"), error = function(e) { NULL })
  if (!is.null(fg)) {
    nms <- names(fg)
    if ("rsids" %in% nms) {
      fg[, rsid := tstrsplit(as.character(rsids), ";", fixed = TRUE, keep = 1)]
    } else if ("variant" %in% nms) {
      fg[, rsid := as.character(variant)]
    } else {
      fg[, rsid := NA_character_]
    }
    xchr <- if ("#chrom" %in% nms) "#chrom" else if ("chrom" %in% nms) "chrom" else NA_character_
    xpos <- if ("pos" %in% nms) "pos" else NA_character_
    xbeta <- if ("beta" %in% nms) "beta" else NA_character_
    xse <- if ("sebeta" %in% nms) "sebeta" else if ("se" %in% nms) "se" else NA_character_
    xp <- if ("pval" %in% nms) "pval" else if ("pvalue" %in% nms) "pvalue" else NA_character_
    fg <- fg[, .(rsid = as.character(rsid), chr = as.character(get(xchr)), pos = as.integer(get(xpos)), beta_gwas = as.numeric(get(xbeta)), se_gwas = as.numeric(get(xse)), p_gwas = if (!is.na(xp)) as.numeric(get(xp)) else NA_real_)]
    fg <- fg[is.finite(beta_gwas) & is.finite(se_gwas)]
    mfg1 <- tryCatch(suppressWarnings(inner_join(eg, fg[, .(rsid, chr, pos, beta_gwas, se_gwas, p_gwas)], by = "rsid")), error = function(e) { NULL })
    merged_fg <- if (is.null(mfg1) || nrow(mfg1) < 10) tryCatch(suppressWarnings(inner_join(eg[, .(chr, pos, beta_eqtl, se_eqtl, p_eqtl, maf)], fg[, .(chr, pos, beta_gwas, se_gwas, p_gwas)], by = c("chr","pos"))), error = function(e) { NULL }) else mfg1
    if (!is.null(merged_fg)) {
      merged_fg <- merged_fg %>% filter(is.finite(beta_eqtl), is.finite(se_eqtl), is.finite(beta_gwas), is.finite(se_gwas))
      merged_fg <- merged_fg %>% filter(se_eqtl > 0, se_gwas > 0)
      if (!("snp" %in% names(merged_fg))) {
        if ("rsid" %in% names(merged_fg) && any(!is.na(merged_fg$rsid))) {
          merged_fg$snp <- as.character(merged_fg$rsid)
        } else if (all(c("chr","pos") %in% names(merged_fg))) {
          merged_fg$snp <- paste0(as.character(merged_fg$chr), ":", as.integer(merged_fg$pos))
        }
      }
      if ("p_eqtl" %in% names(merged_fg) || "p_gwas" %in% names(merged_fg)) {
        merged_fg <- merged_fg %>% mutate(p_min = pmin(p_eqtl, p_gwas, na.rm = TRUE)) %>% arrange(snp, p_min) %>% distinct(snp, .keep_all = TRUE)
      }
      if ("p_eqtl" %in% names(merged_fg)) {
        merged_fg <- merged_fg %>% arrange(p_eqtl)
        if (nrow(merged_fg) > 3000) merged_fg <- merged_fg[1:3000,]
      }
      if (nrow(merged_fg) > 10) {
        d1 <- list(beta = merged_fg$beta_eqtl, varbeta = (merged_fg$se_eqtl)^2, type = "quant", N = 982, sdY = 1, snp = merged_fg$snp)
        d2 <- list(beta = merged_fg$beta_gwas, varbeta = (merged_fg$se_gwas)^2, type = "cc", s = 0.3, N = 50000, snp = merged_fg$snp)
        fg_res <- tryCatch(suppressWarnings(coloc.abf(dataset1 = d1, dataset2 = d2)), error = function(e) { NULL })
        if (!is.null(fg_res)) {
          smfg <- as.data.frame(t(fg_res$summary))
          smfg$n_snps <- nrow(merged_fg)
          write.csv(smfg, fg_out_csv, row.names = FALSE)
          writeLines(as.character(smfg[["PP.H4.abf"]]), fg_num_txt)
        } else {
          dtfg <- data.table(status = "failed", n_snps = nrow(merged_fg))
          write.csv(dtfg, fg_out_csv, row.names = FALSE)
          writeLines("Failed", fg_num_txt)
        }
      } else {
        dtfg <- data.table(status = "failed", n_snps = nrow(merged_fg))
        write.csv(dtfg, fg_out_csv, row.names = FALSE)
        writeLines(paste("subset_too_small", nrow(merged_fg)), fg_num_txt)
      }
    }
  }
}

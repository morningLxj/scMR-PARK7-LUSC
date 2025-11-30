suppressWarnings({ library(data.table); library(dplyr); library(coloc) })

plink_bin <- "My_MR_Project/plink.exe"
if (!file.exists(plink_bin)) plink_bin <- "plink.exe"
bfile <- "My_MR_Project/Reference/g1000_eur"

log_txt <- "My_MR_Project/batch_full_log.txt"
try({ if (file.exists(log_txt)) unlink(log_txt) }, silent = TRUE)
try({ file.create(log_txt) }, silent = TRUE)

ilcco_tsv <- "My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv"
finngen_gz <- "My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz"
mr_csv <- "My_MR_Project/AllCells_MR_Results.csv"
out_csv <- "My_MR_Project/Complete_Analysis_Summary.csv"

# Initialize header if file doesn't exist
header_dt <- data.table(
  gene=character(), cell=character(), beta_ivw=double(), se_ivw=double(), p_ivw=double(),
  pp4_ilcco=character(), n_snps_ilcco=integer(), status_ilcco=character(),
  pp4_finngen=character(), n_snps_finngen=integer(), status_finngen=character()
)

# Only write header if file doesn't exist
if (!file.exists(out_csv)) {
  fwrite(header_dt[0], out_csv)
  completed_tasks <- character(0)
} else {
  # Read existing file to find completed tasks
  existing <- tryCatch(fread(out_csv), error = function(e) data.table())
  if (nrow(existing) > 0 && all(c("gene", "cell") %in% names(existing))) {
    completed_tasks <- paste(existing$cell, existing$gene, sep = "_")
  } else {
    completed_tasks <- character(0)
  }
}

cells <- c(
  "bin", "bmem", "cd4et", "cd4nc", "cd4sox4", 
  "cd8et", "cd8nc", "cd8s100b", "dc", 
  "monoc", "mononc", "nk", "nkr", "plasma"
)
cat(paste0("cells_n=", length(cells), "\n"))

eqtl_files <- list(
  bin = "My_MR_Project/Exposure/bin_eqtl_table.tsv.gz",
  bmem = "My_MR_Project/Exposure/bmem_eqtl_table.tsv.gz",
  cd4et = "My_MR_Project/Exposure/cd4et_eqtl_table.tsv.gz",
  cd4nc = "My_MR_Project/Exposure/cd4nc_eqtl_table.tsv.gz",
  cd4sox4 = "My_MR_Project/Exposure/cd4sox4_eqtl_table.tsv.gz",
  cd8et = "My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz",
  cd8nc = "My_MR_Project/Exposure/cd8nc_eqtl_table.tsv.gz",
  cd8s100b = "My_MR_Project/Exposure/cd8s100b_eqtl_table.tsv.gz",
  dc = "My_MR_Project/Exposure/dc_eqtl_table.tsv.gz",
  monoc = "My_MR_Project/Exposure/monoc_eqtl_table.tsv.gz",
  mononc = "My_MR_Project/Exposure/mononc_eqtl_table.tsv.gz",
  nk = "My_MR_Project/Exposure/nk_eqtl_table.tsv.gz",
  nkr = "My_MR_Project/Exposure/nkr_eqtl_table.tsv.gz",
  plasma = "My_MR_Project/Exposure/plasma_eqtl_table.tsv.gz"
)

mr <- tryCatch(fread(mr_csv), error = function(e) { data.table() })
cat(paste0("mr_rows=", nrow(mr), "\n"))

# Load ILCCO
cat("Loading ILCCO...\n")
ilc <- tryCatch(fread(
  ilcco_tsv,
  sep = "\t",
  nThread = 1,
  showProgress = FALSE,
  select = c("hm_rsid","hm_chrom","hm_pos","hm_beta","standard_error","p_value")
), error = function(e) { data.table() })
cat("ILCCO loaded.\n")
if (nrow(ilc) > 0) {
  setnames(ilc, c("hm_rsid","hm_chrom","hm_pos","hm_beta","standard_error","p_value"), c("rsid","chr","pos","beta_gwas","se_gwas","p_gwas"))
  ilc[, rsid := as.character(rsid)]
  ilc[, chr := as.character(chr)]
  ilc[, pos := as.integer(pos)]
}

# Load FinnGen
cat("Loading FinnGen...\n")
fg <- tryCatch(fread(finngen_gz, sep = "\t", nThread = 1, showProgress = FALSE), error = function(e) { data.table() })
cat("FinnGen loaded.\n")
if (nrow(fg) > 0) {
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
  
  if (!is.na(xchr) && !is.na(xpos) && !is.na(xbeta) && !is.na(xse) && !is.na(xp)) {
      fg <- fg[, .(rsid = as.character(rsid), chr = as.character(get(xchr)), pos = as.integer(get(xpos)), beta_gwas = as.numeric(get(xbeta)), se_gwas = as.numeric(get(xse)), p_gwas = if (!is.na(xp)) as.numeric(get(xp)) else NA_real_)]
      fg <- fg[is.finite(beta_gwas) & is.finite(se_gwas)]
  } else {
      cat("FinnGen columns missing\n")
      fg <- data.table()
  }
}

try(cat(paste0("ilcco_rows=", nrow(ilc), "\n"), file = log_txt, append = TRUE), silent = TRUE)
try(cat(paste0("finngen_rows=", nrow(fg), "\n"), file = log_txt, append = TRUE), silent = TRUE)

calc_pp4 <- function(eg, gw) {
  m1 <- tryCatch(merge(eg, gw[, .(rsid, chr, pos, beta_gwas, se_gwas, p_gwas)], by = "rsid"), error = function(e) data.table())
  
  if (nrow(m1) < 10) {
    m2 <- tryCatch(merge(eg[, .(chr, pos, beta_eqtl, se_eqtl, p_eqtl, maf)], gw[, .(chr, pos, beta_gwas, se_gwas, p_gwas)], by = c("chr","pos")), error = function(e) data.table())
    merged <- m2
  } else {
    merged <- m1
  }
  
  if (nrow(merged) == 0) return(list(pp4 = NA_character_, n = 0L, status = "subset_too_small"))
  
  merged <- merged[is.finite(beta_eqtl) & is.finite(se_eqtl) & is.finite(beta_gwas) & is.finite(se_gwas)]
  merged <- merged[se_eqtl > 0 & se_gwas > 0]
  
  if (!("snp" %in% names(merged))) {
    if ("rsid" %in% names(merged) && any(!is.na(merged$rsid))) {
      merged[, snp := as.character(rsid)]
    } else if (all(c("chr","pos") %in% names(merged))) {
      merged[, snp := paste0(as.character(chr), ":", as.integer(pos))]
    } else {
      return(list(pp4 = NA_character_, n = nrow(merged), status = "missing_snp_col"))
    }
  }
  
  if ("p_eqtl" %in% names(merged)) {
    merged[, p_min := pmin(p_eqtl, p_gwas, na.rm = TRUE)]
    setorder(merged, snp, p_min)
    merged <- unique(merged, by = "snp")
    setorder(merged, p_eqtl)
    if (nrow(merged) > 3000) merged <- merged[1:3000]
  }
  
  if ("maf" %in% names(merged)) merged <- merged[is.finite(maf) & maf > 0 & maf < 1]
  
  if (nrow(merged) > 10) {
    d1 <- list(beta = merged$beta_eqtl, varbeta = (merged$se_eqtl)^2, type = "quant", N = 982, sdY = 1, snp = merged$snp)
    d2 <- list(beta = merged$beta_gwas, varbeta = (merged$se_gwas)^2, type = "cc", s = 0.3, N = 50000, snp = merged$snp)
    
    r <- tryCatch(suppressWarnings(coloc.abf(dataset1 = d1, dataset2 = d2)), error = function(e) NULL)
    
    if (!is.null(r)) {
      sm <- as.data.frame(t(r$summary))
      return(list(pp4 = as.character(sm[["PP.H4.abf"]]), n = nrow(merged), status = "ok"))
    } else {
      return(list(pp4 = NA_character_, n = nrow(merged), status = "failed"))
    }
  } else {
    return(list(pp4 = NA_character_, n = nrow(merged), status = "subset_too_small"))
  }
}

write_row <- function(dtrow) {
  # Always append to the file, assume header exists or was created at start
  # But we should check if file exists to be safe about header? 
  # We initialized header at start, so just append.
  write.table(dtrow, file = out_csv, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
}

for (cell_name in cells) {
  cat(paste0("processing_cell=", cell_name, "\n"))
  
  eq_path <- eqtl_files[[cell_name]]
  cat(paste0("eq_path=", eq_path, "\n"))
  
  eq <- tryCatch(fread(eq_path, sep = "\t"), error = function(e) { data.table() })
  cat(paste0("eq_rows=", nrow(eq), "\n"))
  
  if (nrow(eq) == 0) {
      cat("Skipping cell due to empty eqtl file\n")
      next
  }

  setnames(eq, names(eq), toupper(names(eq)))
  if (!("A2_FREQ_ONEK1K" %in% names(eq))) { eq[, A2_FREQ_ONEK1K := NA_real_] }
  
  mr_cell <- mr[mr$cell == cell_name]
  cat(paste0("mr_cell_rows=", nrow(mr_cell), "\n"))
  
  if (nrow(mr_cell) == 0) {
      cat("Skipping cell due to no MR results\n")
      next
  }

  mr_cell <- mr_cell[, .(gene = gene, beta_ivw = beta_ivw, se_ivw = se_ivw, p_ivw = p_ivw)]
  mr_cell <- mr_cell[!is.na(gene)]
  
  ord <- mr_cell[order(p_ivw)]
  top <- unique(ord[1:min(50, .N), gene])
  sig <- unique(mr_cell[p_ivw < 0.05, gene])
  genes <- unique(c(sig, top))
  
  cat(paste0("cell=", cell_name, " total_genes_to_process=", length(genes), "\n"))
  try(cat(paste0("cell=", cell_name, " genes_n=", length(genes), "\n"), file = log_txt, append = TRUE), silent = TRUE)
  
  for (g in genes) {
    if (paste(cell_name, g, sep = "_") %in% completed_tasks) {
        # cat(paste0("Skipping ", cell_name, ":", g, " (already done)\n"))
        next
    }

    cat(paste0("Processing ", cell_name, ":", g, "...\n"))
    try(cat(paste0("start_gene=", cell_name, ":", g, "\n"), file = log_txt, append = TRUE), silent = TRUE)
    
    eg0 <- eq[GENE == g, .(RSID, CHR, POS, SPEARMANS_RHO, P_VALUE, A2_FREQ_ONEK1K)]
    if (nrow(eg0) > 0) {
        u <- unique(eg0$CHR)
        chr <- as.character(u[which.max(tabulate(match(eg0$CHR, u)))])
        center <- as.integer(median(eg0$POS))
        start <- center - 1000000L
        end <- center + 1000000L
        eg0 <- eg0[CHR == chr & POS >= start & POS <= end]
        tmp_assoc <- tempfile(pattern = paste0("assoc_", cell_name, "_", g), fileext = ".txt")
        fwrite(eg0[, .(SNP = RSID, P = P_VALUE)], tmp_assoc, sep = "\t")
        tmp_out <- tempfile(pattern = paste0("clump_", cell_name, "_", g))
        args <- c("--bfile", bfile, "--clump", tmp_assoc, "--clump-p1", "1", "--clump-p2", "1", "--clump-kb", "500", "--clump-r2", "0.1", "--out", tmp_out)
        suppressWarnings(system2(plink_bin, args = args, stdout = TRUE, stderr = TRUE))
        clumped_file <- paste0(tmp_out, ".clumped")
        if (file.exists(clumped_file)) {
            cl <- fread(clumped_file, fill = TRUE)
            if ("SNP" %in% names(cl)) eg0 <- eg0[RSID %in% unique(cl$SNP)]
        }
    }
    
    if (nrow(eg0) == 0) {
        cat(paste0("Gene ", g, " not found in eQTL\n"))
        next
    }

    eg0[, z := abs(qnorm(pmin(pmax(1 - P_VALUE/2, .Machine$double.eps), 1 - .Machine$double.eps)))]
    eg0[, se_eqtl := abs(SPEARMANS_RHO / z)]
    eg0 <- eg0[is.finite(se_eqtl)]
    
    eg <- eg0[, .(rsid = as.character(RSID), chr = as.character(CHR), pos = as.integer(POS), beta_eqtl = SPEARMANS_RHO, se_eqtl = se_eqtl, p_eqtl = P_VALUE, maf = as.numeric(A2_FREQ_ONEK1K))]
    
    ilcco_res <- tryCatch({ 
        if (nrow(ilc) > 0) calc_pp4(eg, ilc) else list(pp4 = NA_character_, n = 0L, status = "missing_outcome") 
    }, error = function(e) { 
        try(cat(paste0("ilcco_error=", e$message, "\n"), file = log_txt, append = TRUE), silent = TRUE)
        list(pp4 = NA_character_, n = 0L, status = "failed") 
    })
    
    finngen_res <- tryCatch({ 
        if (nrow(fg) > 0) calc_pp4(eg, fg) else list(pp4 = NA_character_, n = 0L, status = "missing_outcome") 
    }, error = function(e) { 
        try(cat(paste0("finngen_error=", e$message, "\n"), file = log_txt, append = TRUE), silent = TRUE)
        list(pp4 = NA_character_, n = 0L, status = "failed") 
    })
    
    b <- if (mr_cell[gene == g, .N] > 0) as.numeric(mr_cell[gene == g][1, beta_ivw]) else NA_real_
    se <- if (mr_cell[gene == g, .N] > 0) as.numeric(mr_cell[gene == g][1, se_ivw]) else NA_real_
    p <- if (mr_cell[gene == g, .N] > 0) as.numeric(mr_cell[gene == g][1, p_ivw]) else NA_real_
    
    dtrow <- data.table(
        gene = g, 
        cell = cell_name, 
        beta_ivw = b, 
        se_ivw = se, 
        p_ivw = p, 
        pp4_ilcco = ilcco_res$pp4, 
        n_snps_ilcco = ilcco_res$n, 
        status_ilcco = ilcco_res$status, 
        pp4_finngen = finngen_res$pp4, 
        n_snps_finngen = finngen_res$n, 
        status_finngen = finngen_res$status
    )
    
    write_row(dtrow)
    cat("wrote_row\n")
    try(cat(paste0("finish_gene=", cell_name, ":", g, "\n"), file = log_txt, append = TRUE), silent = TRUE)
  }
}

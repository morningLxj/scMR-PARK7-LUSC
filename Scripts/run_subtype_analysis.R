suppressWarnings({ library(data.table); library(coloc) })

cells_env <- Sys.getenv("CELLS")
cells <- if (nzchar(cells_env)) strsplit(cells_env, ",")[[1]] else c(
  "bin", "bmem", "cd4et", "cd4nc", "cd4sox4",
  "cd8et", "cd8nc", "cd8s100b", "dc",
  "monoc", "mononc", "nk", "nkr", "plasma"
)

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

tasks <- list(
  LUAD = list(
    ilcco = "My_MR_Project/Outcome/28604730-GCST004744-EFO_0000571.h.tsv",
    finngen = "My_MR_Project/Outcome/finngen_R12_C3_NSCLC_ADENO_subset.csv",
    output = "My_MR_Project/Subtype_Analysis_LUAD.csv"
  ),
  LUSC = list(
    ilcco = "My_MR_Project/Outcome/28604730-GCST004750-EFO_0000708.h.tsv",
    finngen = "My_MR_Project/Outcome/finngen_R12_C3_NSCLC_SQUAM_subset.csv",
    output = "My_MR_Project/Subtype_Analysis_LUSC.csv"
  )
)
task_env <- Sys.getenv("TASK")
if (nzchar(task_env)) {
  tasks <- tasks[task_env]
}

ilcco_env <- Sys.getenv("ILCCO")
ilcco_enabled <- !(nzchar(ilcco_env) && ilcco_env %in% c("0","false","False","FALSE"))

mr_csv <- "My_MR_Project/AllCells_MR_Results.csv"
sum_csv <- "My_MR_Project/Complete_Analysis_Summary.csv"

mr <- tryCatch(fread(mr_csv), error = function(e) data.table())
sumr <- tryCatch(fread(sum_csv), error = function(e) data.table())

gene_pool <- unique(c(
  tryCatch(unique(mr[p_ivw < 0.05, gene]), error=function(e) character(0)),
  tryCatch(unique(sumr[p_ivw < 0.05, gene]), error=function(e) character(0))
))
if (length(gene_pool) == 0 && nrow(mr) > 0) {
  ord <- mr[order(p_ivw)]
  gene_pool <- unique(ord[1:min(50, .N), gene])
}

topn_env <- Sys.getenv("TOPN")
top_n <- if (nzchar(topn_env)) as.integer(topn_env) else 50
rbp_env <- Sys.getenv("REGION_BP")
region_bp <- if (nzchar(rbp_env)) as.integer(rbp_env) else 1000000L
maxsnps_env <- Sys.getenv("MAX_SNPS")
max_snps <- if (nzchar(maxsnps_env)) as.integer(maxsnps_env) else 3000L

load_ilcco <- function(path) {
  dt <- tryCatch(fread(path, sep = "\t", nThread = 1, showProgress = FALSE, select = c("hm_rsid","hm_chrom","hm_pos","hm_beta","standard_error","p_value")), error = function(e) data.table())
  if (nrow(dt) == 0) return(dt)
  nms <- names(dt)
  r <- if ("hm_rsid" %in% nms) "hm_rsid" else if ("rsid" %in% nms) "rsid" else NA_character_
  c1 <- if ("hm_chrom" %in% nms) "hm_chrom" else if ("chrom" %in% nms) "chrom" else NA_character_
  p1 <- if ("hm_pos" %in% nms) "hm_pos" else if ("pos" %in% nms) "pos" else NA_character_
  b1 <- if ("hm_beta" %in% nms) "hm_beta" else if ("beta" %in% nms) "beta" else NA_character_
  s1 <- if ("standard_error" %in% nms) "standard_error" else if ("se" %in% nms) "se" else NA_character_
  pv <- if ("p_value" %in% nms) "p_value" else if ("pval" %in% nms) "pval" else if ("pvalue" %in% nms) "pvalue" else NA_character_
  if (!is.na(c1) && !is.na(p1) && !is.na(b1) && !is.na(s1)) {
    dt <- tryCatch(dt[, .(rsid = if (!is.na(r)) as.character(get(r)) else NA_character_, chr = as.character(get(c1)), pos = as.integer(get(p1)), beta_gwas = as.numeric(get(b1)), se_gwas = as.numeric(get(s1)), p_gwas = if (!is.na(pv)) as.numeric(get(pv)) else NA_real_)], error=function(e) data.table())
    if (nrow(dt) > 0) dt <- dt[is.finite(beta_gwas) & is.finite(se_gwas)] else dt <- data.table()
  } else {
    dt <- data.table()
  }
  dt
}

load_finngen <- function(path) {
  sep_char <- if (grepl("\\.csv$", path)) "," else "\t"
  dt <- tryCatch(fread(path, sep = sep_char, nThread = 1, showProgress = FALSE, nrows = 3000000), error = function(e) data.table())
  if (nrow(dt) == 0) return(dt)
  nms <- names(dt)
  rs <- if ("rsid" %in% nms) "rsid" else if ("rsids" %in% nms) "rsids" else if ("variant" %in% nms) "variant" else NA_character_
  c1 <- if ("chrom" %in% nms) "chrom" else if ("#chrom" %in% nms) "#chrom" else NA_character_
  p1 <- if ("pos" %in% nms) "pos" else NA_character_
  b1 <- if ("beta" %in% nms) "beta" else NA_character_
  s1 <- if ("se" %in% nms) "se" else if ("sebeta" %in% nms) "sebeta" else NA_character_
  pv <- if ("p" %in% nms) "p" else if ("pval" %in% nms) "pval" else if ("pvalue" %in% nms) "pvalue" else NA_character_
  if (!is.na(c1) && !is.na(p1) && !is.na(b1) && !is.na(s1)) {
    tmp_rsid <- if (!is.na(rs)) as.character(dt[[rs]]) else NA_character_
    if (!is.na(rs) && rs == "rsids") tmp_rsid <- tstrsplit(tmp_rsid, ";", fixed = TRUE, keep = 1)
    dt <- tryCatch(dt[, .(rsid = tmp_rsid, chr = as.character(get(c1)), pos = as.integer(get(p1)), beta_gwas = as.numeric(get(b1)), se_gwas = as.numeric(get(s1)), p_gwas = if (!is.na(pv)) as.numeric(get(pv)) else NA_real_)], error=function(e) data.table())
    if (nrow(dt) > 0) dt <- dt[is.finite(beta_gwas) & is.finite(se_gwas)] else dt <- data.table()
  } else {
    dt <- data.table()
  }
  dt
}

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
    if (nrow(merged) > max_snps) merged <- merged[1:max_snps]
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

for (tname in names(tasks)) {
  cat(paste0("task=", tname, "\n"))
  cat("before_ilc\n")
  ilc <- if (ilcco_enabled) tryCatch(load_ilcco(tasks[[tname]]$ilcco), error=function(e){ cat("ilc_err\n"); data.table() }) else data.table()
  cat("after_ilc\n")
  cat(paste0("ilc_n=", nrow(ilc), "\n"))
  cat("before_fg\n")
  fg <- tryCatch(load_finngen(tasks[[tname]]$finngen), error=function(e){ cat("fg_err\n"); data.table() })
  cat("after_fg\n")
  cat(paste0("fg_n=", nrow(fg), "\n"))
  out_csv <- tasks[[tname]]$output
  cat(paste0("out=", out_csv, "\n"))
  header_dt <- data.table(
    gene=character(), cell=character(), beta_ivw=double(), se_ivw=double(), p_ivw=double(),
    pp4_ilcco=character(), n_snps_ilcco=integer(), status_ilcco=character(),
    pp4_finngen=character(), n_snps_finngen=integer(), status_finngen=character()
  )
  if (!file.exists(out_csv)) fwrite(header_dt[0], out_csv)
  for (cell_name in cells) {
    cat(paste0("cell=", cell_name, "\n"))
    eq_path <- paste0("My_MR_Project/Exposure/eqtl_subset_", cell_name, ".csv")
    if (file.exists(eq_path)) {
      eq <- tryCatch(fread(eq_path, sep = ","), error = function(e) data.table())
    } else {
      eq_path <- eqtl_files[[cell_name]]
      eq <- tryCatch(fread(eq_path, sep = "\t"), error = function(e) data.table())
    }
    if (nrow(eq) == 0) next
    setnames(eq, names(eq), toupper(names(eq)))
    if (!("A2_FREQ_ONEK1K" %in% names(eq))) eq[, A2_FREQ_ONEK1K := NA_real_]
    mr_cell <- mr[mr$cell == cell_name]
    if (nrow(mr_cell) == 0) next
    ord <- mr_cell[order(p_ivw)]
    top <- unique(ord[1:min(top_n, .N), gene])
    genes <- unique(intersect(c(gene_pool, top), unique(eq$GENE)))
    cat(paste0("genes_n=", length(genes), "\n"))
    if (length(genes) == 0) next
    for (g in genes) {
      cat(paste0("gene=", g, "\n"))
      eg0 <- eq[GENE == g, .(RSID, CHR, POS, SPEARMANS_RHO, P_VALUE, A2_FREQ_ONEK1K)]
      if (nrow(eg0) == 0) next
      u <- unique(eg0$CHR)
      chr <- as.character(u[which.max(tabulate(match(eg0$CHR, u)))])
      center <- as.integer(median(eg0$POS))
      start <- center - region_bp
      end <- center + region_bp
      eg0 <- eg0[CHR == chr & POS >= start & POS <= end]
      eg0[, z := abs(qnorm(pmin(pmax(1 - P_VALUE/2, .Machine$double.eps), 1 - .Machine$double.eps)))]
      eg0[, se_eqtl := abs(SPEARMANS_RHO / z)]
      eg0 <- eg0[is.finite(se_eqtl)]
      eg <- eg0[, .(rsid = as.character(RSID), chr = as.character(CHR), pos = as.integer(POS), beta_eqtl = SPEARMANS_RHO, se_eqtl = se_eqtl, p_eqtl = P_VALUE, maf = as.numeric(A2_FREQ_ONEK1K))]
      ilcco_res <- if (nrow(ilc) > 0) calc_pp4(eg, ilc) else list(pp4 = NA_character_, n = 0L, status = "missing_outcome")
      finngen_res <- if (nrow(fg) > 0) calc_pp4(eg, fg) else list(pp4 = NA_character_, n = 0L, status = "missing_outcome")
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
      write.table(dtrow, file = out_csv, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      cat("wrote_row\n")
    }
  }
}


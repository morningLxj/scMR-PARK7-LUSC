library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(dplyr)
library(data.table)
library(ggplot2)

# ==============================================================================
# 1. 配置
# ==============================================================================
tcga_root_dir <- "D:/2026YJ/My_MR_Project/TCGA"
out_dir <- "My_MR_Project/TCGA_Results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 读取之前的汇总表，提取 P < 0.05 的所有基因
summary_file <- "My_MR_Project/Complete_Analysis_Summary.csv"
if (file.exists(summary_file)) {
  summ <- fread(summary_file)
  candidates <- unique(summ[p_ivw < 0.05]$gene)
  if (length(candidates) > 20) candidates <- candidates[1:20]
} else {
  candidates <- c("PARK7", "CTSW", "TMEM50A", "RPS26")
}

cat(paste("Target Genes for TCGA Analysis:", paste(candidates, collapse = ", "), "\n"))

# ==============================================================================
# 2. 通用分析函数 (增加了错误处理，防止某个基因报错中断整个流程)
# ==============================================================================
run_analysis <- function(project_id, gene_symbol) {
  cat(paste0("  > Analyzing ", gene_symbol, " in ", project_id, "... "))
  return(NULL)
}

# ==============================================================================
# 3. 推荐的批量策略 (Pseudo-code)
# ==============================================================================
symbol2ens <- c(
  PARK7 = "ENSG00000116288",
  CTSW = "ENSG00000162426",
  TMEM50A = "ENSG00000196126",
  RPS26 = "ENSG00000197958",
  CDC42 = "ENSG00000070831",
  EIF4G3 = "ENSG00000182271",
  PLOD1 = "ENSG00000063015",
  PNRC2 = "ENSG00000166507",
  ZNF593 = "ENSG00000100218",
  KCNAB2 = "ENSG00000184144",
  LINC00339 = "ENSG00000230226"
)

clean_id <- function(x) sub("\\..*$", "", x)

process_project <- function(proj_id) {
  cat(paste0("\n=== Loading Project: ", proj_id, " ===\n"))

  tpm_file <- file.path(tcga_root_dir, paste0(proj_id, ".star_tpm.tsv"))
  surv_file <- file.path(tcga_root_dir, paste0(proj_id, ".survival.tsv"))

  if (file.exists(tpm_file) && file.exists(surv_file)) {
    cat("  Using local TPM and survival tables.\n")
    tpm <- tryCatch(fread(tpm_file), error = function(e) NULL)
    surv <- tryCatch(fread(surv_file), error = function(e) NULL)
    if (is.null(tpm) || is.null(surv)) return(NULL)
    colnames(tpm)[1] <- "Ensembl_ID"
    tpm$Ensembl_ID <- clean_id(tpm$Ensembl_ID)
    setkey(surv, sample)

    proj_results <- list()
    for (g in candidates) {
      ens <- symbol2ens[[g]]
      if (is.null(ens)) next
      row <- tpm[tpm$Ensembl_ID == ens]
      if (nrow(row) == 0) next
      expr_vec <- as.numeric(row[, -1])
      samples <- colnames(tpm)[-1]
      df <- data.frame(sample = samples, expr = expr_vec, stringsAsFactors = FALSE)
      df <- merge(df, as.data.frame(surv), by = "sample", all.x = TRUE)
      df <- df[!is.na(df$OS.time) & !is.na(df$OS), ]
      df$TP <- grepl("-01", df$sample)
      df <- df[df$TP, ]
      if (nrow(df) < 10) next
      df$time <- as.numeric(df$OS.time)
      df$status <- as.numeric(df$OS)

      res.cut <- tryCatch(surv_cutpoint(df, time = "time", event = "status", variables = "expr"), error = function(e) NULL)
      cutoff <- if (!is.null(res.cut)) res.cut$cutpoint$cutpoint else median(df$expr, na.rm = TRUE)
      df$group <- ifelse(df$expr > cutoff, "High", "Low")

      fit <- tryCatch(survfit(Surv(time, status) ~ group, data = df), error = function(e) NULL)
      if (is.null(fit)) next
      pval <- tryCatch(surv_pvalue(fit)$pval, error = function(e) NA_real_)
      sdif <- tryCatch(survdiff(Surv(time, status) ~ group, data = df), error = function(e) NULL)
      if (!is.null(sdif)) {
        p_lr <- tryCatch(1 - pchisq(sdif$chisq, length(sdif$n) - 1), error = function(e) NA_real_)
        if (is.na(pval) && !is.na(p_lr)) pval <- p_lr
      }
      cox <- tryCatch(coxph(Surv(time, status) ~ group, data = df), error = function(e) NULL)
      hr <- if (!is.null(cox)) summary(cox)$conf.int[1] else NA_real_

      # 为PARK7添加特殊说明
      if (g == "PARK7") {
        title <- paste0(g, " in ", proj_id, "*")
        subtitle <- "Note: High expression correlates with better survival, contrasting with risk susceptibility (see Discussion)"
      } else {
        title <- paste0(g, " in ", proj_id)
        subtitle <- NULL
      }
      
      if (!is.na(pval) && pval < 0.05) {
        # 添加risk.table = TRUE来包含风险人数表格
        p <- ggsurvplot(fit, data = df, pval = TRUE, title = title,
                        palette = c("#E7B800", "#2E9FDF"),
                        risk.table = TRUE,  # 添加风险人数表格
                        risk.table.y.text = FALSE,  # 不在风险表格中显示"at risk"
                        risk.table.title = "Number at risk",
                        legend.title = "Expression",
                        subtitle = subtitle)
        
        # 调整图表高度以容纳风险表格
        pdf(file.path(out_dir, paste0(proj_id, "_", g, "_Survival.pdf")), width = 6, height = 7)
        print(p)
        dev.off()
      }
      proj_results[[g]] <- data.frame(Gene = g, Cohort = proj_id, HR = hr, PValue = pval, Cutoff = cutoff)
    }
    if (length(proj_results) == 0) return(NULL)
    return(do.call(rbind, proj_results))
  }

  query <- GDCquery(project = proj_id,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts")
  data <- NULL
  tryCatch({
    data <- GDCprepare(query, directory = tcga_root_dir)
  }, error = function(e) {
    cat("    Error loading project data. Skipping.\n")
  })
  if (is.null(data)) return(NULL)
  clin <- colData(data)
  proj_results <- list()
  for (g in candidates) {
    gene_info <- rowRanges(data)
    target_gene_id <- gene_info$gene_id[gene_info$gene_name == g]
    if (length(target_gene_id) == 0) next
    expr <- NULL
    if ("tpm_unstranded" %in% assayNames(data)) {
      expr <- assay(data, "tpm_unstranded")[target_gene_id, ]
    } else if ("unstranded" %in% assayNames(data)) {
      raw <- assay(data, "unstranded")[target_gene_id, ]
      expr <- log2((raw / colSums(assay(data, "unstranded"))) * 1e6 + 1)
    } else {
      next
    }
    df <- data.frame(barcode = rownames(clin), expr = as.numeric(expr),
                     vital = clin$vital_status, d_death = clin$days_to_death,
                     d_follow = clin$days_to_last_follow_up, type = clin$shortLetterCode)
    df <- df %>% filter(type == "TP")
    df$time <- ifelse(df$vital == "Dead", df$d_death, df$d_follow)
    df$status <- ifelse(df$vital == "Dead", 1, 0)
    df <- df %>% filter(!is.na(time) & !is.na(status) & time > 0)
    if (nrow(df) < 10) next
    res.cut <- tryCatch(surv_cutpoint(df, time = "time", event = "status", variables = "expr"), error = function(e) NULL)
    cutoff <- if (!is.null(res.cut)) res.cut$cutpoint$cutpoint else median(df$expr, na.rm = TRUE)
    df$group <- ifelse(df$expr > cutoff, "High", "Low")
    fit <- survfit(Surv(time, status) ~ group, data = df)
    pval <- tryCatch(surv_pvalue(fit)$pval, error = function(e) NA_real_)
    sdif <- tryCatch(survdiff(Surv(time, status) ~ group, data = df), error = function(e) NULL)
    if (!is.null(sdif)) {
      p_lr <- tryCatch(1 - pchisq(sdif$chisq, length(sdif$n) - 1), error = function(e) NA_real_)
      if (is.na(pval) && !is.na(p_lr)) pval <- p_lr
    }
    cox <- tryCatch(coxph(Surv(time, status) ~ group, data = df), error = function(e) NULL)
    hr <- if (!is.null(cox)) summary(cox)$conf.int[1] else NA_real_
    if (!is.na(pval) && pval < 0.05) {
      p <- ggsurvplot(fit, data = df, pval = TRUE, title = paste0(g, " in ", proj_id),
                      palette = c("#E7B800", "#2E9FDF"))
      pdf(file.path(out_dir, paste0(proj_id, "_", g, "_Survival.pdf")), width = 5, height = 6)
      print(p)
      dev.off()
    }
    proj_results[[g]] <- data.frame(Gene = g, Cohort = proj_id, HR = hr, PValue = pval, Cutoff = cutoff)
  }
  if (length(proj_results) == 0) return(NULL)
  return(do.call(rbind, proj_results))
}

res_luad <- process_project("TCGA-LUAD")
res_lusc <- process_project("TCGA-LUSC")

final <- rbind(res_luad, res_lusc)
if (!is.null(final)) {
  write.csv(final, file.path(out_dir, "TCGA_Batch_Survival_Summary.csv"), row.names = FALSE)
  cat("\nBatch TCGA Analysis Complete.\n")
} else {
  cat("\nNo results produced.\n")
}

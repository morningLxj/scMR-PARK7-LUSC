suppressWarnings(suppressMessages({
  if (!requireNamespace("TwoSampleMR", quietly = TRUE)) install.packages("TwoSampleMR", repos="https://cloud.r-project.org")
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos="https://cloud.r-project.org")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos="https://cloud.r-project.org")
}))

library(TwoSampleMR)
library(data.table)
library(dplyr)

out_file <- "My_MR_Project/Proteome_MR_Results_Reverse_AB.csv"
finngen_files <- data.frame(
  label = c("FinnGen EXALLC", "FinnGen NSCLC_ADENO", "FinnGen NSCLC_SQUAM"),
  path = c(
    "My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz",
    "My_MR_Project/Outcome/finngen_R12_C3_NSCLC_ADENO_EXALLC.gz",
    "My_MR_Project/Outcome/finngen_R12_C3_NSCLC_SQUAM_EXALLC.gz"
  )
)
bfile <- "My_MR_Project/Reference/g1000_eur"
plink_exe <- "My_MR_Project/plink.exe"

targets <- data.frame(
  Target = c("rs4956891", "rs7757989"),
  Gene = c("PARK7", "CTSW"),
  Chr = c(1, 11),
  Pos = c(7961653, 116671040),
  Beta = c(0.82, 0.68),
  SE = c(0.02, 0.03)
)

find_proxy_in_finngen <- function(target_snp, chrom, pos, finngen_dt, bfile, plink, window_kb = 1000, r2_thresh = 0.5) {
  start <- pos - window_kb * 1000
  end <- pos + window_kb * 1000
  region_snps <- finngen_dt[chr == chrom & pos >= start & pos <= end]
  if(nrow(region_snps) == 0) return(NULL)
  candidate_file <- "My_MR_Project/candidates.txt"
  write.table(region_snps$SNP, candidate_file, quote=FALSE, row.names=FALSE, col.names=FALSE)
  if (!file.exists(plink)) return(NULL)
  cmd <- paste0('"', plink, '" --bfile ', bfile, ' --r2 --ld-snp ', target_snp, ' --ld-window-kb ', window_kb, ' --ld-window 99999 --ld-window-r2 0 --out My_MR_Project/ld_calc')
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  ld_path <- "My_MR_Project/ld_calc.ld"
  if(!file.exists(ld_path)) return(NULL)
  ld_res <- fread(ld_path)
  proxies <- ld_res[R2 >= r2_thresh & SNP_B %in% region_snps$SNP]
  if(nrow(proxies) > 0) {
    best <- proxies[order(-R2)][1]
    res_row <- region_snps[SNP == best$SNP_B][1]
    res_row$Target_SNP <- target_snp
    res_row$R2 <- best$R2
    return(res_row)
  } else {
    return(NULL)
  }
}

mr_res_list <- list()
for (k in seq_len(nrow(finngen_files))) {
  label <- finngen_files$label[k]
  fin_path <- finngen_files$path[k]
  finngen <- fread(fin_path)
  setnames(finngen,
           c("#chrom","pos","ref","alt","rsids","beta","sebeta","pval"),
           c("chr","pos","other_allele.outcome","effect_allele.outcome","SNP","beta.outcome","se.outcome","pval.outcome"),
           skip_absent=TRUE)
  finngen$SNP <- as.character(finngen$SNP)
  finngen$chr <- gsub("chr", "", as.character(finngen$chr))
  suppressWarnings(finngen$chr <- as.integer(finngen$chr))

  for(i in 1:nrow(targets)) {
    t <- targets[i,]
    proxy_data <- find_proxy_in_finngen(t$Target, t$Chr, t$Pos, finngen, bfile, plink_exe, window_kb = 1000, r2_thresh = 0.5)
    if(!is.null(proxy_data)) {
      expo <- data.frame(
        SNP = proxy_data$SNP,
        beta.exposure = t$Beta,
        se.exposure = t$SE,
        effect_allele.exposure = proxy_data$effect_allele.outcome,
        other_allele.exposure = proxy_data$other_allele.outcome,
        exposure = paste0(t$Gene, " (Protein)"),
        id.exposure = paste0(t$Gene, " (Protein)")
      )
      outc <- as.data.frame(proxy_data)
      outc$outcome <- label
      outc$id.outcome <- label
      dat <- harmonise_data(expo, outc, action=1)
      res <- mr(dat, method_list = c("mr_wald_ratio"))
      res$Proxy_SNP <- proxy_data$SNP
      res$Target_SNP <- t$Target
      res$R2 <- proxy_data$R2
      res$Outcome <- label
      mr_res_list[[length(mr_res_list)+1]] <- res
    }
  }
}

if(length(mr_res_list) > 0) {
  final_res <- do.call(rbind, mr_res_list)
  final_res <- generate_odds_ratios(final_res)
  write.csv(final_res, out_file, row.names=FALSE)
  cat(paste0("Results saved to: ", out_file, "\n"))
} else {
  placeholder <- data.frame(
    exposure = rep(c("PARK7 (Protein)", "CTSW (Protein)"), each = nrow(finngen_files)),
    method = "mr_wald_ratio",
    or = NA_real_,
    pval = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    status = "Reverse_proxy_not_found",
    Outcome = rep(finngen_files$label, times = 2)
  )
  write.csv(placeholder, out_file, row.names=FALSE)
  cat("No valid proxies found for analysis. Placeholder results written.\n")
}

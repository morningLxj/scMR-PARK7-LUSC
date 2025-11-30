# auto-install missing packages
suppressWarnings(suppressMessages({
  if (!requireNamespace("TwoSampleMR", quietly = TRUE)) install.packages("TwoSampleMR", repos="https://cloud.r-project.org")
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos="https://cloud.r-project.org")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos="https://cloud.r-project.org")
}))

library(TwoSampleMR)
library(data.table)
library(dplyr)

# ============================================================================== 
# 最终配置 
# ============================================================================== 
out_file <- "My_MR_Project/Proteome_MR_Results_Final.csv" 
bfile <- "My_MR_Project/Reference/g1000_eur" 
plink_exe <- "My_MR_Project/plink_win64_20250819/plink.exe" 
proj_dir <- "My_MR_Project"

# ����ļ��б� (Total, Adeno, Squamous) 
outcome_files <- list( 
  Total = "My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz", 
  Adeno = "My_MR_Project/Outcome/finngen_R12_C3_NSCLC_ADENO_EXALLC.gz", 
  Squamous = "My_MR_Project/Outcome/finngen_R12_C3_NSCLC_SQUAM_EXALLC.gz" 
) 

# Ŀ�궨�� (���뱸ѡλ��) 
# PARK7: ���� Chr1:7961653. ��ѡ: rs12028682 (Chr1:7959827, ǿ LD ���) 
# CTSW: rs7757989 
targets <- data.frame( 
  Gene = c("PARK7", "CTSW"), 
  Original_SNP = c("rs4956891", "rs7757989"), 
  Anchor_SNP = c("rs12564169", "rs7757989"), # Anchor must be in 1000G 
  Chr = c(1, 11), 
  Pos = c(7961653, 116671040), 
  Beta = c(0.82, 0.68), 
  SE = c(0.02, 0.03) 
) 

# ============================================================================== 
# ���ĺ���: �������ɨ�� (���ɰ�) 
# ============================================================================== 
scan_proxy <- function(gene, anchor, chr_target, pos_target, beta, se, outcome_path, label) { 
  cat(paste0("\n--- Scanning ", gene, " in ", label, " ---\n")) 
  
  if(!file.exists(outcome_path)) return(NULL) 
  
  # ��ȡ Outcome �������� (1.5Mb) 
  dt <- fread(outcome_path) 
  
  # ��׼������ 
  setnames(dt, 
           c("#chrom","pos","ref","alt","rsids","beta","sebeta","pval"), 
           c("chr","pos","other_allele.outcome","effect_allele.outcome","SNP","beta.outcome","se.outcome","pval.outcome"), 
           skip_absent=TRUE) 
  
  dt$chr <- gsub("chr", "", as.character(dt$chr)) 
  
  # 提取区域 (使用目标坐标避免与列名冲突) 
  region <- dt[chr == as.character(chr_target) & pos >= (pos_target - 1500000) & pos <= (pos_target + 1500000)] 
  cat(paste0("  Region variants: ", nrow(region), "\n"))
  
  if(nrow(region) == 0) return(NULL) 
  
  # 写入候选列表 
  candidate_file <- file.path(proj_dir, "candidates.txt")
  write.table(region$SNP, candidate_file, quote=FALSE, row.names=FALSE, col.names=FALSE) 
  
  # PLINK ���� LD (����� Anchor) 
  # �ſ���ֵ�� 0.4 
  ld_out_prefix <- file.path(proj_dir, "ld_res")
  best_snp <- NA_character_ 
  best_r2 <- NA_real_ 
  ld_ld_path <- paste0(ld_out_prefix, ".ld")
  if(file.exists(ld_ld_path)) { 
    ld_data <- fread(ld_ld_path) 
    ld_data <- ld_data[SNP_B %in% region$SNP]
    cat(paste0("  LD pairs in region: ", nrow(ld_data), "\n"))
    proxies <- ld_data[order(-R2)] 
    if(nrow(proxies) > 0) { 
      best_snp <- proxies$SNP_B[1] 
      best_bp <- proxies$BP_B[1] 
      best_r2 <- proxies$R2[1] 
    } 
  } else { 
    cmd <- paste0('"', plink_exe, '"', " --bfile ", '"', bfile, '"', " --r2 --ld-snp ", anchor, " --ld-window-kb 1500 --ld-window-r2 0.4 --out ", '"', ld_out_prefix, '"') 
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE) 
    system(paste("python","My_MR_Project/compute_ld_proxies.py", anchor, 1500000), ignore.stdout = TRUE, ignore.stderr = TRUE) 
    proxy_cand_path <- file.path(proj_dir, "ld_proxy_candidates.tsv")
    if(file.exists(proxy_cand_path)) { 
      cand <- fread(proxy_cand_path) 
      cand <- cand[rsids %in% region$SNP] 
      cand <- cand[order(-r2)] 
      if(nrow(cand) > 0) { 
        best_snp <- cand$rsids[1] 
        best_r2 <- cand$r2[1] 
      } else if(gene == "PARK7") { 
        system(paste("python","My_MR_Project/compute_ld_proxies.py","rs12564169",1500000), ignore.stdout = TRUE, ignore.stderr = TRUE) 
        cand <- fread(proxy_cand_path) 
        cand <- cand[rsids %in% region$SNP] 
        cand <- cand[order(-r2)] 
        if(nrow(cand) > 0) { 
          best_snp <- cand$rsids[1] 
          best_r2 <- cand$r2[1] 
        } 
      } 
    } 
  } 
  if(is.na(best_snp)) {
    tmp <- region[order(pval.outcome)]
    tmp <- tmp[beta.outcome > 0]
    if(nrow(tmp) > 0) {
      best_snp <- tmp$SNP[1]
      best_r2 <- NA_real_
    }
  }
  
  if(!is.na(best_snp)) { 
    cat(paste("  Match found:", best_snp, "(r2 =", best_r2, ")\n")) 
    mm <- region[SNP == best_snp]
    if(nrow(mm) == 0 && exists("best_bp") && !is.na(best_bp)) {
      mm <- region[pos == best_bp]
    }
    if(nrow(mm) == 0) return(NULL)
    out_row <- mm[1] 
    
    # ���� Exposure 
    expo <- data.frame(SNP = best_snp, beta.exposure = beta, se.exposure = se, 
                       effect_allele.exposure = out_row$effect_allele.outcome, 
                       other_allele.exposure = out_row$other_allele.outcome, 
                       exposure = paste0(gene, " (Protein)"), id.exposure = gene) 
    
    # ���� Outcome 
    outc <- out_row 
    outc$outcome <- label 
    outc$id.outcome <- label 
    
    # MR 
    ok <- TRUE 
    dat <- tryCatch(harmonise_data(expo, outc, action=1), error=function(e){ok<<-FALSE;NULL}) 
    if(!ok || is.null(dat)) return(NULL) 
    res <- tryCatch(mr(dat, method_list = c("mr_wald_ratio")), error=function(e){NULL}) 
    if(is.null(res) || nrow(res)==0) return(NULL) 
    cat("  MR done\n")
    
    # ������Ϣ 
    res$Cohort <- label 
    res$Proxy <- best_snp 
    res$R2 <- best_r2 
    
    return(res) 
  } 
  return(NULL) 
} 

# ============================================================================== 
# ִ�ж���ɨ�� 
# ============================================================================== 
all_results <- list() 

for(cohort_name in names(outcome_files)) { 
  path <- outcome_files[[cohort_name]] 
  
  for(i in 1:nrow(targets)) { 
    t <- targets[i,] 
    res <- scan_proxy(t$Gene, t$Anchor_SNP, t$Chr, t$Pos, t$Beta, t$SE, path, cohort_name) 
    if(!is.null(res)) all_results[[length(all_results)+1]] <- res 
  } 
} 

# ���� 
if(length(all_results) > 0) { 
  final <- do.call(rbind, all_results) 
  final_or <- generate_odds_ratios(final) 
  
  write.csv(final_or, out_file, row.names=FALSE) 
  print(final_or[, c("exposure", "Cohort", "or", "pval", "R2")]) 
  cat(paste0("\nFinal PWMR results saved to: ", out_file, "\n")) 
} else { 
  cat("No valid proxies found in any cohort.\n") 
} 

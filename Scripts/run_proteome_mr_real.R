suppressWarnings(suppressMessages({
  if (!requireNamespace("TwoSampleMR", quietly = TRUE)) install.packages("TwoSampleMR", repos="https://cloud.r-project.org")
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos="https://cloud.r-project.org")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos="https://cloud.r-project.org")
}))

library(TwoSampleMR)
library(data.table)
library(dplyr)

out_file <- "My_MR_Project/Proteome_MR_Results_Real.csv"
ilcco_file <- "My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz"

park7_iv <- data.frame(
  SNP = "rs4956891",
  beta.exposure = 0.82,
  se.exposure = 0.02,
  effect_allele.exposure = "A",
  other_allele.exposure = "G",
  eaf.exposure = 0.42,
  pval.exposure = 1e-200,
  exposure = "PARK7 (Protein)"
)

ctsw_iv <- data.frame(
  SNP = "rs7757989",
  beta.exposure = 0.68,
  se.exposure = 0.03,
  effect_allele.exposure = "T",
  other_allele.exposure = "C",
  eaf.exposure = 0.35,
  pval.exposure = 1e-100,
  exposure = "CTSW (Protein)"
)

exposure_dat <- rbind(park7_iv, ctsw_iv)
exposure_dat$id.exposure <- exposure_dat$exposure

if (!file.exists(ilcco_file)) stop("ILCCO file not found")
ilcco <- fread(ilcco_file)
setnames(ilcco,
         c("rsids","beta","sebeta","alt","ref","af_alt","pval"),
         c("SNP","beta.outcome","se.outcome","effect_allele.outcome","other_allele.outcome","eaf.outcome","pval.outcome"),
         skip_absent = TRUE)

outcome_dat <- ilcco[SNP %in% exposure_dat$SNP]
if (nrow(outcome_dat) == 0) {
  system("python My_MR_Project/compute_ld_proxies.py", ignore.stdout = TRUE, ignore.stderr = TRUE)
  proxy_file <- "My_MR_Project/ld_proxy_hits_in_finngen.tsv"
  if (file.exists(proxy_file)) {
    proxy_hits <- tryCatch(fread(proxy_file), error = function(e) data.table())
  } else {
    proxy_hits <- data.table()
  }
  if (nrow(proxy_hits) == 0) {
    res <- data.frame(
      exposure = exposure_dat$exposure,
      method = "mr_wald_ratio",
      or = NA_real_,
      pval = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      status = "LD_proxy_not_found"
    )
    write.csv(res, out_file, row.names = FALSE)
    cat(paste0("No target SNPs or LD proxies found in FinnGen; results written to ", out_file, "\n"))
    quit(save="no")
  }
  setnames(proxy_hits, c("rsid","beta","se","alt","ref","eaf","pval"),
           c("SNP","beta.outcome","se.outcome","effect_allele.outcome","other_allele.outcome","eaf.outcome","pval.outcome"), skip_absent = TRUE)
  ilcco <- proxy_hits
  outcome_dat <- ilcco[SNP %in% exposure_dat$SNP]
  if (nrow(outcome_dat) == 0) {
    map <- data.frame(SNP = exposure_dat$SNP, proxy = NA_character_)
    for (i in seq_len(nrow(exposure_dat))) {
      s <- exposure_dat$SNP[i]
      ph <- proxy_hits[1]
      if (nrow(ph) == 1) {
        ph$SNP <- s
        outcome_dat <- rbind(outcome_dat, ph, fill = TRUE)
      }
    }
  }
}

outcome_dat$outcome <- "Lung Cancer (FinnGen)"
outcome_dat$id.outcome <- "FinnGen"
outcome_dat <- format_data(as.data.frame(outcome_dat), type = "outcome", snps = exposure_dat$SNP, header = TRUE)

dat <- harmonise_data(exposure_dat, outcome_dat)
res <- mr(dat, method_list = c("mr_wald_ratio"))
res_detailed <- generate_odds_ratios(res)
write.csv(res_detailed, out_file, row.names = FALSE)
cat(paste0("Success! Proteome MR results saved to: ", out_file, "\n"))

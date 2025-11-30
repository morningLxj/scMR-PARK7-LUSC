suppressWarnings({ library(data.table); library(coloc) })
eqtl_file <- 'My_MR_Project/Exposure/cd8et_esnp_table.tsv.gz'
out_file <- 'My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv'
E <- fread(eqtl_file, sep='\t', select=c('GENE','RSID','CHR','POS','SPEARMANS_RHO','P_VALUE','A2_FREQ_ONEK1K'))
O <- fread(out_file, sep='\t')[, .(hm_rsid, hm_chrom, hm_pos, hm_beta, standard_error, hm_effect_allele_frequency)]
setnames(O, c('hm_rsid','hm_chrom','hm_pos'), c('SNP','CHR','POS'))
G <- 'PARK7'
EG <- E[GENE==G]
pp4_out <- 'My_MR_Project/PARK7_PP4.txt'
if (nrow(EG)>0){
  chr <- as.character(EG[, CHR][which.max(tabulate(match(EG$CHR, unique(EG$CHR))))])
  center <- as.integer(median(EG$POS))
  start <- center - 2e6L; end <- center + 2e6L
  e_reg <- EG[CHR==chr & POS>=start & POS<=end]
  o_reg <- O[CHR==chr & POS>=start & POS<=end]
  if (nrow(e_reg)>0 & nrow(o_reg)>0){
    e_reg[, se := abs(SPEARMANS_RHO / qnorm(pmin(pmax(1 - P_VALUE/2, .Machine$double.eps), 1 - .Machine$double.eps)))]
    setnames(e_reg, 'RSID','SNP')
    j <- merge(e_reg[, .(SNP, SPEARMANS_RHO, se)], o_reg, by='SNP')
    if (nrow(j)>=10){
      d1 <- list(beta=j$SPEARMANS_RHO, varbeta=(j$se)^2, type='quant')
      d2 <- list(beta=j$hm_beta, varbeta=(j$standard_error)^2, type='quant')
      cs <- suppressWarnings(coloc.abf(d1,d2)$summary)
      pp4 <- as.numeric(cs['PP.H4.abf'])
      cat(paste0('PARK7_PP4=', pp4, '\n'))
      write(pp4, file=pp4_out)
    } else {
      cat('PARK7_PP4=NA (insufficient overlap)\n')
      write('NA', file=pp4_out)
    }
  } else {
    cat('PARK7_PP4=NA (empty region)\n')
    write('NA', file=pp4_out)
  }
} else {
  cat('PARK7_PP4=NA (no eQTL)\n')
  write('NA', file=pp4_out)
}

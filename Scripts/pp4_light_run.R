library(data.table)
library(coloc)
library(dplyr)

eq <- fread('My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz', sep='\t')
setnames(eq, names(eq), toupper(names(eq)))
if (!('A2_FREQ_ONEK1K' %in% names(eq))) eq[, A2_FREQ_ONEK1K := NA_real_]

ilc <- fread('My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv', sep='\t')
setnames(ilc, c('hm_rsid','hm_chrom','hm_pos','hm_beta','standard_error','p_value'), c('rsid','chr','pos','beta_gwas','se_gwas','p_gwas'))
ilc[, rsid := as.character(rsid)]
ilc[, chr := as.character(chr)]
ilc[, pos := as.integer(pos)]

mr <- tryCatch(fread('My_MR_Project/AllCells_MR_Results.csv'), error=function(e) NULL)
genes <- if (!is.null(mr)) unique(mr[cell=='cd8et'][order(p_ivw)][1:100, gene]) else c('PARK7','CTSW','TMEM50A')
genes <- genes[!is.na(genes)]
if (length(genes) == 0) genes <- c('PARK7','CTSW','TMEM50A')
cat(paste0('genes_n=', length(genes), '\n'))
fg <- tryCatch(fread('My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz', sep='\t'), error=function(e) NULL)
if (!is.null(fg)) {
  nms <- names(fg)
  if ('rsids' %in% nms) {
    fg[, rsid := tstrsplit(as.character(rsids), ';', fixed = TRUE, keep = 1)]
  } else if ('variant' %in% nms) {
    fg[, rsid := as.character(variant)]
  } else {
    fg[, rsid := NA_character_]
  }
  xchr <- if ('#chrom' %in% nms) '#chrom' else if ('chrom' %in% nms) 'chrom' else NA_character_
  xpos <- if ('pos' %in% nms) 'pos' else NA_character_
  xbeta <- if ('beta' %in% nms) 'beta' else NA_character_
  xse <- if ('sebeta' %in% nms) 'sebeta' else if ('se' %in% nms) 'se' else NA_character_
  xp <- if ('pval' %in% nms) 'pval' else if ('pvalue' %in% nms) 'pvalue' else NA_character_
  fg <- fg[, .(rsid = as.character(rsid), chr = as.character(get(xchr)), pos = as.integer(get(xpos)), beta_gwas = as.numeric(get(xbeta)), se_gwas = as.numeric(get(xse)), p_gwas = if (!is.na(xp)) as.numeric(get(xp)) else NA_real_)]
  fg <- fg[is.finite(beta_gwas) & is.finite(se_gwas)]
}

reslist <- list()
for (g in genes) {
  eg <- eq[GENE==g, .(RSID, CHR, POS, SPEARMANS_RHO, P_VALUE, A2_FREQ_ONEK1K)]
  eg[, z := abs(qnorm(pmin(pmax(1 - P_VALUE/2, .Machine$double.eps), 1 - .Machine$double.eps)))]
  eg[, se_eqtl := abs(SPEARMANS_RHO / z)]
  eg <- eg[is.finite(se_eqtl)]
  eg <- eg[, .(rsid=as.character(RSID), chr=as.character(CHR), pos=as.integer(POS), beta_eqtl=SPEARMANS_RHO, se_eqtl=se_eqtl, p_eqtl=P_VALUE, maf=as.numeric(A2_FREQ_ONEK1K))]
  m1 <- suppressWarnings(inner_join(eg, ilc[, .(rsid, chr, pos, beta_gwas, se_gwas, p_gwas)], by = 'rsid'))
  merged <- if (nrow(m1) < 10) suppressWarnings(inner_join(eg[, .(chr,pos,beta_eqtl,se_eqtl,p_eqtl,maf)], ilc[, .(chr,pos,beta_gwas,se_gwas,p_gwas)], by=c('chr','pos'))) else m1
  merged <- merged %>% filter(is.finite(beta_eqtl), is.finite(se_eqtl), is.finite(beta_gwas), is.finite(se_gwas))
  merged <- merged %>% filter(se_eqtl > 0, se_gwas > 0)
  if (!('snp' %in% names(merged))) {
    if ('rsid' %in% names(merged) && any(!is.na(merged$rsid))) {
      merged$snp <- as.character(merged$rsid)
    } else if (all(c('chr','pos') %in% names(merged))) {
      merged$snp <- paste0(as.character(merged$chr), ':', as.integer(merged$pos))
    }
  }
  if ('p_eqtl' %in% names(merged)) {
    merged <- merged %>% mutate(p_min = pmin(p_eqtl, p_gwas, na.rm = TRUE)) %>% arrange(snp, p_min) %>% distinct(snp, .keep_all = TRUE)
    merged <- merged %>% arrange(p_eqtl)
    if (nrow(merged) > 3000) merged <- merged[1:3000,]
  }
  if ('maf' %in% names(merged)) merged <- merged %>% filter(is.finite(maf), maf > 0, maf < 1)
  if (nrow(merged) > 10) {
    d1 <- list(beta=merged$beta_eqtl, varbeta=(merged$se_eqtl)^2, type='quant', N=982, sdY=1, snp=merged$snp)
    d2 <- list(beta=merged$beta_gwas, varbeta=(merged$se_gwas)^2, type='cc', s=0.3, N=50000, snp=merged$snp)
    val <- tryCatch({ suppressWarnings(coloc.abf(dataset1=d1, dataset2=d2)$summary['PP.H4.abf']) }, error=function(e) { paste('Failed', e$message) })
    mfg1 <- if (!is.null(fg)) suppressWarnings(inner_join(eg, fg[, .(rsid, chr, pos, beta_gwas, se_gwas, p_gwas)], by = 'rsid')) else data.table()
    merged_fg <- if (!is.null(fg) && nrow(mfg1) < 10) suppressWarnings(inner_join(eg[, .(chr,pos,beta_eqtl,se_eqtl,p_eqtl,maf)], fg[, .(chr,pos,beta_gwas,se_gwas,p_gwas)], by=c('chr','pos'))) else mfg1
    merged_fg <- merged_fg %>% filter(is.finite(beta_eqtl), is.finite(se_eqtl), is.finite(beta_gwas), is.finite(se_gwas))
    merged_fg <- merged_fg %>% filter(se_eqtl > 0, se_gwas > 0)
    if (!('snp' %in% names(merged_fg))) {
      if ('rsid' %in% names(merged_fg) && any(!is.na(merged_fg$rsid))) {
        merged_fg$snp <- as.character(merged_fg$rsid)
      } else if (all(c('chr','pos') %in% names(merged_fg))) {
        merged_fg$snp <- paste0(as.character(merged_fg$chr), ':', as.integer(merged_fg$pos))
      }
    }
    if ('p_eqtl' %in% names(merged_fg) || 'p_gwas' %in% names(merged_fg)) {
      merged_fg <- merged_fg %>% mutate(p_min = pmin(p_eqtl, p_gwas, na.rm = TRUE)) %>% arrange(snp, p_min) %>% distinct(snp, .keep_all = TRUE)
    }
    if ('p_eqtl' %in% names(merged_fg)) {
      merged_fg <- merged_fg %>% arrange(p_eqtl)
      if (nrow(merged_fg) > 3000) merged_fg <- merged_fg[1:3000,]
    }
    val_fg <- if (!is.null(fg) && nrow(merged_fg) > 10) tryCatch({ suppressWarnings(coloc.abf(dataset1=d1, dataset2=list(beta=merged_fg$beta_gwas, varbeta=(merged_fg$se_gwas)^2, type='cc', s=0.3, N=50000, snp=merged_fg$snp))$summary['PP.H4.abf']) }, error=function(e) { paste('Failed', e$message) }) else NA_character_
    reslist[[length(reslist)+1]] <- data.table(gene=g, pp4_ilcco=as.character(val), status_ilcco=if (grepl('Failed', val)) 'failed' else 'ok', n_snps_ilcco=nrow(merged), pp4_finngen=as.character(val_fg), status_finngen=if (!is.na(val_fg) && grepl('Failed', val_fg)) 'failed' else if (!is.na(val_fg)) 'ok' else 'subset_too_small', n_snps_finngen=if (exists('merged_fg')) nrow(merged_fg) else 0L)
  } else {
    mfg1 <- if (!is.null(fg)) suppressWarnings(inner_join(eg, fg[, .(rsid, chr, pos, beta_gwas, se_gwas, p_gwas)], by = 'rsid')) else data.table()
    merged_fg <- if (!is.null(fg) && nrow(mfg1) < 10) suppressWarnings(inner_join(eg[, .(chr,pos,beta_eqtl,se_eqtl,p_eqtl,maf)], fg[, .(chr,pos,beta_gwas,se_gwas,p_gwas)], by=c('chr','pos'))) else mfg1
    merged_fg <- merged_fg %>% filter(is.finite(beta_eqtl), is.finite(se_eqtl), is.finite(beta_gwas), is.finite(se_gwas))
    merged_fg <- merged_fg %>% filter(se_eqtl > 0, se_gwas > 0)
    if (!('snp' %in% names(merged_fg))) {
      if ('rsid' %in% names(merged_fg) && any(!is.na(merged_fg$rsid))) {
        merged_fg$snp <- as.character(merged_fg$rsid)
      } else if (all(c('chr','pos') %in% names(merged_fg))) {
        merged_fg$snp <- paste0(as.character(merged_fg$chr), ':', as.integer(merged_fg$pos))
      }
    }
    if ('p_eqtl' %in% names(merged_fg)) {
      merged_fg <- merged_fg %>% arrange(p_eqtl)
      if (nrow(merged_fg) > 3000) merged_fg <- merged_fg[1:3000,]
    }
    reslist[[length(reslist)+1]] <- data.table(gene=g, pp4_ilcco=NA_character_, status_ilcco='subset_too_small', n_snps_ilcco=nrow(merged), pp4_finngen=NA_character_, status_finngen=if (nrow(merged_fg)>10) 'ok' else 'subset_too_small', n_snps_finngen=nrow(merged_fg))
  }
}
resdt <- rbindlist(reslist, fill=TRUE)
try({ if (file.exists('My_MR_Project/PP4_ByGene.csv')) unlink('My_MR_Project/PP4_ByGene.csv') }, silent=TRUE)
fwrite(resdt, 'My_MR_Project/PP4_ByGene.csv')
cat(paste0('rows=', nrow(resdt), ' cols=', ncol(resdt), '\n'))
fwrite(resdt, 'My_MR_Project/PP4_ByGene_Dual.csv')

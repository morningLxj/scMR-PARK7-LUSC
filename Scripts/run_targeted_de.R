
library(Matrix)
library(data.table)

# Define paths
counts_rds_path <- "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds"
meta_path <- "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
output_csv <- "D:/2026YJ/My_MR_Project/B_cell_PARK7_High_vs_Low_DEGs_Real.csv"

target_genes <- c("PARK7", "GSR", "TXN", "SOD1", "CAT", "GPX1", "ATM", "ATR", "CHEK1", 
                  "BRCA1", "XRCC1", "NFE2L2", "HMOX1", "NQO1", 
                  "CD19", "CD79A", "CD79B", "MS4A1", "CD40")

tryCatch({
  message("Loading Counts...")
  counts <- readRDS(counts_rds_path)
  
  message("Loading Metadata...")
  meta <- fread(meta_path)
  meta <- as.data.frame(meta)
  
  barcode_col <- colnames(meta)[1] 
  celltype_col <- grep("Cell_type", colnames(meta), value=TRUE)[1]
  if(is.na(celltype_col)) celltype_col <- "Cell_type"
  rownames(meta) <- meta[[barcode_col]]
  
  # Subset B cells
  common_cells <- intersect(colnames(counts), rownames(meta))
  meta_common <- meta[common_cells, ]
  
  b_cell_indices <- grep("B cell|B_cell|B-cell", meta_common[[celltype_col]], ignore.case=TRUE)
  b_cells <- rownames(meta_common)[b_cell_indices]
  message("B cells found: ", length(b_cells))
  
  if(length(b_cells) < 10) stop("Not enough B cells")
  
  # Get expression matrix for B cells
  # Note: counts is dgCMatrix, subsetting is efficient
  counts_b <- counts[, b_cells]
  
  # Normalize (CPM-like or LogNormalize manually)
  # Simple Log(CPM + 1)
  lib_size <- colSums(counts_b)
  norm_b <- t(t(counts_b) / lib_size * 10000)
  norm_b@x <- log1p(norm_b@x)
  
  # Stratify by PARK7
  park7_expr <- norm_b["PARK7", ]
  median_val <- median(park7_expr[park7_expr > 0])
  high_cells <- names(park7_expr)[park7_expr > median_val]
  low_cells <- names(park7_expr)[park7_expr <= median_val]
  
  message("High: ", length(high_cells), ", Low: ", length(low_cells))
  
  # Calculate stats for target genes
  results <- data.frame()
  
  available_genes <- intersect(target_genes, rownames(norm_b))
  
  for(gene in available_genes) {
    val_high <- norm_b[gene, high_cells]
    val_low <- norm_b[gene, low_cells]
    
    # Wilcoxon test
    w_test <- wilcox.test(val_high, val_low)
    p_val <- w_test$p.value
    
    # Fold Change
    mean_high <- mean(expm1(val_high)) # convert back to linear for mean
    mean_low <- mean(expm1(val_low))
    log2fc <- log2((mean_high + 0.1) / (mean_low + 0.1))
    
    # Percentages
    pct_1 <- mean(val_high > 0)
    pct_2 <- mean(val_low > 0)
    
    results <- rbind(results, data.frame(
      gene = gene,
      p_val = p_val,
      avg_log2FC = log2fc,
      pct.1 = pct_1,
      pct.2 = pct_2
    ))
  }
  
  # Adjust P-values
  results$p_val_adj <- p.adjust(results$p_val, method = "bonferroni", n = nrow(counts_b)) # Adjust for genome-wide testing
  
  # Format
  results <- results[order(results$p_val), ]
  
  write.csv(results, output_csv, row.names = FALSE)
  message("Saved targeted DEGs to: ", output_csv)
  print(results)

}, error = function(e) {
  message("Error: ", e$message)
  quit(status = 1)
})

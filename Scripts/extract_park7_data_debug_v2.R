
library(Seurat)
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
  message("Loading Counts (Seurat loaded)...")
  counts <- readRDS(counts_rds_path)
  message("Counts loaded. Dim: ", paste(dim(counts), collapse=" x "))
  
  message("Loading Metadata...")
  meta <- fread(meta_path)
  meta <- as.data.frame(meta)
  
  barcode_col <- colnames(meta)[1] 
  celltype_col <- grep("Cell_type", colnames(meta), value=TRUE)[1]
  rownames(meta) <- meta[[barcode_col]]
  
  common_cells <- intersect(colnames(counts), rownames(meta))
  meta_common <- meta[common_cells, ]
  
  b_cell_indices <- grep("B cell|B_cell|B-cell", meta_common[[celltype_col]], ignore.case=TRUE)
  b_cells <- rownames(meta_common)[b_cell_indices]
  message("B cells found: ", length(b_cells))
  
  if(length(b_cells) < 10) stop("Not enough B cells")
  
  # Calculate lib sizes manually to stay light
  message("Calculating lib sizes...")
  # We only need lib sizes for B cells.
  # Subsetting counts[, b_cells] might crash.
  # Let's try to calculate lib sizes for B cells by iterating? No, that's slow.
  # colSums on the whole matrix is fast.
  all_lib_sizes <- colSums(counts)
  b_lib_sizes <- all_lib_sizes[b_cells]
  
  # PARK7 Stratification
  message("Stratifying...")
  park7_raw <- counts["PARK7", b_cells]
  park7_norm <- log1p(park7_raw / b_lib_sizes * 10000)
  
  median_val <- median(park7_norm[park7_norm > 0])
  if(is.na(median_val)) median_val <- 0
  
  high_cells <- names(park7_norm)[park7_norm > median_val]
  low_cells <- names(park7_norm)[park7_norm <= median_val]
  message("High: ", length(high_cells), ", Low: ", length(low_cells))
  
  results <- data.frame()
  available_genes <- intersect(target_genes, rownames(counts))
  
  for(gene in available_genes) {
    gene_raw <- counts[gene, b_cells]
    gene_norm <- log1p(gene_raw / b_lib_sizes * 10000)
    
    val_high <- gene_norm[high_cells]
    val_low <- gene_norm[low_cells]
    
    w_test <- wilcox.test(val_high, val_low)
    p_val <- w_test$p.value
    
    mean_high <- mean(expm1(val_high))
    mean_low <- mean(expm1(val_low))
    log2fc <- log2((mean_high + 0.1) / (mean_low + 0.1))
    
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
  
  results$p_val_adj <- p.adjust(results$p_val, method = "bonferroni", n = 20000)
  results <- results[order(results$p_val), ]
  
  write.csv(results, output_csv, row.names = FALSE)
  message("SUCCESS: Saved DEGs to ", output_csv)
  print(results)

}, error = function(e) {
  message("Error: ", e$message)
  quit(status = 1)
})

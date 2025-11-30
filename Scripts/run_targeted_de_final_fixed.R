
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
  message("Loading Metadata...")
  meta <- fread(meta_path)
  meta <- as.data.frame(meta)
  
  barcode_col <- colnames(meta)[1] 
  celltype_col <- "Cell_type" 
  rownames(meta) <- meta[[barcode_col]]
  
  message("Loading Counts...")
  counts <- readRDS(counts_rds_path)
  
  common_cells <- intersect(colnames(counts), rownames(meta))
  message("Common cells: ", length(common_cells))
  
  meta_common <- meta[common_cells, ]
  b_cells <- rownames(meta_common)[meta_common[[celltype_col]] == "B lymphocytes"]
  message("B lymphocytes found: ", length(b_cells))
  
  if(length(b_cells) < 10) stop("Not enough B cells")
  
  message("Calculating lib sizes (B cells only)...")
  # Use counts[, b_cells] but ensure it returns a matrix-like object we can use colSums on
  # counts is likely a dgCMatrix
  counts_b <- counts[, b_cells]
  b_lib_sizes <- colSums(counts_b)
  
  message("Processing PARK7...")
  gene_name <- "PARK7"
  if(!gene_name %in% rownames(counts_b)) {
      if("DJ-1" %in% rownames(counts_b)) gene_name <- "DJ-1"
      else if("DJ1" %in% rownames(counts_b)) gene_name <- "DJ1"
      else stop("PARK7 not found")
  }
  message("Using gene: ", gene_name)
  
  # Force numeric conversion for the vector
  park7_raw <- as.numeric(counts_b[gene_name, ])
  park7_norm <- log1p(park7_raw / b_lib_sizes * 10000)
  
  median_val <- median(park7_norm[park7_norm > 0])
  if(is.na(median_val)) median_val <- 0
  
  # Define groups using indices or logical vector since names might be lost with as.numeric
  high_mask <- park7_norm > median_val
  low_mask <- park7_norm <= median_val
  
  message("High: ", sum(high_mask), ", Low: ", sum(low_mask))
  
  results <- data.frame()
  available_genes <- intersect(target_genes, rownames(counts_b))
  
  for(gene in available_genes) {
    # Force numeric
    gene_raw <- as.numeric(counts_b[gene, ])
    gene_norm <- log1p(gene_raw / b_lib_sizes * 10000)
    
    val_high <- gene_norm[high_mask]
    val_low <- gene_norm[low_mask]
    
    # Check if we have enough data for test
    if(length(val_high) > 0 && length(val_low) > 0) {
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
  }
  
  if(nrow(results) > 0) {
      results$p_val_adj <- p.adjust(results$p_val, method = "bonferroni", n = 20000)
      results <- results[order(results$p_val), ]
      write.csv(results, output_csv, row.names = FALSE)
      message("Success!")
      print(results)
  } else {
      message("No results generated.")
  }

}, error = function(e) {
  message("Error: ", e$message)
  quit(status = 1)
})

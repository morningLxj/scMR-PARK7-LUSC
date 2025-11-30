
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
  
  # INSPECT METADATA
  message("Metadata loaded. Dim: ", paste(dim(meta), collapse=" x "))
  message("Colnames: ", paste(colnames(meta), collapse=", "))
  
  # Try to find relevant columns
  barcode_col <- colnames(meta)[1] 
  
  # Look for any column that might contain cell type info
  possible_type_cols <- grep("type|cell|cluster|annot", colnames(meta), value=TRUE, ignore.case=TRUE)
  message("Possible cell type columns: ", paste(possible_type_cols, collapse=", "))
  
  celltype_col <- "Cell_type" # Default expectation
  if(!celltype_col %in% colnames(meta)) {
      if(length(possible_type_cols) > 0) {
          celltype_col <- possible_type_cols[1]
          message("Using column: ", celltype_col)
      } else {
          stop("Could not identify cell type column")
      }
  }
  
  # Check unique values in cell type column
  unique_types <- unique(meta[[celltype_col]])
  message("Unique cell types found (head): ", paste(head(unique_types, 10), collapse=", "))
  
  # Check if any contain "B"
  b_matches <- grep("B", unique_types, value=TRUE, ignore.case=TRUE)
  message("Cell types matching 'B': ", paste(b_matches, collapse=", "))
  
  rownames(meta) <- meta[[barcode_col]]
  
  message("Loading Counts...")
  counts <- readRDS(counts_rds_path)
  
  common_cells <- intersect(colnames(counts), rownames(meta))
  message("Common cells: ", length(common_cells))
  
  meta_common <- meta[common_cells, ]
  
  # TRY EXACT MATCHING IF GREP FAILED
  # Or broaden search
  b_cell_indices <- grep("B_cell|B-cell|B cell|Plasma|Lymphocyte", meta_common[[celltype_col]], ignore.case=TRUE)
  b_cells <- rownames(meta_common)[b_cell_indices]
  message("B cells found (broad search): ", length(b_cells))
  
  if(length(b_cells) < 10) {
      # Print table of cell types to help debug
      print(table(meta_common[[celltype_col]]))
      stop("Not enough B cells")
  }
  
  # ... rest of analysis ...
  # (Only proceed if B cells found)
  
  message("Calculating lib sizes...")
  all_lib_sizes <- colSums(counts)
  b_lib_sizes <- all_lib_sizes[b_cells]
  
  message("Processing PARK7...")
  # Check if PARK7 exists
  if(!"PARK7" %in% rownames(counts)) {
      message("PARK7 not found. Trying aliases...")
      if("DJ-1" %in% rownames(counts)) {
          gene_name <- "DJ-1"
      } else if("DJ1" %in% rownames(counts)) {
          gene_name <- "DJ1"
      } else {
          stop("PARK7 gene not found")
      }
  } else {
      gene_name <- "PARK7"
  }
  message("Using gene: ", gene_name)
  
  park7_raw <- counts[gene_name, b_cells]
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
  message("Success!")
  print(results)

}, error = function(e) {
  message("Error: ", e$message)
  quit(status = 1)
})

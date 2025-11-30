
library(Seurat)
library(Matrix)
library(data.table)
library(dplyr)

# Define paths
counts_rds_path <- "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds"
meta_path <- "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
output_csv <- "D:/2026YJ/My_MR_Project/B_cell_PARK7_High_vs_Low_DEGs_Real.csv"

message("Starting DE Analysis v2...")

tryCatch({
  # 1. Load Metadata First (Lightweight)
  message("Loading Metadata...")
  meta <- fread(meta_path)
  meta <- as.data.frame(meta)
  
  barcode_col <- colnames(meta)[1] 
  celltype_col <- grep("Cell_type", colnames(meta), value=TRUE)[1]
  if(is.na(celltype_col)) celltype_col <- "Cell_type"
  rownames(meta) <- meta[[barcode_col]]
  
  # Filter for B cells in metadata FIRST to save memory if we could partial load, 
  # but readRDS loads all. We'll subset after loading.
  
  # 2. Load Counts
  message("Loading Counts RDS...")
  counts <- readRDS(counts_rds_path)
  message("Counts loaded.")

  # 3. Create Seurat Object (Minimal)
  # Only keep cells present in both
  common_cells <- intersect(colnames(counts), rownames(meta))
  message("Common cells: ", length(common_cells))
  
  counts <- counts[, common_cells]
  meta <- meta[common_cells, ]
  
  seu <- CreateSeuratObject(counts = counts, meta.data = meta)
  
  # 4. Subset B cells
  message("Subsetting B cells...")
  cell_types <- unique(seu@meta.data[[celltype_col]])
  b_cell_label <- grep("B cell|B_cell|B-cell", cell_types, value=TRUE, ignore.case=TRUE)
  message("B cell labels: ", paste(b_cell_label, collapse=", "))
  
  if(length(b_cell_label) == 0) stop("No B cells found")
  
  seu_b <- subset(seu, cells = rownames(seu@meta.data)[seu@meta.data[[celltype_col]] %in% b_cell_label])
  message("B cells: ", ncol(seu_b))
  
  # Free up memory
  rm(counts, seu)
  gc()
  
  # 5. Normalize (LogNormalize)
  message("Normalizing...")
  seu_b <- NormalizeData(seu_b, verbose = FALSE)
  
  # 6. Stratify
  gene_name <- "PARK7"
  if(!gene_name %in% rownames(seu_b)) gene_name <- "DJ-1"
  
  expr <- GetAssayData(seu_b, slot = "data")[gene_name, ]
  median_val <- median(expr[expr > 0]) # Median of expressing cells
  if(is.na(median_val)) median_val <- 0
  
  seu_b$PARK7_Group <- ifelse(expr > median_val, "High", "Low")
  Idents(seu_b) <- "PARK7_Group"
  
  message("Groups: ", paste(table(seu_b$PARK7_Group), collapse=" vs "))
  
  # 7. Find Markers (Wilcoxon)
  message("Running FindMarkers...")
  # Use a looser threshold to ensure we get results, then filter
  markers <- FindMarkers(seu_b, ident.1 = "High", ident.2 = "Low", 
                         logfc.threshold = 0.1, 
                         min.pct = 0.1,
                         verbose = TRUE)
  
  markers$gene <- rownames(markers)
  markers <- markers %>% 
             arrange(p_val_adj) %>%
             select(gene, avg_log2FC, p_val_adj, pct.1, pct.2)
  
  # 8. Save
  write.csv(markers, output_csv, row.names = FALSE)
  message("Saved to: ", output_csv)
  
  # Print top results
  print(head(markers, 15))

}, error = function(e) {
  message("CRITICAL ERROR: ", e$message)
  quit(status = 1)
})

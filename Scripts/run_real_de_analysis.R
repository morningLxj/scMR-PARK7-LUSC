
library(Seurat)
library(Matrix)
library(data.table)
library(dplyr)

# Define paths
counts_rds_path <- "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds"
meta_path <- "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
output_csv <- "D:/2026YJ/My_MR_Project/B_cell_PARK7_High_vs_Low_DEGs_Real.csv"

tryCatch({
  message("Step 1: Loading Counts...")
  counts <- readRDS(counts_rds_path)
  message("Counts dim: ", paste(dim(counts), collapse="x"))

  message("Step 2: Loading Metadata...")
  meta <- fread(meta_path)
  meta <- as.data.frame(meta)
  
  # Align Barcodes
  barcode_col <- colnames(meta)[1] 
  celltype_col <- grep("Cell_type", colnames(meta), value=TRUE)[1]
  if(is.na(celltype_col)) celltype_col <- "Cell_type"
  rownames(meta) <- meta[[barcode_col]]
  
  # Create Seurat Object
  message("Step 3: Creating Seurat Object...")
  seu <- CreateSeuratObject(counts = counts, meta.data = meta)
  
  # Normalize
  message("Step 4: Normalizing Data...")
  seu <- NormalizeData(seu)
  
  # Subset B cells
  message("Step 5: Subsetting B cells...")
  # Check available cell types
  cell_types <- unique(seu@meta.data[[celltype_col]])
  message("Cell types found: ", paste(head(cell_types), collapse=", "))
  
  # Try to find B cells (flexible matching)
  b_cell_label <- grep("B cell|B_cell|B-cell", cell_types, value=TRUE, ignore.case=TRUE)
  if(length(b_cell_label) == 0) stop("No B cells found in metadata labels")
  
  seu_b <- subset(seu, cells = rownames(seu@meta.data)[seu@meta.data[[celltype_col]] %in% b_cell_label])
  message("B cells count: ", ncol(seu_b))
  
  if(ncol(seu_b) < 10) stop("Too few B cells for analysis")

  # Step 6: Stratify by PARK7 Expression
  gene_name <- "PARK7"
  if(!gene_name %in% rownames(seu_b)) {
      if("DJ-1" %in% rownames(seu_b)) gene_name <- "DJ-1"
      else stop("PARK7 not found")
  }
  
  expr <- GetAssayData(seu_b, slot = "data")[gene_name, ]
  median_val <- median(expr[expr > 0]) # Median of expressed cells or global median
  
  seu_b$PARK7_Group <- ifelse(expr > median_val, "High", "Low")
  Idents(seu_b) <- "PARK7_Group"
  
  message("High group: ", sum(seu_b$PARK7_Group == "High"))
  message("Low group: ", sum(seu_b$PARK7_Group == "Low"))

  # Step 7: Find Markers
  message("Step 7: Running FindMarkers...")
  markers <- FindMarkers(seu_b, ident.1 = "High", ident.2 = "Low", min.pct = 0.1, logfc.threshold = 0.25)
  
  # Add gene column
  markers$gene <- rownames(markers)
  markers <- markers %>% select(gene, everything())
  
  # Save
  write.csv(markers, output_csv, row.names = FALSE)
  message("Successfully saved REAL DEGs to: ", output_csv)
  
  # Print top 10 for log
  print(head(markers, 10))

}, error = function(e) {
  message("Error: ", e$message)
  quit(status = 1)
})

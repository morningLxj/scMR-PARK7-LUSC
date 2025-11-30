
library(Seurat)
library(Matrix)
library(data.table)

# Define paths
counts_rds_path <- "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds"
meta_path <- "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
output_csv <- "D:/2026YJ/My_MR_Project/PARK7_Expression_Data.csv"

tryCatch({
  # 1. Load Counts
  message("Reading RDS: ", counts_rds_path)
  counts <- readRDS(counts_rds_path)
  message("Counts loaded. Class: ", class(counts), ", Dim: ", paste(dim(counts), collapse=" x "))
  
  # 2. Load Metadata
  message("Reading Metadata: ", meta_path)
  meta <- fread(meta_path)
  meta <- as.data.frame(meta)
  message("Metadata loaded. Dim: ", paste(dim(meta), collapse=" x "))
  
  # 3. Align Barcodes
  barcode_col <- colnames(meta)[1] 
  celltype_col <- grep("Cell_type", colnames(meta), value=TRUE)[1]
  if(is.na(celltype_col)) celltype_col <- "Cell_type"
  
  rownames(meta) <- meta[[barcode_col]]
  
  common_cells <- intersect(colnames(counts), rownames(meta))
  message("Found ", length(common_cells), " common cells")
  
  if (length(common_cells) == 0) {
    message("First 5 counts cols: ", paste(head(colnames(counts)), collapse=", "))
    message("First 5 meta rows: ", paste(head(rownames(meta)), collapse=", "))
    stop("No overlapping cells between counts and metadata")
  }
  
  # 4. Extract PARK7 Expression
  gene_name <- "PARK7"
  if (!gene_name %in% rownames(counts)) {
      if ("DJ-1" %in% rownames(counts)) gene_name <- "DJ-1"
      else if ("DJ1" %in% rownames(counts)) gene_name <- "DJ1"
      else {
          message("Available genes (head): ", paste(head(rownames(counts)), collapse=", "))
          stop("PARK7 gene not found in matrix")
      }
  }
  message("Using gene symbol: ", gene_name)
  
  expr_data <- counts[gene_name, common_cells]
  
  # 5. Create result data frame
  df_out <- data.frame(
    Barcode = common_cells,
    Cell_Type = meta[common_cells, celltype_col],
    Expression = as.numeric(expr_data)
  )
  
  # 6. Save to CSV
  write.csv(df_out, output_csv, row.names = FALSE)
  message("Successfully saved expression data to: ", output_csv)
  
}, error = function(e) {
  message("Error occurred: ", e$message)
  quit(status = 1)
})

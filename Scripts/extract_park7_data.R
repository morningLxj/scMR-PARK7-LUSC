
library(Seurat)
library(Matrix)
library(data.table)

# Define paths
counts_rds_path <- "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds"
counts_gz_path <- "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz"
meta_path <- "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
output_csv <- "D:/2026YJ/My_MR_Project/PARK7_Expression_Data.csv"

# Function to read RDS (handling potential gzip)
read_counts <- function() {
  if (file.exists(counts_rds_path)) {
    message("Reading RDS: ", counts_rds_path)
    con <- file(counts_rds_path, "rb")
    if (summary(con)$class == "gzfile") {
       close(con)
       con <- gzcon(file(counts_rds_path, "rb"))
    }
    x <- readRDS(con)
    close(con)
    return(x)
  } else if (file.exists(counts_gz_path)) {
    message("Reading RDS.gz: ", counts_gz_path)
    con <- gzcon(file(counts_gz_path, "rb"))
    x <- readRDS(con)
    close(con)
    return(x)
  }
  stop("Counts file not found")
}

tryCatch({
  # 1. Load Counts
  counts <- read_counts()
  
  # 2. Load Metadata
  message("Reading Metadata: ", meta_path)
  meta <- fread(meta_path)
  meta <- as.data.frame(meta)
  
  # 3. Align Barcodes
  # Assuming first column is barcode/Index
  barcode_col <- colnames(meta)[1] 
  # Find cell type column
  celltype_col <- grep("Cell_type", colnames(meta), value=TRUE)[1]
  if(is.na(celltype_col)) celltype_col <- "Cell_type"
  
  rownames(meta) <- meta[[barcode_col]]
  
  common_cells <- intersect(colnames(counts), rownames(meta))
  message("Found ", length(common_cells), " common cells")
  
  if (length(common_cells) == 0) stop("No overlapping cells between counts and metadata")
  
  # 4. Extract PARK7 Expression
  # Check gene name variations
  gene_name <- "PARK7"
  if (!gene_name %in% rownames(counts)) {
      if ("DJ-1" %in% rownames(counts)) gene_name <- "DJ-1"
      else if ("DJ1" %in% rownames(counts)) gene_name <- "DJ1"
      else stop("PARK7 gene not found in matrix")
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

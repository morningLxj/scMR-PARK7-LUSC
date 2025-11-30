
message("DEBUG: Script started")
flush.console()

# Check packages
if (!require("Matrix", quietly = TRUE)) stop("Matrix package missing")
if (!require("data.table", quietly = TRUE)) stop("data.table package missing")
message("DEBUG: Packages loaded")
flush.console()

counts_rds_path <- "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds"
output_csv <- "D:/2026YJ/My_MR_Project/B_cell_PARK7_High_vs_Low_DEGs_Real.csv"

if (!file.exists(counts_rds_path)) {
  stop("File not found: ", counts_rds_path)
}
message("DEBUG: File exists, size: ", file.size(counts_rds_path))
flush.console()

tryCatch({
  message("DEBUG: Attempting readRDS...")
  flush.console()
  
  counts <- readRDS(counts_rds_path)
  
  message("DEBUG: readRDS success!")
  message("DEBUG: Class: ", class(counts))
  message("DEBUG: Dimensions: ", paste(dim(counts), collapse=" x "))
  flush.console()
  
  # Create a dummy result to prove we can write
  write.csv(data.frame(Status="Success", Rows=nrow(counts)), output_csv)
  message("DEBUG: Wrote test output")
  
}, error = function(e) {
  message("DEBUG: Error caught: ", e$message)
  quit(status = 1)
})

message("DEBUG: Script finished")
quit(status = 0)

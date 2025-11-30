con_path <- "My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds"
out_path <- "My_MR_Project/scRNA_Data/diagnose_counts.txt"
sink(out_path)
cat("diagnose start\n")
obj <- readRDS(con_path)
cat("class=", paste(class(obj), collapse=","), "\n")
cat("typeof=", typeof(obj), "\n")
if (!is.null(dim(obj))) {
  cat("dim=", paste(dim(obj), collapse=","), "\n")
}
if (is.list(obj)) {
  cat("names=", paste(names(obj), collapse=","), "\n")
}
if (inherits(obj, "dgCMatrix")) {
  cat("dgCMatrix dims=", paste(dim(obj), collapse=","), "\n")
}
str(obj, max.level=1)
cat("diagnose end\n")
sink()

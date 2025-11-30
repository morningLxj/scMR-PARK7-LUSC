library(data.table)
out <- "My_MR_Project/scRNA_Data/inspect.txt"
sink(out)
cat("inspect start\n")
con <- gzcon(file("My_MR_Project/scRNA_Data/GSE131907_counts.rds","rb"))
counts <- readRDS(con)
close(con)
cat("counts_class=", paste(class(counts), collapse=","), "\n")
cat("counts_dim=", paste(dim(counts), collapse=","), "\n")
md <- fread("My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_cell_annotation.txt.gz")
cat("md_cols=", paste(names(md), collapse=","), "\n")
if (!is.null(colnames(counts))) {
  cat("counts_cols_head=", paste(head(colnames(counts),5), collapse=","), "\n")
}
best_col <- NA_character_
best_overlap <- -1L
for (nm in names(md)) {
  ov <- sum(as.character(md[[nm]]) %in% colnames(counts))
  if (ov > best_overlap) { best_overlap <- ov; best_col <- nm }
  if (ov > 0) cat("overlap:", nm, ov, "\n")
}
cat("best_col=", best_col, " best_overlap=", best_overlap, "\n")
sink()

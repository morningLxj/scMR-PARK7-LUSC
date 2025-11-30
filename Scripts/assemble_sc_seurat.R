library(Seurat)
library(data.table)

cat("Loading raw counts\n")
counts_path_gz <- "My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz"
counts_path_rds <- "My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds"
if (file.exists(counts_path_rds)) {
  counts <- tryCatch(readRDS(counts_path_rds), error = function(e) {
    con <- gzcon(file(counts_path_rds, "rb"))
    on.exit(close(con))
    readRDS(con)
  })
} else if (file.exists(counts_path_gz)) {
  counts <- tryCatch(readRDS(counts_path_gz), error = function(e) {
    con <- gzcon(file(counts_path_gz, "rb"))
    on.exit(close(con))
    readRDS(con)
  })
} else {
  stop("Counts file not found")
}

cat("Loading metadata\n")
meta_path <- "My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
metadata <- fread(meta_path)
metadata <- as.data.frame(metadata)
cols_counts <- colnames(counts)
best_col <- NA_character_
best_overlap <- -1L
for (nm in names(metadata)) {
  ov <- sum(as.character(metadata[[nm]]) %in% cols_counts)
  if (ov > best_overlap) { best_overlap <- ov; best_col <- nm }
}
if (is.na(best_col)) best_col <- names(metadata)[1]
rownames(metadata) <- as.character(metadata[[best_col]])

common_cells <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common_cells]
if (is.data.frame(counts)) counts <- as.matrix(counts)
metadata <- metadata[common_cells, , drop = FALSE]

ctype_candidates <- c("Cell_type","cell_type","celltype","CellType","MajorCellType","cell_type_level","annotation","Annotation","cluster_label","cluster","Type","Identity","Cell","cellclass","cell_class","celltypes")
ctype_col <- ctype_candidates[ctype_candidates %in% names(metadata)]
ctype_use <- NA_character_
if (length(ctype_col) > 0) {
  ctype_use <- ctype_col[1]
  metadata$Cell_type <- as.character(metadata[[ctype_use]])
}
sorigin_candidates <- c("Sample_Origin","sample_origin","Origin","origin","Condition","condition","Tumor_Normal","tumor_normal","TumorOrNormal","Group","group","TissueType","tissue","Source","source")
sorigin_col <- sorigin_candidates[sorigin_candidates %in% names(metadata)]
if (length(sorigin_col) > 0) metadata$Sample_Origin <- as.character(metadata[[sorigin_col[1]]])

sc_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)

if ("Cell_type" %in% names(sc_obj@meta.data)) print(table(sc_obj$Cell_type))
if ("Sample_Origin" %in% names(sc_obj@meta.data)) print(table(sc_obj$Sample_Origin))

saveRDS(sc_obj, "My_MR_Project/scRNA_Data/GSE131907_expr.rds")
cat("Saved My_MR_Project/scRNA_Data/GSE131907_expr.rds\n")

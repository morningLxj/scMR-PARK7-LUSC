library(data.table)
library(Matrix)

meta_path <- "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
meta <- data.table::fread(meta_path)
meta <- as.data.frame(meta)
id_candidates <- c("Index","CellID","Barcode","cell_id","cellid","Cell")
id_col <- id_candidates[id_candidates %in% names(meta)]
if (length(id_col) > 0) rownames(meta) <- as.character(meta[[id_col[1]]]) else rownames(meta) <- as.character(meta[[1]])
ctype_candidates <- c("Cell_type","cell_type","celltype","CellType","MajorCellType","annotation","Annotation","cluster_label","cluster","Type","Identity")
cell_type_col <- ctype_candidates[ctype_candidates %in% names(meta)][1]
type_vec <- as.character(meta[[cell_type_col]])
keep_bplasma <- grepl("B|Plasma", type_vec, ignore.case = TRUE)
keep_cd8 <- grepl("CD8", type_vec, ignore.case = TRUE)
cells_bplasma <- rownames(meta)[keep_bplasma]
cells_cd8 <- rownames(meta)[keep_cd8]
if (length(cells_cd8) > 30000) { set.seed(1); cells_cd8 <- sample(cells_cd8, 30000) }
target_cells <- unique(c(cells_bplasma, cells_cd8))
cat("target_cells_n=", length(target_cells), "\n")

counts_path <- "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds"
con <- gzcon(file(counts_path, "rb"))
counts <- readRDS(con)
close(con)
cat("counts_dim=", paste(dim(counts), collapse=","), "\n")
common <- intersect(colnames(counts), target_cells)
cat("common_n=", length(common), "\n")
sub <- counts[, common, drop = FALSE]
if (is.data.frame(sub)) sub <- as.matrix(sub)
sub <- Matrix(sub, sparse = TRUE)
saveRDS(sub, "D:/2026YJ/My_MR_Project/scRNA_Data/subset_bcd8_counts.rds")
cat("saved subset_bcd8_counts.rds\n")

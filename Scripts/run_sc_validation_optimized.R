library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(Matrix)
if ("memory.limit" %in% ls(envir = baseenv())) try(suppressWarnings(memory.limit(16000)), silent = TRUE)

cat("Step 1: Loading Metadata to filter cells...\n")
meta_path <- "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
meta <- fread(meta_path)
meta <- as.data.frame(meta)

id_candidates <- c("Index","CellID","Barcode","cell_id","cellid","Cell")
id_col <- id_candidates[id_candidates %in% names(meta)]
if (length(id_col) > 0) {
  rownames(meta) <- as.character(meta[[id_col[1]]])
} else {
  rownames(meta) <- as.character(meta[[1]])
}

ctype_candidates <- c("Cell_type","cell_type","celltype","CellType","MajorCellType","annotation","Annotation","cluster_label","cluster","Type","Identity")
cell_type_col <- ctype_candidates[ctype_candidates %in% names(meta)]
if (length(cell_type_col) == 0) stop("Cannot find Cell_type column!")
cell_type_col <- cell_type_col[1]

type_vec <- as.character(meta[[cell_type_col]])
keep_bplasma <- grepl("B|Plasma", type_vec, ignore.case = TRUE)
keep_cd8 <- grepl("CD8", type_vec, ignore.case = TRUE)
cells_bplasma <- rownames(meta)[keep_bplasma]
cells_cd8 <- rownames(meta)[keep_cd8]
if (length(cells_cd8) > 30000) {
  set.seed(1)
  cells_cd8 <- sample(cells_cd8, 30000)
}
target_cells <- unique(c(cells_bplasma, cells_cd8))

cat(paste("Selected", length(target_cells), "cells for analysis (out of", nrow(meta), "total).\n"))

cat("Step 2: Loading Expression Matrix...\n")
counts_paths <- c(
  "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_counts.rds",
  "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds",
  "D:/2026YJ/My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz"
)

counts <- NULL
for (p in counts_paths) {
  if (file.exists(p)) {
    counts <- tryCatch(readRDS(p), error = function(e) {
      con <- gzcon(file(p, "rb"))
      on.exit(close(con))
      readRDS(con)
    })
    break
  }
}
if (is.null(counts)) stop("Counts file not found")
cat("counts_class=", paste(class(counts), collapse=","), "\n")
if (!is.null(dim(counts))) cat("counts_dim=", paste(dim(counts), collapse=","), "\n")

cols <- colnames(counts)
common_cells <- intersect(cols, target_cells)
cat(paste("Matched", length(common_cells), "cells in expression matrix.\n"))

sub_counts <- counts[, common_cells, drop = FALSE]
cat("sub_counts_dim=", paste(dim(sub_counts), collapse=","), "\n")
sub_meta <- meta[common_cells, , drop = FALSE]
rm(counts); gc()

if (is.data.frame(sub_counts)) {
  sub_counts <- as.matrix(sub_counts)
}
if (!inherits(sub_counts, "dgCMatrix")) {
  sub_counts <- Matrix(sub_counts, sparse = TRUE)
}

cat("Step 3: Creating Seurat Object...\n")
sc <- CreateSeuratObject(counts = sub_counts, meta.data = sub_meta)
sc <- NormalizeData(sc)
Idents(sc) <- sc[[cell_type_col]]

cat("Step 4: Generating Validation Plots...\n")
pdf_dir <- "D:/2026YJ/My_MR_Project/scRNA_Validation_Plots"
dir.create(pdf_dir, showWarnings = FALSE)

genes_try <- c("PARK7","DJ1","DJ-1")
gene_use <- genes_try[genes_try %in% rownames(sc)]
if (length(gene_use) == 0) gene_use <- "PARK7"

p1 <- VlnPlot(sc, features = gene_use[1], pt.size = 0, sort = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ggtitle("PARK7 Expression across Immune Cells")
ggsave(file.path(pdf_dir, "Fig3A_PARK7_CellType_Violin.pdf"), p1, width = 8, height = 6)

origin_candidates <- c("Sample_Origin","sample_origin","Origin","origin","Tissue","tissue","Condition","Tumor_Normal","Group")
origin_col <- origin_candidates[origin_candidates %in% names(meta)]
if (length(origin_col) > 0) {
  origin_col <- origin_col[1]
  b_levels <- levels(Idents(sc))
  b_sel <- grep("B|Plasma", b_levels, value = TRUE, ignore.case = TRUE)
  if (length(b_sel) > 0) {
    b_cells <- subset(sc, idents = b_sel)
    p2 <- VlnPlot(b_cells, features = gene_use[1], group.by = origin_col, pt.size = 0.1) +
      ggtitle("PARK7 in B Cells: Tumor vs Normal")
    ggsave(file.path(pdf_dir, "Fig3B_PARK7_TumorVsNormal.pdf"), p2, width = 6, height = 6)
  }
}

cat("Validation Complete! Plots saved to:", pdf_dir, "\n")

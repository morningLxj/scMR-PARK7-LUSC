library(Seurat)
library(ggplot2)
library(dplyr)

obj_path <- "My_MR_Project/scRNA_Data/GSE131907_expr.rds"
if (file.exists(obj_path)) {
  sc_obj <- readRDS(obj_path)
} else {
  c_gz <- "My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz"
  c_rds <- "My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds"
  if (file.exists(c_rds)) {
    con <- gzcon(file(c_rds, "rb"))
    counts <- readRDS(con)
    close(con)
  } else if (file.exists(c_gz)) {
    con <- gzcon(file(c_gz, "rb"))
    counts <- readRDS(con)
    close(con)
  } else {
    stop("counts file missing")
  }
  md_path <- "My_MR_Project/scRNA_Data/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
  metadata <- data.table::fread(md_path)
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
  sc_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
}

DefaultAssay(sc_obj) <- "RNA"

if (!("Cell_type" %in% names(sc_obj@meta.data))) {
  cands <- c("cell_type","celltype","CellType","MajorCellType","annotation","Annotation","cluster_label","cluster","Type","Identity","Cell")
  hit <- cands[cands %in% names(sc_obj@meta.data)]
  if (length(hit) > 0) sc_obj$Cell_type <- as.character(sc_obj@meta.data[[hit[1]]])
}
if (!("Sample_Origin" %in% names(sc_obj@meta.data))) {
  cands2 <- c("sample_origin","Origin","origin","Condition","Tumor_Normal","TumorOrNormal","Group","TissueType","Source")
  hit2 <- cands2[cands2 %in% names(sc_obj@meta.data)]
  if (length(hit2) > 0) sc_obj$Sample_Origin <- as.character(sc_obj@meta.data[[hit2[1]]])
}

if ("Cell_type" %in% names(sc_obj@meta.data)) print(table(sc_obj$Cell_type))

genes_try <- c("PARK7","DJ1","DJ-1")
gene_use <- genes_try[genes_try %in% rownames(sc_obj)]
if (length(gene_use) > 0) {
  p1 <- VlnPlot(sc_obj, features = gene_use[1], group.by = "Cell_type", pt.size = 0) +
    geom_boxplot(width = 0.1, fill = "white") +
    ggtitle("PARK7 Expression in Lung Adenocarcinoma TME") +
    theme(legend.position = "none")
  ggsave("TME_PARK7_Expression_Violin.pdf", p1, width = 8, height = 6)
}

b_cells <- NULL
if ("Cell_type" %in% names(sc_obj@meta.data)) {
  patt <- "(^B$|B cells|B_cell|B-cell|B_lineage|^Plasma$|Plasma cells|Plasmablast)"
  b_cells <- subset(sc_obj, subset = grepl(patt, Cell_type, ignore.case = TRUE))
}

if (!is.null(b_cells) && ("Sample_Origin" %in% names(sc_obj@meta.data))) {
  if (length(gene_use) > 0) {
    p2 <- VlnPlot(b_cells, features = gene_use[1], group.by = "Sample_Origin") +
      ggtitle("PARK7 in B Cells: Tumor vs Normal")
    ggsave("PARK7_Bcell_Tumor_vs_Normal.pdf", p2, width = 5, height = 5)
  }
}

if (length(gene_use) > 0) {
  expr_val <- GetAssayData(sc_obj, slot = "data")[gene_use[1], ]
  med_pos <- median(expr_val[expr_val > 0])
  if (!is.finite(med_pos)) med_pos <- median(expr_val)
  sc_obj$PARK7_Status <- ifelse(expr_val > med_pos, "PARK7_High", "PARK7_Low")
}

if (!is.null(b_cells) && ("PARK7_Status" %in% names(sc_obj@meta.data))) {
  b_cells <- SetIdent(b_cells, value = "PARK7_Status")
  de_markers <- FindMarkers(b_cells, ident.1 = "PARK7_High", ident.2 = "PARK7_Low")
  write.csv(de_markers, "B_cell_PARK7_High_vs_Low_DEGs.csv")
}

if (!("umap" %in% names(Reductions(sc_obj)))) {
  sc_obj <- NormalizeData(sc_obj)
  sc_obj <- FindVariableFeatures(sc_obj)
  sc_obj <- ScaleData(sc_obj)
  sc_obj <- RunPCA(sc_obj)
  sc_obj <- RunUMAP(sc_obj, dims = 1:20)
}

if (length(gene_use) == 0) genes_try <- c("PARK7")
feature_genes <- c(genes_try[1], "CD19", "MS4A1", "CD274", "HLA-DRA", "NFE2L2")
feature_genes <- feature_genes[feature_genes %in% rownames(sc_obj)]
if (length(feature_genes) > 0) {
  p3 <- FeaturePlot(sc_obj, features = feature_genes, ncol = 3, reduction = "umap")
  ggsave("PARK7_CoExpression_UMAP.pdf", p3, width = 12, height = 8)
}

cat("Validation done\n")

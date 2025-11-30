library(ggplot2)
library(dplyr)
library(data.table)

res_file <- "My_MR_Project/TCGA_Results/TCGA_Batch_Survival_Summary.csv"
if (!file.exists(res_file)) {
  stop("Result file not found. Run run_tcga_local.R first.")
}

data <- fread(res_file)
data$Significant <- ifelse(data$PValue < 0.05, "P < 0.05", "ns")
data$Cohort_Label <- ifelse(data$Cohort == "TCGA-LUAD", "Adenocarcinoma (LUAD)", "Squamous (LUSC)")

p <- ggplot(data, aes(x = HR, y = reorder(Gene, HR), color = Significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(size = 3.5) +
  geom_errorbarh(aes(xmin = HR - 0.1, xmax = HR + 0.1), height = 0.2, alpha = 0.5) +
  facet_wrap(~Cohort_Label, scales = "free_y") +
  scale_color_manual(values = c("ns" = "grey70", "P < 0.05" = "firebrick")) +
  labs(
    title = "Prognostic Landscape of sc-MR Candidates",
    subtitle = "Hazard Ratio > 1 implies Poor Prognosis (Risk Factor)",
    x = "Hazard Ratio (HR)",
    y = NULL,
    color = "Significance"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

out_file <- "My_MR_Project/TCGA_Results/Supplementary_Figure_Survival_Forest.pdf"
ggsave(out_file, p, width = 10, height = 7)
cat(paste("Forest plot saved to:", out_file, "\n"))

# 0. R
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
  library(readxl)
  library(pheatmap)
  library(ggrepel)
  library(xgboost)        
  library(SHAPforxgboost) 
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
})

# 1.
response_data <- read_excel("MSD/data/info/clinical_info.xlsx") %>%
  filter(pathological_response %in% c("MPR","non-MPR", "PD")) %>%
  filter(histology %in% c("LUSC")) %>%
  mutate(
    pathological_response = case_when(
      patientID %in% c("Pa169", "Pa160", "Pa83") ~ "MPR",
      pathological_response == "MPR" ~ "MPR",
      TRUE ~ "non-MPR"
    )
  ) %>% as.data.frame()

# 2. data.table
cat("\n===== Starting gene-level analysis =====\n")

# Pseudobulk
score <- read.csv("MSD/data/h5ad2/all_pseudobulk.csv") %>% as.data.frame()
rownames(score) <- score$X
score <- score[, -1, drop = FALSE]
score <- as.data.frame(t(score))

score$sampleID <- rownames(score)
score$patientID <- sub("S$", "", score$sampleID)
score$state <- ifelse(grepl("S$", score$sampleID), "post", "pre")

# prepost
paired_patients_gene <- score %>%
  group_by(patientID) %>%
  summarise(has_pre = "pre" %in% state, has_post = "post" %in% state, .groups = "drop") %>%
  filter(has_pre & has_post) %>%
  pull(patientID)

valid_patients <- intersect(paired_patients_gene, response_data$patientID)
cat("Gene-level analysis: ", length(valid_patients), " eligible paired patients selected\n")

score <- score %>% filter(patientID %in% valid_patients)

"IL4I1" %in% colnames(score)

df_wide_gene <- score %>%
  select(patientID, state, everything()) %>%
  select(-sampleID)

df_pre_gene <- df_wide_gene %>% filter(state == "pre") %>% arrange(patientID) %>% select(-state)
df_post_gene <- df_wide_gene %>% filter(state == "post") %>% arrange(patientID) %>% select(-state)

genes <- setdiff(colnames(df_pre_gene), "patientID")

eps <- 1e-4

# Log2FC
df_delta_gene <- df_pre_gene %>% select(patientID)
pre_mat <- as.matrix(df_pre_gene[, genes])
post_mat <- as.matrix(df_post_gene[, genes])
delta_mat <- log2((post_mat + eps) / (pre_mat + eps))
df_delta_gene <- cbind(df_delta_gene, as.data.frame(delta_mat))

df_pre_gene <- df_pre_gene %>% `rownames<-`(NULL)
df_post_gene <- df_post_gene %>% `rownames<-`(NULL)
df_delta_gene <- df_delta_gene %>% `rownames<-`(NULL)

resp <- response_data %>% filter(patientID %in% valid_patients) %>% arrange(patientID)
y_response <- ifelse(resp$pathological_response == "MPR", 1, 0)

mat_pre <- df_pre_gene %>% select(-patientID) %>% as.matrix() %>% `rownames<-`(df_pre_gene$patientID)
mat_delta <- df_delta_gene %>% select(-patientID) %>% as.matrix() %>% `rownames<-`(df_delta_gene$patientID)

# 1
mat_pre_rank <- apply(mat_pre, 2, rank)
rownames(mat_pre_rank) <- rownames(mat_pre)

mat_delta_rank <- apply(mat_delta, 2, rank)
rownames(mat_delta_rank) <- rownames(mat_delta)

# 2
xgb_params <- list(
  max_depth = 2,
  eta = 0.05,
  subsample = 0.8,
  colsample_bytree = 0.8,
  alpha = 0,
  lambda = 1,
  objective = "binary:logistic"
)

set.seed(123)
xgb_model_pre_gene <- xgboost(
  data = mat_pre_rank, 
  label = y_response, 
  params = xgb_params, 
  nrounds = 100,
  verbose = 0
)

set.seed(123)
xgb_model_delta_gene <- xgboost(
  data = mat_delta_rank, 
  label = y_response, 
  params = xgb_params, 
  nrounds = 100, 
  verbose = 0
)

# SHAP
cat("Computing SHAP values...\n")

# data.table
shap_mat_pre <- as.matrix(shap.values(xgb_model_pre_gene, mat_pre_rank)$shap_score)
shap_mat_delta <- as.matrix(shap.values(xgb_model_delta_gene, mat_delta_rank)$shap_score)

importance_pre <- colMeans(abs(shap_mat_pre))
importance_delta <- colMeans(abs(shap_mat_delta))

cor_pre <- sapply(1:ncol(mat_pre), function(i) cor(mat_pre[,i], shap_mat_pre[,i], method="spearman"))
cor_delta <- sapply(1:ncol(mat_delta), function(i) cor(mat_delta[,i], shap_mat_delta[,i], method="spearman"))

cor_pre[is.na(cor_pre)] <- 0
cor_delta[is.na(cor_delta)] <- 0
cor_pre[colSds(mat_pre) == 0] <- 0
cor_delta[colSds(mat_delta) == 0] <- 0

ml_results_gene <- tibble(
  feature = genes,
  Effect_Pre = importance_pre * sign(cor_pre),
  Effect_Delta = importance_delta * sign(cor_delta)
)

rm(shap_mat_pre, shap_mat_delta, mat_pre, mat_delta, mat_pre_rank, mat_delta_rank, pre_mat, post_mat, delta_mat)
gc()

# 3SHAP
shap_threshold <- 0.01

ml_results_gene <- ml_results_gene %>%
  mutate(
    Pre_Status = case_when(
      Effect_Pre > shap_threshold ~ "Pre: Good",
      Effect_Pre < -shap_threshold ~ "Pre: Bad",
      TRUE ~ "Pre: None"
    ),
    Delta_Status = case_when(
      Effect_Delta > shap_threshold ~ "Delta: Good",
      Effect_Delta < -shap_threshold ~ "Delta: Bad",
      TRUE ~ "Delta: None"
    ),
    Group = paste(Pre_Status, Delta_Status, sep = " | ")
  )

ml_results_filtered_gene <- ml_results_gene %>%
  filter(Group != "Pre: None | Delta: None") %>%
  arrange(desc(abs(Effect_Pre) + abs(Effect_Delta)))

cat("Gene-level analysis: ", nrow(ml_results_filtered_gene), " genes retained after filtering\n")
print(table(ml_results_filtered_gene$Group))

write.csv(ml_results_filtered_gene, "MSD/data/h5ad2/Gene_SHAP_Groups.csv", row.names = FALSE)

# 3.
p_scatter_gene <- ggplot(ml_results_filtered_gene, aes(x = Effect_Pre, y = Effect_Delta, color = Group)) +
  geom_hline(yintercept = c(-shap_threshold, shap_threshold), linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = c(-shap_threshold, shap_threshold), linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey90", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey90", size = 0.5) +
  
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = feature), size = 3, max.overlaps = 20, show.legend = FALSE) +
  
  scale_color_manual(values = c(
    "Pre: Good | Delta: Good" = "#D73027",
    "Pre: Bad | Delta: Bad"   = "#4575B4",
    "Pre: Good | Delta: Bad"  = "#F46D43",
    "Pre: Bad | Delta: Good"  = "#74ADD1",
    "Pre: Good | Delta: None" = "#FDAE61",
    "Pre: Bad | Delta: None"  = "#ABD9E9",
    "Pre: None | Delta: Good" = "#1A9850",
    "Pre: None | Delta: Bad"  = "#006837"
  )) +
  theme_bw(base_size = 14) +
  labs(
    title = "Gene Biomarker Quadrants (Decoupled XGBoost)",
    subtitle = paste("Threshold = ±", shap_threshold, "| Good = Promotes MPR, Bad = Inhibits MPR"),
    x = "SHAP Effect of Baseline Expression (Pre Model)",
    y = "SHAP Effect of Log2 Fold Change (Delta Model)",
    color = "Biomarker Role"
  ) +
  theme(legend.position = "right")

print(p_scatter_gene)
ggsave("MSD/figure/BEYOND/Gene_SHAP_Scatter_9Grid.pdf", plot = p_scatter_gene, width = 11, height = 7)

top_genes <- ml_results_filtered_gene %>% 
  top_n(50, abs(Effect_Pre) + abs(Effect_Delta)) %>%
  arrange(Group)

heatmap_mat_gene <- top_genes %>% 
  select(feature, Effect_Pre, Effect_Delta) %>% 
  column_to_rownames("feature") %>% 
  as.matrix()

colnames(heatmap_mat_gene) <- c("Baseline (Pre) Effect", "Change (Delta) Effect")

row_groups_gene <- top_genes$Group

row_ha_gene <- rowAnnotation(
  Biomarker_Role = row_groups_gene,
  col = list(Biomarker_Role = c(
    "Pre: Good | Delta: Good" = "#D73027",
    "Pre: Bad | Delta: Bad"   = "#4575B4",
    "Pre: Good | Delta: Bad"  = "#F46D43",
    "Pre: Bad | Delta: Good"  = "#74ADD1",
    "Pre: Good | Delta: None" = "#FDAE61",
    "Pre: Bad | Delta: None"  = "#ABD9E9",
    "Pre: None | Delta: Good" = "#1A9850",
    "Pre: None | Delta: Bad"  = "#006837"
  )),
  show_annotation_name = FALSE
)

max_val_gene <- max(abs(heatmap_mat_gene))
col_fun_gene <- colorRamp2(c(-max_val_gene, 0, max_val_gene), c("#4575B4", "white", "#D73027"))

ht_gene <- Heatmap(
  heatmap_mat_gene,
  name = "SHAP Effect",
  col = col_fun_gene,
  
  row_split = row_groups_gene,
  cluster_row_slices = FALSE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 8),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  
  cluster_columns = FALSE,
  show_column_names = TRUE,
  column_names_side = "top",
  column_names_rot = 0,
  column_names_centered = TRUE,
  
  left_annotation = row_ha_gene,
  border = TRUE,
  row_gap = unit(3, "mm"),
  
  column_title = "Gene Biomarker Roles (Decoupled XGBoost)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)

pdf("MSD/figure/BEYOND/Gene_SHAP_Heatmap_Grouped.pdf", width = 8, height = 15)
draw(ht_gene, padding = unit(c(2, 15, 2, 2), "mm"))
dev.off()

cat("Gene-level heatmap saved to MSD/figure/BEYOND/Gene_SHAP_Heatmap_Grouped.pdf\n")

# 4.
top3_genes <- ml_results_filtered_gene %>% 
  top_n(3, abs(Effect_Pre) + abs(Effect_Delta)) %>%
  pull(feature)

for(target_gene in top3_genes) {
  plot_df <- score_filtered %>%
    filter(patientID %in% valid_patients) %>%
    select(patientID, state, all_of(target_gene)) %>%
    left_join(resp %>% select(patientID, pathological_response), by = "patientID") %>%
    mutate(state = factor(state, levels = c("pre", "post")))
  
  colnames(plot_df)[3] <- "expression"
  
  p_verify <- ggplot(plot_df, aes(x = state, y = expression, color = pathological_response)) +
    geom_boxplot(aes(fill = pathological_response), alpha = 0.2, outlier.shape = NA, width = 0.5) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(aes(group = patientID), alpha = 0.5, size = 1) +
    facet_wrap(~ pathological_response) +
    scale_color_manual(values = c("MPR" = "#3BA0B0", "non-MPR" = "#D75A5A")) +
    scale_fill_manual(values = c("MPR" = "#3BA0B0", "non-MPR" = "#D75A5A")) +
    theme_bw(base_size = 14) +
    labs(
      title = paste("Expression Changes of", target_gene),
      x = "Treatment State",
      y = "Normalized Expression"
    ) +
    theme(legend.position = "none")
  
  print(p_verify)
  ggsave(paste0("MSD/figure/BEYOND/Gene_", target_gene, "_Changes.pdf"), 
         plot = p_verify, width = 8, height = 6)
}

cat("\n===== All gene-level analyses completed =====\n")
cat("Gene result table: MSD/data/h5ad2/Gene_SHAP_Groups.csv\n")

target_gene <- "IL4I1" # "CXCL13", "GZMK", "TNFRSF9"

plot_df <- score_filtered %>%
  filter(patientID %in% valid_patients) %>%
  select(patientID, state, all_of(target_gene)) %>%
  left_join(resp %>% select(patientID, pathological_response), by = "patientID") %>%
  mutate(state = factor(state, levels = c("pre", "post")))

p_verify_gene <- ggplot(plot_df, aes(x = state, y = .data[[target_gene]], color = pathological_response)) +
  geom_boxplot(aes(fill = pathological_response), alpha = 0.2, outlier.shape = NA, width = 0.5) +
  geom_point(size = 3, alpha = 0.8) +
  geom_line(aes(group = patientID), alpha = 0.5, size = 1) +
  facet_wrap(~ pathological_response) +
  scale_color_manual(values = c("MPR" = "#3BA0B0", "non-MPR" = "#D75A5A")) +
  scale_fill_manual(values = c("MPR" = "#3BA0B0", "non-MPR" = "#D75A5A")) +
  theme_bw(base_size = 14) +
  labs(
    title = paste("Real Expression Changes of", target_gene),
    x = "Treatment State",
    y = "Normalized Expression"
  ) +
  theme(legend.position = "none")

print(p_verify_gene)

ggsave(paste0("MSD/figure/BEYOND/Gene_", target_gene, "_Changes.pdf"), 
       plot = p_verify_gene, width = 8, height = 6)

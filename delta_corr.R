# CD4T_Treg_TNFRSF9high Δ Δ

# 0.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readxl)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
})

# 1. &
OUTDIR <- "MSD/figure/Treg_TNFRSF9high_GeneCor"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

target_cell <- "CD4T_Treg_TNFRSF9high"
eps         <- 1e-4 # log2

gene_list <- c(
  "MT1G", "IL4I1", "PTGS1", "AREG", "SYAP1", "MAN2B1", "OTUD1", "C4orf48", "CD4", "TMEM176A",
  "PLAT", "APOD", "TMEM45A", "CERCAM", "EVA1B", "RARRES2", "CHI3L1", "SELENOM", "RAB9A", "FUCA2",
  "GORASP2", "TFPI", "AKR1A1", "HLA-DPA1", "HLA-DQB1", "GPR183", "NANS", "NPDC1", "EGLN3", "USP12"
)

# 2.
response_data <- read_excel("MSD/data/info/clinical_info.xlsx") %>%
  filter(pathological_response %in% c("MPR", "non-MPR", "PD"),
         histology %in% c("LUSC")) %>%
  mutate(
    pathological_response = ifelse(pathological_response == "PD", "non-MPR", pathological_response),
    pathological_response = factor(pathological_response, levels = c("MPR", "non-MPR"))
  ) %>% as.data.frame()

cat("Clinical information: ", nrow(response_data), " patients\n")
print(table(response_data$pathological_response))

# 3. T
process_cell_freq <- function(file_path, cell_type_col) {
  df <- read.csv(file_path) %>% dplyr::select(sampleID, all_of(cell_type_col))
  table(df$sampleID, df[[cell_type_col]]) %>% as.data.frame() %>%
    dplyr::mutate(Freq = Freq / ave(Freq, Var1, FUN = sum)) %>%
    dplyr::rename(sampleID = Var1, cell.type = Var2)
}

all_df <- process_cell_freq("MSD/data/h5ad2/T_harmony.csv", "sub_cell_type")

all_df$patientID <- sub("S$", "", all_df$sampleID)
all_df$state     <- ifelse(grepl("S$", all_df$sampleID), "post", "pre")

# 4. pre/post &
paired_patients <- all_df %>%
  group_by(patientID) %>%
  summarise(has_pre  = "pre"  %in% state,
            has_post = "post" %in% state, .groups = "drop") %>%
  filter(has_pre & has_post) %>% pull(patientID)

valid_patients <- intersect(paired_patients, response_data$patientID)
cat("Paired patients with clinical information: ", length(valid_patients), "\n")

# 5.
df_wide <- all_df %>% filter(patientID %in% valid_patients) %>%
  dplyr::select(patientID, state, cell.type, Freq) %>%
  pivot_wider(names_from = cell.type, values_from = Freq, values_fill = 0) %>%
  as.data.frame()

df_pre  <- df_wide %>% filter(state == "pre")  %>% arrange(patientID) %>% dplyr::select(-state)
df_post <- df_wide %>% filter(state == "post") %>% arrange(patientID) %>% dplyr::select(-state)
stopifnot(identical(df_pre$patientID, df_post$patientID))

# 6. all_pseudobulk
score <- read.csv("MSD/data/h5ad2/Myeloid_pseudobulk.csv")
rownames(score) <- score$X
score <- score[, -1, drop = FALSE]
score <- as.data.frame(t(score))
score$sampleID <- rownames(score)

# 7. CD4T_Treg_TNFRSF9high Δ
if (!target_cell %in% colnames(df_pre)) {
  stop("Cell type not found: ", target_cell,
       "\nAvailable cell types: ", paste(setdiff(colnames(df_pre), "patientID"), collapse = ", "))
}

freq_pre   <- df_pre[[target_cell]]
freq_post  <- df_post[[target_cell]]
delta_cell <- log2((freq_post + eps) / (freq_pre + eps)) # freq_post - freq_pre

cell_df <- data.frame(
  patientID  = df_pre$patientID,
  freq_pre   = freq_pre,
  freq_post  = freq_post,
  delta_cell = delta_cell
)

# 8. Δ
score_pb <- score %>%
  mutate(
    state     = ifelse(grepl("S$", sampleID), "post", "pre"),
    patientID = sub("S$", "", sampleID)
  ) %>%
  filter(patientID %in% cell_df$patientID)

# pre / post
paired_pat_gene <- score_pb %>% group_by(patientID) %>%
  filter(all(c("pre", "post") %in% state)) %>% pull(patientID) %>% unique()

common_pats <- intersect(cell_df$patientID, paired_pat_gene)
cat("\nPaired patients (with both cell frequency and gene expression): ", length(common_pats), "\n")

gene_in_data <- intersect(gene_list, colnames(score_pb))
cat("Genes found: ", length(gene_in_data), "/", length(gene_list), "\n")
missing_genes <- setdiff(gene_list, gene_in_data)
if (length(missing_genes) > 0) cat("Missing genes: ", paste(missing_genes, collapse = ", "), "\n")

gene_pre  <- score_pb %>% filter(state == "pre",  patientID %in% common_pats) %>%
  arrange(patientID) %>% dplyr::select(patientID, all_of(gene_in_data))
gene_post <- score_pb %>% filter(state == "post", patientID %in% common_pats) %>%
  arrange(patientID) %>% dplyr::select(patientID, all_of(gene_in_data))
stopifnot(identical(gene_pre$patientID, gene_post$patientID))

gene_pre_mat  <- as.matrix(gene_pre[,  -1])
gene_post_mat <- as.matrix(gene_post[, -1])
delta_gene    <- log2((gene_post_mat + eps) / (gene_pre_mat + eps))
rownames(delta_gene) <- gene_pre$patientID

# 9. + response
cell_sub <- cell_df %>% filter(patientID %in% common_pats) %>% arrange(patientID)
stopifnot(identical(cell_sub$patientID, rownames(delta_gene)))

resp_sub <- response_data %>% filter(patientID %in% common_pats) %>% arrange(patientID)
stopifnot(identical(as.character(resp_sub$patientID), cell_sub$patientID))

cat("\nFinal number of patients used for correlation: ", nrow(cell_sub), "\n")
print(table(resp_sub$pathological_response))

# 10. Pearson + Spearman
cor_results <- lapply(colnames(delta_gene), function(g) {
  x <- cell_sub$delta_cell
  y <- delta_gene[, g]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4) return(NULL)
  pe <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
  sp <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  data.frame(
    gene         = g,
    n            = sum(ok),
    pearson_r    = unname(pe$estimate),
    pearson_p    = pe$p.value,
    spearman_rho = unname(sp$estimate),
    spearman_p   = sp$p.value
  )
}) %>% bind_rows() %>%
  mutate(
    pearson_fdr  = p.adjust(pearson_p,  method = "BH"),
    spearman_fdr = p.adjust(spearman_p, method = "BH")
  ) %>%
  arrange(pearson_p)

write.csv(cor_results,
          file.path(OUTDIR, paste0(target_cell, "_GeneDelta_Cor.csv")),
          row.names = FALSE)

cat("\n--- Top correlated genes (by Pearson p) ---\n")
print(head(cor_results, 15))

# 11. ①
lol_df <- cor_results %>%
  mutate(
    sig = case_when(
      pearson_p < 0.01 ~ "**",
      pearson_p < 0.05 ~ "*",
      pearson_p < 0.1  ~ ".",
      TRUE             ~ ""
    ),
    direction = ifelse(pearson_r > 0, "Positive", "Negative")
  ) %>%
  arrange(pearson_r) %>%
  mutate(gene = factor(gene, levels = gene))

p_lol <- ggplot(lol_df, aes(x = pearson_r, y = gene, color = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_segment(aes(x = 0, xend = pearson_r, yend = gene), linewidth = 0.5) +
  geom_point(aes(size = -log10(pearson_p + 1e-10)), alpha = 0.85) +
  geom_text(aes(label = sig, x = pearson_r + sign(pearson_r) * 0.04),
            color = "black", size = 4.5, fontface = "bold") +
  scale_color_manual(values = c("Positive" = "#B2182B", "Negative" = "#2166AC")) +
  scale_size_continuous(range = c(1.5, 5)) +
  labs(
    title    = paste0("Correlation: Δ", target_cell, " vs Δ Gene"),
    subtitle = "log2(Post/Pre) — Pearson; * p<0.05  ** p<0.01",
    x = "Pearson r", y = NULL,
    color = "Direction", size = expression(-log[10](p))
  ) +
  theme_pubr() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        axis.text.y   = element_text(face = "italic"))

ggsave(file.path(OUTDIR, paste0(target_cell, "_Lollipop.pdf")),
       p_lol, width = 7, height = 0.25 * nrow(lol_df) + 2)

# 12. ② +
plot_scatter <- function(g) {
  df_p <- data.frame(
    patientID  = cell_sub$patientID,
    delta_cell = cell_sub$delta_cell,
    delta_gene = delta_gene[, g],
    Response   = resp_sub$pathological_response
  ) %>% filter(is.finite(delta_cell) & is.finite(delta_gene))
  
  stat <- cor_results %>% filter(gene == g)
  sub  <- paste0("r = ", signif(stat$pearson_r, 2),
                 ",  p = ", signif(stat$pearson_p, 2),
                 "   |   rho = ", signif(stat$spearman_rho, 2),
                 ",  p = ", signif(stat$spearman_p, 2))
  
  ggplot(df_p, aes(x = delta_cell, y = delta_gene)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    geom_smooth(method = "lm", se = TRUE, color = "grey30",
                fill = "grey80", alpha = 0.4, linewidth = 0.6) +
    geom_point(color = "#6B7FA6", size = 2.4, alpha = 0.9) +
    labs(title    = g,
         subtitle = sub,
         x = paste0("Δ ", target_cell, "  (log2 Post/Pre)"),
         y = paste0("Δ ", g, "  (log2 Post/Pre)")) +
    theme_pubr() +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold.italic"),
          plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey30"),
          legend.position = "right")
}

scatter_list <- lapply(colnames(delta_gene), plot_scatter)
names(scatter_list) <- colnames(delta_gene)

for (g in names(scatter_list)) {
  ggsave(file.path(OUTDIR, paste0("Scatter_", target_cell, "_vs_", g, ".pdf")),
         scatter_list[[g]], width = 4, height = 4.2)
}

comb <- wrap_plots(scatter_list, ncol = 5) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(OUTDIR, paste0("ALL_Scatter_", target_cell, "_vs_Genes.pdf")),
       comb, width = 22, height = 4.2 * ceiling(length(scatter_list)/5),
       limitsize = FALSE)

# 13. ③ p<0.05
sig_genes <- cor_results %>% filter(pearson_p < 0.05) %>% pull(gene)
cat("\nSignificant genes at Pearson p<0.05: ", length(sig_genes), " -> ",
    paste(sig_genes, collapse = ", "), "\n")

if (length(sig_genes) > 0) {
  comb_sig <- wrap_plots(scatter_list[sig_genes],
                         ncol = min(3, length(sig_genes))) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path(OUTDIR, paste0("SIG_Scatter_", target_cell, ".pdf")),
         comb_sig, width = 5 * min(3, length(sig_genes)),
         height = 4.2 * ceiling(length(sig_genes)/3), limitsize = FALSE)
}

# 14. ④
heat_mat <- cor_results %>%
  arrange(desc(pearson_r)) %>%
  dplyr::select(gene, pearson_r) %>%
  tibble::column_to_rownames("gene") %>% as.matrix()

sig_mat <- cor_results %>%
  arrange(desc(pearson_r)) %>%
  mutate(label = case_when(
    pearson_p < 0.01 ~ "**",
    pearson_p < 0.05 ~ "*",
    pearson_p < 0.1  ~ ".",
    TRUE             ~ ""
  )) %>% dplyr::select(gene, label) %>%
  tibble::column_to_rownames("gene") %>% as.matrix()

pheatmap(
  heat_mat,
  cluster_rows = FALSE, cluster_cols = FALSE,
  color  = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  display_numbers = sig_mat,
  number_color = "black", fontsize_number = 10,
  border_color = "white", fontsize_row = 9,
  main = paste0("Δ", target_cell, " vs Δ Gene (Pearson r)"),
  filename = file.path(OUTDIR, paste0(target_cell, "_CorHeatmap.pdf")),
  width = 4.5, height = 0.28 * nrow(heat_mat) + 2
)

# 15.
cat("\n========================================\n")
cat("Correlation analysis completed: ", target_cell, " vs ", length(gene_in_data), " genes\n", sep = "")
cat("- Number of paired samples: ", length(common_pats), "\n", sep = "")
cat("- Pearson  p<0.05: ", sum(cor_results$pearson_p  < 0.05, na.rm = TRUE), "\n")
cat("- Spearman p<0.05: ", sum(cor_results$spearman_p < 0.05, na.rm = TRUE), "\n")
cat("- Output directory: ", OUTDIR, "\n", sep = "")
cat("========================================\n")

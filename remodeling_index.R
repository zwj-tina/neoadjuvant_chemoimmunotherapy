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
})

 # 1-7. pre-process
process_cell_freq <- function(file_path, cell_type_col) {
  df <- read.csv(file_path) %>% dplyr::select(sampleID, all_of(cell_type_col))
  table(df$sampleID, df[[cell_type_col]]) %>% as.data.frame() %>%
    dplyr::mutate(Freq = Freq / ave(Freq, Var1, FUN = sum)) %>%
    dplyr::rename(sampleID = Var1, cell.type = Var2)
}

T_df       <- process_cell_freq("./T_harmony.csv", "sub_cell_type2")
B_df       <- process_cell_freq("./B_harmony.csv", "sub_cell_type")
Myeloid_df <- process_cell_freq("./Myeloid_harmony.csv", "sub_cell_type")
Endo_df    <- process_cell_freq("./Endothelial_harmony.csv", "sub_cell_type")
Fibro_df   <- process_cell_freq("./Fibroblast_harmony.csv", "sub_cell_type2")
Epi_df     <- process_cell_freq("./Epithelial_harmony.csv", "sub_cell_type")

all_df <- T_df
all_df$patientID <- sub("S$", "", all_df$sampleID)
all_df$state <- ifelse(grepl("S$", all_df$sampleID), "post", "pre")

response_data <- read_excel("./clinical_info.xlsx") %>%
  filter(histology %in% c("LUSC")) %>%
  filter(!is.na(pathological_response_rate)) %>%
  mutate(
    pathological_response_rate = as.numeric(gsub("%", "", pathological_response_rate))
  ) %>%
  filter(pathological_response_rate >= 0 & pathological_response_rate <= 100) %>%
  as.data.frame()

PATHE <- read.csv("./PATHE.csv")
rownames(PATHE) <- PATHE$sampleID
PATHE$clusters <- as.character(PATHE$clusters)
table(PATHE$clusters, PATHE$pathological_response)
PATHE <- PATHE[PATHE$clusters %in% c("1"),]
PATHE <- PATHE[PATHE$pathological_response %in% c("non-MPR"),]
head(PATHE)

all_df <- all_df[all_df$patientID %in% PATHE$patientID,]

paired_patients <- all_df %>%
  group_by(patientID) %>%
  summarise(has_pre = "pre" %in% state, has_post = "post" %in% state, .groups = "drop") %>%
  filter(has_pre & has_post) %>% pull(patientID)

valid_patients <- intersect(paired_patients, response_data$patientID)

 # step1
df_wide <- all_df %>%
  filter(patientID %in% valid_patients) %>%
  select(patientID, state, cell.type, Freq) %>%
  pivot_wider(names_from = cell.type, values_from = Freq, values_fill = 0)

df_pre <- df_wide %>%
  filter(state == "pre") %>%
  arrange(patientID) %>%
  select(-state)

df_post <- df_wide %>%
  filter(state == "post") %>%
  arrange(patientID) %>%
  select(-state)

cell_types <- setdiff(colnames(df_pre), "patientID")

 # step2
patient_remodeling <- lapply(valid_patients, function(pid) {
  pre_vec <- as.numeric(df_pre[df_pre$patientID == pid, cell_types])
  post_vec <- as.numeric(df_post[df_post$patientID == pid, cell_types])
  
  dist <- sqrt(sum((post_vec - pre_vec)^2))
  
  data.frame(patientID = pid, remodeling_index = dist)
}) %>% bind_rows()

head(patient_remodeling)

 # step3
delta_list <- lapply(cell_types, function(ct) {
  pre_v <- as.numeric(df_pre[[ct]])
  post_v <- as.numeric(df_post[[ct]])
  data.frame(
    patientID = df_pre$patientID,
    cell.type = ct,
    delta = post_v - pre_v  
  )
}) %>% bind_rows()

delta_wide <- delta_list %>%
  pivot_wider(names_from = cell.type, values_from = delta)

 # step4
cell_contribution <- delta_list %>%
  left_join(patient_remodeling, by = "patientID") %>%
  group_by(cell.type) %>%
  summarise(
    cor = cor.test(delta, remodeling_index, method = "spearman")$estimate,
    pval = cor.test(delta, remodeling_index, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(abs_cor = abs(cor)) %>%
  arrange(desc(abs_cor))

head(cell_contribution, 20)

 #cell_contribution <- cell_contribution[cell_contribution$pval < 0.05,]

library(ComplexHeatmap)
library(circlize)

 # step4
mat4 <- cell_contribution %>%
  select(cell.type, cor) %>%
  column_to_rownames("cell.type") %>%
  as.matrix()

p_mat4 <- cell_contribution %>%
  select(cell.type, pval) %>%
  column_to_rownames("cell.type") %>%
  as.matrix()

max_val <- max(abs(mat4))
col_fun <- colorRamp2(c(-max_val, 0, max_val), c("#4575B4", "white", "#D73027"))

ht4 <- Heatmap(
  mat4,
  name = "Spearman r",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  cell_fun = function(j, i, x, y, width, height, fill) {
    p <- p_mat4[i, j]
    if(p < 0.001) {
      grid.text("***", x, y, gp = gpar(fontsize = 10, col = "black"))
    } else if(p < 0.01) {
      grid.text("**", x, y, gp = gpar(fontsize = 10, col = "black"))
    } else if(p < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 10, col = "black"))
    }
  },
  column_title = "Cell Delta vs Remodeling Index (Continuous)"
)

pdf("/home/zwj/raid1/MSD/figure/BEYOND/Treg_Step4_Cell_Remodeling_Cor_Heatmap_sig.pdf", width = 4, height = 6)
draw(ht4)
dev.off()

 # step5
remodel_response <- patient_remodeling %>%
  left_join(response_data %>% select(patientID, pathological_response_rate), by = "patientID")

cor_test <- cor.test(
  remodel_response$remodeling_index,
  remodel_response$pathological_response_rate,
  method = "spearman"
)

cor_test

ggplot(remodel_response,
       aes(x = remodeling_index, y = pathological_response_rate)) +
  geom_point(size=3, alpha=0.7) +
  geom_smooth(method="lm", color="red") +
  labs(
    x = "Continuous Tumor Microenvironment Remodeling Index",
    y = "Pathological Response Rate (%)",
    title = paste("rho =", round(cor_test$estimate, 2),
                  "p =", signif(cor_test$p.value, 2))
  ) +
  theme_bw()
ggsave("/home/zwj/raid1/MSD/figure/BEYOND/Treg_remodeling_index_PRR.pdf",width = 4, height = 4)

 # step6
cell_response_corr <- delta_list %>%
  left_join(response_data %>% select(patientID, pathological_response_rate), by = "patientID") %>%
  group_by(cell.type) %>%
  summarise(
    cor = cor.test(delta, pathological_response_rate, method = "spearman")$estimate,
    pval = cor.test(delta, pathological_response_rate, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  arrange(pval)

head(cell_response_corr, 30)

mat6 <- cell_response_corr %>%
  select(cell.type, cor) %>%
  column_to_rownames("cell.type") %>%
  as.matrix()

p_mat6 <- cell_response_corr %>%
  select(cell.type, pval) %>%
  column_to_rownames("cell.type") %>%
  as.matrix()

max_val <- max(abs(mat6))
col_fun <- colorRamp2(c(-max_val, 0, max_val), c("#4575B4", "white", "#D73027"))

ht6 <- Heatmap(
  mat6,
  name = "Spearman r",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  cell_fun = function(j, i, x, y, width, height, fill) {
    p <- p_mat6[i, j]
    if(p < 0.001) {
      grid.text("***", x, y, gp = gpar(fontsize = 10, col = "black"))
    } else if(p < 0.01) {
      grid.text("**", x, y, gp = gpar(fontsize = 10, col = "black"))
    } else if(p < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 10, col = "black"))
    }
  },
  column_title = "Cell Delta vs Pathological Response Rate"
)

pdf("./Treg_Step6_Cell_Response_Cor_Heatmap.pdf", width = 4, height = 12)
draw(ht6)
dev.off()

 # delta
score <- read.csv("./Myeloid_pseudobulk.csv") %>% as.data.frame()
rownames(score) <- score$X
score <- score[, -1, drop = FALSE]

 #score <- scale(score)
score[1:5,1:5]
clean_gene_names <- rownames(score)

 # delete
remove_patterns <- c(
  "^MT-", "^MT\\.",
  "^RPL|^RPS",
  "^AC[0-9]|^AL[0-9]", # AC00... AL00...
  "^LINC", "^LINC\\.", # LINC
  "^RN7SL", "^SNOR", # RNA
  "^NEAT1$", "^MALAT1$",
  "^KRTAP", "^NME\\.",
  "AS1$|AS2$",
  "\\."
)

for(p in remove_patterns){
  clean_gene_names <- clean_gene_names[!grepl(p, clean_gene_names, ignore.case = F)]
}

score <- score[rownames(score) %in% clean_gene_names, ]

score[1:5,1:5]

score <- as.data.frame(t(score))

score$sampleID <- rownames(score)
score$patientID <- sub("S$", "", score$sampleID)
score$state <- ifelse(grepl("S$", score$sampleID), "post", "pre")

paired_patients_gene <- score %>%
  group_by(patientID) %>%
  summarise(has_pre = "pre" %in% state, has_post = "post" %in% state, .groups = "drop") %>%
  filter(has_pre & has_post) %>%
  pull(patientID)

valid_patients <- intersect(paired_patients_gene, response_data$patientID)

score <- score %>% filter(patientID %in% valid_patients)

df_wide_gene <- score %>%
  select(patientID, state, everything()) %>%
  select(-sampleID)

df_pre_gene <- df_wide_gene %>%
  filter(state == "pre") %>%
  arrange(patientID) %>%
  select(-state)

df_post_gene <- df_wide_gene %>%
  filter(state == "post") %>%
  arrange(patientID) %>%
  select(-state)

genes <- setdiff(colnames(df_pre_gene), "patientID")

df_delta_gene <- df_pre_gene %>% select(patientID)
for (g in genes) {
  df_delta_gene[[g]] <- df_post_gene[[g]] - df_pre_gene[[g]]
}

df_delta_gene <- df_delta_gene %>%
  filter(patientID %in% patient_remodeling$patientID)

df_delta_gene[1:5,1:5]

remodel_vec <- patient_remodeling %>%
  filter(patientID %in% df_delta_gene$patientID) %>%
  arrange(patientID) %>%
  pull(remodeling_index)

gene_remodel_corr <- lapply(genes, function(g) {
  vals <- df_delta_gene[[g]]
  ct <- cor.test(vals, remodel_vec, method = "spearman")
  tibble(
    gene = g,
    cor = ct$estimate,
    pval = ct$p.value
  )
}) %>% bind_rows()

gene_remodel_corr[gene_remodel_corr$gene %in% c("IL4I1"),]

gene_sig <- gene_remodel_corr %>%
  filter(pval < 0.01) %>%
  arrange(desc(abs(cor)))

library(ComplexHeatmap)
library(circlize)

mat_gene <- gene_sig %>%
  select(gene, cor) %>%
  column_to_rownames("gene") %>%
  as.matrix()

p_mat_gene <- gene_sig %>%
  select(gene, pval) %>%
  column_to_rownames("gene") %>%
  as.matrix()

max_val <- max(abs(mat_gene))
col_fun <- colorRamp2(c(-max_val, 0, max_val), c("#4575B4", "white", "#D73027"))

ht_gene <- Heatmap(
  mat_gene,
  name = "Spearman r",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, width, height, fill) {
    p <- p_mat_gene[i, j]
    if(p < 0.001) {
      grid.text("***", x, y, gp = gpar(col="black", fontsize=10))
    } else if(p < 0.01) {
      grid.text("**", x, y, gp = gpar(col="black", fontsize=10))
    } else if(p < 0.05) {
      grid.text("*", x, y, gp = gpar(col="black", fontsize=10))
    }
  },
  column_title = "Gene Delta vs Remodeling Index (sig < 0.05)"
)

pdf("./Treg_Myeloid_Gene_Delta_vs_Remodeling_Heatmap.pdf",
    width = 4, height = 20)
draw(ht_gene)
dev.off()

plot_single_gene_remodel <- function(target_gene){
  
 # 1.
  plot_df <- data.frame(
    patientID = df_delta_gene$patientID,
    logFC = df_delta_gene[[target_gene]],
    remodeling_index = remodel_vec
  )
  
 # 2.
  ct <- cor.test(plot_df$logFC, plot_df$remodeling_index, method = "spearman")
  rho <- round(ct$estimate, 3)
  pval <- signif(ct$p.value, 3)
  
 # 3.
  p <- ggplot(plot_df, aes(x = logFC, y = remodeling_index)) +
    geom_point(size = 3, alpha = 0.7, color = "#3A7BB5") +
    geom_smooth(method = "lm", color = "#E65656", fill = "#FFE5E5") +
    labs(
      x = paste0(target_gene, " Post vs Pre"),
      y = "Tumor Microenvironment Remodeling Index",
      title = paste0("rho = ", rho, " | p = ", pval)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 11)
    )
  
  ggsave(
    paste0("./Gene_", target_gene, "_logFC_vs_Remodeling.pdf"),
    plot = p, width = 5, height = 4
  )
  
}

# ==========================================
# 0. 加载所有必要的R包
# ==========================================
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

# ==========================================
# 1-7. 数据读取与预处理
# ==========================================
process_cell_freq <- function(file_path, cell_type_col) {
  df <- read.csv(file_path) %>% dplyr::select(sampleID, all_of(cell_type_col))
  table(df$sampleID, df[[cell_type_col]]) %>% as.data.frame() %>%
    dplyr::mutate(Freq = Freq / ave(Freq, Var1, FUN = sum)) %>%
    dplyr::rename(sampleID = Var1, cell.type = Var2)
}

T_df       <- process_cell_freq("/home/zwj/raid1/MSD/data/h5ad2/T_harmony.csv", "sub_cell_type2")
B_df       <- process_cell_freq("/home/zwj/raid1/MSD/data/h5ad2/B_harmony.csv", "sub_cell_type")
Myeloid_df <- process_cell_freq("/home/zwj/raid1/MSD/data/h5ad2/Myeloid_harmony.csv", "sub_cell_type")
Endo_df    <- process_cell_freq("/home/zwj/raid1/MSD/data/h5ad2/Endothelial_harmony.csv", "sub_cell_type")
Fibro_df   <- process_cell_freq("/home/zwj/raid1/MSD/data/h5ad2/Fibroblast_harmony.csv", "sub_cell_type2")
Epi_df     <- process_cell_freq("/home/zwj/raid1/MSD/data/h5ad2/Epithelial_harmony.csv", "sub_cell_type")

all_df <- T_df
all_df$patientID <- sub("S$", "", all_df$sampleID)
all_df$state <- ifelse(grepl("S$", all_df$sampleID), "post", "pre")

# ==============================
# 【关键修改】使用连续病理缓解率，不再二分类
# ==============================
response_data <- read_excel("/home/zwj/raid1/MSD/data/info/clinical_info.xlsx") %>%
  filter(histology %in% c("LUSC")) %>%
  filter(!is.na(pathological_response_rate)) %>%
  mutate(
    pathological_response_rate = as.numeric(gsub("%", "", pathological_response_rate))
  ) %>%
  filter(pathological_response_rate >= 0 & pathological_response_rate <= 100) %>%
  as.data.frame()


PATHE <- read.csv("/home/zwj/raid1/MSD/data/h5ad2/PATHE.csv")
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

# ------------------------------------------------------
# 步骤1：筛选有效患者 → 变成宽格式（pre/post分开）
# ------------------------------------------------------
df_wide <- all_df %>%
  filter(patientID %in% valid_patients) %>%
  select(patientID, state, cell.type, Freq) %>%
  pivot_wider(names_from = cell.type, values_from = Freq, values_fill = 0)

# 拆分 pre 和 post
df_pre <- df_wide %>%
  filter(state == "pre") %>%
  arrange(patientID) %>%
  select(-state)

df_post <- df_wide %>%
  filter(state == "post") %>%
  arrange(patientID) %>%
  select(-state)

# 统一细胞类型顺序
cell_types <- setdiff(colnames(df_pre), "patientID")

# ------------------------------------------------------
# 步骤2：计算每个患者的连续重塑指数
# ------------------------------------------------------
patient_remodeling <- lapply(valid_patients, function(pid) {
  pre_vec <- as.numeric(df_pre[df_pre$patientID == pid, cell_types])
  post_vec <- as.numeric(df_post[df_post$patientID == pid, cell_types])
  
  # 欧氏距离 = 整体微环境变化大小（连续值）
  dist <- sqrt(sum((post_vec - pre_vec)^2))
  
  data.frame(patientID = pid, remodeling_index = dist)
}) %>% bind_rows()

# 查看每个患者的“重塑程度”（完全连续）
head(patient_remodeling)

# ------------------------------------------------------
# 步骤3：计算每个细胞在每个患者的变化量 Δ = Post - Pre
# ------------------------------------------------------
delta_list <- lapply(cell_types, function(ct) {
  pre_v <- as.numeric(df_pre[[ct]])
  post_v <- as.numeric(df_post[[ct]])
  data.frame(
    patientID = df_pre$patientID,
    cell.type = ct,
    delta = post_v - pre_v  # 连续变化量
  )
}) %>% bind_rows()

# 变成宽格式（方便后续建模）
delta_wide <- delta_list %>%
  pivot_wider(names_from = cell.type, values_from = delta)

# ------------------------------------------------------
# 步骤4：细胞对重塑的贡献度（连续变量相关性）
# ------------------------------------------------------
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

# 输出：驱动重塑的细胞排名
head(cell_contribution, 20)

#cell_contribution <- cell_contribution[cell_contribution$pval < 0.05,]

# ==========================================
# 加载热图包
# ==========================================
library(ComplexHeatmap)
library(circlize)

# ==========================================
# 【步骤4 可视化】细胞 → 重塑指数 相关性热图
# ==========================================
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




# ------------------------------------------------------
# 步骤5：重塑指数 ↔ 病理缓解率（连续相关）
# ------------------------------------------------------
remodel_response <- patient_remodeling %>%
  left_join(response_data %>% select(patientID, pathological_response_rate), by = "patientID")

# 相关性检验
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

# ------------------------------------------------------
# 步骤6：细胞变化 ↔ 病理缓解率（连续）
# ------------------------------------------------------
cell_response_corr <- delta_list %>%
  left_join(response_data %>% select(patientID, pathological_response_rate), by = "patientID") %>%
  group_by(cell.type) %>%
  summarise(
    cor = cor.test(delta, pathological_response_rate, method = "spearman")$estimate,
    pval = cor.test(delta, pathological_response_rate, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  arrange(pval)

# 输出：与缓解率最相关的细胞变化
head(cell_response_corr, 30)

# ==========================================
# 【步骤6 可视化】细胞变化 → 病理缓解率 相关性热图
# ==========================================
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

pdf("/home/zwj/raid1/MSD/figure/BEYOND/Treg_Step6_Cell_Response_Cor_Heatmap.pdf", width = 4, height = 12)
draw(ht6)
dev.off()

# ==========================================
# 【基因水平】计算治疗前后变化 delta
# ==========================================
score <- read.csv("/home/zwj/raid1/MSD/data/h5ad2/Myeloid_pseudobulk.csv") %>% as.data.frame()
rownames(score) <- score$X
score <- score[, -1, drop = FALSE]

#score <- scale(score)
score[1:5,1:5]
# ==============================
clean_gene_names <- rownames(score)

# 黑名单模式：去掉以下开头/包含的基因
remove_patterns <- c(
  "^MT-", "^MT\\.",       # 线粒体
  "^RPL|^RPS",           # 核糖体
  "^AC[0-9]|^AL[0-9]",   # AC00... AL00... 乱码基因
  "^LINC", "^LINC\\.",    # LINC 长链非编码
  "^RN7SL", "^SNOR",     # RNA 家族
  "^NEAT1$", "^MALAT1$",  # 经典长非编码
  "^KRTAP", "^NME\\.",    # 乱七八糟
  "AS1$|AS2$",           # 反义链
  "\\."                   # 带点的杂基因
)

# 批量过滤
for(p in remove_patterns){
  clean_gene_names <- clean_gene_names[!grepl(p, clean_gene_names, ignore.case = F)]
}

# 只保留干净基因
score <- score[rownames(score) %in% clean_gene_names, ]

score[1:5,1:5]

score <- as.data.frame(t(score))

score$sampleID <- rownames(score)
score$patientID <- sub("S$", "", score$sampleID)
score$state <- ifelse(grepl("S$", score$sampleID), "post", "pre")

# 筛选同时拥有pre和post配对数据且符合临床条件的患者
paired_patients_gene <- score %>%
  group_by(patientID) %>%
  summarise(has_pre = "pre" %in% state, has_post = "post" %in% state, .groups = "drop") %>%
  filter(has_pre & has_post) %>%
  pull(patientID)

valid_patients <- intersect(paired_patients_gene, response_data$patientID)
cat("基因层面分析：筛选出", length(valid_patients), "位符合条件的配对患者\n")

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

# 计算 delta = post - pre（连续变化量）
df_delta_gene <- df_pre_gene %>% select(patientID)
for (g in genes) {
  df_delta_gene[[g]] <- df_post_gene[[g]] - df_pre_gene[[g]]
}
# ==========================================
# 【核心】基因 delta ~ 重塑指数 相关性（连续不分组）
# ==========================================
# 只保留有重塑指数的患者
df_delta_gene <- df_delta_gene %>%
  filter(patientID %in% patient_remodeling$patientID)

df_delta_gene[1:5,1:5]

remodel_vec <- patient_remodeling %>%
  filter(patientID %in% df_delta_gene$patientID) %>%
  arrange(patientID) %>%
  pull(remodeling_index)

# 批量计算相关性
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

# 筛选显著基因
gene_sig <- gene_remodel_corr %>%
  filter(pval < 0.01) %>%
  arrange(desc(abs(cor)))

cat("显著关联重塑指数的基因数量：", nrow(gene_sig), "\n")

# ==========================================
# 可视化：显著基因热图（你指定的配色）
# ==========================================
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

pdf("/home/zwj/raid1/MSD/figure/BEYOND/Treg_Myeloid_Gene_Delta_vs_Remodeling_Heatmap.pdf",
    width = 4, height = 20)
draw(ht_gene)
dev.off()

cat("✅ 基因 delta ~ 重塑指数 热图已保存！\n")

# ==========================================
# 单个基因 logFC ~ remodeling index 散点图
# ==========================================
plot_single_gene_remodel <- function(target_gene){
  
  # 1. 准备数据
  plot_df <- data.frame(
    patientID = df_delta_gene$patientID,
    logFC = df_delta_gene[[target_gene]],
    remodeling_index = remodel_vec
  )
  
  # 2. 相关性检验
  ct <- cor.test(plot_df$logFC, plot_df$remodeling_index, method = "spearman")
  rho <- round(ct$estimate, 3)
  pval <- signif(ct$p.value, 3)
  
  # 3. 绘图
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
  
  # 显示图片
  print(p)
  
  # 保存
  ggsave(
    paste0("/home/zwj/raid1/MSD/figure/BEYOND/Gene_", target_gene, "_logFC_vs_Remodeling.pdf"),
    plot = p, width = 5, height = 4
  )
  
  cat("✅ 已保存：", target_gene, "散点图\n")
}

# ======================
# 在这里替换你想看的基因
# ======================
plot_single_gene_remodel("CCL18")

# 你可以继续画任意基因，例如：
plot_single_gene_remodel("LAG3")
# plot_single_gene_remodel("CHI3L1")

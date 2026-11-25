suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(mgcv)
})

# ---------- 前半部分：数据读取与对齐（与你之前相同） ----------
Lin1 <- read.csv("/home/zwj/raid1/MSD/data/h5ad2/Lin1.csv") |> as.data.frame()
#Lin1 <- Lin1[Lin1$pseudotime > 0.25,]
score <- read.csv("/home/zwj/raid1/MSD/data/h5ad2/all_pseudobulk.csv") |> as.data.frame()
rownames(score) <- score$X
score <- score[, -1, drop = FALSE]
score[1:10,1:10]

ig_genes <- grep("^IG", rownames(score), value = TRUE)
score <- score[!rownames(score) %in% ig_genes,]

expr_mat <- t(score) |> as.data.frame(check.names = FALSE)

samples <- intersect(Lin1$sampleID, rownames(expr_mat))
stopifnot(length(samples) > 0)
Lin1 <- Lin1[match(samples, Lin1$sampleID), , drop = FALSE]
rownames(Lin1) <- Lin1$sampleID
expr_mat <- expr_mat[samples, , drop = FALSE]
stopifnot(identical(rownames(expr_mat), rownames(Lin1)))

drop_ids <- c("Pa164")
keep_rows <- !rownames(expr_mat) %in% drop_ids
expr_mat <- expr_mat[keep_rows, , drop = FALSE]
Lin1 <- Lin1[keep_rows, , drop = FALSE]

# --------- 先构造目标细胞中的表达强度/覆盖度指标 ---------
# expr_mat 行=样本, 列=基因；Lin1 里应有目标细胞的比例列
target_cluster <- "CD8T_Temra_CX3CR1"

Lin1$CD8T_Tem_GZMK <- Lin1$CD8T_Tem_NR4A1high + Lin1$CD8T_Tem_NR4A1low

# 1) 基因“在目标细胞中表达”的证据：用该细胞比例作为权重做加权均值与加权非零率
w <- as.numeric(Lin1[[target_cluster]])
w <- ifelse(is.na(w), 0, w)
w <- w / (sum(w) + 1e-12)

# 如果你的 score 来自多细胞伪bulk，总表达可能包含其他细胞的贡献。
# 可用加权相关偏向于目标细胞占比高的样本。
weighted_mean <- colSums(sweep(expr_mat, 1, w, `*`), na.rm = TRUE)

# 加权非零率：在目标细胞占比高的样本中该基因为非零的权重和
nz <- expr_mat > 0
weighted_nz_rate <- colSums(sweep(nz, 1, w, `*`), na.rm = TRUE)

# 2) 过滤标准（可调）
min_w_mean <- quantile(weighted_mean, 0.5, na.rm = TRUE)   # 至少高于加权均值的中位数
min_w_nz   <- 0.1                                          # 在权重意义上 ≥10% 的“非零覆盖”
genes_keep <- names(which((weighted_mean >= min_w_mean) & (weighted_nz_rate >= min_w_nz)))

message("基因库：原始 ", ncol(expr_mat), " 个；通过目标细胞表达过滤后剩余 ", length(genes_keep), " 个。")

# --------- 预处理：对极端值做 winsorize，并在样本维度做简单标准化 ---------
winsorize <- function(v, p = 0.01) {
  q <- quantile(v, probs = c(p, 1-p), na.rm = TRUE)
  v[v < q[1]] <- q[1]; v[v > q[2]] <- q[2]; v
}
expr_clean <- expr_mat


# --------- 拟合与相似度（加入权重） ----------
x <- as.numeric(Lin1$pseudotime)
grid_x <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 200)

predict_smooth <- function(x, y, newx, span = 0.75) {
  fit <- try(loess(y ~ x, span = span, degree = 2,
                   control = loess.control(surface = "direct")), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, length(newx)))
  as.numeric(predict(fit, newdata = data.frame(x = newx)))
}
z <- function(v) if (sd(v, na.rm = TRUE) > 0) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE) else v
cosine_sim <- function(a, b) {
  if (anyNA(a) || anyNA(b)) return(NA_real_)
  den <- sqrt(sum(a^2)) * sqrt(sum(b^2)); if (den == 0) return(NA_real_)
  sum(a*b)/den
}

# 目标亚群曲线（对 y 也可做 winsorize，避免极端点）
y_target <- winsorize(as.numeric(Lin1[[target_cluster]]), p = 0.01)
target_curve <- predict_smooth(x, y_target, grid_x, span = 0.75)
target_curve_z <- z(target_curve)

# 加权相关函数：在样本时域的拟合前，先将原始点做权重相关也可以，这里对拟合后的网格做相关，无需权重；
# 但为进一步保守，我们在样本层面再算一个“点相关（加权）”，作为附加排序参考。
w_norm <- w / (sum(w) + 1e-12)

weighted_cor <- function(a, b, w) {
  # 计算加权皮尔逊相关（样本点）
  ok <- is.finite(a) & is.finite(b) & is.finite(w)
  a <- a[ok]; b <- b[ok]; w <- w[ok]
  if (length(a) < 3 || sum(w) == 0) return(NA_real_)
  a_mean <- sum(w*a) / sum(w); b_mean <- sum(w*b) / sum(w)
  num <- sum(w*(a - a_mean)*(b - b_mean))
  den <- sqrt(sum(w*(a - a_mean)^2) * sum(w*(b - b_mean)^2))
  if (den == 0) return(NA_real_) else num/den
}

res <- purrr::map_dfr(colnames(expr_clean), function(g) {
  yg <- winsorize(as.numeric(expr_clean[[g]]), p = 0.01)
  curve_g <- predict_smooth(x, yg, grid_x, span = 0.75)
  if (all(is.na(curve_g))) return(NULL)
  curve_g_z <- z(curve_g)
  
  tibble::tibble(
    gene = g,
    pearson_curve = suppressWarnings(cor(target_curve_z, curve_g_z, use = "pairwise.complete.obs")),
    cosine_curve  = cosine_sim(target_curve_z, curve_g_z),
    w_cor_points  = weighted_cor(y_target, yg, w_norm),   # 样本点的加权相关
    w_mean_expr   = weighted_mean[g],
    w_nz_rate     = weighted_nz_rate[g]
  )
}) |>
  dplyr::filter(is.finite(pearson_curve), is.finite(cosine_curve)) |>
  dplyr::arrange(dplyr::desc(pearson_curve))

# 进一步保守：要求样本点的加权相关也为正且显著
res <- res |>
  dplyr::filter(w_cor_points > 0)

# 查看Top
head(res, 40)


res_plot <- res |>
  dplyr::filter(w_mean_expr >= quantile(w_mean_expr, 0.3, na.rm = TRUE)) |>
  dplyr::slice_max(order_by = pearson_curve, n = 30)

res_plot


#res_plot <- res[res$gene %in% c("ERO1A","IFNGR2","IFNGR1","HIF1A","ID2","ID3"),]

plot_top <- function(gene_vec) {
  df_plot <- purrr::map_dfr(c("TARGET" = target_cluster, gene_vec), function(name) {
    yv <- if (name == target_cluster) y_target else winsorize(as.numeric(expr_clean[[name]]), 0.01)
    yc <- predict_smooth(x, yv, grid_x)
    tibble::tibble(gene = name, pseudotime = grid_x, value = z(yc))
  }, .id = "src")
  df_plot$gene <- factor(df_plot$gene, levels = c(target_cluster, gene_vec))
  ggplot2::ggplot(df_plot, ggplot2::aes(pseudotime, value, color = gene, size = gene == target_cluster)) +
    ggplot2::geom_line() +
    ggplot2::scale_size_manual(values = c(`TRUE` = 1.2, `FALSE` = 0.5), guide = "none") +
    ggplot2::theme_bw() + ggplot2::labs(y = "Fitted (z)", title = paste0("Similarity to ", target_cluster))
}

print(plot_top(res_plot$gene))
print(res_plot$gene)



suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(viridis)     # 更好的配色
  library(gridExtra)   # 拼图用
  library(tibble)
})

# ==============================================================================
# 准备绘图数据 (基于你已有的 res_plot 和 predict_smooth 函数)
# ==============================================================================

# 1. 提取 Top 基因列表
top_genes <- res_plot$gene
# 确保目标细胞也在列表里，方便对比
all_features <- unique(c(target_cluster, top_genes))

# 2. 构建矩阵数据 (Long Format)
# 计算所有基因在 grid_x 上的平滑曲线
plot_data <- purrr::map_dfr(all_features, function(feat) {
  # 获取原始数据
  if (feat == target_cluster) {
    raw_val <- y_target # 使用之前计算好的 y_target
  } else {
    raw_val <- winsorize(as.numeric(expr_clean[[feat]]), 0.01)
  }
  
  # 平滑拟合
  smooth_val <- predict_smooth(x, raw_val, grid_x)
  
  # Z-score 标准化 (为了在热图中可比)
  z_val <- z(smooth_val)
  
  # 归一化到 0-1 (为了在影子图中比较形态)
  min_v <- min(smooth_val, na.rm = TRUE)
  max_v <- max(smooth_val, na.rm = TRUE)
  norm_val <- (smooth_val - min_v) / (max_v - min_v)
  
  tibble::tibble(
    feature = feat,
    pseudotime = grid_x,
    z_score = z_val,
    norm_value = norm_val,
    is_target = (feat == target_cluster)
  )
})

# ==============================================================================
# 方案一：伪时间热图 (Pseudotime Heatmap)
# 优点：一目了然地看到基因随时间的表达波峰，不拥挤
# ==============================================================================

# 对基因进行排序：按照最大表达峰出现的时间点排序
# 1. 确保 plot_data 已经生成 (如果之前的步骤没报错的话)
# 如果 plot_data 也是空的，说明上一步就有问题，可以用 head(plot_data) 检查一下

# 2. 修正后的排序代码
gene_order <- plot_data %>%
  dplyr::filter(!is_target) %>%
  dplyr::group_by(feature) %>%
  # 关键修改：显式调用 dplyr::summarise 并强制解除分组
  dplyr::summarise(peak_time = pseudotime[which.max(z_score)], .groups = "drop") %>%
  dplyr::arrange(peak_time) %>%
  # 关键修改：使用引号 "feature" 避免对象查找错误
  dplyr::pull("feature")

# 检查一下是否成功
print(head(gene_order))

# 3. 重新运行绘图部分 (确保 gene_order 存在后)
plot_data_heatmap <- plot_data %>%
  mutate(feature = factor(feature, levels = c(target_cluster, rev(gene_order))))

p1_heatmap <- ggplot(plot_data_heatmap, aes(x = pseudotime, y = feature, fill = z_score)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", name = "Z-score") +
  scale_x_continuous(expand = c(0, 0), name = "Pseudotime") +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Gene Expression Dynamics along Pseudotime",
       subtitle = paste("Top genes correlated with", target_cluster))

print(p1_heatmap)

# ==============================================================================
# 方案二：分面“影子”图 (Small Multiples with Reference Shadow)
# 优点：可以清晰地看到每个基因与 Target 的具体差异（是早一点、晚一点还是完全同步）
# ==============================================================================

# 提取目标细胞的曲线作为背景
target_curve_data <- plot_data %>% 
  filter(is_target) %>% 
  select(pseudotime, target_val = norm_value)

# 提取基因数据
gene_curve_data <- plot_data %>% 
  filter(!is_target) %>%
  left_join(target_curve_data, by = "pseudotime")

p2_shadow <- ggplot(gene_curve_data, aes(x = pseudotime)) +
  # 1. 画影子 (目标细胞趋势)，用灰色虚线
  geom_line(aes(y = target_val), color = "grey70", linetype = "dashed", size = 0.8) +
  # 2. 画基因趋势，用彩色实线
  geom_line(aes(y = norm_value, color = feature), size = 1) +
  # 分面展示
  facet_wrap(~feature, ncol = 5) + 
  scale_color_viridis_d(option = "turbo", guide = "none") + # 不显示图例，颜色仅用于区分
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 9),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Gene Trends (Color) vs Target Treg (Grey Dashed)",
       y = "Normalized Expression (0-1)", x = "Pseudotime")
p2_shadow
# ==============================================================================
# 方案三：棒棒糖图 (Lollipop Chart) - 展示相关性排序与表达量
# 优点：清晰展示谁的相关性最强，以及这个基因是否高表达
# ==============================================================================

p3_rank <- ggplot(res_plot, aes(x = reorder(gene, pearson_curve), y = pearson_curve)) +
  geom_segment(aes(xend = gene, yend = 0), color = "grey") +
  geom_point(aes(size = w_mean_expr, color = pearson_curve)) +
  coord_flip() +
  scale_color_viridis_c(option = "plasma", name = "Correlation") +
  scale_size_continuous(name = "Weighted\nMean Expr") +
  theme_bw() +
  labs(x = "Gene", y = "Pearson Correlation with Target",
       title = "Top Correlated Genes Ranking")
p3_rank
# ==============================================================================
# 输出图像
# ==============================================================================

# 打印热图
print(p1_heatmap)

# 打印分面图 (建议保存为大尺寸 PDF)
print(p2_shadow)

# 打印排序图
print(p3_rank)

# 如果你想把它们拼在一起 (可选)
# grid.arrange(p1_heatmap, p3_rank, ncol = 2, widths = c(1.5, 1))
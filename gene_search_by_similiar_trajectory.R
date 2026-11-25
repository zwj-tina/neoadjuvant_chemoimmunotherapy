suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(mgcv)
})

Lin1 <- read.csv("./Lin1.csv") |> as.data.frame()
 #Lin1 <- Lin1[Lin1$pseudotime > 0.25,]
score <- read.csv("./all_pseudobulk.csv") |> as.data.frame()
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

target_cluster <- "CD8T_Temra_CX3CR1"

Lin1$CD8T_Tem_GZMK <- Lin1$CD8T_Tem_NR4A1high + Lin1$CD8T_Tem_NR4A1low

w <- as.numeric(Lin1[[target_cluster]])
w <- ifelse(is.na(w), 0, w)
w <- w / (sum(w) + 1e-12)

weighted_mean <- colSums(sweep(expr_mat, 1, w, `*`), na.rm = TRUE)

nz <- expr_mat > 0
weighted_nz_rate <- colSums(sweep(nz, 1, w, `*`), na.rm = TRUE)

min_w_mean <- quantile(weighted_mean, 0.5, na.rm = TRUE)   
min_w_nz   <- 0.1                                          
genes_keep <- names(which((weighted_mean >= min_w_mean) & (weighted_nz_rate >= min_w_nz)))

winsorize <- function(v, p = 0.01) {
  q <- quantile(v, probs = c(p, 1-p), na.rm = TRUE)
  v[v < q[1]] <- q[1]; v[v > q[2]] <- q[2]; v
}
expr_clean <- expr_mat

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

y_target <- winsorize(as.numeric(Lin1[[target_cluster]]), p = 0.01)
target_curve <- predict_smooth(x, y_target, grid_x, span = 0.75)
target_curve_z <- z(target_curve)

w_norm <- w / (sum(w) + 1e-12)

weighted_cor <- function(a, b, w) {

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
    w_cor_points  = weighted_cor(y_target, yg, w_norm),   
    w_mean_expr   = weighted_mean[g],
    w_nz_rate     = weighted_nz_rate[g]
  )
}) |>
  dplyr::filter(is.finite(pearson_curve), is.finite(cosine_curve)) |>
  dplyr::arrange(dplyr::desc(pearson_curve))

res <- res |>
  dplyr::filter(w_cor_points > 0)

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
  library(viridis)
  library(gridExtra)
  library(tibble)
})

top_genes <- res_plot$gene

all_features <- unique(c(target_cluster, top_genes))

plot_data <- purrr::map_dfr(all_features, function(feat) {

  if (feat == target_cluster) {
    raw_val <- y_target 
  } else {
    raw_val <- winsorize(as.numeric(expr_clean[[feat]]), 0.01)
  }
  
  smooth_val <- predict_smooth(x, raw_val, grid_x)
  
  z_val <- z(smooth_val)
  
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

gene_order <- plot_data %>%
  dplyr::filter(!is_target) %>%
  dplyr::group_by(feature) %>%

  dplyr::summarise(peak_time = pseudotime[which.max(z_score)], .groups = "drop") %>%
  dplyr::arrange(peak_time) %>%

  dplyr::pull("feature")

print(head(gene_order))

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

target_curve_data <- plot_data %>% 
  filter(is_target) %>% 
  select(pseudotime, target_val = norm_value)

gene_curve_data <- plot_data %>% 
  filter(!is_target) %>%
  left_join(target_curve_data, by = "pseudotime")

p2_shadow <- ggplot(gene_curve_data, aes(x = pseudotime)) +

  geom_line(aes(y = target_val), color = "grey70", linetype = "dashed", size = 0.8) +

  geom_line(aes(y = norm_value, color = feature), size = 1) +

  facet_wrap(~feature, ncol = 5) + 
  scale_color_viridis_d(option = "turbo", guide = "none") + 
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

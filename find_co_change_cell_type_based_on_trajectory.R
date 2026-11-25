suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(mgcv)
  library(tidyr)
  library(tibble)
})

Lin1 <- read.csv("/home/zwj/raid1/MSD/data/h5ad2/Lin1.csv") |> as.data.frame()

T_cell <- read.csv("/home/zwj/raid1/MSD/data/h5ad2/T_harmony.csv") |>
  dplyr::select(sampleID, sub_cell_type2)
T_df <- table(T_cell$sampleID, T_cell$sub_cell_type2) |>
  as.data.frame() |>
  dplyr::mutate(Freq = Freq / ave(Freq, Var1, FUN = sum)) |>
  dplyr::rename(sampleID = Var1, cell.type = Var2)

B_cell <- read.csv("/home/zwj/raid1/MSD/data/h5ad2/B_harmony.csv") |>
  dplyr::select(sampleID, sub_cell_type)
B_df <- table(B_cell$sampleID, B_cell$sub_cell_type) |>
  as.data.frame() |>
  dplyr::mutate(Freq = Freq / ave(Freq, Var1, FUN = sum)) |>
  dplyr::rename(sampleID = Var1, cell.type = Var2)

Myeloid_cell <- read.csv("/home/zwj/raid1/MSD/data/h5ad2/Myeloid_harmony.csv") |>
  dplyr::select(sampleID, sub_cell_type)
Myeloid_df <- table(Myeloid_cell$sampleID, Myeloid_cell$sub_cell_type) |>
  as.data.frame() |>
  dplyr::mutate(Freq = Freq / ave(Freq, Var1, FUN = sum)) |>
  dplyr::rename(sampleID = Var1, cell.type = Var2)

Endo_cell <- read.csv("/home/zwj/raid1/MSD/data/h5ad2/Endothelial_harmony.csv") |>
  dplyr::select(sampleID, sub_cell_type)
Endo_df <- table(Endo_cell$sampleID, Endo_cell$sub_cell_type) |>
  as.data.frame() |>
  dplyr::mutate(Freq = Freq / ave(Freq, Var1, FUN = sum)) |>
  dplyr::rename(sampleID = Var1, cell.type = Var2)

Fibro_cell <- read.csv("/home/zwj/raid1/MSD/data/h5ad2/Fibroblast_harmony.csv") |>
  dplyr::select(sampleID, sub_cell_type2)
Fibro_df <- table(Fibro_cell$sampleID, Fibro_cell$sub_cell_type2) |>
  as.data.frame() |>
  dplyr::mutate(Freq = Freq / ave(Freq, Var1, FUN = sum)) |>
  dplyr::rename(sampleID = Var1, cell.type = Var2)

all_df <- dplyr::bind_rows(T_df, B_df, Myeloid_df, Endo_df, Fibro_df)

ratio <- all_df |>
  tidyr::pivot_wider(names_from = cell.type, values_from = Freq, values_fill = 0) |>
  tibble::column_to_rownames("sampleID")

 # S
ratio <- ratio[!grepl("S", rownames(ratio)), , drop = FALSE]
ratio

ratio$CD8T_Tem_GZMK <- ratio$CD8T_Tem_NR4A1high + ratio$CD8T_Tem_NR4A1low
ratio <- ratio[, !names(ratio) %in% c("CD8T_Tem_NR4A1high", "CD8T_Tem_NR4A1low")]

 # Lin1
expr_mat <- ratio
samples <- intersect(Lin1$sampleID, rownames(expr_mat))
stopifnot(length(samples) > 0)
Lin1 <- Lin1[match(samples, Lin1$sampleID), , drop = FALSE]
rownames(Lin1) <- Lin1$sampleID
expr_mat <- expr_mat[samples, , drop = FALSE]
stopifnot(identical(rownames(expr_mat), rownames(Lin1)))

head(Lin1)

Lin1$CD8T_Tem_GZMK <- Lin1$CD8T_Tem_NR4A1high + Lin1$CD8T_Tem_NR4A1low
drop_ids <- c("Pa164")
keep_rows <- !rownames(expr_mat) %in% drop_ids
expr_mat <- expr_mat[keep_rows, , drop = FALSE]
Lin1 <- Lin1[keep_rows, , drop = FALSE]

# --------- Target cell-type metric column (from Lin1) ---------
# e.g. Lin1 stores CD8Tex/CD4Treg proportions or a combined metric
target_cluster <- "CD4T_Treg_TNFRSF9high"

winsorize <- function(v, p = 0.01) {
  q <- stats::quantile(v, probs = c(p, 1 - p), na.rm = TRUE)
  v[v < q[1]] <- q[1]; v[v > q[2]] <- q[2]; v
}
z <- function(v) if (stats::sd(v, na.rm = TRUE) > 0) (v - mean(v, na.rm = TRUE)) / stats::sd(v, na.rm = TRUE) else v
cosine_sim <- function(a, b) {
  if (anyNA(a) || anyNA(b)) return(NA_real_)
  den <- sqrt(sum(a^2)) * sqrt(sum(b^2)); if (den == 0) return(NA_real_)
  sum(a * b) / den
}
predict_smooth <- function(x, y, newx, span = 0.75) {
  fit <- try(stats::loess(y ~ x, span = span, degree = 2,
                          control = stats::loess.control(surface = "direct")),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, length(newx)))
  as.numeric(stats::predict(fit, newdata = data.frame(x = newx)))
}
weighted_cor <- function(a, b, w) {
  ok <- is.finite(a) & is.finite(b) & is.finite(w)
  a <- a[ok]; b <- b[ok]; w <- w[ok]
  if (length(a) < 3 || sum(w) == 0) return(NA_real_)
  a_mean <- sum(w * a) / sum(w); b_mean <- sum(w * b) / sum(w)
  num <- sum(w * (a - a_mean) * (b - b_mean))
  den <- sqrt(sum(w * (a - a_mean)^2) * sum(w * (b - b_mean)^2))
  if (den == 0) return(NA_real_) else num / den
}

x <- as.numeric(Lin1$pseudotime)
grid_x <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 200)

w <- as.numeric(Lin1[[target_cluster]])
w <- ifelse(is.na(w), 0, w)
w <- w / (sum(w) + 1e-12)
w_norm <- w

y_target <- winsorize(as.numeric(Lin1[[target_cluster]]), p = 0.01)
target_curve <- predict_smooth(x, y_target, grid_x, span = 0.75)
target_curve_z <- z(target_curve)

 # --------- ( P ) ----------
res <- purrr::map_dfr(colnames(expr_mat), function(ct) {
  y_ct <- winsorize(as.numeric(expr_mat[[ct]]), p = 0.01)
  curve_ct <- predict_smooth(x, y_ct, grid_x, span = 0.75)
  
  if (all(is.na(curve_ct))) return(NULL)
  
  curve_ct_z <- z(curve_ct)
  
 # cor.test p
  test_res <- suppressWarnings(
    stats::cor.test(target_curve_z, curve_ct_z, use = "pairwise.complete.obs", method = "pearson")
  )
  
  tibble::tibble(
    cell_type = ct,
    pearson_curve = as.numeric(test_res$estimate),
    p_value = as.numeric(test_res$p.value), # P
    cosine_curve  = cosine_sim(target_curve_z, curve_ct_z),
    w_cor_points  = weighted_cor(y_target, y_ct, w_norm)
  )
}) |>
  dplyr::filter(is.finite(pearson_curve), is.finite(cosine_curve)) |>
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    sig_label = ifelse(p_adj < 0.001, "***", 
                       ifelse(p_adj < 0.01, "**", 
                              ifelse(p_adj < 0.05, "*", "ns")))
  ) |>
  dplyr::arrange(dplyr::desc(pearson_curve))

 # P
head(res, 30)

 # res <- res |> dplyr::filter(w_cor_points > 0)

write.csv(res, "/home/zwj/raid1/MSD/data/h5ad2/CD4T_Treg_corr.csv")

 # --------- Top ----------
res_plot <- res |>
  dplyr::slice_max(order_by = pearson_curve, n = 10)
head(res_plot)

plot_top <- function(ct_vec) {
  ct_vec <- unique(ct_vec[ct_vec != target_cluster])
  
  df_plot <- purrr::map_dfr(c("TARGET" = target_cluster, ct_vec), function(name) {
    yv <- if (name == target_cluster) y_target else winsorize(as.numeric(expr_mat[[name]]), 0.01)
    yc <- predict_smooth(x, yv, grid_x)
    tibble::tibble(feature = name, pseudotime = grid_x, value = z(yc))
  }, .id = "src")
  df_plot$feature <- factor(df_plot$feature, levels = c(target_cluster, ct_vec))
  ggplot2::ggplot(df_plot, ggplot2::aes(pseudotime, value, color = feature, size = feature == target_cluster)) +
    ggplot2::geom_line() +
    ggplot2::scale_size_manual(values = c(`TRUE` = 1.2, `FALSE` = 0.5), guide = "none") +
    ggplot2::theme_bw() +
    ggplot2::labs(y = "Fitted (z)", title = paste0("Similarity to ", target_cluster)) 
}

print(plot_top(res_plot$cell_type))
print(res_plot$cell_type)

ggsave("/home/zwj/raid1/MSD/figure/BEYOND/Lin2_CD4.pdf",width = 6, height = 4.5)


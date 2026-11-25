library(FNN)
library(ggrepel)
library(dplyr)
library(stringr)
library(tidyverse)
library(readxl)
library(ggsci) # ggsci scale_color_npg()

# 1.
T_cell <- read.csv("MSD/data/h5ad2/T_harmony.csv")

info <- as.data.frame(read_excel("MSD/data/info/clinical_info.xlsx"))
info <- info[info$pathological_response %in% c("MPR","non-MPR","PD"),]

T_cell <- T_cell[,c("sampleID","sub_cell_type2")]
T_df <- table(T_cell$sampleID, T_cell$sub_cell_type2)
T_df <- as.data.frame(prop.table(T_df, margin = 1)) # prop.table
colnames(T_df) <- c("sampleID","cell.type","Freq")

# "S$" S S
T_df$patientID <- sub("S$", "", T_df$sampleID)

# PHATE/UMAP
phate_ref <- read.csv("MSD/data/h5ad2/PATHE.csv")
phate_ref$patientID <- phate_ref$sampleID
rownames(phate_ref) <- phate_ref$sampleID
phate_ref <- phate_ref[,-1]

T_df <- T_df[T_df$patientID %in% phate_ref$patientID,]

# pre post
all_df <- T_df %>%
  mutate(
    state     = if_else(str_detect(sampleID, "S$"), "post", "pre"),
    patientID = sub("S$", "", sampleID)
  )

# pre_df (Reference)
pre_df <- all_df[all_df$state %in% c("pre"),]
ratio_ref <- pre_df %>%
  select(sampleID, cell.type, Freq) %>%
  tidyr::pivot_wider(names_from = cell.type, values_from = Freq, values_fill = 0) %>%
  column_to_rownames("sampleID") %>%
  as.data.frame()

# ratio_ref phate_ref phate_ref
intersect_samples <- intersect(rownames(ratio_ref), rownames(phate_ref))
ratio_ref <- ratio_ref[intersect_samples, ]
phate_ref_aligned <- phate_ref[intersect_samples, ]
head(phate_ref_aligned)

# post_df (Query)
post_df <- all_df[all_df$state %in% c("post"),]
ratio_new <- post_df %>%
  select(sampleID, cell.type, Freq) %>%
  tidyr::pivot_wider(names_from = cell.type, values_from = Freq, values_fill = 0) %>%
  column_to_rownames("sampleID") %>%
  as.data.frame()

common_cols <- intersect(colnames(ratio_ref), colnames(ratio_new))
ratio_ref <- ratio_ref[, common_cols, drop = FALSE]
ratio_new <- ratio_new[, common_cols, drop = FALSE]

# 2.
knn_res <- get.knnx(data = ratio_ref, query = ratio_new, k = 3)

# 3. as.matrix
# matrixcluster
phate_coords <- as.matrix(phate_ref_aligned[, c("PHATE1", "PHATE2")])
head(phate_coords)

phate_clusters <- phate_ref_aligned$clusters

# 4. KNN
# 1 ncol = 3 ncol = 2 phate_coords UMAP1 UMAP2
mapped_phate <- matrix(NA, nrow = nrow(ratio_new), ncol = 2)
predicted_clusters <- rep(NA, nrow(ratio_new))

for (i in 1:nrow(ratio_new)) {
  neighbor_idx <- knn_res$nn.index[i, ]
  
# 1. ( 2)
  mapped_phate[i, ] <- colMeans(phate_coords[neighbor_idx, , drop = FALSE])
  
# 2. cluster cluster mode
  cluster_vals <- phate_clusters[neighbor_idx]
  predicted_clusters[i] <- names(sort(table(cluster_vals), decreasing = TRUE))[1]
}

# 2 PATHE -> PHATE
mapped_df <- data.frame(
  PHATE1 = mapped_phate[, 1],
  PHATE2 = mapped_phate[, 2],
  predicted_cluster = as.factor(predicted_clusters),
  row.names = rownames(ratio_new)
)

mapped_df$sampleID <- rownames(mapped_df)
mapped_df$patientID <- sub("S$", "", mapped_df$sampleID)

# 5.
info_csv <- as.data.frame(read.csv("MSD/data/info/clinical_info.csv"))
info_csv <- info_csv[info_csv$pathological_response %in% c("MPR","non-MPR","PD"),]

mapped_df <- merge(mapped_df, info_csv, by = "patientID")

write.csv(mapped_df, "MSD/data/h5ad2/LUSC_post_PHATE.csv", row.names = FALSE)
table(mapped_df$predicted_cluster, mapped_df$pathological_response)

# 6. (Pre + Post )
# reference cluster factor
phate_ref_aligned$clusters <- as.factor(phate_ref_aligned$clusters)

p_combined <- ggplot() +
# Pre (Reference)
# alpha = 0.3 shape = 16
  geom_point(data = phate_ref_aligned, 
             aes(x = PHATE1, y = PHATE2, color = clusters), 
             size = 2, alpha = 0.3, shape = 16) +
  
# Post (Mapped)
# shape = 21 () fill color="black"
  geom_point(data = mapped_df, 
             aes(x = PHATE1, y = PHATE2, fill = predicted_cluster), 
             size = 3, shape = 21, color = "black", stroke = 0.6) +
  
  theme_classic(base_size = 12) +
  scale_color_npg(name = "Pre (Reference)") + 
  scale_fill_npg(name = "Post (Mapped)") +
  
  labs(
    title = "PHATE Projection: Pre vs Post",
    subtitle = "Faded points: Pre (Reference) | Solid points with borders: Post (Mapped)",
    x = "PHATE 1", 
    y = "PHATE 2"
  )

print(p_combined)

ggsave("MSD/figure/BEYOND/LUSC_combined_PHATE.pdf", plot = p_combined, width = 7, height = 5)

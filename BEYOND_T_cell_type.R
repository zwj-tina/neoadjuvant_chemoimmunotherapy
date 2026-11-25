 # 0)
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggrepel)
  library(ggsci)
  library(ComplexHeatmap)
  library(reticulate)
  library(slingshot)
  library(mgcv)
  library(ggalluvial)
  library(limma)
})

 # embedding.similarity / palantir
source("/home/zwj/raid1/MSD/BEYOND/BEYOND_DLPFC-main/2. Cell-type analysis/utils/load.code.env.R")
source("/home/zwj/raid1/MSD/BEYOND/BEYOND_DLPFC-main/4. BEYOND/utils/utils.R")

 # Python
sc      <- reticulate::import("scanpy")
anndata <- reticulate::import("anndata")
pd      <- reticulate::import("pandas")
np      <- reticulate::import("numpy")
 # reticulate::py_install("phate", pip = TRUE)

set.seed(1); np$random$seed(1L)

# User option: keep only the pre state; set to FALSE for a pre-to-post alluvial plot
keep_only_state <- NULL # set to "pre" or "post"; NULL keeps all states

# 1) Convert sample x subcluster counts into sample x subcluster proportions
T_cell <- read.csv("/home/zwj/raid1/MSD/data/h5ad2/T_harmony.csv")
number <- table(T_cell$sampleID)
T_cell <- T_cell[T_cell$sampleID %in% names(which(number > 200)),]

T_cell <- T_cell[,c("sampleID","sub_cell_type2")]
T_df <- table(T_cell$sampleID,T_cell$sub_cell_type2)
T_df <- as.data.frame(T_df / rowSums(T_df))
colnames(T_df) <- c("sampleID","cell.type","Freq")
head(T_df)

all_df <- T_df
head(all_df)

 # T
all_df <- all_df %>%
  mutate(
    state     = if_else(str_detect(sampleID, "S$"), "post", "pre"),
    patientID = sub("S$", "", sampleID)
  )

all_df <- all_df[all_df$state %in% c("pre"),]
length(unique(all_df$patientID))

if (!is.null(keep_only_state)) {
  all_df <- dplyr::filter(all_df, state == keep_only_state)
}

 # 2) ratio==
ratio <- all_df %>%
  select(sampleID, cell.type, Freq) %>%
  tidyr::pivot_wider(names_from = cell.type, values_from = Freq, values_fill = 0) %>%
  column_to_rownames("sampleID") %>%
  as.data.frame()
ratio[1:5,1:5]

 # 3)
info <- read.csv("/home/zwj/raid1/MSD/data/info/clinical_info.csv", check.names = FALSE) %>%
  mutate(across(everything(), ~ replace(., is.na(.), "unknown"))) %>%
  filter(pathological_response != "unknown")
info <- info[info$histology %in% c("LUSC"),]

meta_smpl <- tibble(sampleID = rownames(ratio)) %>%
  mutate(
    state     = if_else(str_detect(sampleID, "S$"), "post", "pre"),
    patientID = sub("S$", "", sampleID)
  ) %>%
  left_join(info, by = "patientID")

keep_ids <- meta_smpl %>% filter(!is.na(pathological_response), pathological_response %in% c("MPR","non-MPR","PD")) %>% pull(sampleID)
ratio    <- ratio[keep_ids, , drop = FALSE]
meta_smpl <- meta_smpl %>% filter(sampleID %in% keep_ids) %>% column_to_rownames("sampleID")

stopifnot(identical(rownames(ratio), rownames(meta_smpl)))
dim(ratio)

 # 4) AnnData + / + UMAP/TSNE/PHATE
 # AnnData
obs_pd = pd$DataFrame(meta_smpl, index=rownames(meta_smpl))
ad = anndata$AnnData(X = as.matrix(ratio), obs = obs_pd)
sc$pp$pca(ad, n_comps = as.integer(10))
sc$pp$neighbors(ad, n_neighbors = as.integer(5), metric="cosine")

sc$tl$leiden(ad, resolution=0.6)
ad$obs$clusters <- ad$obs$leiden

sc$tl$tsne(ad, use_rep = "X", learning_rate = as.integer(100))
sc$tl$umap(ad, maxiter = as.integer(3000), spread = 5)

sc$pl$umap(ad, color = "clusters",size = 70)
 # PHATE
sc$external$tl$phate(ad, n_components = as.integer(4), k = as.integer(3), a = as.integer(4))

umap_2d <- as.data.frame(ad$obsm$X_umap)
colnames(umap_2d)[1:2] <- c("UMAP1","UMAP2")
umap_2d$clusters <- as.vector(ad$obs$clusters)
head(umap_2d)
sample_ids <- rownames(ad)
umap_2d$sampleID <- sample_ids

ggplot(umap_2d, aes(UMAP1, UMAP2, color = clusters)) +
  geom_point(alpha = 1, size = 2) +
  theme_classic(base_size = 12) +
  stat_ellipse(aes(fill = clusters), 
               geom = "polygon", 
               alpha = 0.2, 
               level = 0.95, 
               linetype = 2)+
  scale_color_npg() + scale_fill_npg()

ggsave("/home/zwj/raid1/MSD/figure/pre_LUSC_cluster_UMAP_T_subtypes.pdf", width = 5.5, height = 5)

umap_2d
write.csv(umap_2d,"/home/zwj/raid1/MSD/data/h5ad2/UMAP_cluster.csv")

 # 5) ggplot
df_2d <- as.data.frame(ad$obsm$X_phate)
colnames(df_2d)[1:2] <- c("PHATE1","PHATE2")
 # pandas index
sample_ids <- rownames(ad)
df_2d$sampleID <- sample_ids

obs_df <- as.data.frame(ad$obs)
obs_df$sampleID <- rownames(obs_df)

df_2d <- df_2d %>%
  left_join(obs_df, by = "sampleID") %>%
  mutate(clusters = as.character(clusters)) %>%
  arrange(match(sampleID, rownames(ratio)))

ratio_df <- ratio %>% as.data.frame() %>% rownames_to_column("sampleID")
df_2d <- df_2d %>% left_join(ratio_df, by = "sampleID")

stopifnot(!any(is.na(df_2d$PHATE1)), !any(is.na(df_2d$PHATE2)))

 # 6) PHATE
p1 <- ggplot(df_2d, aes(PHATE1, PHATE2, color = clusters)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(aes(label = sampleID), size = 3, max.overlaps = 20) +
  theme_classic(base_size = 12) +
  labs(title = "2D PHATE with Clusters")
print(p1)
ggsave("/home/zwj/raid1/MSD/figure/BEYOND/pre_LUSC_cluster_PHATE_T_subtypes.pdf", p1, width = 7, height = 6)

p1 <- ggplot(df_2d, aes(PHATE1, PHATE2, color = clusters)) +
  geom_point(size = 3) +
  theme_classic(base_size = 12) + scale_color_npg() + scale_fill_npg()
print(p1)
ggsave("/home/zwj/raid1/MSD/figure/BEYOND/pre_LUSC_cluster_PHATE_T_subtypes.pdf", p1, width = 5.5, height = 5)

p2 <- ggplot(df_2d, aes(PHATE1, PHATE2, color = histology)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(aes(label = sampleID), size = 3, max.overlaps = 20) +
  theme_classic(base_size = 12) +
  scale_color_manual(values = c("LUAD"="#004784","LUSC"="#FF7F0E")) +
  labs(title = "PHATE colored by histology (post only)")
print(p2)
ggsave("/home/zwj/raid1/MSD/figure/BEYOND/pre_all_histology_PHATE_T_subtypes.pdf", p2, width = 6.5, height = 6)

p2 <- ggplot(df_2d, aes(PHATE1, PHATE2, color = histology)) +
  geom_point(size = 3) +
  theme_classic(base_size = 12) +
  scale_color_manual(values = c("LUAD"="#004784","LUSC"="#FF7F0E")) +
  labs(title = "PHATE colored by histology (post only)")
print(p2)
ggsave("/home/zwj/raid1/MSD/figure/BEYOND/pre_all_histology_PHATE_T_subtypes.pdf", p2, width = 6.5, height = 6)

df_resp <- df_2d
df_resp <- df_2d %>% filter(pathological_response %in% c("MPR","PD","non-MPR"))
pathological_response <- c()
for(each in df_resp$pathological_response){
  if(each %in% c("PD","non-MPR")){
    pathological_response <- c(pathological_response,"non-MPR")
  }else{
    pathological_response <- c(pathological_response,"MPR")
  }
}
df_resp$pathological_response <- pathological_response

p3 <- ggplot(df_resp, aes(PHATE1, PHATE2, color = pathological_response)) +
  geom_point(size = 3) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  geom_text_repel(aes(label = sampleID), size = 3, max.overlaps = 20) + 
  scale_color_manual(values = c("MPR"="#3BA0B0","non-MPR"="#D75A5A")) +
  labs(title = "PHATE colored by pathological response")
print(p3)
ggsave("/home/zwj/raid1/MSD/figure/BEYOND/pre_PHATE_by_response_T_subtypes.pdf", p3, width = 6.3, height = 5)

target_feature <- "CD8T_Tem_NR4A1low"
if (target_feature %in% colnames(df_2d)) {
  p_feat <- ggplot(df_2d, aes(PHATE1, PHATE2, color = .data[[target_feature]])) +
    geom_point(alpha = 0.85, size = 2) +
    geom_text_repel(aes(label = sampleID), size = 3, max.overlaps = 20) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black")) +
    scale_color_gradient(low = "#D3D3D3", high = "#B22222") +
    labs(title = paste("PHATE colored by", target_feature))
  print(p_feat)
}
ggsave("/home/zwj/raid1/MSD/figure/BEYOND/pre_CD4T_Treg_TNFRSF9high.pdf", p_feat, width = 7, height = 6)

 # 7) Palantir all
suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(reticulate)
  suppressWarnings({
    if (!requireNamespace("compositions", quietly = TRUE)) message("Installing the compositions package is recommended for CLR transformation")
  })
})

 # ========== 1) Palantir ==========
 # Python
palantir_py <- import("palantir", convert = FALSE)
pandas_py   <- import("pandas",   convert = FALSE)
numpy_py    <- import("numpy",    convert = FALSE)

 # ========== 2) ==========
 # ratio: cell
stopifnot(!anyDuplicated(rownames(ratio)))

X_mat <- as.matrix(ratio)
storage.mode(X_mat) <- "double"

 # Python
X_pd <- pandas_py$DataFrame(
  data = r_to_py(X_mat, convert = TRUE),
  index = r_to_py(rownames(ratio), convert = TRUE),
  columns = r_to_py(colnames(ratio), convert = TRUE),
  dtype = numpy_py$float64
)

 # ========== 3) ==========
root_id <- "Pa59"
message("Palantir early cell (root): ", root_id)

 # ========== 4) Palantir ==========
 # The 'knn' argument is correctly used here in run_diffusion_maps
dm_res <- palantir_py$utils$run_diffusion_maps(X_pd, n_components = as.integer(18),knn = as.integer(30))
ms_res <- palantir_py$utils$determine_multiscale_space(dm_res)

pr_res <- palantir_py$core$run_palantir(
  ms_res,
  early_cell = root_id,
  num_waypoints = as.integer(40),
  knn = as.integer(30) # <--- Added this argument
)
 # ========== 5) pseudotime fate bias ==========
 # R
pseudotime <- pandas_py$DataFrame(pr_res$pseudotime)
pseudotime

fate_probs <- pandas_py$DataFrame(pr_res$branch_probs)

 # data frame
palantir_py_df <- pandas_py$concat(
  list(pseudotime, fate_probs),
  axis = as.integer(1)
)

palantir_df <- py_to_r(palantir_py_df)
palantir_df <- palantir_df %>%
  rownames_to_column("sampleID") 
head(palantir_df)
colnames(palantir_df) <- c("sampleID","pseudotime")
head(palantir_df)

df_2d$pseudotime <- palantir_df$pseudotime

p3 <- ggplot(df_2d, aes(PHATE1, PHATE2, color = pseudotime)) +
  geom_point(size = 3) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_color_gradient(low = "#D3D3D3", high = "#B22222") +
  geom_text_repel(aes(label = sampleID), size = 3, max.overlaps = 20)

print(p3)

df_2d[order(df_2d$pseudotime),c("sampleID","pseudotime")]

write.csv(df_2d,"/home/zwj/raid1/MSD/data/h5ad2/PATHE.csv")

 # 7) Palantir part1
suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(reticulate)
  suppressWarnings({
    if (!requireNamespace("compositions", quietly = TRUE)) message("Installing the compositions package is recommended for CLR transformation")
  })
})

 # ========== 1) Palantir ==========
 # Python
palantir_py <- import("palantir", convert = FALSE)
pandas_py   <- import("pandas",   convert = FALSE)
numpy_py    <- import("numpy",    convert = FALSE)

 # ========== 2) ==========
 # ratio: cell
stopifnot(!anyDuplicated(rownames(ratio)))

Lin1 <- df_2d[df_2d$clusters %in% c("3","1","2"),]
Lin1_ratio <- ratio[Lin1$sampleID,]
Lin1_ratio

X_mat <- as.matrix(Lin1_ratio)
storage.mode(X_mat) <- "double"

 # Python
X_pd <- pandas_py$DataFrame(
  data = r_to_py(X_mat, convert = TRUE),
  index = r_to_py(rownames(Lin1_ratio), convert = TRUE),
  columns = r_to_py(colnames(Lin1_ratio), convert = TRUE),
  dtype = numpy_py$float64
)

 # ========== 3) ==========
root_id <- "Pa59"
message("Palantir early cell (root): ", root_id)

 # ========== 4) Palantir ==========
 # The 'knn' argument is correctly used here in run_diffusion_maps
dm_res <- palantir_py$utils$run_diffusion_maps(X_pd, n_components = as.integer(20),knn = as.integer(10))
ms_res <- palantir_py$utils$determine_multiscale_space(dm_res)

pr_res <- palantir_py$core$run_palantir(
  ms_res,
  early_cell = root_id,
  num_waypoints = as.integer(40),
  knn = as.integer(10) # <--- Added this argument
)
 # ========== 5) pseudotime fate bias ==========
 # R
pseudotime <- pandas_py$DataFrame(pr_res$pseudotime)
pseudotime

fate_probs <- pandas_py$DataFrame(pr_res$branch_probs)

 # data frame
palantir_py_df <- pandas_py$concat(
  list(pseudotime, fate_probs),
  axis = as.integer(1)
)

palantir_df <- py_to_r(palantir_py_df)
palantir_df <- palantir_df %>%
  rownames_to_column("sampleID") 
head(palantir_df)
colnames(palantir_df) <- c("sampleID","pseudotime")
head(palantir_df)

Lin1$pseudotime <- palantir_df$pseudotime
head(Lin1)

p3 <- ggplot(Lin1, aes(PHATE1, PHATE2, color = pseudotime)) +
  geom_point(size = 3) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_color_gradient(low = "#D3D3D3", high = "#B22222") +
  geom_text_repel(aes(label = sampleID), size = 3, max.overlaps = 20)

print(p3)

Lin1[order(Lin1$pseudotime),c("sampleID","pseudotime")]

write.csv(Lin1,"/home/zwj/raid1/MSD/data/h5ad2/Lin1.csv")

 # 7) Palantir part2
suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(reticulate)
  suppressWarnings({
    if (!requireNamespace("compositions", quietly = TRUE)) message("Installing the compositions package is recommended for CLR transformation")
  })
})

 # ========== 1) Palantir ==========
 # Python
palantir_py <- import("palantir", convert = FALSE)
pandas_py   <- import("pandas",   convert = FALSE)
numpy_py    <- import("numpy",    convert = FALSE)

 # ========== 2) ==========
 # ratio: cell
stopifnot(!anyDuplicated(rownames(ratio)))

Lin2 <- df_2d[df_2d$clusters %in% c("3","0"),]
Lin2_ratio <- ratio[Lin2$sampleID,]
Lin2_ratio

X_mat <- as.matrix(Lin2_ratio)
storage.mode(X_mat) <- "double"

 # Python
X_pd <- pandas_py$DataFrame(
  data = r_to_py(X_mat, convert = TRUE),
  index = r_to_py(rownames(Lin2_ratio), convert = TRUE),
  columns = r_to_py(colnames(Lin2_ratio), convert = TRUE),
  dtype = numpy_py$float64
)

 # ========== 3) ==========
root_id <- "Pa59"
message("Palantir early cell (root): ", root_id)

 # ========== 4) Palantir ==========
 # The 'knn' argument is correctly used here in run_diffusion_maps
dm_res <- palantir_py$utils$run_diffusion_maps(X_pd, n_components = as.integer(18),knn = as.integer(13))
ms_res <- palantir_py$utils$determine_multiscale_space(dm_res)

pr_res <- palantir_py$core$run_palantir(
  ms_res,
  early_cell = root_id,
  num_waypoints = as.integer(40),
  knn = as.integer(13) # <--- Added this argument
)
 # ========== 5) pseudotime fate bias ==========
 # R
pseudotime <- pandas_py$DataFrame(pr_res$pseudotime)
pseudotime

fate_probs <- pandas_py$DataFrame(pr_res$branch_probs)

 # data frame
palantir_py_df <- pandas_py$concat(
  list(pseudotime, fate_probs),
  axis = as.integer(1)
)

palantir_df <- py_to_r(palantir_py_df)
palantir_df <- palantir_df %>%
  rownames_to_column("sampleID") 
head(palantir_df)
colnames(palantir_df) <- c("sampleID","pseudotime")
head(palantir_df)

Lin2$pseudotime <- palantir_df$pseudotime
head(Lin2)

p3 <- ggplot(Lin2, aes(PHATE1, PHATE2, color = pseudotime)) +
  geom_point(size = 3) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_color_gradient(low = "#D3D3D3", high = "#B22222") +
  geom_text_repel(aes(label = sampleID), size = 3, max.overlaps = 20)

print(p3)

Lin2[order(Lin2$pseudotime),c("sampleID","pseudotime")]

write.csv(Lin2,"/home/zwj/raid1/MSD/data/h5ad2/Lin2.csv")


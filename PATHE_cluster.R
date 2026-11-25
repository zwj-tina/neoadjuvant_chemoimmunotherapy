library(ggpubr)
library(ggplot2)
library(tidyr) # gather & spread
library(reshape2) # melt & dcast
library(readxl)
library(ggrepel)
library(ggsci)
PATHE <- read.csv("MSD/data/h5ad2/PATHE.csv")
rownames(PATHE) <- PATHE$sampleID
PATHE$clusters <- as.character(PATHE$clusters)
table(PATHE$clusters, PATHE$pathological_response)
head(PATHE)

a <- PATHE[,c("sampleID","clusters")]

ggplot(PATHE, aes(x = CD8T_Tem_NR4A1high, y = CD8T_Tem_NR4A1low)) +
  geom_point(size = 2, color = "#87B9D9") +
  stat_cor(method = "pearson") +
  theme_classic() +
  geom_smooth(method = "lm", se = FALSE, color = "#87B9D9")

PATHE$CD8T_Tem_GZMK <- PATHE$CD8T_Tem_NR4A1high + PATHE$CD8T_Tem_NR4A1low

PATHE$clusters <- factor(PATHE$clusters, c("0","2","1","3"))
library(tidyverse)

compaired <- list(c("0","2"),
                  c("0","1"),
                  c("0","3"))
ggplot(PATHE, 
       aes(x = clusters, y = CD8T_Tem_NR4A1low, color = clusters,
           fill = clusters)) +
  geom_boxplot(outlier.shape = NA,width = 0.4,alpha = 0.5) +
  geom_jitter(size = 0.4, alpha = 1, width = 0.1) +
  theme_pubr() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0)) + 
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test,color = "black") +
  scale_color_manual(values = c("0" = "#E64B35FF",
                                "1" = "#4DBBD5FF",
                                "2" = "#00A087FF",
                                "3" = "#3C5488FF")) + 
  scale_fill_manual(values = c("0" = "#E64B35FF",
                               "1" = "#4DBBD5FF",
                               "2" = "#00A087FF",
                               "3" = "#3C5488FF"))

ggsave("MSD/figure/BEYOND/cluster_CD8Tem.pdf", width = 3.5, height = 3.5)

PATHE$CD8T <- PATHE$CD8T_FOXP3 + PATHE$CD8T_ISG15 + PATHE$CD8T_MAIT_SLC4A10 + PATHE$CD8T_NKlike_KIR3DL2 + 
  PATHE$CD8T_prf_MKI67 + PATHE$CD8T_Tem_NR4A1high + PATHE$CD8T_Tem_NR4A1low + PATHE$CD8T_Temra_CX3CR1 + 
  PATHE$CD8T_Tex_CXCL13 + PATHE$CD8T_Tm_IL7R + PATHE$CD8T_Tm_TCF7 + PATHE$CD8T_Trm_ZNF683
  
PATHE$CD4T <-PATHE$CD4T_prf_MKI67 + PATHE$CD4T_Tem_GZMK + PATHE$CD4T_Tfh_CXCL13 + PATHE$CD4T_Th17_IL17A + 
  PATHE$CD4T_Th1like_CXCL13 + PATHE$CD4T_Tm_ANXA1 + PATHE$CD4T_Tm_FOSB + PATHE$CD4T_Tm_GZMA + 
  PATHE$CD4T_Tm_TCF7 + PATHE$CD4T_Treg_TNFRSF9high + PATHE$CD4T_Treg_TNFRSF9low
  
PATHE$NK <- PATHE$NK_CD16high_CX3CR1 + PATHE$NK_CD16low_GZMK + PATHE$NK_CD16low_GZMK + 
  PATHE$NK_CD16low_IL7R + PATHE$NK_ISG15 + PATHE$NK_prf_MKI67

compaired <- list(c("MPR","non-MPR"),
                  c("MPR","PD"),
                  c("non-MPR","PD"))

ggplot(PATHE, 
       aes(x = pathological_response, y = CD8T, color = pathological_response,
           fill = pathological_response)) +
  geom_boxplot(outlier.shape = NA,width = 0.4,alpha = 0.5) +
  geom_jitter(size = 0.4, alpha = 0.8, width = 0.1) +
  theme_pubr() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0)) + 
  scale_color_manual(values = c('MPR'=' #2C7C75',
                                'non-MPR'=' #B35C1E',
                                'PD'=' #5F5AA2')) +
  scale_fill_manual(values = c('MPR'=' #2C7C75',
                               'non-MPR'=' #B35C1E',
                               'PD'=' #5F5AA2')) +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test,color = "black")

library(ggridges)
library(tidyverse)
PATHE$clusters <- factor(PATHE$clusters,levels = c("3","1","0","2"))
ggplot(PATHE, aes(x = pseudotime, y = clusters, fill = clusters)) +
  geom_density_ridges(scale = 1.2, rel_min_height = 0.01,
                      color = "grey20", alpha = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Pseudotime", y = "Group") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  scale_color_npg() + scale_fill_npg() 
  
ggsave("MSD/figure/BEYOND/cluster_pseudotime.pdf", width = 5, height = 4.5)

mm <- table(PATHE$clusters, PATHE$pathological_response)
mm <- as.data.frame(mm / rowSums(mm))
colnames(mm) <- c("cluster","pathological_response","Freq")
head(mm)

ggbarplot(mm,
          x = "cluster",
          y = "Freq",
          fill = "pathological_response") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_manual(values = c('MPR'=' #6B7FA6',
                                'non-MPR'=' #CDAA59',
                                'PD'=' #8E3B2E')) +
  scale_fill_manual(values = c('MPR'=' #6B7FA6',
                               'non-MPR'=' #CDAA59',
                               'PD'=' #8E3B2E'))

ggsave("MSD/figure/BEYOND/response_bar.pdf",width = 4, height = 5)

mat <- matrix(c(
  1, 5,
  4, 10,
  11, 1,
  13, 1
), ncol = 2, byrow = TRUE)

rownames(mat) <- c("3", "1", "2", "0")
colnames(mat) <- c("MPR", "non-MPR")

# Fisher test state vs
library(dplyr)
library(ggplot2)

result_list <- lapply(1:nrow(mat), function(i) {
  state <- rownames(mat)[i]
  this <- mat[i, ]
  others <- colSums(mat[-i, ])
  
# 2x2
  mat2 <- rbind(this, others)
  rownames(mat2) <- c("this", "others")
  
  ft <- fisher.test(mat2)
  
  data.frame(
    state = state,
    OR = as.numeric(ft$estimate),
    OR_low = ft$conf.int[1],
    OR_high = ft$conf.int[2],
    p_value = ft$p.value
  )
})

plot_dat <- bind_rows(result_list)

plot_dat <- plot_dat %>%
  mutate(
    FDR = p.adjust(p_value, method = "fdr"),
    sig = FDR < 0.05,
    p_label = ifelse(sig, paste0("FDR=", signif(FDR, 2)), ""),
    y_f = factor(state, levels = rev(state)),
    label_x = max(OR_high, na.rm = TRUE) * 1.5
  )

ggplot(plot_dat, aes(x = OR, y = y_f, color = sig)) +
  geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.2, size = 0.8, alpha = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = 1, linetype = 2, color = "gray40") +
  geom_text(aes(x = label_x, label = p_label), hjust = 0, size = 3.3, color = "black") +
  scale_x_log10(limits = c(0.01, 1000)) +
  scale_color_manual(
    values = c("TRUE" = "#D55E00", "FALSE" = "#4D4D4D"),
    labels = c("TRUE" = "FDR < 0.05", "FALSE" = "ns / FDR ≥ 0.05"),
    name = ""
  ) +
  labs(x = "Odds Ratio (log scale)", y = "State") +
  theme_classic(base_size = 13) +
  theme(legend.position = "top") + 
  scale_x_log10(limits = c(0.001, 1000), breaks = c(0.001, 0.1, 1, 10, 100, 1000))

ggsave("MSD/figure/BEYOND/state_response_V2.pdf",width = 4, height = 4.5)

all_immune <- as.data.frame(read.csv("MSD/data/h5ad2/T_harmony_obs.csv"))

all_immune$state <- ifelse(grepl("S$", all_immune$sampleID), "post", "pre")
all_immune$patientID <- sub("S","",all_immune$sampleID)
all_immune <- all_immune[all_immune$state %in% c("pre"),]
a <- table(all_immune$sampleID,all_immune$sub_cell_type)
a <- as.data.frame(a / rowSums(a))
colnames(a) <- c("sampleID","cell.type","Freq")
a$patientID <- sub("S","",a$sampleID)
info <- PATHE[,c("patientID","pathological_response","histology","clusters")]
a <- a[a$patientID %in% PATHE$patientID,]
a <- merge(a, info, by = "patientID")

a$clusters <- factor(a$clusters, levels = c("0","2","1","3"))

compaired <- list(c("0","1"),
                  c("0","2"))
ggplot(a[a$cell.type %in% c("CD4T_Treg_TNFRSF9high"),], 
       aes(x = clusters, y = Freq, color = clusters,
           fill = clusters)) +
  geom_boxplot(outlier.shape = NA,width = 0.4,alpha = 0.5) +
  geom_jitter(size = 0.4, alpha = 1, width = 0.1) +
  theme_pubr() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0)) + 
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test,color = "black") +
  scale_color_manual(values = c("0" = "#E64B35FF",
                                "1" = "#4DBBD5FF",
                                "2" = "#00A087FF",
                                "3" = "#3C5488FF")) + 
  scale_fill_manual(values = c("0" = "#E64B35FF",
                               "1" = "#4DBBD5FF",
                               "2" = "#00A087FF",
                               "3" = "#3C5488FF"))+
  facet_wrap(~ cell.type,nrow = 1) 
ggsave("MSD/figure/sub_cell_type/cluster_post.pdf",width = 3.5, height = 3.5)


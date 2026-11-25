library(plyr)
library(dplyr)
library(ggpubr)
library(ggsci)
library(reshape2)
library(survival)
library(survminer)
library(GSVA)
library(Seurat)

# ====================== 1. TCGA & LUAD ======================
sce = readRDS('neoadjuvant/epi/TCGA/TCGA_Toil_Recomp_moregene.rds')
sce = sce[, sce$sampleType == "Primary.Solid.Tumor"]
# SeuratTPM
se = as.Seurat(sce, data = 'tpm')
# LUAD
se.luad = subset(se, cancer.type.abbreviation %in% c("LUAD"))

expr_tpm <- GetAssayData(se.luad, assay = "originalexp", slot = "data")
meta_tcga <- se.luad@meta.data
head(meta_tcga)

# ID
sample_ids <- colnames(expr_tpm)
cat("Number of primary TCGA LUAD tumor samples: ", length(sample_ids), "\n")

# ====================== 2. GSVACAF ======================
# CAF
geneSets = list(
  CAF_Score =  c("SPHK1","POSTN","MMP11","SULF1","COL11A1","TMEM158","TNFAIP6",
                 "ADAMTS12","SLC16A3","LOXL2","MICAL2","SDC1","NREP","CYTOR",
                 "EPSTI1","EGFL6","KIAA1217","IL7R","TDO2","CREB3L1","CHST11",
                 "ADAMTS2","C12orf75","COL27A1","PLAU","GJB2","UNC5B","LRRC15",
                 "LIMD2","TES","VCAM1","C1QTNF6","ISG15","GPM6B","INHBA","PLPP4",
                 "ADAM12","HAPLN3","EVA1A","ARL4C","MSX2","MYO1B","NTM","APBA2",
                 "RUNX2","MMP1","SLC39A14","COL10A1","MX1","COL7A1","TPST2","OAS2",
                 "PHLDA2","PPIF","SPECC1","DERL3","SERINC2","GBP1","COL5A3",
                 "HOPX","AGRN","CDKL5","POU2F2","GBP5","STRA6","DYRK2","C17orf49",
                 "SNTB1","AFAP1L2","SLAMF8","CRNDE","C3orf80","ZNF469","MKI67",
                 "PTTG1","BIRC5","HEY1","CDCP1","HMGB3","DSG2","BHLHE41","OPN3",
                 "MMP9","SMC4","TRIM59","TK1","CDC20","NPTX2","CENPP","CXCL9",
                 "BATF3","ASPM","ANLN","HOXB-AS1","LMNB2","RSAD2","IL4R",
                 "CXCL10","PCLAF","PRC1","MMP28","IFIT3","IL4I1","FOXM1",
                 "CKS2","CENPK","GRP","IDO1","MSMO1","CDKN3","CENPF","KIF20A",
                 "CEP55","KIF4A","LAMP5","ITGA4","RGS19","CCNB2","HMMR","TPX2",
                 "KIF20B","VSNL1","CNN1","CENPE","UPP1","ARNTL2","TNFSF11","RNF19B",
                 "UBE2C","GUCA1A","CXCL5","PSEN2","PTPRU","TOP2A","CA12","CXCL13",
                 "KIF23","ATP6AP1L","LAMC3","PBK","LZTS1","EEF1AKMT4","SPOCD1",
                 "LFNG","KNL1","LMNB1","AK5","CXCL11","LGR4","CCDC74A","GPR68",
                 "OSBPL3","TIMELESS","TMSB15B","GGH","NUSAP1","BMP8A","ETV4",
                 "TMEM26","USP31","GJC1","PRAG1","MAD2L1","CKAP2","MXD3","HERC6",
                 "NRARP","RNASEH2A","CCNB1","TRPA1","MLKL","RACGAP1","S100A3","GAL",
                 "MIR193BHG","CENPH","TYMS","GTPBP2","CST1","RAB40C","IFFO2",
                 "SMIM13","STMN1","MMP13","OASL","FAM174B","SGO2","ATP10A","IL24",
                 "PTP4A3","PPP2R3B","COL24A1","ZWILCH","HIST1H2AG","PARPBP",
                 "UBE2T","CHI3L2","IGKV3-11"
  )
)

# GSVA
gsva_results <- gsva(
  expr_tpm,
  geneSets,
  method = "ssgsea",
  kcdf = "Gaussian",
  abs.ranking = TRUE
)

# CAF
score_df <- as.data.frame(t(gsva_results))
score_df$sample <- rownames(score_df)
meta_tcga$sample <- rownames(meta_tcga)
meta_tcga <- left_join(meta_tcga, score_df, by = "sample")
rownames(meta_tcga) <- meta_tcga$sample
head(meta_tcga)
# ====================== 3. IFNG & / ======================
cytokine <- read.csv("MSD/data/cytosig/TCGA_LUAD_cytosig.csv")
cytokine[1:10,1:10]
rownames(cytokine) <- cytokine$X
cytokine <- cytokine[, -1]

common_pat <- intersect(rownames(meta_tcga), rownames(cytokine))
meta_tcga <- meta_tcga[common_pat, ]
cytokine <- cytokine[common_pat, ]
meta_tcga$IFNG <- cytokine$IFNG
head(meta_tcga)

median_CAF <- median(meta_tcga$CAF_Score, na.rm = TRUE)
median_IFNG <- median(meta_tcga$IFNG, na.rm = TRUE)

head(meta_tcga)

meta_tcga$CAF_Status <- ifelse(meta_tcga$CAF_Score >= median_CAF, "high", "low")
meta_tcga$IFNG_Status <- ifelse(meta_tcga$IFNG >= median_IFNG, "high", "low")

meta_tcga$Final_Group <- case_when(
  meta_tcga$CAF_Status == "high" & meta_tcga$IFNG_Status == "high" ~ "Double_High",
  TRUE ~ "Other_Groups"
)
meta_tcga$Final_Group <- factor(meta_tcga$Final_Group, levels = c("Double_High", "Other_Groups"))

cat("Sample counts by group:\n")
print(table(meta_tcga$Final_Group))

head(meta_tcga)

part <- meta_tcga[meta_tcga$ajcc_pathologic_tumor_stage %in% c("Stage IV","Stage IIIA",
                                                               "Stage IIIB"),]
# ====================== 4. TCGAOS ======================
# TCGAOS.time()OS(0/1)
meta_surv <- part[!is.na(part$OS.time) & !is.na(part$OS), ]

surv_obj <- Surv(time = meta_surv$OS.time, event = meta_surv$OS)
km_fit <- survfit(surv_obj ~ Final_Group, data = meta_surv)

# ====================== 5. KM ======================
p <- ggsurvplot(
  km_fit,
  data = meta_surv,
  conf.int = FALSE,
  pval = TRUE,
  pval.coord = c(20, 0.8),
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.col = "strata",
  censor = TRUE,
  censor.shape = 124,
  censor.size = 4,
  surv.median.line = "hv",
  legend.title = "Group",
  legend.labs = c("Double High (CAF+IFNG high)", "Other Groups"),
  ggtheme = theme_minimal(),
  title = "OS: Double-High vs Other Groups (TCGA-LUAD)",
  xlab = "Time (months)",
  ylab = "Overall survival",
  palette = c("#DC0000B2", "#3C5488B2")
)

p$plot <- p$plot +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(size = 0.7),
    axis.line.y = element_line(size = 0.7)
  )
p$table <- p$table +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(size = 0.7),
    axis.line.y = element_line(size = 0.7)
  )

print(p)

# TCGA
ggsave(
  "neoadjuvant/epi/figure/TCGA_LUAD_KM_TwoGroups_DoubleHigh_vs_Other.pdf",
  width = 8, height = 8
)

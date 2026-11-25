library(matrixStats)
library(tidyverse)
library(readxl)
library(ggpubr)
library(survival)
library(ggsurvfit)
library(survminer)
library(GSVA)
library(ggsci)

OAK.meta <- read.csv("./EGAF00005797822/go28915_anon_subsetted_BYN_n699.csv")
dim(OAK.meta)
head(OAK.meta)
table(OAK.meta$HIST)

OAK <- read.csv("./EGAF00005797824/anonymized_OAK-TPMs2.csv")
dim(OAK)
OAK[1:5,1:5]

colnames(OAK) <- gsub("^.*\\.", "", colnames(OAK))
OAK[1:10,1:10]

head(OAK.meta)

dim(OAK.meta)
dim(OAK)

OAK.meta$sample <- colnames(OAK)
head(OAK.meta)

OAK.meta <- OAK.meta[OAK.meta$HIST %in% c("SQUAMOUS"),]
OAK <- OAK[,OAK.meta$sample]

OAK <- t(t(OAK) / rowSums(t(OAK)) * 1000000)

OAK.meta.ICI <- OAK.meta %>% filter(ACTARM=="MPDL3280A")
ICI.sampleID <- OAK.meta.ICI$sample
OAK.ICI <- OAK[, ICI.sampleID]
dim(OAK.ICI)

OAK.meta.chemo <- OAK.meta %>% filter(ACTARM=="Docetaxel")
chemo.sampleID <- OAK.meta.chemo$sample
OAK.chemo <- OAK[, chemo.sampleID]
dim(OAK.chemo)

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

gsva_results <- gsva(as.matrix(OAK.ICI), 
                     geneSets, 
                     method="ssgsea",
                     kcdf="Gaussian", 
                     abs.ranking = TRUE)

gsva_t <- as.data.frame(t(gsva_results))
gsva_t$sample <- rownames(gsva_t)

fit.data <- OAK.meta.ICI %>%
  left_join(gsva_t, by = "sample")

head(fit.data)

cytokine <- read.csv("./OAK_cytosig.csv")
rownames(cytokine) <- cytokine$X
cytokine <- cytokine[,-1]
rownames(cytokine) <- gsub("^.*-", "", rownames(cytokine))
head(cytokine)

cytokine <- cytokine[fit.data$sample,]
head(cytokine)

fit.data$IFNG <- cytokine$IFNG
head(fit.data)

event_col <- "OS_CENSOR" 
time_col <- "OS_MONTHS"

fit.data <- fit.data[!is.na(fit.data[[time_col]]) & !is.na(fit.data[[event_col]]), ]

median_CAF <- median(fit.data$CAF_Score, na.rm = TRUE)
median_IFNG <- median(fit.data$IFNG, na.rm = TRUE)

fit.data$CAF_Status <- ifelse(fit.data$CAF_Score >= median_CAF, "high", "low")
fit.data$IFNG_Status <- ifelse(fit.data$IFNG >= median_IFNG, "high", "low")

Q3_CAF <- quantile(fit.data$CAF_Score, 2/3, na.rm = TRUE)
Q3_IFNG <- quantile(fit.data$IFNG, 2/3, na.rm = TRUE)

fit.data$Final_Group <- case_when(
  fit.data$CAF_Score >= Q3_CAF & fit.data$IFNG >= Q3_IFNG ~ "Double_High",
  TRUE ~ "Others"
)

print(table(fit.data$Final_Group, fit.data$BCOR))

fit <- survfit(Surv(OS_MONTHS, OS_CENSOR) ~ Final_Group, data = fit.data)

my_palette <- c("#DC0000B2", "#3C5488B2") 

p <- ggsurvplot(fit,
                data = fit.data,
                pval = TRUE,
                pval.coord = c(5, 0.1),
                conf.int = FALSE,
                risk.table = TRUE,
                risk.table.col = "strata",
                linetype = "strata",
                ggtheme = theme_classic(),
                xlab = "Months from surgery",
                title = "OS by CAF Score and IFNG", 
                legend.title = "Combined Group",
                legend.labs = c("Both High", "Others"), 
                palette = my_palette
)

print(p)

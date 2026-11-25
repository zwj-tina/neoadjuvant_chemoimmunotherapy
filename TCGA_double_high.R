library(plyr)
library(dplyr)
library(ggpubr)
library(ggsci)
library(reshape2)
library(survival)
library(survminer)
library(Seurat)
library(GSVA) 

# 1.
sce = readRDS('neoadjuvant/epi/TCGA/TCGA_Toil_Recomp_moregene.rds')
sce = sce[,sce$sampleType=="Primary.Solid.Tumor"]
se = as.Seurat(sce, data='tpm')

metadata <- as.data.frame(sce@colData@listData)

# 2. ( cDC2 Treg)
# cDC2_CD207: CD207
# Treg: Treg (FOXP3, CCR8, IL2RA, CTLA4)
geneSets = list(
  cDC2_Score = c("IL4I1"), 
  Treg_Score =  c("FOXP3","TNFRSF4","TNFRSF18","IL2RA","TBC1D4","BATF","ICOS","MIR4435-2HG","ARID5B","DNPH1", "PMAIP1","TNFRSF1B","CARD16","PIM2", "LTB","BIRC3","PRDM1","NAMPT","GBP2","CD27","SPOCK2", "SYNGR2","SAT1","PBXIP1","GBP5","CORO1B","CTSC","PKM", "UCP2","SH2D2A","ENO1","REL","NEAT1","IL32")
)

# 3. (LUSC)
se.c = subset(se, cancer.type.abbreviation %in% c("LUSC"))
dim(se.c)
# 4. (GSVA/ssGSEA)
gsva_matrix <- gsva(as.matrix(se.c@assays$originalexp@data), geneSets,
                    method="ssgsea",
                    kcdf="Gaussian", abs.ranking = TRUE)

ssGsva.score <- t(as.data.frame(gsva_matrix))
ssGsva.score <- as.data.frame(ssGsva.score)
ssGsva.score$sample <- rownames(ssGsva.score)

# 5.
dat = metadata[,c("sample","X_PATIENT","age_at_initial_pathologic_diagnosis",
                  "gender","race","ajcc_pathologic_tumor_stage",
                  "OS","OS.time","sampleType")]

rownames(dat) = dat$sample
use = intersect(rownames(dat), ssGsva.score$sample)
dat = dat[use,]
ssGsva.score = ssGsva.score[use,]

dat = cbind(dat, ssGsva.score[, c("cDC2_Score", "Treg_Score")])

dat = dat[!is.na(dat$OS.time),]
dat = dat[!is.na(dat$age_at_initial_pathologic_diagnosis),]

dat$OS.time = round(dat$OS.time/365, 2) 
dat$OS.use = dat$OS
dat$OS.time.use = dat$OS.time

cutoff_time = 7
dat$OS.use[dat$OS.time.use > cutoff_time] <- 0
dat$OS.time.use[dat$OS.time.use > cutoff_time] <- cutoff_time

# 6. (HighLow)
MP = 0.3
res.cut <- surv_cutpoint(dat, time = "OS.time.use",
                         event = "OS.use",
                         variables = c("cDC2_Score", "Treg_Score"), 
                         minprop = MP)

res.cat <- surv_categorize(res.cut)

dat$cDC2_Group = res.cat$cDC2_Score # high / low
dat$Treg_Group = res.cat$Treg_Score # high / low

# 7. ()

# cDC2 Treg "high" "Double High"
# "Others"
dat$Final_Group <- ifelse(dat$cDC2_Group == "high" & dat$Treg_Group == "high", 
                          "Double High", 
                          "Others")

# Double High ()
dat$Final_Group <- factor(dat$Final_Group, levels = c("Double High", "Others"))

print("Group size summary:")
print(table(dat$Final_Group))

# 8.
fit = survfit(Surv(OS.time.use, OS.use) ~ Final_Group, data=dat)

# 9.
# Double High/ Others
my_palette <- c("#DC0000B2", "#4DBBD5B2") 

p = ggsurvplot(fit, dat, size=0.8, 
               legend.labs = levels(dat$Final_Group),
               surv.median.line="none", 
               pval=T, 
               conf.int=F,
               palette = my_palette, 
               xlab="Years", 
               xlim = c(0, 7), 
               break.time.by = 1,
               risk.table = TRUE,
               risk.table.height = 0.28) 

print(p)

ggsave(paste0("MSD/figure/LUSC_TCGA_DoubleHigh_vs_Others.pdf"), 
       width = 6, height = 6.5)

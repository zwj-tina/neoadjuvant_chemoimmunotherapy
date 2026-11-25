library(matrixStats)
library(tidyverse)
library(readxl)
library(ggpubr)
library(survival)
library(ggsurvfit)
library(survminer)
library(GSVA)
library(ggsci)

# 1.
OAK.meta <- read.csv("neoadjuvant/epi/ICB_bulk/Genetech/EGAF00005797822/go28915_anon_subsetted_BYN_n699.csv")
dim(OAK.meta)
head(OAK.meta)
table(OAK.meta$HIST)

OAK <- read.csv("neoadjuvant/epi/ICB_bulk/Genetech/EGAF00005797824/anonymized_OAK-TPMs2.csv")
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
dim(OAK)
OAK <- t(t(OAK) / rowSums(t(OAK)) * 1000000)

OAK.meta.ICI <- OAK.meta %>% filter(ACTARM=="MPDL3280A")
ICI.sampleID <- OAK.meta.ICI$sample
OAK.ICI <- OAK[, ICI.sampleID]
dim(OAK.ICI)

OAK.meta.chemo <- OAK.meta %>% filter(ACTARM=="Docetaxel")
chemo.sampleID <- OAK.meta.chemo$sample
OAK.chemo <- OAK[, chemo.sampleID]
dim(OAK.chemo)

OAK.meta.ICI$IL4I1 <- as.data.frame(t(OAK.ICI))$IL4I1
head(OAK.meta.ICI)

# 3. ( CAF_Score IFNG vs )
event_col <- "PFS_CENSOR" 
time_col <- "PFS_MONTHS"

# NA
OAK.meta.ICI <- OAK.meta.ICI[!is.na(OAK.meta.ICI[[time_col]]) & !is.na(OAK.meta.ICI[[event_col]]), ]

median_IL4I1 <- median(OAK.meta.ICI$IL4I1, na.rm = TRUE)

# high low ( high low)
OAK.meta.ICI$IL4I1_Status <- ifelse(OAK.meta.ICI$IL4I1 >= median_IL4I1, "high", "low")

head(OAK.meta.ICI)

print("Joint group size summary:")
print(table(OAK.meta.ICI$IL4I1_Status, OAK.meta.ICI$BCOR))

# 4.
# survfit ggsurvplot (symbol)bug
fit <- survfit(Surv(PFS_MONTHS, PFS_CENSOR) ~ IL4I1_Status, data = OAK.meta.ICI)

# 22
my_palette <- c("#DC0000B2", "#3C5488B2") 

p <- ggsurvplot(fit,
                data = OAK.meta.ICI,
                pval = TRUE,
                pval.coord = c(5, 0.1),
                conf.int = FALSE,
                risk.table = TRUE,
                risk.table.col = "strata",
                linetype = "strata",
                ggtheme = theme_classic(),
                xlab = "Months from surgery",
                title = "PFS by CAF Score and IFNG", 
                legend.title = "Combined Group",
                legend.labs = c("Both High", "Others"),
                palette = my_palette
)

print(p)

# (4)
pdf("neoadjuvant/epi/figure/OAK_CAF_IFNG_Combined_TwoGroups.pdf", width = 8, height = 7)
print(p)
dev.off()

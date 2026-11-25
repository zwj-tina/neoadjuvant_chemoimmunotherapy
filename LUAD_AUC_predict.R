library(tidyverse)
library(pROC)
library(ggplot2)
library(ggpubr)
library(caret)
library(readxl)

Fibro_cell <- read.csv("MSD/data/h5ad2/Fibroblast_harmony.csv")
Fibro_cell <- Fibro_cell[,c("sampleID","sub_cell_type2")]

group <- ifelse(Fibro_cell$sub_cell_type2 %in% c("SMC_MYH11","Pericyte_RGS5","Pericyte_PTN","Fibro_ADH1B"),
                "CM02", "CM01")
Fibro_cell$group <- group

info <- as.data.frame(read_excel("MSD/data/info/clinical_info.xlsx"))
info <- info[info$pathological_response %in% c("MPR","non-MPR","PD"),]
info$pathological_response <- ifelse(info$pathological_response == "MPR", "MPR", "non-MPR")

a <- table(Fibro_cell$sampleID, Fibro_cell$group)
a <- as.data.frame(a / rowSums(a))
colnames(a) <- c("sampleID","cell.type","Freq")
a$patientID <- sub("S","",a$sampleID)
a$state <- ifelse(grepl("S$", a$sampleID), "post", "pre")

a <- a[a$patientID %in% info$patientID,]
a <- merge(a, info, by = "patientID")
a$state <- factor(a$state, levels = c("pre","post"))
a$pathological_response <- factor(a$pathological_response, levels = c("non-MPR", "MPR"))

a <- a[a$state %in% c("pre"),]
a <- a[a$histology %in% c("LUAD"),]

CM01 <- a[a$cell.type %in% c("CM01"),]

# 1. PD-L1
CM01 <- CM01 %>%
  mutate(
# MPR0/1
    MPR_binary = ifelse(pathological_response == "MPR", 1, 0),
    
# PD-L1
    PDL1_category = factor(`PD-L1`, 
                           levels = c("<1%", "1-49%", ">=50%"),
                           ordered = TRUE),
    
# PD-L1
    PDL1_numeric = case_when(
      `PD-L1` == "<1%" ~ 0,
      `PD-L1` == "1-49%" ~ 1,
      `PD-L1` == ">=50%" ~ 2
    ),
    
# PD-L1
    PDL1_high = ifelse(`PD-L1` == ">=50%", 1, 0)
  )

# 2.
# Model 1: Activated stroma
model1 <- glm(MPR_binary ~ Freq, 
              data = CM01, 
              family = binomial())

# Model 2: PD-L1 ()
model2_numeric <- glm(MPR_binary ~ PDL1_numeric, 
                      data = CM01, 
                      family = binomial())

model2_factor <- glm(MPR_binary ~ PDL1_category, 
                     data = CM01, 
                     family = binomial())

# Model 3: Combined
model3 <- glm(MPR_binary ~ Freq + PDL1_numeric, 
              data = CM01, 
              family = binomial())

# 3.
CM01$pred_module <- predict(model1, type = "response")
CM01$pred_pdl1 <- predict(model2_numeric, type = "response")
CM01$pred_combined <- predict(model3, type = "response")

# 4. ROC (MPR_binary)
roc_module <- roc(CM01$MPR_binary, CM01$pred_module, levels = c(0, 1), direction = "<")
roc_pdl1 <- roc(CM01$MPR_binary, CM01$pred_pdl1, levels = c(0, 1), direction = "<")
roc_combined <- roc(CM01$MPR_binary, CM01$pred_combined, levels = c(0, 1), direction = "<")

# AUC
auc_module <- auc(roc_module)
ci_module <- ci.auc(roc_module)

auc_pdl1 <- auc(roc_pdl1)
ci_pdl1 <- ci.auc(roc_pdl1)

auc_combined <- auc(roc_combined)
ci_combined <- ci.auc(roc_combined)

# 5.
test_module_vs_pdl1 <- roc.test(roc_module, roc_pdl1, method = "delong")
test_combined_vs_pdl1 <- roc.test(roc_combined, roc_pdl1, method = "delong")
test_combined_vs_module <- roc.test(roc_combined, roc_module, method = "delong")

# 6. ROC
# ggroc
p_roc <- ggroc(list(
  "Activated stroma" = roc_module,
  "PD-L1" = roc_pdl1,
  "Combined" = roc_combined
), legacy.axes = TRUE, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = c("Activated stroma" = "#E64B35",
               "PD-L1" = "#4DBBD5",
               "Combined" = "#00A087"),
    labels = c(
      sprintf("Activated stroma (AUC = %.3f)", auc_module),
      sprintf("PD-L1 (AUC = %.3f)", auc_pdl1),
      sprintf("Combined (AUC = %.3f)", auc_combined)
    )
  ) +
  labs(
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    title = "ROC Curves for MPR Prediction in LUAD",
    color = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.65, 0.25),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1))

print(p_roc)

ggsave("MSD/figure/ROC_curve_LUAD.pdf", p_roc, width = 5, height = 5)

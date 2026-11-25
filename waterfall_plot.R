library(ggplot2)
library(readxl)
library(tidyverse)
library(ggpubr)
data <- as.data.frame(read_excel("MSD/data/info/clinical_info.xlsx"))
head(data)
data[is.na(data)] <- "unknown"
unique(data$pathological_response_rate)
data <- data[!data$pathological_response_rate %in% c("unknown"),]
data$pathological_response_rate <- as.numeric(data$pathological_response_rate)
dim(data)
data <- data[data$histology %in% c("LUAD"),]

data <- data[data$pathological_response_detail %in% c("PD","nPR","pPR","MPR","pCR"),]

data <- data %>% 
  arrange(pathological_response_rate) %>%
  mutate(cumulative_response_rate = cumsum(pathological_response_rate) / sum(pathological_response_rate))
data$pathological_response_detail <- factor(data$pathological_response_detail,
                                            levels = c("PD","nPR","pPR","MPR","pCR"))

data <- data[order(data$pathological_response_rate), ]
data$patientID <- factor(data$patientID, levels = data$patientID)

ggplot(data, aes(x = patientID, 
                 y = pathological_response_rate, 
                 fill = pathological_response_detail,
                 color = pathological_response_detail)) +
  geom_col() +
  theme_bw() +
  theme(
    axis.text.x = element_blank(), # x
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_line(color = "black"), # y
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = c(
    "nPR" = "white",
    "pPR" = "#40B3CD",
    "MPR" = "#2695C7",
    "pCR" = "#09649D",
    "PD"  = "grey"
  )) +
  scale_color_manual(values = c(
    "nPR" = "black",
    "pPR" = "black",
    "MPR" = "black",
    "pCR" = "black",
    "PD"  = "black"
  )) +
  scale_y_reverse() +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "gray", linewidth = 1)

ggsave("MSD/figure/LUAD_water_fall_plot.pdf",width = 6, height = 4)
#PDL1
data <- as.data.frame(read_excel("MSD/data/info/clinical_info.xlsx"))
head(data)
data[is.na(data)] <- "unknown"
unique(data$`PD-L1`)
data <- data[!data$pathological_response_rate %in% c("unknown"),]
data <- data[!data$`PD-L1` %in% c("unknown"),]

PDL1.group <- c()
for(each in data$`PD-L1`){
  if(each %in% "0"){
    PDL1.group <- c(PDL1.group, "negative")
  }else{
    PDL1.group <- c(PDL1.group, "positive")
  }
}
data$PDL1.group <- PDL1.group
head(data)
#data$PDL1.group <- factor(PDL1.group, levels = c("<1%","1%~49%",">=50%"))
LUSC <- data[data$histology %in% c("LUAD"),]
a <- table(LUSC$PDL1.group,LUSC$pathological_response)
a <- as.data.frame(a / rowSums(a))
colnames(a) <- c("PDL1.group","pathological_response","Freq")
a
ggbarplot(a,x = "PDL1.group",y = "Freq",fill = "pathological_response") + 
  scale_fill_manual(values = c("MPR" = "#B2D3A4",
                               "non-MPR"="#E89DA0",
                               "PD" = "#B696B6"))
ggsave("MSD/figure/LUAD_PDL1.pdf",width = 3.5, height = 4.5)
#stage
data <- as.data.frame(read_excel("MSD/data/info/clinical_info.xlsx"))
head(data)
data[is.na(data)] <- "unknown"
unique(data$pathological_response_rate)
data <- data[!data$pathological_response_rate %in% c("unknown"),]
data <- data[!data$TNM %in% c("unknown"),]
data$pathological_response_rate <- as.numeric(data$pathological_response_rate)
head(data)
data$TNM <- factor(data$TNM, levels = c("IIA","IIB","IIIA","IIIB"))
LUSC <- data[data$histology %in% c("LUSC"),]
a <- as.data.frame(table(LUSC$TNM,LUSC$pathological_response))
colnames(a) <- c("TNM","pathological_response","number")
a
ggbarplot(a,x = "TNM",y = "number",fill = "pathological_response") + 
  scale_fill_manual(values = c("MPR" = "#9FBA95",
                               "non-MPR"="#E89DA0"))
ggsave("MSD/figure/LUSC_TNM.pdf",width = 4, height = 4.5)
#barplot
a <- data.frame(histology = c("LUSC","LUAD"),
                   MPR = c(0.5517,0.4615))
a
library(ggpubr)
ggbarplot(a,x = "histology",y = "MPR",fill = "histology") + 
  theme(legend.position = 'none') + 
  scale_fill_manual(values = c("LUSC" = "#92A5D1",
                               "LUAD" = "#D9B9D4"))

ggsave("MSD/figure/barplot.pdf",width = 2.5, height = 3.5)

data <- matrix(c(16, 13, 6,7), nrow = 2, byrow = TRUE)
colnames(data) <- c("MPR", "Non-MPR")
rownames(data) <- c("SCC", "ADC")
data <- as.table(data)
data
chi_test <- chisq.test(data)
print(chi_test)


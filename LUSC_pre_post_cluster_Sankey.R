library(dplyr)
library(ggplot2)
library(ggalluvial)
library(ggsci) 

# ---------------- 1. ----------------
pre <- read.csv("MSD/data/h5ad2/PATHE.csv")
pre <- pre[,c("patientID","sampleID","clusters","pathological_response","pathological_response_rate","PD.L1")]

# pathological_response
pre$pathological_response <- ifelse(pre$pathological_response == "MPR", "MPR", "non-MPR")

post <- read.csv("MSD/data/h5ad2/LUSC_post_PATHE.csv")
post <- post[,c("patientID","sampleID","predicted_cluster")]

# ---------------- 2. ----------------
# pre post
merged_data <- inner_join(pre, post, by = "patientID", suffix = c("_pre", "_post"))

merged_data <- merged_data %>%
# patientID pathological_response "MPR"
  mutate(
    pathological_response = if_else(
      patientID %in% c("Pa169", "Pa160", "Pa83"), 
      "MPR", 
      pathological_response
    )
  ) 

flow_data <- merged_data %>%
# patientID pathological_response "MPR"
  mutate(
    pathological_response = if_else(
      patientID %in% c("Pa169", "Pa160", "Pa83"), 
      "MPR", 
      pathological_response
    )
  ) %>%
  filter(!is.na(pathological_response) & pathological_response != "") %>%
  mutate(
    Pre_Cluster = paste0("Cluster ", clusters),
    Post_Cluster = paste0("Cluster ", predicted_cluster)
  ) %>%
  group_by(pathological_response, Pre_Cluster, Post_Cluster) %>%
  summarise(Freq = n(), .groups = "drop")

flow_data$Pre_Cluster <- factor(flow_data$Pre_Cluster, levels = c("Cluster 0","Cluster 2","Cluster 1",
                                                                  "Cluster 3"))
flow_data$Post_Cluster <- factor(flow_data$Post_Cluster, levels = c("Cluster 0","Cluster 2","Cluster 1",
                                                                  "Cluster 3"))
# ---------------- 3. ----------------
p <- ggplot(data = flow_data,
            aes(axis1 = Pre_Cluster,   
                axis2 = Post_Cluster,  
                y = Freq)) +           
  
# () - curve_type
  geom_alluvium(aes(fill = Pre_Cluster), width = 1/4, alpha = 0.7) +
  
  geom_stratum(width = 1/4, fill = "grey90", color = "black") +
  
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
  
# X ( limits bug)
  scale_x_discrete(
    limits = c("axis1", "axis2"), 
    labels = c("Pre-treatment", "Post-treatment"), 
    expand = c(0.15, 0.05)
  ) +
  scale_color_manual(values = c("Cluster 0" = "#E64B35FF",
                                "Cluster 1" = "#4DBBD5FF",
                                "Cluster 2" = "#00A087FF",
                                "Cluster 3" = "#3C5488FF")) + 
  scale_fill_manual(values = c("Cluster 0" = "#E64B35FF",
                               "Cluster 1" = "#4DBBD5FF",
                               "Cluster 2" = "#00A087FF",
                               "Cluster 3" = "#3C5488FF")) +
  theme_minimal() +
  facet_wrap(~ pathological_response, scales = "free_y") +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    strip.background = element_rect(fill = "grey95", color = NA)
  ) +
  labs(title = "Patient Cluster Transition by Pathological Response")

print(p)
# 4.
ggsave("MSD/figure/Filtered_Cluster_Sankey_by_Response.pdf", plot = p, width = 10, height = 5)

mm <- table(merged_data$predicted_cluster, merged_data$pathological_response)
mm <- as.data.frame(mm / rowSums(mm))
colnames(mm) <- c("cluster","pathological_response","Freq")
head(mm)

ggbarplot(mm, x = "cluster", y = "Freq",fill = "pathological_response") + 
  scale_color_manual(values = c('MPR'=' #6B7FA6',
                                'non-MPR'=' #CDAA59',
                                'PD'=' #8E3B2E')) +
  scale_fill_manual(values = c('MPR'=' #6B7FA6',
                               'non-MPR'=' #CDAA59',
                               'PD'=' #8E3B2E'))
ggsave("MSD/figure/post_cluster.pdf", width = 4, height = 5)

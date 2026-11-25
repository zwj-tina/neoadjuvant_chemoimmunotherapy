library(readxl)
library(dplyr)
library(stringr)
library(broom)
library(ggplot2)
library(forcats)
library(ggpubr)

##################################################################
# data preparation
##################################################################
info <-as.data.frame(read_excel("/home/zwj/raid1/MSD/data/info/clinical_info.xlsx"))
head(info)

info <- info[info$histology %in% c("LUSC"),]

info <- info[info$pathological_response %in% c("MPR","non-MPR","PD"),]

# 把 PD 归为 non-MPR
info$pathological_response <- ifelse(info$pathological_response == "PD",
                                     "non-MPR", info$pathological_response)
head(info)

# 【修改点】：保留原本的三个分组，同时正确处理 NA，并设为因子指定参照水平
info <- info %>%
  mutate(
    PD.L1 = case_when(
      is.na(`PD-L1`) | `PD-L1` == "" | `PD-L1` == "NA" ~ NA_character_, # 遇到空值保留为 NA
      `PD-L1` == "<1%"   ~ "<1%",
      `PD-L1` == "1-49%" ~ "1-49%",
      `PD-L1` == "≥50%"  ~ "≥50%",
      TRUE ~ NA_character_
    ),
    # 转换为因子，指定 levels 使得 <1% 成为参照组 (Ref)
    PD.L1 = factor(PD.L1, levels = c("<1%", "1-49%", "≥50%"))
  )

# 只保留需要的列
info <- info[, c("patientID","histology","pathological_response","PD.L1","TNM","smoking_history","ECOG")]

# 读入基础信息
base <- as.data.frame(read_excel("/home/zwj/raid1/MSD/data/info/basic_info.xlsx"))

base <- base[, c("patientID","性别","年龄")]

base <- base[base$patientID %in% info$patientID,]
colnames(base) <- c("patientID","sex","age")
base <- base %>%
  mutate(
    sex = recode(sex, "男" = "male", "女" = "female"),
    age = as.numeric(gsub("岁", "", age))
  )

median(base$age, na.rm = TRUE)

# 合并
info <- merge(info, base, by = "patientID")
head(info)
dim(info)

info <- info %>%
  mutate(
    pathological_response = factor(pathological_response,
                                   levels = c("non-MPR","MPR")), # 参照水平=non-MPR
    MPR_bin = as.integer(pathological_response == "MPR")        # 0/1
  )
head(info)

# ==== 依赖包 ====
suppressPackageStartupMessages({
  library(broom)
  library(ggplot2)
  library(forcats)
})

# ==== 载入到分析数据框 ====
df_baseline <- info

# ==== 配置 ====
# 以0/1的列作为结局，更直观
outcome <- "MPR_bin"

# 连续变量
cont_vars <- c("age")

# 分类变量
cat_vars  <- c("TNM", "PD.L1","smoking_history","ECOG")

# 是否z-score
standardize <- TRUE

# 多变量模型自变量
multi_vars <- c("age", "PD.L1", "TNM","smoking_history","ECOG")

# ==== 预处理 ====
df <- df_baseline %>%
  dplyr::select(all_of(c(outcome, cont_vars, cat_vars))) %>%  
  dplyr::mutate(across(all_of(cat_vars), ~ fct_drop(factor(.x))))  # 因子化并去掉空水平

# 不再把 outcome 强制 as.integer（已经是0/1）
# 连续变量 z-score（可选）
if (standardize) {
  df <- df %>%
    dplyr::mutate(across(all_of(cont_vars), ~ as.numeric(scale(.x)), .names = "z_{.col}"))
  cont_vars_z <- paste0("z_", cont_vars)
} else {
  cont_vars_z <- cont_vars
}

# ==== 单变量逻辑回归函数 ====
fit_uni <- function(var){
  fml <- as.formula(paste(outcome, "~", var))
  # na.action = na.exclude 保证了计算当前变量时，仅剔除当前变量或结局为 NA 的样本，不影响其他变量
  fit <- glm(fml, data = df, family = binomial(), na.action = na.exclude)
  tab <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  tab %>%
    filter(term != "(Intercept)") %>%
    transmute(
      variable = var,
      term     = term,
      OR       = estimate,
      OR_low   = conf.low,
      OR_high  = conf.high,
      p        = p.value
    )
}

# ==== 1) 单变量分析 ====
uni_list <- lapply(c(cont_vars_z, cat_vars), fit_uni)
uni_res <- bind_rows(uni_list) %>%
  group_by(variable) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  ungroup()

# ==== 从这里开始：生成“分面分组森林图”的数据并绘图 ====
# ==== 从这里开始：生成“分面分组森林图”的数据并绘图 ====
library(scales)
# 变量名美化
pretty_var <- function(x){
  x <- as.character(x)
  x <- str_replace(x, "^z_", "")
  x <- str_replace_all(x, "_", " ")
  str_to_title(x)
}

# 1) 连续变量（per SD） + 中位数[IQR]
cont_plot <- uni_res %>%
  filter(variable %in% c(cont_vars_z)) %>%
  mutate(variable = str_replace(variable, "^z_", ""),
         variable_label = pretty_var(variable),
         level = "per SD",
         level_label = "per SD",
         ref = FALSE)

cont_summ <- lapply(cont_vars, function(v){
  x <- df[[v]]
  tibble(variable = v,
         med = median(x, na.rm = TRUE),
         p25 = quantile(x, 0.25, na.rm = TRUE),
         p75 = quantile(x, 0.75, na.rm = TRUE))
}) %>% bind_rows() %>%
  mutate(dist_label = paste0("median ", round(med,1),
                             " [", round(p25,1), ", ", round(p75,1), "]"))

cont_plot <- cont_plot %>%
  left_join(cont_summ, by = "variable") %>%
  mutate(y_label = paste0(level_label, "  (", dist_label, ")"))

# 2) 分类变量：参照行 + 各水平的 n 与 MPR率
cat_plot_list <- lapply(cat_vars, function(v){
  ref_level <- levels(df[[v]])[1]
  sub <- uni_res %>% filter(variable == v)
  if (nrow(sub) > 0) {
    sub <- sub %>%
      mutate(level = gsub(paste0("^", v), "", term),
             ref = FALSE)
  } else {
    sub <- tibble(variable = v, term = character(), OR = numeric(),
                  OR_low = numeric(), OR_high = numeric(),
                  p = numeric(), q = numeric(), level = character(), ref = logical())
  }
  ref_row <- tibble(variable = v, term = paste0(v, ref_level),
                    OR = 1, OR_low = NA_real_, OR_high = NA_real_,
                    p = NA_real_, q = NA_real_,   # q 为 NA，后面会替换为1
                    level = ref_level, ref = TRUE)
  
  # 这里的 filter(!is.na(.data[[v]])) 保证了计算各组 N 和 Rate 时，只剔除当前变量为 NA 的样本
  cnt <- df %>%
    filter(!is.na(.data[[v]])) %>%
    group_by(level = .data[[v]]) %>%
    summarise(n = n(),
              events = sum(.data[[outcome]] == 1, na.rm = TRUE),
              rate = events / n, .groups = "drop")
  
  bind_rows(ref_row, sub) %>%
    left_join(cnt, by = "level") %>%
    mutate(
      variable_label = pretty_var(v),
      level_label = ifelse(ref, paste0(level, " (Ref)"),
                           paste0(level, " vs ", ref_level)),
      y_label = paste0(level_label, "  (n=", n, ", MPR ", percent(rate, accuracy = 1), ")")
    )
})

cat_plot <- bind_rows(cat_plot_list)

# 3) 合并、处理显著性与顺序
plot_dat <- bind_rows(
  cont_plot %>% dplyr::select(variable, variable_label, y_label, OR, OR_low, OR_high, p, q, ref),
  cat_plot  %>% dplyr::select(variable, variable_label, y_label, OR, OR_low, OR_high, p, q, ref)
) %>%
  mutate(
    q = ifelse(is.na(q), 1, q),          # 参照行/缺失设为1（不显著）
    sig = q < 0.05
  )

plot_dat <- plot_dat %>%
  group_by(variable_label) %>%
  mutate(y_f = fct_rev(fct_inorder(y_label))) %>%
  ungroup()

# 【修改点 1】：将原先的 format_p 改为 format_q，并修改输出文本为 FDR
format_q <- function(q_val) {
  # 使用 ifelse 处理 NA 值
  ifelse(is.na(q_val), "", 
         ifelse(q_val < 0.05,
                "FDR < 0.05",
                paste0("FDR = ", signif(q_val, 2))))
}

# 4) 横轴范围
x_min <- max(0.1, quantile(plot_dat$OR_low, 0.05, na.rm = TRUE))
x_max <- quantile(pmax(plot_dat$OR_high, plot_dat$OR), 0.95, na.rm = TRUE) * 1.5
plot_dat <- plot_dat %>%
  mutate(
    # 【修改点 2】：这里传入 q 值，并生成 q_label
    q_label = ifelse(ref, "", format_q(q)),
    # 文本放在误差条右侧一些；若没区间则基于OR放置
    label_x = ifelse(!ref,
                     pmin(ifelse(is.finite(OR_high), OR_high, OR) * 1.12, x_max * 0.98),
                     NA_real_)
  )

plot_dat

# 5) 绘图：非参照行一层 + 参照行单独叠加空心点（x=1）
p <- ggplot() +
  # 非参照行（有点有误差线）
  geom_errorbarh(data = subset(plot_dat, !ref),
                 aes(xmin = OR_low, xmax = OR_high, y = y_f, color = sig),
                 height = 0.18, alpha = 0.8) +
  geom_point(data = subset(plot_dat, !ref),
             aes(x = OR, y = y_f, color = sig), size = 2) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
  scale_x_log10(limits = c(x_min, x_max)) +
  scale_color_manual(values = c("FALSE" = "#2C3E50", "TRUE" = "#D55E00"),
                     labels = c("FALSE" = "ns / FDR≥0.05", "TRUE" = "FDR<0.05"),
                     na.value = "#4D4D4D") +
  facet_grid(variable_label ~ ., scales = "free_y", space = "free_y") +
  # 【修改点 3】：将 geom_text 中映射的 label 改为 q_label
  geom_text(data = subset(plot_dat, !ref & q_label != ""),
            aes(x = label_x, y = y_f, label = q_label),
            hjust = 0, size = 3.2, color = "black") +
  labs(x = "Odds ratio (log scale)", y = "", color = "") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    strip.text.y = element_text(face = "bold"),
    legend.position = "top"
  )

print(p)

ggsave("/home/zwj/raid1/MSD/figure/clinical_info.pdf",width = 6, height = 5)

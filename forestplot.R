library(readxl)
library(dplyr)
library(stringr)
library(broom)
library(ggplot2)
library(forcats)
library(ggpubr)

 # data preparation
info <-as.data.frame(read_excel("/home/zwj/raid1/MSD/data/info/clinical_info.xlsx"))
head(info)

info <- info[info$histology %in% c("LUSC"),]

info <- info[info$pathological_response %in% c("MPR","non-MPR","PD"),]

 # PD non-MPR
info$pathological_response <- ifelse(info$pathological_response == "PD",
                                     "non-MPR", info$pathological_response)
head(info)

 # NA
info <- info %>%
  mutate(
    PD.L1 = case_when(
      is.na(`PD-L1`) | `PD-L1` == "" | `PD-L1` == "NA" ~ NA_character_, # NA
      `PD-L1` == "<1%"   ~ "<1%",
      `PD-L1` == "1-49%" ~ "1-49%",
      `PD-L1` == "≥50%"  ~ "≥50%",
      TRUE ~ NA_character_
    ),
 # levels <1% (Ref)
    PD.L1 = factor(PD.L1, levels = c("<1%", "1-49%", "≥50%"))
  )

info <- info[, c("patientID","histology","pathological_response","PD.L1","TNM","smoking_history","ECOG")]

base <- as.data.frame(read_excel("/home/zwj/raid1/MSD/data/info/basic_info.xlsx"))

base <- base[, c("patientID","sex","age")]

base <- base[base$patientID %in% info$patientID,]
colnames(base) <- c("patientID","sex","age")
base <- base %>%
  mutate(
    sex = as.character(sex),
    age = as.numeric(age)
  )

median(base$age, na.rm = TRUE)

info <- merge(info, base, by = "patientID")
head(info)
dim(info)

info <- info %>%
  mutate(
    pathological_response = factor(pathological_response,
                                   levels = c("non-MPR","MPR")), # =non-MPR
    MPR_bin = as.integer(pathological_response == "MPR") # 0/1
  )
head(info)

suppressPackageStartupMessages({
  library(broom)
  library(ggplot2)
  library(forcats)
})

df_baseline <- info

 # 0/1
outcome <- "MPR_bin"

cont_vars <- c("age")

cat_vars  <- c("TNM", "PD.L1","smoking_history","ECOG")

 # z-score
standardize <- TRUE

multi_vars <- c("age", "PD.L1", "TNM","smoking_history","ECOG")

df <- df_baseline %>%
  dplyr::select(all_of(c(outcome, cont_vars, cat_vars))) %>%  
  dplyr::mutate(across(all_of(cat_vars), ~ fct_drop(factor(.x))))

 # outcome as.integer0/1
 # z-score
if (standardize) {
  df <- df %>%
    dplyr::mutate(across(all_of(cont_vars), ~ as.numeric(scale(.x)), .names = "z_{.col}"))
  cont_vars_z <- paste0("z_", cont_vars)
} else {
  cont_vars_z <- cont_vars
}

fit_uni <- function(var){
  fml <- as.formula(paste(outcome, "~", var))
 # na.action = na.exclude NA
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

 # ==== 1) ====
uni_list <- lapply(c(cont_vars_z, cat_vars), fit_uni)
uni_res <- bind_rows(uni_list) %>%
  group_by(variable) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  ungroup()

library(scales)
pretty_var <- function(x){
  x <- as.character(x)
  x <- str_replace(x, "^z_", "")
  x <- str_replace_all(x, "_", " ")
  str_to_title(x)
}

 # 1) per SD + [IQR]
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

 # 2) + n MPR
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
                    p = NA_real_, q = NA_real_, # q NA1
                    level = ref_level, ref = TRUE)
  
 # filter(!is.na(.data[[v]])) N Rate NA
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

 # 3)
plot_dat <- bind_rows(
  cont_plot %>% dplyr::select(variable, variable_label, y_label, OR, OR_low, OR_high, p, q, ref),
  cat_plot  %>% dplyr::select(variable, variable_label, y_label, OR, OR_low, OR_high, p, q, ref)
) %>%
  mutate(
    q = ifelse(is.na(q), 1, q), # /1
    sig = q < 0.05
  )

plot_dat <- plot_dat %>%
  group_by(variable_label) %>%
  mutate(y_f = fct_rev(fct_inorder(y_label))) %>%
  ungroup()

 # 1 format_p format_q FDR
format_q <- function(q_val) {
 # ifelse NA
  ifelse(is.na(q_val), "", 
         ifelse(q_val < 0.05,
                "FDR < 0.05",
                paste0("FDR = ", signif(q_val, 2))))
}

 # 4)
x_min <- max(0.1, quantile(plot_dat$OR_low, 0.05, na.rm = TRUE))
x_max <- quantile(pmax(plot_dat$OR_high, plot_dat$OR), 0.95, na.rm = TRUE) * 1.5
plot_dat <- plot_dat %>%
  mutate(
 # 2 q q_label
    q_label = ifelse(ref, "", format_q(q)),
 # OR
    label_x = ifelse(!ref,
                     pmin(ifelse(is.finite(OR_high), OR_high, OR) * 1.12, x_max * 0.98),
                     NA_real_)
  )

plot_dat

 # 5) + x=1
p <- ggplot() +
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
 # 3 geom_text label q_label
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

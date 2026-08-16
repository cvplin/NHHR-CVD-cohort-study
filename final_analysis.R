# ============================================================
# final_analysis.R
# NHHR and cardiovascular mortality in US adults with diabetes
# NHANES 1999-2018
# ============================================================

# ============================================================
# 01. 加载R包
# ============================================================

library(nhanesdata)
library(dplyr)
library(survey)
library(survival)
library(ggplot2)
library(rms)
library(cmprsk)
library(foreign)
library(readr)

options(timeout = 600)
options(download.file.method = "libcurl")

# 设置重试函数
read_nhanes_retry <- function(dataset, max_attempts = 5) {
  for (attempt in seq_len(max_attempts)) {
    message("正在下载 ", dataset, "：第 ", attempt, "/", max_attempts, " 次")
    result <- try(nhanesdata::read_nhanes(dataset), silent = TRUE)
    if (!inherits(result, "try-error")) {
      message(dataset, " 下载成功")
      return(result)
    }
    message(dataset, " 下载失败，准备重试")
  }
  stop(dataset, " 连续下载5次仍失败。")
}

# ============================================================
# 02. 读取数据
# ============================================================

demo <- read_nhanes_retry("demo")
bmx <- read_nhanes_retry("bmx")
bpq <- read_nhanes_retry("bpq")
diq <- read_nhanes_retry("diq")
hdl <- read_nhanes_retry("hdl")
tchol <- read_nhanes_retry("tchol")
trigly <- read_nhanes_retry("trigly")
glu <- read_nhanes_retry("glu")
ghb <- read_nhanes_retry("ghb")
smq <- read_nhanes_retry("smq")
mcq <- read_nhanes_retry("mcq")
mort <- read_nhanes_retry("mortality")

# ============================================================
# 03. 数据合并
# ============================================================

analysis_data <- demo %>%
  inner_join(bmx, by = c("seqn", "year")) %>%
  inner_join(bpq, by = c("seqn", "year")) %>%
  inner_join(diq, by = c("seqn", "year")) %>%
  inner_join(hdl, by = c("seqn", "year")) %>%
  inner_join(tchol, by = c("seqn", "year")) %>%
  inner_join(trigly, by = c("seqn", "year")) %>%
  inner_join(glu, by = c("seqn", "year")) %>%
  inner_join(ghb, by = c("seqn", "year")) %>%
  inner_join(smq, by = c("seqn", "year")) %>%
  inner_join(mcq, by = c("seqn", "year")) %>%
  left_join(mort, by = c("seqn", "year"))

# ============================================================
# 04. 数据清洗与糖尿病定义
# ============================================================

analysis_data <- analysis_data %>%
  filter(
    ridageyr >= 18,
    eligstat == 1,
    !is.na(permth_int),
    permth_int > 0,
    !is.na(lbxtc),
    !is.na(lbdhdd),
    lbdhdd > 0
  ) %>%
  mutate(
    nhhr = (lbxtc - lbdhdd) / lbdhdd,
    # 糖尿病定义
    self_reported_dm = ifelse(diq010 == 1, 1, 0),
    medication_dm = ifelse(diq050 == 1 | diq070 == 1, 1, 0),
    hba1c_dm = ifelse(lbxgh >= 6.5, 1, 0),
    fbg_dm = ifelse(lbxglu >= 7.0, 1, 0),
    diabetes = ifelse(
      self_reported_dm == 1 | medication_dm == 1 | 
        hba1c_dm == 1 | fbg_dm == 1, 1, 0
    )
  ) %>%
  filter(
    diabetes == 1,
    nhhr > 0,
    !is.na(nhhr)
  )

# ============================================================
# 05. 协变量编码
# ============================================================

analysis_data <- analysis_data %>%
  mutate(
    # 教育程度
    edu_cat = case_when(
      dmdeduc2 %in% c(1, 2) ~ "Below High School",
      dmdeduc2 == 3 ~ "High School Graduate",
      dmdeduc2 %in% c(4, 5) ~ "College or Above",
      TRUE ~ NA_character_
    ),
    # 吸烟
    ever_smoker = ifelse(smq020 == 1, "Yes", "No"),
    # 高血压
    hypertension = ifelse(bpq020 == 1 | bpq040a == 1, "Yes", "No"),
    # 性别
    sex = ifelse(riagendr == 1, "Male", "Female")
  ) %>%
  filter(
    !is.na(edu_cat),
    !is.na(ever_smoker),
    !is.na(hypertension),
    !is.na(sex),
    !is.na(indfmpir),
    !is.na(bmxbmi)
  )

# ============================================================
# 06. 最终分析数据集
# ============================================================

analysis_clean <- analysis_data %>%
  mutate(
    # 心血管死亡事件（竞争风险模型用）
    compete_event = case_when(
      mortstat == 0 ~ 0,
      mortstat == 1 & ucod_leading %in% c(6, 7) ~ 1,
      mortstat == 1 & !ucod_leading %in% c(6, 7) ~ 2,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(compete_event))

cat("最终样本量:", nrow(analysis_clean), "\n")
cat("心血管死亡数:", sum(analysis_clean$compete_event == 1), "\n")

# ============================================================
# 07. 计算NHHR四分位数
# ============================================================

nhhr_cutpoints <- quantile(analysis_clean$nhhr, 
                           probs = c(0.25, 0.5, 0.75), 
                           na.rm = TRUE)

analysis_clean <- analysis_clean %>%
  mutate(
    nhhr_q = cut(nhhr, 
                 breaks = c(-Inf, nhhr_cutpoints, Inf),
                 labels = c("Q1", "Q2", "Q3", "Q4"),
                 include.lowest = TRUE)
  )

# ============================================================
# 08. Table 1 基线特征（按NHHR四分位数）
# ============================================================

table(analysis_clean$nhhr_q)

# 连续变量汇总
analysis_clean %>%
  group_by(nhhr_q) %>%
  summarise(
    n = n(),
    age_mean = mean(ridageyr, na.rm = TRUE),
    age_sd = sd(ridageyr, na.rm = TRUE),
    nhhr_mean = mean(nhhr, na.rm = TRUE),
    nhhr_sd = sd(nhhr, na.rm = TRUE)
  ) %>%
  print()

# ============================================================
# 09. Cox回归分析
# ============================================================

# Model 1: 未调整
model1 <- coxph(Surv(permth_int, compete_event == 1) ~ nhhr_q, 
                data = analysis_clean)

# Model 2: 调整人口学因素
model2 <- coxph(Surv(permth_int, compete_event == 1) ~ nhhr_q + 
                  ridageyr + sex + ridreth1 + edu_cat + indfmpir, 
                data = analysis_clean)

# Model 3: 完全调整
model3 <- coxph(Surv(permth_int, compete_event == 1) ~ nhhr_q + 
                  ridageyr + sex + ridreth1 + edu_cat + indfmpir + 
                  bmxbmi + ever_smoker + hypertension, 
                data = analysis_clean)

summary(model3)

# ============================================================
# 10. Fine-Gray竞争风险模型
# ============================================================

analysis_clean$nhhr_num <- as.numeric(analysis_clean$nhhr_q)

fg_data <- analysis_clean %>%
  select(permth_int, compete_event, nhhr_num, ridageyr, sex) %>%
  filter(
    !is.na(permth_int),
    !is.na(compete_event),
    !is.na(nhhr_num),
    !is.na(ridageyr),
    !is.na(sex)
  )

fg_data$sex_num <- ifelse(fg_data$sex == "Male", 1, 0)
covariates <- as.matrix(fg_data[, c("nhhr_num", "ridageyr", "sex_num")])

fg_model <- crr(
  ftime = fg_data$permth_int,
  fstatus = fg_data$compete_event,
  cov1 = covariates,
  failcode = 1,
  cencode = 0
)

summary(fg_model)

# ============================================================
# 11. 比例风险假设检验
# ============================================================

ph_test <- cox.zph(model3)
print(ph_test)

# ============================================================
# 12. 随访人年数计算
# ============================================================

person_years <- analysis_clean %>%
  group_by(nhhr_q) %>%
  summarise(
    n = n(),
    total_followup_months = sum(permth_int, na.rm = TRUE),
    total_followup_years = total_followup_months / 12,
    cvd_deaths = sum(compete_event == 1, na.rm = TRUE),
    cvd_rate_per_1000_py = (cvd_deaths / total_followup_years) * 1000
  )

print(person_years)

# ============================================================
# 13. 限制性立方样条（RCS）
# ============================================================

dd <- datadist(analysis_clean[, c("nhhr", "ridageyr", "sex", 
                                  "ridreth1", "edu_cat", "indfmpir",
                                  "bmxbmi", "ever_smoker", "hypertension",
                                  "permth_int", "compete_event")])
options(datadist = "dd")

rcs_model <- cph(Surv(permth_int, compete_event == 1) ~ 
                   rcs(nhhr, 5) + ridageyr + sex + ridreth1 + 
                   edu_cat + indfmpir + bmxbmi + ever_smoker + hypertension,
                 data = analysis_clean, x = TRUE, y = TRUE)

anova(rcs_model)

# ============================================================
# 14. 生成Figure 2 (RCS曲线)
# ============================================================

p <- Predict(rcs_model, nhhr = seq(0.5, 10, by = 0.1),
             ridageyr = median(analysis_clean$ridageyr, na.rm = TRUE),
             sex = "Male",
             ridreth1 = "Non-Hispanic White",
             edu_cat = "High School Graduate",
             indfmpir = median(analysis_clean$indfmpir, na.rm = TRUE),
             bmxbmi = median(analysis_clean$bmxbmi, na.rm = TRUE),
             ever_smoker = "No",
             hypertension = "No",
             fun = exp)

ggplot(p, aes(x = nhhr, y = yhat)) +
  geom_line(color = "#2E86C1", size = 1.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#2E86C1") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 0.8) +
  labs(x = "NHHR", y = "Hazard Ratio (HR)") +
  theme_minimal()

ggsave("Figure_2_RCS_English.tiff", width = 6, height = 5, dpi = 300, compression = "lzw")

# ============================================================
# 15. 生成Figure 3 (KM生存曲线)
# ============================================================

library(survminer)

km_fit <- survfit(Surv(permth_int, compete_event == 1) ~ nhhr_q, 
                  data = analysis_clean)

ggsurvplot(km_fit, data = analysis_clean,
           palette = c("#2E86C1", "#28B463", "#F39C12", "#E74C3C"),
           xlab = "Follow-up Time (months)",
           ylab = "Survival Probability",
           legend.title = "NHHR Quartile",
           legend.labs = c("Q1 (0.20–1.85)", "Q2 (1.85–2.54)", 
                           "Q3 (2.54–3.45)", "Q4 (3.45–26.67)"),
           risk.table = TRUE)

ggsave("Figure_3_KM_English.tiff", width = 8, height = 6, dpi = 300, compression = "lzw")
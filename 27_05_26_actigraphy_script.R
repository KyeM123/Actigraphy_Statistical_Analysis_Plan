
### ACTIGRAPHY DATA ANALYSES

#Authors: R.V.R, K.M.
#Date: 27/05/2026
#Purpose: To address wristwatch actigraphy data collected pre and post each session

# --------------------------------------------------------------------------------------------
# CODING DATA FOR ANALYSIS 
# --------------------------------------------------------------------------------------------

#  Import data
actigraphy = read.csv("Z:/Analysis/Kye/Actigraphy Analysis/Actigraphy_data_sheet_for_analysis.csv")
setwd("Z:/Analysis/Kye/Actigraphy Analysis/R Scripts/Output")

# label NAs
if (exists("actigraphy")) {
  actigraphy[actigraphy == "na"] <- NA
}

# install.packages
# install.packages(c(
#   "tidyverse", "dplyr", "lubridate", "summarytools", "xfun", "ggplot2",
#   "naniar", "MuMIn", "car", "lme4", "lmerTest", "emmeans", "glmmTMB",
#   "performance", "DHARMa", "broom.mixed", "tibble", "readr", "forcats", "see", "nortest"
# ))

# loading packages
library(tidyverse)
library(dplyr)
library(lubridate)
library(summarytools)
library(xfun)
library(ggplot2)
library(naniar)
library(MuMIn)
library(car)
library(lme4)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(performance)
library(DHARMa)
library(broom.mixed)
library(tibble)
library(readr)
library(forcats)
library(see)
library(nortest)

# Creating (nominal) categorical variables

actigraphy$id <- as.factor(actigraphy$id)
actigraphy$arm <- as.factor(actigraphy$arm)
actigraphy$treatment <- as.factor(actigraphy$treatment)
actigraphy$preorpost <- as.factor(actigraphy$preorpost)
actigraphy <- actigraphy %>%
  group_by(id) %>%
  mutate(order = if_else(first(treatment[arm == 1]) == "A", "A_first", "B_first")) %>%
  ungroup()
actigraphy$order <- as.factor(actigraphy$order)
actigraphy$mean_SPT <- as.numeric(actigraphy$mean_SPT)
actigraphy$SD_SPT <- as.numeric(actigraphy$SD_SPT)
actigraphy$mean_sleep_duration <- as.numeric(actigraphy$mean_sleep_duration)
actigraphy$SD_sleep_duration <- as.numeric(actigraphy$SD_sleep_duration)
actigraphy$sleep_efficiency <- as.numeric(actigraphy$sleep_efficiency)
actigraphy$mean_WASO <- as.numeric(actigraphy$mean_WASO)
actigraphy$SD_WASO <- as.numeric(actigraphy$SD_WASO)

class(actigraphy$mean_SPT)
class(actigraphy$mean_sleep_duration)
class(actigraphy$sleep_efficiency)
class(actigraphy$mean_WASO)


## Check the structure of data

str(actigraphy) 
View(actigraphy) #checking that the data looks fine at first glance
View(dfSummary(actigraphy)) #and make sure it looks good in more detail 


# -----------------------------------------------------------------------------------------
# CLEANING DATA - OUTLIER DETECTION AND REMOVAL
# -----------------------------------------------------------------------------------------

# Flag outliers using the Median Absolute Deviation (MAD) method across all outcomes
flag_mad <- function(x, threshold = 3) {
  med <- median(x, na.rm = TRUE)
  m   <- mad(x, na.rm = TRUE)
  abs(x - med) > threshold * m
}

actigraphy <- actigraphy %>%
  mutate(
    outlier_MAD_SPT = flag_mad(mean_SPT),
    outlier_MAD_duration = flag_mad(mean_sleep_duration),
    outlier_MAD_eff = flag_mad(sleep_efficiency),
    outlier_MAD_WASO = flag_mad(mean_WASO),
    # Heuristic plausibility flags
    flag_implausibly_low_SPT = mean_SPT < 2,
    flag_suspicious_perfect = sleep_efficiency >= 99 & mean_WASO <= 0.1
  )

# Inspect the flagged rows
actigraphy %>% 
  filter(outlier_MAD_SPT | flag_implausibly_low_SPT | flag_suspicious_perfect) %>%
  arrange(id, arm, preorpost) %>%
  select(id, arm, treatment, preorpost, mean_SPT, sleep_efficiency, mean_WASO, SD_SPT)
  # R01_2pre: SPT = 1.667hr, sleep efficiency = 100%, WASO = 0. Strongly suggestive of non-wear or device failure
  # R33_1post: SPT = 0.496hr, sleep efficiency = 100%, WASO = 0. Strongly suggestive of non-wear or device failure
  # R27_1post: SPT = 1.29hr, sleep efficiency = 95%, WASO = 0.059. Likely non-wear
  # Other values are less clearly invalid - could be severe sleep restriction or non-wear/device failure

# Visualise flagged SPT values - spaghetti plot with flags
ggplot(actigraphy, aes(preorpost, mean_SPT, group = id, colour = outlier_MAD_SPT)) +
  geom_line(alpha = 0.35) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "#d62728")) +
  facet_wrap(~ treatment + arm, ncol = 2) +
  labs(title = "SPT by period with MAD-flagged outliers",
       x = "Period", y = "Sleep Period Time (h)", colour = "MAD outlier")


# Spaghetti plot looking at specific flagged participants
flag_ids <- actigraphy %>% filter(outlier_MAD_SPT) %>% pull(id) %>% unique()

ggplot(filter(actigraphy, id %in% flag_ids),
       aes(preorpost, mean_SPT, group = id, colour = id)) +
  geom_line() + geom_point(size = 2) +
  theme(legend.position = "bottom") +
  labs(title = "Flagged IDs only", x = "", y = "SPT (h)")
  # R01 and R27 show sudden implausible changes from pre to post


# R01 session 3 (2_post), R33 session 2 (1_post) and R27 session 2 (1_post) flagged due to highly likely non-wear or device failure
# Comparison with self-reported PSQI SPT data:
  # R01 session 3: 6.5 hrs reported sleep vs. 1.667 hrs recorded.
  # R33 session 2: PSQI not completed. However, SPT of 0.496 hrs with perfect sleep efficiency remains physiologically implausible
  # R27 session 2: PSQI not completed. However, SPT of 1.29hr with 95% SE and WASO 0.059 remains physiologically implausible. There is also demonstrated above an implausible change from pre to post
# Decision made to exclude the following rows due to likely non-wear or device error:
  # R01 session 3 (arm 2, pre): SPT = 1.7 hr vs PSQI 6.5 hr; SE = 100%; WASO = 0
  # R33 session 2 (arm 1, post): SPT = 0.496 h; SE = 100%; WASO = 0
  # R27 session 2 (arm 1, post): SPT = 1.29 h; SE = 95%; implausible pre-to-post change

# Standardise key columns to character so joins match exactly
actigraphy2 <- actigraphy %>%
  mutate(id = as.character(id), arm = as.character(arm), preorpost = as.character(preorpost))

to_drop <- tribble(
  ~id,  ~arm, ~preorpost,
  "R33","1",  "post",
  "R27","1",  "post",
  "R01","2",  "pre"
)

actigraphy2 %>%
  semi_join(to_drop, by = c("id","arm","preorpost")) %>%
  arrange(id, arm, preorpost)

# Remove the rows
actigraphy_clean <- actigraphy2 %>%
  anti_join(to_drop, by = c("id","arm","preorpost"))

# Verify count removed (should be 3)
nrow(actigraphy2) - nrow(actigraphy_clean)


# Re-check histograms before and after removal
par(mfrow = c(1,2))
hist(actigraphy$mean_SPT, breaks = 20, col = "grey85",
     main = "SPT — BEFORE", xlab = "Hours")
hist(actigraphy_clean$mean_SPT, breaks = 20, col = "grey85",
     main = "SPT — AFTER", xlab = "Hours")
par(mfrow = c(1,1))
# Improved left skew


# Re-code all variables on clean dataset (factors and numerics)
# Explicit re-casting of numeric columns is necessary here because anti_join() can return columns as character depending on how the join resolves types.
actigraphy_clean <- actigraphy_clean %>%
  mutate(
    # Factors
    id              = as.factor(id),
    arm             = as.factor(arm),
    treatment       = as.factor(treatment),
    preorpost       = as.factor(preorpost),
    order           = as.factor(order),
    preorpost       = fct_relevel(preorpost, "pre", "post"),
    treatment       = fct_relevel(treatment, "A"),
    # Numeric outcomes
    mean_SPT            = as.numeric(mean_SPT),
    SD_SPT              = as.numeric(SD_SPT),
    mean_sleep_duration = as.numeric(mean_sleep_duration),
    SD_sleep_duration   = as.numeric(SD_sleep_duration),
    sleep_efficiency    = as.numeric(sleep_efficiency),
    mean_WASO           = as.numeric(mean_WASO),
    SD_WASO             = as.numeric(SD_WASO)
  )

# ----------------------------------------------------------------------------------------
# DATA CHECKS AFTER CLEANING
# ----------------------------------------------------------------------------------------

# Structure check - to confirm all variables are coded correctly after cleaning
str(actigraphy_clean, list.len=ncol(actigraphy_clean)) 
View(actigraphy_clean) #checking that the data looks fine at first glance
view(dfSummary(actigraphy_clean)) #and make sure it looks good in more detail 

# Visually inspect frequency and unique value for each categorical variable
cat("\n--- Categorical variable checks ---\n")

cat("\nID (n participants):\n"); print(table(actigraphy_clean$id))
cat("\nUnique IDs:\n");   print(unique(actigraphy_clean$id))

cat("\nArm:\n");  print(table(actigraphy_clean$arm))
cat("\nUnique arm values:\n");  print(unique(actigraphy_clean$arm))

cat("\nTreatment:\n");  print(table(actigraphy_clean$treatment))
cat("\nUnique treatment values:\n"); print(unique(actigraphy_clean$treatment))

cat("\nPre or post:\n");  print(table(actigraphy_clean$preorpost))
cat("\nUnique preorpost values:\n"); print(unique(actigraphy_clean$preorpost))

cat("\nOrder:\n");  print(table(actigraphy_clean$order))
cat("\nUnique order values:\n");  print(unique(actigraphy_clean$order))

# Inspect missing data
cat("\n--- Missing data ---\n")
vis_miss(actigraphy_clean) # Visual heatmap of missingness
print(colSums(is.na(actigraphy_clean)))  # Numeric count per column
# ~9.4% missing due to insufficient or invalid wear data

# Duplicate record check — each participant should appear once per condition
# Any n > 1 indicates a data entry or merging error
cat("\n--- Duplicate record check ---\n")
dupe_check <- actigraphy_clean %>%
  mutate(condition = paste0(str_to_title(as.character(preorpost)), treatment)) %>%
  group_by(id, condition) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

if (nrow(dupe_check) == 0) {
  cat("No duplicate records found. Each participant appears once per condition.\n")
} else {
  cat("WARNING: Duplicate records detected for the following id/condition combinations:\n")
  print(dupe_check)
}
# No duplicate records found

# ------------------------------------------------------------------------------------------
# EXPORT AND RE-IMPORT CLEANED DATA
# ------------------------------------------------------------------------------------------
readr::write_csv(actigraphy_clean, "actigraphy_clean.csv")
cat("Cleaned dataset exported to actigraphy_clean.csv\n")
cat("Rows:", nrow(actigraphy_clean), "| Columns:", ncol(actigraphy_clean), "\n")

# Re-import the cleaned dataset 
actigraphy_clean <- read.csv("actigraphy_clean.csv")

# Re-code factors and numerics after re-import
actigraphy_clean <- actigraphy_clean %>%
  mutate(
    id              = as.factor(id),
    arm             = as.factor(arm),
    treatment       = as.factor(treatment),
    preorpost       = as.factor(preorpost),
    order           = as.factor(order),
    preorpost       = fct_relevel(preorpost, "pre", "post"),
    treatment       = fct_relevel(treatment, "A"),
    mean_SPT            = as.numeric(mean_SPT),
    SD_SPT              = as.numeric(SD_SPT),
    mean_sleep_duration = as.numeric(mean_sleep_duration),
    SD_sleep_duration   = as.numeric(SD_sleep_duration),
    sleep_efficiency    = as.numeric(sleep_efficiency),
    mean_WASO           = as.numeric(mean_WASO),
    SD_WASO             = as.numeric(SD_WASO)
  )

cat("Cleaned dataset re-imported and recoded successfully.\n")
str(actigraphy_clean)


# -----------------------------------------------------------------------------------------
# DESCRIPTIVE STATISTICS
# -----------------------------------------------------------------------------------------
outcomes <- c("mean_SPT", "mean_sleep_duration", "sleep_efficiency", "mean_WASO")

desc_list <- lapply(outcomes, function(v) {
  actigraphy_clean %>%
    group_by(treatment, preorpost) %>%
    summarise(
      n = sum(!is.na(.data[[v]])),
      mean = mean(.data[[v]], na.rm = TRUE),
      sd = sd(.data[[v]],   na.rm = TRUE),
      median = median(.data[[v]], na.rm = TRUE),
      q1 = quantile(.data[[v]], 0.25, na.rm = TRUE),
      q3 = quantile(.data[[v]], 0.75, na.rm = TRUE),
      min = min(.data[[v]], na.rm = TRUE),
      max = max(.data[[v]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(outcome = v) # labels the outcomes
})

desc_all <- bind_rows(desc_list) %>%
  select(outcome, treatment, preorpost, everything())

print(desc_all)
readr::write_csv(desc_all, "descriptive_statistics_all_outcomes.csv") 

#Normality tests on raw outcome distributions 
# Anderson-Darling test used as supplementary normality assessment to Shapiro-Wilk as Shapiro-Wilk can be oversensitive in small samples
sw_results <- lapply(outcomes, function(v) {
  x <- actigraphy_clean[[v]]
  x <- x[!is.na(x)]
  sw  <- shapiro.test(x)
  ad  <- ad.test(x)
  data.frame(outcome = v,
             SW_W = sw$statistic, SW_p = sw$p.value,
             AD_A = ad$statistic, AD_p = ad$p.value)
})
sw_table <- bind_rows(sw_results)
print(sw_table)
readr::write_csv(sw_table, "normality_tests_raw_outcomes.csv")
# All outcomes demonstrate p < 0.05 on Anderson-Darling and Shapiro-Wilk test, suggesting departure of the raw data from normality
# However, Shapiro-Wilk test can be overly sensitive and sleep data is prone to non-normality
# Normality of residuals of the fitted model will be assessed later

# ------------------------------------------------------------------------------------------
# CARRYOVER TESTING
# ------------------------------------------------------------------------------------------
period2_pre <- actigraphy_clean %>%
  filter(arm == "2", preorpost == "pre")

outcomes_carryover <- c("mean_SPT", "mean_sleep_duration", "sleep_efficiency",
                        "mean_WASO")

outcome_labels_carryover <- c(
  mean_SPT = "Sleep Period Time (h)",
  mean_sleep_duration = "Sleep Duration (h)",
  sleep_efficiency = "Sleep Efficiency (%)",
  mean_WASO = "Wake After Sleep Onset (h)"
)

cat("\n--- Carryover tests: Period 2 pre-treatment values by treatment order ---\n")

carryover_results <- list()

for (v in outcomes_carryover) {
  lbl <- outcome_labels_carryover[v]
  
  # Boxplot: Period 2 pre-treatment values by order group
  p <- ggplot(period2_pre, aes(x = order, y = .data[[v]], fill = order)) +
    geom_boxplot(alpha = 0.6, outlier.shape = 21) +
    scale_fill_manual(values = c("A_first" = "#2196F3", "B_first" = "#F44336")) +
    labs(title = paste("Carryover check —", lbl),
         subtitle = "Period 2 pre-treatment values by treatment order",
         x = "Treatment order", y = lbl) +
    theme_bw() +
    theme(legend.position = "none")
  print(p)
  
  # Independent samples t-test: A_first vs B_first at Period 2 pre
  tt <- t.test(period2_pre[[v]] ~ period2_pre$order)
  
  cat("\nOutcome:", lbl, "\n")
  cat("  t =", round(tt$statistic, 3),
      "| df =", round(tt$parameter, 1),
      "| p =", round(tt$p.value, 4), "\n")
  if (tt$p.value > 0.05) {
    cat("  → No evidence of carryover.\n")
  } else {
    cat("  → WARNING: Significant difference detected. Consider carryover effects in interpretation.\n")
  }
  
  # Store a tidy summary row
  carryover_results[[v]] <- data.frame(
    outcome     = v,
    label       = lbl,
    t           = round(tt$statistic, 3),
    df          = round(tt$parameter, 1),
    p.value     = round(tt$p.value, 4),
    mean_Afirst = round(tt$estimate[1], 3),
    mean_Bfirst = round(tt$estimate[2], 3)
  )
}

carryover_tbl <- bind_rows(carryover_results)
print(carryover_tbl)
readr::write_csv(carryover_tbl, "carryover_tests_all_outcomes.csv")
# No evidence of carryover in any outcome


# -------------------------------------------------------------------------------------------
# HELPER FUNCTIONS FOR STATISTICAL ANALYSIS PIPELINE
# -------------------------------------------------------------------------------------------

# In our secondary outcome analysis we will be looking at sleep parameters (sleep period time, sleep duration, sleep efficiency, wake after sleep onset).
# The aim is to observe whether the parameters are affected by the treatment phase (i.e. active or inactive) and type of treatment (i.e. CBD or placebo).

# (1) We have the following factors:
# Treatment = 2-Level (WS) Categorical: CBD, Placebo - treatment
# Time = 2-Level (WS) Categorical: Pre-Treatment, Post-Treatment - preorpost
# Treatment Order = 2-Level (BS) Categorical: CBD First, Placebo First - order
# Treatment Period = 2-Level (WS) Categorical: Period 1, Period 2 - arm

# (2) We are primarily looking for one of the following effects (and we will select the model with the lowest AIC):
# a) Treatment * Time
# b) Treatment * Time + Treatment Order - This will identify and control for systematic differences between the CBD-PLA and PLA-CBD groups (e.g., if randomisation is not 100% effective)
# c) Treatment * Time + Treatment Period - This will identify and control for systematic differences between the first and second treatment period
# d) Treatment * Time * Treatment Order - This will identify and control for carryover effects (i.e., time- and/or treatment-dependent differences between the CBD-PLA and PLA-CBD groups)

# All data is continuous in nature and will be measured using a linear model. 
# The model with the lowest AIC will be identified and undergo assumption testing (i.e. collinearity and normality (Anderson-Darling test, p>0.05) and homoscedasticity (Breusch-Pagan test, p>0.05))
# If the continuous model violates normality or homoscedasticity, the data will undergo SQRT transformation, then will progress to a LOG transformation
# If the LOG and SQRT structure does not pass, a Gamma generalised linear model will be used
# The models will be manually checked through visual inspection as the formal statistical tests for normality and homoscedasticity can be overly sensitive
# There will be the option to override the automatic model selection based on visual inspection and manual checks

# Through the above we have the capacity to identify situations where participants who receive CBD First differ from those who receive Placebo First (particularly at the Pre-Treatment time point) during Treatment Period 2.

# RANDOM EFFECTS STRUCTURE SELECTION
# Using mean_SPT. The random effects structure is a feature of the study design rather than the outcome, so the winning structure can be applied across outcomes

re_m1 <- lmer(mean_SPT ~ + (1 + treatment | id), data = actigraphy_clean) # model with random intercept + random slope for treatment per participant
re_m2 <- lmer(mean_SPT ~ + (1 | id), data = actigraphy_clean) # model with random intercept only

cat("\n--- Random effects structure selection (AICc) ---\n")
re_aic <- MuMIn::AICc(re_m1, re_m2)
print(re_aic)

# Inspect which model is preferred
best_re <- rownames(re_aic)[which.min(re_aic$AICc)]
cat("\nBest random effects structure:", best_re, "\n")
cat("→ (1|id) selected and applied consistently across all outcomes.\n") # re_m2 appears to be the better model with the lowest AICc
# Note: if re_m1 had been preferred, the random effects formula in select_fixed_effects() below would need to be updated to (1 + treatment | id) for all outcomes.

# Clean up temporary models
rm(re_m1, re_m2, re_aic, best_re)


# FIXED-EFFECTS MODEL SELECTION USING AICc
select_fixed_effects <- function(outcome_var, data) {
  
  f_list <- list(
    a = as.formula(paste(outcome_var, "~ treatment*preorpost + (1|id)")),
    b = as.formula(paste(outcome_var, "~ treatment*preorpost + order + (1|id)")),
    c = as.formula(paste(outcome_var, "~ treatment*preorpost + arm + (1|id)")),
    d = as.formula(paste(outcome_var, "~ treatment*preorpost*order + (1|id)"))
  )
  
  fits <- lapply(f_list, function(f) {
    tryCatch(
      lmer(f, data = data, REML = FALSE, na.action = na.exclude),
      error = function(e) NULL
    )
  })
  
  valid  <- Filter(Negate(is.null), fits)
  aic_df <- MuMIn::AICc(valid[[1]], valid[[2]], valid[[3]], valid[[4]])
  best   <- which.min(aic_df$AICc)
  
  list(best_model = valid[[best]], aicc_table = aic_df, selected = names(valid)[best])
}

# Fixed-effects model selection via AICc for Gamma GLMM (glmmTMB)
# Parallel to select_fixed_effects() but uses glmmTMB with Gamma(link = "log")
# Called in the Gamma GLMM fallback branch of check_and_transform() to ensure AICc-based fixed-effects selection is applied consistently across all model types.
# Only positive outcome values are used (Gamma requires y > 0).
select_fixed_effects_gamma <- function(outcome_var, data) {
  
  # Restrict to positive values as required by Gamma family
  dat_pos <- data %>% filter(.data[[outcome_var]] > 0)
  
  f_list <- list(
    a = as.formula(paste(outcome_var, "~ treatment*preorpost + (1|id)")),
    b = as.formula(paste(outcome_var, "~ treatment*preorpost + order + (1|id)")),
    c = as.formula(paste(outcome_var, "~ treatment*preorpost + arm + (1|id)")),
    d = as.formula(paste(outcome_var, "~ treatment*preorpost*order + (1|id)"))
  )
  
  fits <- lapply(f_list, function(f) {
    tryCatch(
      glmmTMB(f, data = dat_pos, family = Gamma(link = "log")),
      error = function(e) NULL
    )
  })
  
  valid <- Filter(Negate(is.null), fits)
  
  if (length(valid) == 0) {
    message("All Gamma GLMM models failed to fit.")
    return(NULL)
  }
  
  # MuMIn::AICc does not support glmmTMB — compute AICc manually
  aic_df <- do.call(rbind, lapply(seq_along(valid), function(i) {
    m <- valid[[i]]
    k <- base::attr(logLik(m), "df")
    n <- nrow(dat_pos)
    aic <- AIC(m)
    aicc <- aic + (2 * k^2 + 2 * k) / (n - k - 1)
    data.frame(row.names = names(valid)[i], df = k, AICc = aicc)
  }))
  best <- which.min(aic_df$AICc) 
  
  cat("Selected Gamma GLMM fixed-effects structure:", names(valid)[best], "\n")
  print(aic_df)
  
  list(best_model = valid[[best]], aicc_table = aic_df, selected = names(valid)[best])
}

# ASSUMPTION CHECKS AND TRANSFORMATIONS
  # This function fits the best-fitting fixed-effects model for a given outcome, checks normality (Shapiro-Wilk) and homoscedasticity (Breusch-Pagan) formally and generates visual inspection plots (histogram of residuals, QQ plot, fitted vs residuals) at each stage (raw, SQRT, LOG) to allow human review alongside the automated decision.
# Overriding the formal test:
  #   The function selects a transform (or no transform) based on p-values alone by default. However, the Shapiro-Wilk test can be overly sensitive
  # If the QQ plot and residual histogram look approximately normal despite a significant Shapiro-Wilk p-value, can override the automatic decision by passing the `override` argument.
# The `override` argument accepts one of:
  # NULL   — default; the function decides automatically based on p-values
  # "raw"  — keep the raw (untransformed) model regardless of test results
  # "sqrt" — force the SQRT-transformed model regardless of test results
  # "log"  — force the LOG-transformed model regardless of test results
  # "gamma"— force the Gamma GLMM regardless of test results

# Example usage to override for a specific outcome:
  # Run normally for all outcomes first to review plots, then re-run any outcome where you want to override, e.g.: result <- check_and_transform("mean_SPT", actigraphy_clean, override = "raw")

# Normality check: returns TRUE if residuals are considered normal
# Uses Anderson-Darling as primary test; Shapiro-Wilk reported alongside for transparency
check_normality_ad <- function(model) {
  r <- resid(model)
  r <- r[!is.na(r)]
  ad  <- ad.test(r)
  sw  <- shapiro.test(r)
  cat("    Anderson-Darling p =", round(ad$p.value, 4),
      "| Shapiro-Wilk p =", round(sw$p.value, 4), "\n")
  ad$p.value > 0.05  
}


check_and_transform <- function(outcome_var, data, override = NULL) {
  # Internal helper: produce visual inspection plots for a fitted lmer model
  plot_assumptions <- function(mod, label) {
    r <- resid(mod)
    f <- fitted(mod)
    
    dev.new(width = 10, height = 4)
    par(mfrow = c(1, 3))
    
    # Histogram of residuals
    hist(r, breaks = 20, col = "grey85",
         main = paste0("Residuals — ", label, "\n(", outcome_var, ")"),
         xlab = "Residual")
    
    # QQ plot
    qqnorm(r, main = paste0("QQ plot — ", label, "\n(", outcome_var, ")"))
    qqline(r, col = "red", lty = 2)
    
    # Fitted vs residuals
    plot(f, r,
         main = paste0("Fitted vs Residuals — ", label, "\n(", outcome_var, ")"),
         xlab = "Fitted values", ylab = "Residuals")
    abline(h = 0, lty = 2, col = "red")
    
    par(mfrow = c(1, 1))
  }
  
  message("\n============================================================")
  message("Outcome: ", outcome_var)
  message("============================================================")
  
  sel <- select_fixed_effects(outcome_var, data)
  model <- sel$best_model
  cat("Selected fixed-effects structure:", sel$selected, "\n")
  print(sel$aicc_table)
  
  # Raw scale assumption checks
  bp <- check_heteroscedasticity(model)
  p_bp <- as.numeric(bp)
  homosc <- p_bp > 0.05
  
  cat("\n--- Assumption checks (raw scale) ---\n")
  cat("Breusch-Pagan p =", p_bp, "| Homoscedastic:", homosc, "\n")
  normal <- check_normality_ad(model)
  cat("→ Inspect plots before accepting or overriding the automated decision.\n")
  
  plot_assumptions(model, "Raw scale")
  
  transform_label <- "raw"
  
  # Override: return raw model immediately if requested
  if (!is.null(override) && override == "raw") {
    cat("\n→ OVERRIDE applied: retaining raw scale model as instructed.\n")
    cat("→ Final model transform: raw\n\n")
    return(list(model = model, transform = "raw", data = data))
  }
  
  if (!normal || !homosc) {
    message("  → Attempting SQRT transform...")
    
    data[[paste0(outcome_var, "_sqrt")]] <- sqrt(data[[outcome_var]])
    sel2 <- select_fixed_effects(paste0(outcome_var, "_sqrt"), data)
    model2  <- sel2$best_model
    bp2 <- check_heteroscedasticity(model2)
    p_bp2 <- as.numeric(bp2)
    homosc2 <- p_bp2 > 0.05
    
    cat("\n--- Assumption checks (SQRT transform) ---\n")
    cat("Breusch-Pagan p =", p_bp2, "| Homoscedastic:", homosc2, "\n")
    normal2 <- check_normality_ad(model2)
    cat("→ Inspect plots before accepting or overriding the automated decision.\n")
    
    plot_assumptions(model2, "SQRT transform")
    
    # Override: return SQRT model immediately if requested
    if (!is.null(override) && override == "sqrt") {
      cat("\n→ OVERRIDE applied: retaining SQRT model as instructed.\n")
      cat("→ Final model transform: sqrt\n\n")
      return(list(model = model2, transform = "sqrt", data = data))
    }
    
    if (normal2 & homosc2) {
      model <- model2
      transform_label <- "sqrt"
    } else {
      message("  → SQRT failed. Attempting LOG transform...")
      
      # Shift if non-positive values exist (LOG requires y > 0)
      min_val <- min(data[[outcome_var]], na.rm = TRUE)
      shift <- if (min_val <= 0) abs(min_val) + 0.001 else 0
      data[[paste0(outcome_var, "_log")]] <- log(data[[outcome_var]] + shift)
      
      sel3 <- select_fixed_effects(paste0(outcome_var, "_log"), data)
      model3 <- sel3$best_model
      bp3 <- check_heteroscedasticity(model3)
      p_bp3 <- as.numeric(bp3)
      homosc3 <- p_bp3 > 0.05
      
      cat("\n--- Assumption checks (LOG transform) ---\n")
      cat("Breusch-Pagan p =", p_bp3, "| Homoscedastic:", homosc3, "\n")
      normal3 <- check_normality_ad(model3)
      cat("→ Inspect plots before accepting or overriding the automated decision.\n")
      
      plot_assumptions(model3, "LOG transform")
      
      # Override: return LOG model immediately if requested
      if (!is.null(override) && override == "log") {
        cat("\n→ OVERRIDE applied: retaining LOG model as instructed.\n")
        cat("→ Final model transform: log\n\n")
        return(list(model = model3, transform = "log", data = data))
      }
      
      if (normal3 & homosc3) {
        model <- model3
        transform_label <- "log"
      } else {
        message("  → LOG failed. Falling back to Gamma GLMM...")
        transform_label <- "gamma"
        
        gamma_sel <- select_fixed_effects_gamma(outcome_var, data)
        
        if (!is.null(gamma_sel)) {
          model <- gamma_sel$best_model
        } else {
          glmm_formula <- as.formula(paste(outcome_var, "~ treatment*preorpost + (1|id)"))
          dat_pos <- data %>% filter(.data[[outcome_var]] > 0)
          model <- tryCatch(
            glmmTMB(glmm_formula, data = dat_pos, family = Gamma(link = "log")),
            error = function(e) {
              message("Gamma GLMM failed: ", e$message)
              NULL
            }
          )
        }
      }
    }
  }
  
  # Override: return Gamma GLMM immediately if requested
  if (!is.null(override) && override == "gamma") {
    cat("\n→ OVERRIDE applied: fitting Gamma GLMM as instructed.\n")
    gamma_sel <- select_fixed_effects_gamma(outcome_var, data)
    if (!is.null(gamma_sel)) {
      model <- gamma_sel$best_model
    } else {
      glmm_formula <- as.formula(paste(outcome_var, "~ treatment*preorpost + (1|id)"))
      dat_pos <- data %>% filter(.data[[outcome_var]] > 0)
      model <- tryCatch(
        glmmTMB(glmm_formula, data = dat_pos, family = Gamma(link = "log")),
        error = function(e) { message("Gamma GLMM failed: ", e$message); NULL }
      )
    }
    transform_label <- "gamma"
  }

  # Collinearity check (manual inspection)
  # Runs on the final selected model, regardless of transform or override.
  # check_model() produces a full diagnostic panel including VIF.
  # Review the VIF plot: values < 5 are acceptable; values > 10 indicate a serious collinearity problem that would need to be addressed before interpreting fixed effects
  if (!is.null(model) && !inherits(model, "glmmTMB")) {
    cat("\n--- Collinearity check (check_model panel) —", outcome_var, "---\n")
    cat("Review the VIF panel: values < 5 are acceptable.\n")
    dev.new(width = 12, height = 10)
    print(check_model(model))
  } else if (inherits(model, "glmmTMB")) {
    cat("\n--- Collinearity check (VIF only — Gamma GLMM) —", outcome_var, "---\n")
    print(performance::check_collinearity(model))
  }
  
  cat("\n→ Final model transform:", transform_label, "\n\n")
  list(model = model, transform = transform_label, data = data)
}


# EXTRACT AND EXPORT RESULTS TABLE
export_results <- function(outcome_var, model, transform, data, prefix) {
  
  is_glmm <- inherits(model, "glmmTMB")
  
  # Fixed effects table
  coef_tbl <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(outcome = outcome_var, transform = transform) %>%
    select(outcome, transform, everything())
  
  readr::write_csv(coef_tbl, paste0(prefix, "_fixed_effects.csv"))
  
  # Type-III ANOVA
  if (is_glmm) {
    anova_tbl <- car::Anova(model, type = "III") %>%
      as.data.frame() %>%
      rownames_to_column("Effect") %>%
      mutate(outcome = outcome_var, transform = transform)
  } else {
    anova_tbl <- car::Anova(model, type = "III") %>%
      as.data.frame() %>%
      rownames_to_column("Effect") %>%
      rename_with(~ sub("^Pr\\(>F\\)$", "p.value", .x)) %>%
      rename_with(~ sub("^Pr\\(>Chisq\\)$", "p.value", .x)) %>%
      mutate(outcome = outcome_var, transform = transform)
  }
  
  readr::write_csv(anova_tbl, paste0(prefix, "_typeIII.csv"))
  cat("\n--- Type III ANOVA:", outcome_var, "---\n")
  print(anova_tbl)
  
  # Difference-in-differences (DID) contrast (primary estimand)
  # DID = (B_post − B_pre) − (A_post − A_pre)
  # emmeans cell order: A_post, A_pre, B_post, B_pre (alphabetical by treatment then preorpost)
  # Verify order before computing contrast:
  emm <- emmeans(model, ~ treatment * preorpost)
  cat("\nemmeans cell order check:\n")
  print(emm)
  
  did_tbl <- tryCatch({
    did_raw <- contrast(emm, list(DID = c(1, -1, -1, 1))) %>%
      summary(infer = TRUE, adjust = "none") %>%
      as.data.frame()
      # Gamma GLMM emmeans use asymp.LCL/asymp.UCL; lmer uses lower.CL/upper.CL
      # Rename to a consistent pair so transmute works for both model types
    if ("asymp.LCL" %in% names(did_raw)) {
      did_raw <- did_raw %>%
        rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL)
    }
    
    did_raw %>%
      transmute(
        outcome   = outcome_var,
        transform = transform,
        contrast  = "DID: (B_post−B_pre) − (A_post−A_pre)",
        estimate, SE, df, lower.CL, upper.CL, p.value
      )
  }, error = function(e) {
    message("DID contrast failed: ", e$message)
    NULL
  })
  
  if (!is.null(did_tbl)) {
    readr::write_csv(did_tbl, paste0(prefix, "_DID.csv"))
    cat("\n--- DID:", outcome_var, "---\n")
    print(did_tbl)
  }
  
  # Pairwise marginal means (provides descriptive context)
  emm_summary <- emmeans(model, ~ treatment * preorpost) %>%
    as.data.frame() %>%
    mutate(outcome = outcome_var, transform = transform)
  
  readr::write_csv(emm_summary, paste0(prefix, "_emmeans.csv"))
  
  invisible(list(coef = coef_tbl, anova = anova_tbl, did = did_tbl, emm = emm_summary))
}

# CHECKING FOR INFLUENTIAL CASES (only for lmer models)
# influence.ME is incompatible with dplyr being loaded (dplyr::filter masks stats::filter which influence.ME uses internally), and influence.ME itself imports dplyr so a clean subprocess does not resolve this either
# This function replicates the Cook's distance calculation manually:neach participant is removed in turn, the model is refit, and the change in fixed effect estimates is used to compute a Cook's distance analogue.
# Results are equivalent to influence.ME for this use case.
check_influential <- function(model, outcome_var, data) {
  if (inherits(model, "glmmTMB")) {
    message("Influential case check skipped for Gamma GLMM.")
    return(invisible(NULL))
  }
  
  explicit_formula <- as.formula(paste(outcome_var, "~ treatment * preorpost + (1 | id)"))
  
  # Fit full model explicitly (avoids stored formula issue from select_fixed_effects)
  m_full <- tryCatch(
    lmer(explicit_formula, data = data, REML = FALSE),
    error = function(e) { message("Full model refit failed: ", e$message); NULL }
  )
  if (is.null(m_full)) return(invisible(NULL))
  
  beta_full <- fixef(m_full)
  vcov_full <- as.matrix(vcov(m_full))
  ids       <- unique(data$id[!is.na(data[[outcome_var]])])
  p         <- length(beta_full)
  
  cd <- sapply(ids, function(i) {
    dat_i <- data[data$id != i, ]
    m_i   <- tryCatch(
      lmer(explicit_formula, data = dat_i, REML = FALSE),
      error = function(e) NULL
    )
    if (is.null(m_i)) return(NA_real_)
    diff  <- fixef(m_i) - beta_full
    # Cook's distance: scaled squared Mahalanobis distance of coefficient change
    tryCatch(
      as.numeric(t(diff) %*% solve(vcov_full) %*% diff) / p,
      error = function(e) NA_real_
    )
  })
  
  names(cd) <- as.character(ids)
  cd  <- cd[!is.na(cd)]
  cut <- 4 / length(cd)
  
  cat("\n--- Influential cases (Cook's D > 4/n =", round(cut, 3), ") for", outcome_var, "---\n")
  flagged <- cd[cd > cut]
  if (length(flagged) == 0) {
    cat("  None flagged.\n")
  } else {
    print(sort(flagged, decreasing = TRUE))
  }
  
  invisible(list(cooks_d = cd, cutoff = cut, flagged = flagged))
}

# LOIO ANALYSIS (drops top Cook's-D case if flagged) - only for lmer models
loio_sensitivity <- function(outcome_var, influential_result, data, prefix) {
  if (is.null(influential_result)) return(invisible(NULL))
  
  flagged_ids <- names(influential_result$flagged)
  if (length(flagged_ids) == 0) {
    cat("No influential cases to drop for LOIO analysis.\n")
    return(invisible(NULL))
  }
  
  top_id <- names(sort(influential_result$flagged, decreasing = TRUE))[1]
  cat("\nLOIO: dropping", top_id, "for", outcome_var, "\n")
  
  dat_drop <- filter(data, id != top_id)
  
  sel_drop  <- select_fixed_effects(outcome_var, dat_drop)
  m_drop    <- sel_drop$best_model
  
  # Helper function to extract DID
  did_from  <- function(fit) { 
    emm <- emmeans(fit, ~ treatment * preorpost)
    contrast(emm, list(DID = c(1, -1, -1, 1))) %>%
      summary(infer = TRUE, adjust = "none") %>%
      as.data.frame() %>%
      transmute(estimate, SE, df, lower.CL, upper.CL, p.value)
  }
  # Compute DID for full and reduced data
  m_all_fit <- select_fixed_effects(outcome_var, data)$best_model
  
  did_all  <- did_from(m_all_fit)  %>% mutate(spec = "All participants")
  did_drop <- did_from(m_drop)     %>% mutate(spec = paste0("Drop ", top_id))
  
  # Compare estimates
  did_comp <- bind_rows(did_all, did_drop) %>%
    mutate(outcome = outcome_var, delta_est = estimate - first(estimate)) %>%
    select(outcome, spec, estimate, SE, df, lower.CL, upper.CL, p.value, delta_est)
  
  readr::write_csv(did_comp, paste0(prefix, "_LOIO.csv"))
  cat("\n--- LOIO comparison:", outcome_var, "---\n")
  print(did_comp)
  
  invisible(did_comp)
}

# DHARMA DIAGNOSTICS (for lmer models only)
dharma_diagnostics <- function(outcome_var, data) {
  mf  <- model.frame(
    as.formula(paste(outcome_var, "~ treatment * preorpost")),
    data = data, na.action = na.omit
  )
  dat_cc <- data[as.integer(rownames(mf)), ]
  
  m_cc <- tryCatch(
    lmer(as.formula(paste(outcome_var, "~ treatment * preorpost + (1|id)")),
         data = dat_cc, REML = FALSE, na.action = na.omit),
    error = function(e) NULL
  )
  if (is.null(m_cc)) { message("DHARMa model failed to fit."); return(invisible(NULL)) }
  
  set.seed(123)
  simres <- simulateResiduals(fittedModel = m_cc, n = 1000, refit = FALSE)
  plot(simres, main = paste("DHARMa —", outcome_var))
  testUniformity(simres)
  testOutliers(simres)
  
  invisible(simres)
}

# -------------------------------------------------------------------------------------------------
# RUN FULL PIPELINE FOR EACH OUTCOME
# ------------------------------------------------------------------------------------------------
# NOTE: results for mean_SPT, mean_sleep_duration and sleep_efficiency are overwritten in the manual checks section below

outcomes_to_run <- list(
  list(var = "mean_SPT", label = "Sleep Period Time (h)"),
  list(var = "mean_sleep_duration", label = "Sleep Duration (h)"),
  list(var = "sleep_efficiency", label = "Sleep Efficiency (%)"),
  list(var = "mean_WASO", label = "Wake After Sleep Onset (h)")
)

# Store all results for a combined summary at the end
all_did_results <- list()
all_anova_results <- list()

for (o in outcomes_to_run) {
  v   <- o$var
  lbl <- o$label
  pfx <- paste0("results_", v)   # prefix for output CSV files
  
  cat("\n\n##############################################################\n")
  cat("OUTCOME:", lbl, "\n")
  cat("##############################################################\n")
  
  # Histogram of raw data
  hist(actigraphy_clean[[v]], breaks = 20, col = "grey85",
       main = paste("Raw distribution —", lbl), xlab = lbl)
  
  # Assumption checking & model selection with transform fallback
  # By default, the function decides the transform automatically based on Shapiro-Wilk and Breusch-Pagan p-values. Review the plots it generates before accepting its decision.
  # To override for a specific outcome after reviewing the plots, stop the loop and re-run check_and_transform() manually with the override argument
    # e.g.: result <- check_and_transform("mean_SPT", actigraphy_clean, override = "raw")
    # Then pass result$model and result$transform to export_results() manually.
    # Valid override values: "raw", "sqrt", "log", "gamma"
  result <- check_and_transform(v, actigraphy_clean)
  model  <- result$model
  transform  <- result$transform
  data_used  <- result$data   # may include derived transform columns
  
  if (is.null(model)) {
    warning("Model fitting failed for outcome: ", v, ". Skipping.")
    next
  }

  # Model summary 
  cat("\n--- Model summary:", v, "---\n")
  print(summary(model))
  
  # Residual plots
  if (!inherits(model, "glmmTMB")) {
    par(mfrow = c(1, 2))
    hist(resid(model), main = paste("Residuals —", v), xlab = "Residual")
    plot(fitted(model), resid(model), main = paste("Fitted vs Residuals —", v))
    abline(h = 0, lty = 2)
    par(mfrow = c(1, 1))
    
    # DHARMa on raw outcome variable (complete cases)
    dharma_diagnostics(v, actigraphy_clean)
  }
  
  # Export results
  res <- export_results(v, model, transform, data_used, pfx)
  
  # Accumulate DID and ANOVA results
  if (!is.null(res$did))   all_did_results[[v]]   <- res$did
  if (!is.null(res$anova)) all_anova_results[[v]] <- res$anova
  
  # Influential case check
  if (!inherits(model, "glmmTMB")) {
    infl_res <- check_influential(model, v, data_used)
    
    # LOIO sensitivity
    loio_sensitivity(v, infl_res, actigraphy_clean, pfx)
  }
}

# ------------------------------------------------------------------------------------------
# MANUAL MODEL CHECKS AND OVERRIDE
# ------------------------------------------------------------------------------------------
# MEAN_SPT

# Manually review plots
result_SPT <- check_and_transform("mean_SPT", actigraphy_clean)
# Gamma GLMM selected but raw model appears mostly normal on visual inspection. Slight left skew on histogram and tail on QQ plot
# Anderson-Darling and Shapiro-Wilk test significant but can be overly sensitive
# Homoscedastic on Breusch-Pagan test and visual inspection

# Assess collinearity
result_SPT_raw <- check_and_transform("mean_SPT", actigraphy_clean, override = "raw")
dev.new(width = 12, height = 10)
check_model(result_SPT_raw$model)
# Nil issues of collinearity in raw model (VIF < 5)
# Mild nonlinearity

# Decision made to override Gamma GLMM with raw model
v_override         <- "mean_SPT"
lbl_override       <- "Sleep Period Time (h)"
pfx_override       <- paste0("results_", v_override)

result_override    <- check_and_transform(v_override, actigraphy_clean, override = "raw")
model_override     <- result_override$model
transform_override <- result_override$transform
data_override      <- result_override$data

print(summary(model_override))

if (!inherits(model_override, "glmmTMB")) {
  par(mfrow = c(1, 2))
  hist(resid(model_override), main = paste("Residuals —", v_override), xlab = "Residual")
  plot(fitted(model_override), resid(model_override),
       main = paste("Fitted vs Residuals —", v_override))
  abline(h = 0, lty = 2)
  par(mfrow = c(1, 1))
  dharma_diagnostics(v_override, actigraphy_clean)
}

res_override <- export_results(v_override, model_override, transform_override,
                               data_override, pfx_override)

# Add to accumulator lists (overwrites any earlier automated entry for this outcome)
all_did_results[[v_override]]   <- res_override$did
all_anova_results[[v_override]] <- res_override$anova

# Influential case check and LOIO
# check_influential() now handles the refit internally using data_override,
# so no separate refit step is needed here.
if (!inherits(model_override, "glmmTMB")) {
  infl_res_override <- check_influential(model_override, v_override, data_override)
  loio_sensitivity(v_override, infl_res_override, actigraphy_clean, pfx_override)
}
# R06 was flagged with a Cook's distance of 0.179 (>4/n)
# LOIO analysis: slight increase in DID estimate when R06 dropped but CIs overlap and remain non-significant, with CIs crossing zero
# Thus, the results are robust to exclusion of R06

#----------------------------------------------------------
# MEAN_SLEEP_DURATION

# Manually review plots
result_duration <- check_and_transform("mean_sleep_duration", actigraphy_clean)
# Gamma GLMM selected but raw model appears mostly normal on visual inspection. Slight left skew on histogram and tail on QQ plot
# Anderson-Darling and Shapiro-Wilk test significant
# Homoscedastic on Breusch-Pagan test and visual inspection

# Assess collinearity
result_duration_raw <- check_and_transform("mean_sleep_duration", actigraphy_clean, override = "raw")
dev.new(width = 12, height = 10)
check_model(result_duration_raw$model)
# Nil issues of collinearity in raw model (VIF < 5)

# Decision made to override Gamma GLMM with raw model
v_override         <- "mean_sleep_duration"
lbl_override       <- "Sleep Duration (h)"
pfx_override       <- paste0("results_", v_override)

result_override    <- check_and_transform(v_override, actigraphy_clean, override = "raw")
model_override     <- result_override$model
transform_override <- result_override$transform
data_override      <- result_override$data

print(summary(model_override))

if (!inherits(model_override, "glmmTMB")) {
  par(mfrow = c(1, 2))
  hist(resid(model_override), main = paste("Residuals —", v_override), xlab = "Residual")
  plot(fitted(model_override), resid(model_override),
       main = paste("Fitted vs Residuals —", v_override))
  abline(h = 0, lty = 2)
  par(mfrow = c(1, 1))
  dharma_diagnostics(v_override, actigraphy_clean)
}

res_override <- export_results(v_override, model_override, transform_override,
                               data_override, pfx_override)

# Add to accumulator lists (overwrites any earlier automated entry for this outcome)
all_did_results[[v_override]]   <- res_override$did
all_anova_results[[v_override]] <- res_override$anova

# Influential case check and LOIO
# check_influential() now handles the refit internally using data_override,
# so no separate refit step is needed here.
if (!inherits(model_override, "glmmTMB")) {
  infl_res_override <- check_influential(model_override, v_override, data_override)
  loio_sensitivity(v_override, infl_res_override, actigraphy_clean, pfx_override)
}
# R06 was flagged with a Cook's distance of 0.126 (>4/n)
# LOIO analysis: slight increase in DID estimate when R06 dropped but CIs overlap and remain non-significant, with CIs crossing zero
# Thus, the results are robust to exclusion of R06

#---------------------------------------------------
# SLEEP EFFICIENCY
# Manually review plots
result_efficiency <- check_and_transform("sleep_efficiency", actigraphy_clean)
# Raw, SQRT and LOG models all failed formal normality testing - Gamma GLMM was automatically selected based on this
# Raw model appears largely normal on histogram and has mild heavy-tailedness on QQ plot. No clear evidence of heteroscedasticity

# Decision made to run both models and compare key inferential outputs (DID estimate, SE, CI, p-value, Type-III interaction term) between both models
# Decision rule: 
  # If the DID estimates are in the same direction, CIs overlap substantially and both models agree on statistical significance/non-significance => the material are materially consistent, use the raw LMM
  # If the models diverge in direction, significance or CI overlap is poor, retain Gamma GLMM

# Step 1: Fit both models
# Raw LMM (override the automated Gamma selection)
result_eff_raw   <- check_and_transform("sleep_efficiency", actigraphy_clean, override = "raw")
model_eff_raw    <- result_eff_raw$model
# Gamma GLMM (automated selection; re-run to ensure it is stored cleanly)
result_eff_gamma <- check_and_transform("sleep_efficiency", actigraphy_clean, override = "gamma")
model_eff_gamma  <- result_eff_gamma$model

# Step 2: Assumption diagnostics for the raw model

cat("\n--- Assumption checks: Sleep Efficiency RAW model ---\n")

# Normality of residuals (Anderson-Darling primary; Shapiro-Wilk supplementary)
r_raw <- resid(model_eff_raw)
ad_raw <- nortest::ad.test(r_raw[!is.na(r_raw)])
sw_raw <- shapiro.test(r_raw[!is.na(r_raw)])
cat("Anderson-Darling p =", round(ad_raw$p.value, 4),
    "| Shapiro-Wilk p =",   round(sw_raw$p.value, 4), "\n")

# Homoscedasticity
bp_raw   <- performance::check_heteroscedasticity(model_eff_raw)
p_bp_raw <- as.numeric(bp_raw)
cat("Breusch-Pagan p =", round(p_bp_raw, 4),
    "| Homoscedastic:", p_bp_raw > 0.05, "\n")

# Visual inspection — residual plots for raw model
dev.new(width = 10, height = 4)
par(mfrow = c(1, 3))
hist(r_raw, breaks = 20, col = "grey85",
     main = "Residuals — Sleep Efficiency\n(Raw LMM)", xlab = "Residual")
qqnorm(r_raw, main = "QQ plot — Sleep Efficiency\n(Raw LMM)")
qqline(r_raw, col = "red", lty = 2)
plot(fitted(model_eff_raw), r_raw,
     main = "Fitted vs Residuals — Sleep Efficiency\n(Raw LMM)",
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, lty = 2, col = "red")
par(mfrow = c(1, 1))

# check_model panel for collinearity (raw LMM)
dev.new(width = 12, height = 10)
print(check_model(model_eff_raw))
cat("Review VIF panel: values < 5 acceptable.\n")
# No problematic collinearity

# DHARMa diagnostics for raw model (complete cases)
dharma_diagnostics("sleep_efficiency", actigraphy_clean)
# No outliers flagged; KS, dispersion and Levene test non-significant

# DHARMa diagnostics for Gamma GLMM
set.seed(123)
simres_eff_gamma <- simulateResiduals(fittedModel = model_eff_gamma, n = 1000)
plot(simres_eff_gamma, main = "DHARMa diagnostics — Sleep Efficiency (Gamma GLMM)")
testUniformity(simres_eff_gamma)
testOutliers(simres_eff_gamma)
testDispersion(simres_eff_gamma)
# The Gamma GLMM appears to fit sleep_efficiency well
# One observation is flagged in red but is not statistically significant (p = 0.196)
# No evidence of over/underdispersion

# Step 3: Extract inferential outputs from both models 
# Helper: extract DID contrast from any supported model type
extract_did <- function(mod, model_label) {
  emm <- emmeans(mod, ~ treatment * preorpost)
  
  did_raw <- tryCatch(
    contrast(emm, list(DID = c(1, -1, -1, 1))) %>%
      summary(infer = TRUE, adjust = "none") %>%
      as.data.frame(),
    error = function(e) { message("DID contrast failed: ", e$message); NULL }
  )
  if (is.null(did_raw)) return(NULL)
  
  # Gamma GLMM emmeans use asymp.LCL/UCL; lmer uses lower.CL/upper.CL
  if ("asymp.LCL" %in% names(did_raw)) {
    did_raw <- did_raw %>% rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL)
  }
  
  did_raw %>%
    transmute(
      model = model_label,
      contrast = "DID: (B_post−B_pre) − (A_post−A_pre)",
      estimate, SE, df, lower.CL, upper.CL, p.value
    )
}

did_raw_tbl <- extract_did(model_eff_raw,   "Raw LMM")
did_gamma_tbl <- extract_did(model_eff_gamma, "Gamma GLMM")

# Helper: extract the treatment × preorpost interaction row from Type-III ANOVA
extract_interaction <- function(mod, model_label) {
  car::Anova(mod, type = "III") %>%
    as.data.frame() %>%
    rownames_to_column("Effect") %>%
    rename_with(~ sub("^Pr\\(>F\\)$",    "p.value", .x)) %>%
    rename_with(~ sub("^Pr\\(>Chisq\\)$","p.value", .x)) %>%
    filter(grepl("treatment:preorpost|treatment\\.preorpost", Effect, ignore.case = TRUE)) %>%
    mutate(model = model_label) %>%
    select(model, Effect, everything())
}

anova_raw_tbl <- extract_interaction(model_eff_raw,   "Raw LMM")
anova_gamma_tbl <- extract_interaction(model_eff_gamma, "Gamma GLMM")

# Step 4: Side-by-side comparison tables

cat("\n\n=======================================================\n")
cat("SLEEP EFFICIENCY — DID CONTRAST COMPARISON\n")
cat("=======================================================\n")
did_comparison <- bind_rows(did_raw_tbl, did_gamma_tbl)
print(did_comparison)

cat("\n=======================================================\n")
cat("SLEEP EFFICIENCY — TYPE-III INTERACTION COMPARISON\n")
cat("=======================================================\n")
anova_comparison <- bind_rows(anova_raw_tbl, anova_gamma_tbl)
print(anova_comparison)

# Export both comparison tables
readr::write_csv(did_comparison,   "sleep_efficiency_DID_raw_vs_gamma.csv")
readr::write_csv(anova_comparison, "sleep_efficiency_typeIII_raw_vs_gamma.csv")
cat("\nComparison tables written to:\n",
    "  sleep_efficiency_DID_raw_vs_gamma.csv\n",
    "  sleep_efficiency_typeIII_raw_vs_gamma.csv\n")

# Step 5: Decision and pipeline entry
# Review the comparison tables above and make decision based on decision rule
  # Both models show a trivial DID (raw = 0.53, gamma = 0.007 on log scale) with CIs overlapping and crossing zero. p-values are similar and non-significant (raw = 0.799, gamma = 0.766)
  # Near-identical conclusions in Type III interaction (Chi-sq: raw = 0.068, gamma = 0.088)


# NOTE: would have changed to "gamma" if models diverged and the Gamma GLMM was to be retained
eff_model_choice <- "raw" 

if (eff_model_choice == "raw") {
  cat("\n→ DECISION: Raw LMM selected for sleep efficiency (consistent with Gamma GLMM).\n")
  model_eff_final <- model_eff_raw
  transform_eff_final <- "raw"
  data_eff_final <- result_eff_raw$data
} else {
  cat("\n→ DECISION: Gamma GLMM retained for sleep efficiency (models diverged).\n")
  model_eff_final <- model_eff_gamma
  transform_eff_final <- "gamma"
  data_eff_final <- result_eff_gamma$data
}

# Export results for the chosen model and add to accumulator lists
v_eff <- "sleep_efficiency"
pfx_eff <- paste0("results_", v_eff)

cat("\n--- Model summary: sleep_efficiency (", transform_eff_final, ") ---\n")
print(summary(model_eff_final))

if (!inherits(model_eff_final, "glmmTMB")) {
  par(mfrow = c(1, 2))
  hist(resid(model_eff_final), main = paste("Residuals — sleep_efficiency"), xlab = "Residual")
  plot(fitted(model_eff_final), resid(model_eff_final),
       main = "Fitted vs Residuals — sleep_efficiency")
  abline(h = 0, lty = 2)
  par(mfrow = c(1, 1))
}

res_eff <- export_results(v_eff, model_eff_final, transform_eff_final,
                          data_eff_final, pfx_eff)

# Overwrite (or create) the accumulator entries for this outcome
all_did_results[[v_eff]] <- res_eff$did
all_anova_results[[v_eff]] <- res_eff$anova

# Influential case check and LOIO (lmer only)
if (!inherits(model_eff_final, "glmmTMB")) {
  infl_res_eff <- check_influential(model_eff_final, v_eff, data_eff_final)
  loio_sensitivity(v_eff, infl_res_eff, actigraphy_clean, pfx_eff)
}

# R36 was flagged as the most influential case with a Cook's distance of 0.154 (>4/n)
# LOIO analysis: reversal in direction of DID estimate when R36 dropped but CIs overlap and cross zero, p-values remain non-significant
# Thus, the results are robust to exclusion of R36

#----------------------------------------
# MEAN WASO
# Manually review plots
result_WASO <- check_and_transform("mean_WASO", actigraphy_clean)
# SQRT transformation model appears to be the appropriate model for mean WASO. Appears normal and homoscedastic on visual inspection and passed normality and homoscedasticity testing
# Nil evidence of problematic collinearity (VIF < 5)
# AUtomatic SQRT model accepted - no override required

# R06 was flagged as the most influential case with a Cook's distance of 0.292 (>4/n)
# LOIO analysis: slight increase in DID estimate when R06 dropped but CIs overlap and remain non-significant, with CIs crossing zero
# Thus, the results are robust to exclusion of R06


# ------------------------------------------------------------------------------------------
# GENERATE COMBINED SUMMARY TABLES
# -------------------------------------------------------------------------------------------

# Combined DID results across all outcomes
did_combined <- bind_rows(all_did_results)
readr::write_csv(did_combined, "SUMMARY_DID_all_outcomes.csv")

cat("\n\n=== COMBINED DID SUMMARY (all outcomes) ===\n")
print(did_combined %>% select(outcome, transform, estimate, SE, lower.CL, upper.CL, p.value))

# Combined Type-III ANOVA — interaction row only (most relevant - "is the pre-to-post change different between CBD and placebo?")
anova_combined <- bind_rows(all_anova_results) %>%
  filter(grepl("treatment:preorpost|treatment.preorpost", Effect, ignore.case = TRUE))

readr::write_csv(anova_combined, "SUMMARY_typeIII_interaction_all_outcomes.csv")

cat("\n=== TYPE-III INTERACTION TERMS (all outcomes) ===\n")
print(anova_combined %>% select(outcome, transform, Effect, everything()))


# ----------------------------------------------------------------------------------------
# GENERATE VISUALISATIONS
# ----------------------------------------------------------------------------------------

# Spaghetti plots for all outcomes (pre/post by treatment)
for (o in outcomes_to_run) {
  v   <- o$var
  lbl <- o$label
  
  p <- ggplot(actigraphy_clean, aes(x = preorpost, y = .data[[v]],
                                           group = id, colour = id)) +
    geom_line(alpha = 0.5) +
    geom_point(size = 1.8) +
    facet_wrap(~ treatment, labeller = labeller(treatment = c("A" = "Treatment A", "B" = "Treatment B"))) +
    labs(title = paste("Individual trajectories —", lbl),
         x = "Period (pre/post)", y = lbl) +
    theme_bw() +
    theme(legend.position = "none")
  
  print(p)
  ggsave(paste0("spaghetti_", v, ".png"), plot = p, width = 8, height = 5, dpi = 150)
}

# Mean ± SE summary plot for all outcomes
for (o in outcomes_to_run) {
  v   <- o$var
  lbl <- o$label
  
  summ <- actigraphy_clean %>%
    group_by(treatment, preorpost) %>%
    summarise(
      mean_val = mean(.data[[v]], na.rm = TRUE),
      se_val   = sd(.data[[v]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[v]]))),
      .groups  = "drop"
    )
  
  p <- ggplot(summ, aes(x = preorpost, y = mean_val,
                        colour = treatment, group = treatment)) +
    geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                  width = 0.1, position = position_dodge(0.2)) +
    geom_line(position = position_dodge(0.2), linewidth = 0.9) +
    geom_point(size = 3, position = position_dodge(0.2)) +
    scale_colour_manual(values = c("A" = "#2196F3", "B" = "#F44336"),
                        labels = c("A" = "Treatment A", "B" = "Treatment B")) +
    labs(title = paste("Mean ± SE —", lbl),
         x = "Period (pre/post)", y = lbl, colour = "") +
    theme_bw()
  
  print(p)
  ggsave(paste0("mean_SE_", v, ".png"), plot = p, width = 6, height = 4, dpi = 150)
}

cat("\n\nAnalysis complete. All CSV and PNG files written to the working directory.\n")






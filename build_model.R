###############################################################################
## Delayed graft function (DGF) after DCD kidney transplant -- main analysis
##
## Pipeline (run in this order):
##   1. build_model.R  <- THIS FILE
##        preprocessing, 5-fold CV (RF + XGBoost), final XGBoost fit,
##        isotonic calibration. Writes the cached inputs every other
##        script depends on:
##          data/df1_DCD_vital.csv            (analysis dataset)
##          data/plot_objects.RData           (model + matrices for figures)
##          data/center_calibration_inputs.RData
##   2. compute_feature_importance.R  -> data/xgb_feature_importance_top30.csv
##      compute_pdp.R                 -> data/pdp_12_features_data.csv
##      compute_permutation_importance.R -> permutation_importance.csv
##   3. plot_feature_importance.R / plot_pdp.R / plot_permutation_importance.R
##        cheap, read the CSVs above, never touch the model -- rerun freely.
##
## Inputs (not distributed with this code; SRTR/OPTN data use agreement
## required):
##   KIDPAN_DATA_CLEAN_MERGED.dta        recipient/donor/transplant records
##   DECEASED_DONOR_DCD_MEASURES.DTA     donor hemodynamic measurements
##   Diagnosis Categorization.xlsx       DIAG_KI code -> category map
##
## Outcome: FIRST_WK_DIAL (dialysis in the first week post-transplant).
## Cohort:  DCD donors (identified by presence of hemodynamic records),
##          adult donor and adult recipient, excluding donors who
##          progressed to brain death.
###############################################################################

## ---- paths -----------------------------------------------------------------
## Defaults to the repository root, so running this script from the repo
## folder needs no editing. Override it if the data lives elsewhere.
## PROJECT_DIR/data must contain the three input files listed above; all
## intermediate CSV/RData files are written there too.
PROJECT_DIR <- "."

PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)
DATA_DIR <- normalizePath(file.path(PROJECT_DIR, "data"), mustWork = TRUE)
OUT_DIR  <- PROJECT_DIR
setwd(DATA_DIR)

library(haven)
library(tidyverse)
library(recipes)
library(ranger)
library(pROC)
library(readxl)
library(ResourceSelection)   # hoslem.test()
library(xgboost)
library(ggpubr)              # ggarrange()

# load full dataset
data = read_dta("KIDPAN_DATA_CLEAN_MERGED.dta")

# load longitudinal donor hemodynamic dataset
df_long1 = read_dta("DECEASED_DONOR_DCD_MEASURES.DTA")

#####################################################
############## helper function ######################
#####################################################
prepare_matrix_data <- function(train_df, test_df = NULL, calib_df = NULL,
                                outcome = "FIRST_WK_DIAL") {
  train_df <- droplevels(train_df)
  if (!is.null(test_df))  test_df  <- droplevels(test_df)
  if (!is.null(calib_df)) calib_df <- droplevels(calib_df)
  
  # factor columns in training
  factor_cols <- names(train_df)[sapply(train_df, is.factor)]
  
  # drop factor columns that have < 2 levels in training
  bad_factor_cols <- factor_cols[
    sapply(train_df[factor_cols], function(x) nlevels(droplevels(x)) < 2)
  ]
  
  if (length(bad_factor_cols) > 0) {
    cat("\nDropping factor columns with <2 levels in training:\n")
    print(bad_factor_cols)
  }
  
  train_df <- train_df %>% select(-any_of(bad_factor_cols))
  if (!is.null(test_df))  test_df  <- test_df  %>% select(-any_of(bad_factor_cols))
  if (!is.null(calib_df)) calib_df <- calib_df %>% select(-any_of(bad_factor_cols))
  
  # align factor levels in test/calib to training
  factor_cols2 <- names(train_df)[sapply(train_df, is.factor)]
  for (col in factor_cols2) {
    # add explicit missing level to training
    train_chr <- as.character(train_df[[col]])
    train_chr[is.na(train_chr)] <- "(Missing)"
    if (col != outcome) {
      train_df[[col]] <- factor(train_chr,
                                 levels = sort(union(unique(train_chr), "(Missing)")))
    } else {
      train_df[[col]] <- factor(train_chr)
    }
    
    train_levels <- levels(train_df[[col]])
    
    if (!is.null(test_df)) {
      test_chr <- as.character(test_df[[col]])
      test_chr[is.na(test_chr)] <- "(Missing)"
      test_chr[!(test_chr %in% train_levels)] <- "(Missing)"
      test_df[[col]] <- factor(test_chr, levels = train_levels)
    }
    
    if (!is.null(calib_df)) {
      calib_chr <- as.character(calib_df[[col]])
      calib_chr[is.na(calib_chr)] <- "(Missing)"
      calib_chr[!(calib_chr %in% train_levels)] <- "(Missing)"
      calib_df[[col]] <- factor(calib_chr, levels = train_levels)
    }
  }
  
  form <- as.formula(paste(outcome, "~ . -1"))
  
  # preserve all rows; do not silently drop rows with NA
  train_mat <- model.matrix(form, data = train_df, na.action = na.pass)
  test_mat  <- if (!is.null(test_df))  model.matrix(form, data = test_df,  na.action = na.pass) else NULL
  calib_mat <- if (!is.null(calib_df)) model.matrix(form, data = calib_df, na.action = na.pass) else NULL
  
  # quick checks
  if (nrow(train_mat) != nrow(train_df)) {
    stop("train_mat row count does not match train_df row count")
  }
  if (!is.null(test_df) && nrow(test_mat) != nrow(test_df)) {
    stop("test_mat row count does not match test_df row count")
  }
  if (!is.null(calib_df) && nrow(calib_mat) != nrow(calib_df)) {
    stop("calib_mat row count does not match calib_df row count")
  }
  
  list(
    train_df = train_df,
    test_df = test_df,
    calib_df = calib_df,
    train_mat = train_mat,
    test_mat = test_mat,
    calib_mat = calib_mat
  )
}


#####################################################
############## data preprocess ######################
#####################################################

# select out the rows with all NA values
df_long_DCD_noinfo = df_long1 %>% 
  select(1:7) %>% mutate(across(where(is.numeric), ~na_if(., 0))) %>%
  filter(if_all(-c(DONOR_ID, MEASURE_DATE, MEASURE_TIME), is.na)) 


df_long2 = df_long1 %>% 
  select(1:7) %>% mutate(across(where(is.numeric), ~na_if(., 0))) %>%
  filter(!if_all(-c(DONOR_ID, MEASURE_DATE, MEASURE_TIME), is.na)) 

# find out the intersected DONOR ID
length(unique(data$DONOR_ID))
length(unique(df_long2$DONOR_ID))

length(intersect(unique(data$DONOR_ID), unique(df_long2$DONOR_ID)))

# only keep the donors who appear in main dataset
df_long3 = df_long2 %>%
  filter(DONOR_ID %in% intersect(unique(data$DONOR_ID), unique(df_long2$DONOR_ID)))

df_long4 = df_long3 %>% 
  mutate(
    MEASURE_TIME = as.POSIXct(
      paste(MEASURE_DATE, format(MEASURE_TIME, "%H:%M:%S")),
      format = "%Y-%m-%d %H:%M:%S"
    )
  ) %>% 
  select(-MEASURE_DATE) %>% 
  arrange(DONOR_ID, MEASURE_TIME)

# impute missing values(forward filling within each DONOR_ID)
df_long5 = df_long4 %>% 
  group_by(DONOR_ID) %>% 
  fill(-c(DONOR_ID, MEASURE_TIME), .direction = "down") %>%
  ungroup()

# return last non-NA
last_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_) else return(tail(x, 1))
}

safe_stat <- function(f) {
  function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_) else return(f(x))
  }
}

safe_range <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0) NA_real_ else max(v) - min(v)
}

safe_slope <- function(y, t) {
  if (inherits(t, "POSIXt")) t_num <- as.numeric(t) / 60 else
    if (inherits(t, "Date")) t_num <- as.numeric(t) else t_num <- as.numeric(t)
  ok <- is.finite(y) & is.finite(t_num)
  y2 <- y[ok]
  t2 <- t_num[ok]
  if (length(y2) < 2 || sd(t2) == 0) return(NA_real_)
  t2c <- t2 - mean(t2)
  coef(lm(y2 ~ t2c))[2]
}

# summary statistics: max, range, sd, last computed on forward-filled data (df_long5)
df_long_flat <- df_long5 %>%
  select(-MEAN_PRESSURE) %>%
  group_by(DONOR_ID) %>%
  summarise(
    across(
      -c(MEASURE_TIME),
      list(
        max   = safe_stat(max),
        range = ~safe_range(.),
        sd    = safe_stat(sd),
        last  = ~last_non_na(.)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# slopes computed on pre-forward-fill data (df_long4) to avoid bias from repeated imputed values
df_long_flat_slopes <- df_long4 %>%
  select(-MEAN_PRESSURE) %>%
  group_by(DONOR_ID) %>%
  summarise(
    across(
      -c(MEASURE_TIME),
      list(slope = ~safe_slope(., MEASURE_TIME)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

df_long_flat <- df_long_flat %>%
  left_join(df_long_flat_slopes, by = "DONOR_ID")

# winsorize slope columns at 1st/99th percentile to remove extreme outliers
# caused by donors with only 2 measurements very close in time
slope_cols_all <- grep("_slope$", names(df_long_flat), value = TRUE)
df_long_flat <- df_long_flat %>%
  mutate(across(
    all_of(slope_cols_all),
    ~{
      lo <- quantile(., 0.01, na.rm = TRUE)
      hi <- quantile(., 0.99, na.rm = TRUE)
      pmin(pmax(., lo), hi)
    }
  ))

dcd_donors_no_info <- setdiff(unique(df_long_DCD_noinfo$DONOR_ID), unique(df_long2$DONOR_ID))
dcd_donors_no_info_final <- intersect(dcd_donors_no_info, unique(data$DONOR_ID))

df_long_flat_add <- df_long_flat %>%
  bind_rows(tibble(DONOR_ID = dcd_donors_no_info_final))

# missing indicators: 1 if slope could not be computed (too few obs, all-NA vitals,
# or no monitoring record at all), 0 otherwise. Created BEFORE mean imputation so
# they are never NA when the model sees them.
df_long_flat_add <- df_long_flat_add %>%
  mutate(
    SYSTOLIC_PRESSURE_slope_missing = as.integer(is.na(SYSTOLIC_PRESSURE_slope)),
    O2_SATURATION_slope_missing     = as.integer(is.na(O2_SATURATION_slope))
  )

# NOTE: vitals NAs are intentionally left unimputed here. They are filled in
# later by the per-split recipe (step_impute_mean(all_numeric())), which is
# fit on train_df_DCD only and applied to calib/test via bake() -- fixing a
# leakage bug where they were previously mean-imputed over the full
# train+calib+test population before the split even happened.

### add the variable DIAG_KI ####
diag_KI_decoder <- read_excel("Diagnosis Categorization.xlsx") %>%
  rename(DIAG_KI = `Diagnosis Code`) %>% 
  select(DIAG_KI, Category)

data2 <- data %>%
  left_join(diag_KI_decoder, by = "DIAG_KI") %>%
  mutate(DIAG_KI = Category) %>%  
  select(-Category) %>%
  filter(DCD_PROGRESSED_BRAIN_DEATH != "Y")

df_subset2 = data2 %>%
  select(
    DONOR_ID, CTR_CODE, FIRST_WK_DIAL, GENDER, PREV_TX, DIABETES, INIT_CPRA,
    HLAMIS, BMI_CALC, DIALYSIS_DURATION,
    CREAT_TRR, NUM_PREV_TX, ABO, ABO_MAT, PREV_PREGNANCY,
    GENDER_DON, AGE_DON, BMI_DON_CALC, NON_HRT_DON,
    BUN_DON, SGPT_DON, ALCOHOL_HEAVY_DON,
    HIST_CIG_DON, DIABETES_DON, HIST_CAD, HIST_HYPERTENS_DON,
    ABO_DON, DCD_PROGRESSED_BRAIN_DEATH,
    AGE, COLD_ISCH_KI,
    CREAT_DON, COD_CAD_DON, PUMP_KI, DIAG_KI
  )

df = df_subset2 %>%  
  mutate(across(where(is.character), as.factor)) %>%
  mutate(
    PREV_PREGNANCY = as.factor(PREV_PREGNANCY),
    PUMP_KI = as.factor(PUMP_KI)
  ) %>%
  mutate(across(where(is.character) | where(is.factor), ~ {
    x <- trimws(as.character(.x))
    x[x == ""] <- NA
    factor(x)
  }))

df = df %>%
  mutate(
    COD_CAD_DON = na_if(COD_CAD_DON, 999)
  ) %>%
  mutate(
    DIABETES = case_when(
      DIABETES == 0 ~ "N",
      DIABETES == 1 ~ "Y",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(DIABETES = factor(DIABETES, levels = c("N", "Y")))

# Row-level filter: keep only transplant records where BOTH donor and
# recipient are adults. This drops a pediatric donor's row regardless of
# recipient age, and a pediatric recipient's row regardless of donor age,
# without discarding a donor's OTHER (adult-recipient) rows in the process.
# Also drops rows with unknown outcome (FIRST_WK_DIAL is NA) -- previously
# these survived into the model and were silently mode-imputed by the
# per-split recipe, since recipe(~ ., data = ...) gives the outcome column
# the generic "predictor" role with no outcome exclusion.
df1 = df %>% filter(AGE_DON >= 18, AGE >= 18, !is.na(FIRST_WK_DIAL))
df_vital = df_long_flat_add

df1_DCD = df1 %>% filter(DONOR_ID %in% unique(df_vital$DONOR_ID))
df1_DBD = df1 %>% filter(!DONOR_ID %in% unique(df_vital$DONOR_ID))

table(df1_DCD$FIRST_WK_DIAL)
table(df1_DBD$FIRST_WK_DIAL)

##########################################################
#### Analysis dataset: DCD cohort + hemodynamic features #
##########################################################
print("Random forest Modeling DCD group with vital signs added...")

df1_DCD_vital = df1_DCD %>% inner_join(df_long_flat_add, by = "DONOR_ID")

cat("DCD with vital signs: rows =", nrow(df1_DCD_vital),
    ", cols =", ncol(df1_DCD_vital), "\n")
write.csv(df1_DCD_vital, "df1_DCD_vital.csv", row.names = FALSE)

##########################################################
#### 5-fold cross-validation (RF + XGBoost) ##############
##########################################################

# 1) Re-derive same training donors as the final model (70/30, seed 333)
set.seed(333)
all_donors_cv   <- unique(df1_DCD_vital$DONOR_ID)
test_ids_cv     <- sample(all_donors_cv, size = round(0.3 * length(all_donors_cv)))
train_df_full   <- df1_DCD_vital %>% filter(!(DONOR_ID %in% test_ids_cv))

# 2) Assign fold IDs at the donor level
K <- 5
set.seed(333)
donor_fold_map <- train_df_full %>%
  distinct(DONOR_ID) %>%
  mutate(fold = sample(rep(1:K, length.out = n())))

write.csv(
  train_df_full %>%
    select(DONOR_ID, FIRST_WK_DIAL) %>%
    left_join(donor_fold_map, by = "DONOR_ID"),
  "cv_fold_assignments.csv",
  row.names = FALSE
)

compute_auprc <- function(scores, labels) {
  ord  <- order(scores, decreasing = TRUE)
  labs <- labels[ord]
  tp   <- cumsum(labs)
  fp   <- cumsum(1 - labs)
  precision <- c(1, tp / (tp + fp))
  recall    <- c(0, tp / sum(labs))
  sum(diff(recall) * (head(precision, -1) + tail(precision, -1)) / 2)
}

# 3) Per-fold fitting function
run_one_fold <- function(fold_k, train_df_full, donor_fold_map,
                         xgb_params, n_rounds = 500) {
  val_donors   <- donor_fold_map$DONOR_ID[donor_fold_map$fold == fold_k]
  fold_train   <- train_df_full %>% filter(!(DONOR_ID %in% val_donors))
  fold_val     <- train_df_full %>% filter(DONOR_ID %in% val_donors)

  # Imputation recipe fit on fold training only (no leakage)
  fold_rec <- recipe(~ ., data = fold_train) %>%
    step_impute_mean(all_numeric()) %>%
    step_impute_mode(all_nominal()) %>%
    prep()

  fold_train_imp <- bake(fold_rec, new_data = NULL)
  fold_val_imp   <- bake(fold_rec, new_data = fold_val)

  fold_train_mod <- fold_train_imp %>%
    select(-DONOR_ID, -CTR_CODE) %>%
    mutate(across(where(is.numeric), ~ifelse(is.nan(.), 0, .))) %>%
    mutate(FIRST_WK_DIAL = factor(FIRST_WK_DIAL, levels = c("N", "Y"),
                                  labels = c("No", "Yes"))) %>%
    mutate(across(where(is.character), as.factor))
  fold_val_mod <- fold_val_imp %>%
    select(-DONOR_ID, -CTR_CODE) %>%
    mutate(across(where(is.numeric), ~ifelse(is.nan(.), 0, .))) %>%
    mutate(FIRST_WK_DIAL = factor(FIRST_WK_DIAL, levels = c("N", "Y"),
                                  labels = c("No", "Yes"))) %>%
    mutate(across(where(is.character), as.factor))

  # Random forest
  rf_fold <- ranger(
    FIRST_WK_DIAL ~ .,
    data      = fold_train_mod,
    num.trees = 1000,
    probability = TRUE,
    seed      = 333
  )
  rf_probs <- predict(rf_fold, data = fold_val_mod)$predictions[, "Yes"]
  rf_auc   <- as.numeric(suppressMessages(
    pROC::auc(pROC::roc(fold_val_mod$FIRST_WK_DIAL, rf_probs))
  ))

  # XGBoost
  mm_fold <- prepare_matrix_data(fold_train_mod, fold_val_mod,
                                 outcome = "FIRST_WK_DIAL")
  label_tr  <- as.numeric(mm_fold$train_df$FIRST_WK_DIAL == "Yes")
  label_val <- as.numeric(mm_fold$test_df$FIRST_WK_DIAL  == "Yes")
  spw       <- sum(label_tr == 0) / sum(label_tr == 1)

  xgb_params_fold        <- xgb_params
  xgb_params_fold$scale_pos_weight <- spw

  dtrain_fold <- xgb.DMatrix(mm_fold$train_mat, label = label_tr)
  dval_fold   <- xgb.DMatrix(mm_fold$test_mat,  label = label_val)

  xgb_fold <- xgb.train(
    params        = xgb_params_fold,
    data          = dtrain_fold,
    nrounds       = n_rounds,
    verbose       = 0
  )
  xgb_probs   <- predict(xgb_fold, newdata = dval_fold)
  xgb_roc_obj <- suppressMessages(pROC::roc(mm_fold$test_df$FIRST_WK_DIAL, xgb_probs))
  xgb_auc     <- as.numeric(pROC::auc(xgb_roc_obj))

  # Optimal threshold: point on ROC curve closest to (0, 1)
  sens_vec   <- xgb_roc_obj$sensitivities
  spec_vec   <- xgb_roc_obj$specificities
  opt_idx    <- which.min(sqrt((1 - sens_vec)^2 + (1 - spec_vec)^2))
  opt_thresh <- xgb_roc_obj$thresholds[opt_idx]
  opt_sens   <- sens_vec[opt_idx]
  opt_spec   <- spec_vec[opt_idx]

  pred_bin <- as.integer(xgb_probs >= opt_thresh)
  tp       <- sum(pred_bin == 1 & label_val == 1)
  fp       <- sum(pred_bin == 1 & label_val == 0)
  opt_prec <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)

  xgb_auprc <- compute_auprc(xgb_probs, label_val)

  tibble(
    fold       = fold_k,
    rf_auc     = rf_auc,
    xgb_auc    = xgb_auc,
    xgb_auprc  = xgb_auprc,
    xgb_sens   = opt_sens,
    xgb_spec   = opt_spec,
    xgb_prec   = opt_prec,
    xgb_thresh = opt_thresh,
    xgb_roc    = list(xgb_roc_obj)
  )
}

# 4) XGBoost hyperparameters matching the final model
xgb_params_cv <- list(
  objective        = "binary:logistic",
  eval_metric      = "auc",
  eta              = 0.03,
  max_depth        = 5,
  subsample        = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 5,
  lambda           = 1,
  alpha            = 0.5,
  tree_method      = "hist"
)

# 5) Run all 5 folds
cv_results <- map_dfr(1:K, run_one_fold,
                      train_df_full  = train_df_full,
                      donor_fold_map = donor_fold_map,
                      xgb_params     = xgb_params_cv,
                      n_rounds       = 500)

write.csv(cv_results %>% select(-xgb_roc), "cv_fold_auc.csv", row.names = FALSE)

# 6) Summarise with mean, SD, 95% CI (t-distribution)
t_crit <- qt(0.975, df = K - 1)

cv_summary <- cv_results %>%
  summarise(
    rf_mean   = mean(rf_auc),
    rf_sd     = sd(rf_auc),
    rf_ci_lo  = mean(rf_auc) - t_crit * sd(rf_auc) / sqrt(K),
    rf_ci_hi  = mean(rf_auc) + t_crit * sd(rf_auc) / sqrt(K),
    xgb_mean  = mean(xgb_auc),
    xgb_sd    = sd(xgb_auc),
    xgb_ci_lo = mean(xgb_auc) - t_crit * sd(xgb_auc) / sqrt(K),
    xgb_ci_hi = mean(xgb_auc) + t_crit * sd(xgb_auc) / sqrt(K)
  )

print(cv_summary)
write.csv(cv_summary, "cv_summary_auc.csv", row.names = FALSE)

# XGBoost threshold metrics summary (mean ± SD across folds)
cv_thresh_summary <- cv_results %>%
  summarise(
    auroc_mean  = mean(xgb_auc),   auroc_sd  = sd(xgb_auc),
    auroc_ci_lo = mean(xgb_auc)  - t_crit * sd(xgb_auc)  / sqrt(K),
    auroc_ci_hi = mean(xgb_auc)  + t_crit * sd(xgb_auc)  / sqrt(K),
    auprc_mean  = mean(xgb_auprc), auprc_sd  = sd(xgb_auprc),
    auprc_ci_lo = mean(xgb_auprc) - t_crit * sd(xgb_auprc) / sqrt(K),
    auprc_ci_hi = mean(xgb_auprc) + t_crit * sd(xgb_auprc) / sqrt(K),
    sens_mean   = mean(xgb_sens),  sens_sd   = sd(xgb_sens),
    spec_mean   = mean(xgb_spec),  spec_sd   = sd(xgb_spec),
    prec_mean   = mean(xgb_prec, na.rm = TRUE), prec_sd = sd(xgb_prec, na.rm = TRUE)
  )
cat("\nXGBoost CV threshold metrics (mean ± SD, 95% CI across 5 folds):\n")
cat(sprintf("  AUROC:       %.3f ± %.3f  (%.3f, %.3f)\n",
            cv_thresh_summary$auroc_mean, cv_thresh_summary$auroc_sd,
            cv_thresh_summary$auroc_ci_lo, cv_thresh_summary$auroc_ci_hi))
cat(sprintf("  AUPRC:       %.3f ± %.3f  (%.3f, %.3f)\n",
            cv_thresh_summary$auprc_mean, cv_thresh_summary$auprc_sd,
            cv_thresh_summary$auprc_ci_lo, cv_thresh_summary$auprc_ci_hi))
cat(sprintf("  Sensitivity: %.3f ± %.3f\n", cv_thresh_summary$sens_mean,  cv_thresh_summary$sens_sd))
cat(sprintf("  Specificity: %.3f ± %.3f\n", cv_thresh_summary$spec_mean,  cv_thresh_summary$spec_sd))
cat(sprintf("  Precision:   %.3f ± %.3f\n", cv_thresh_summary$prec_mean,  cv_thresh_summary$prec_sd))

# 7) Plot AUC per fold
cv_long <- cv_results %>%
  pivot_longer(cols = c(rf_auc, xgb_auc),
               names_to  = "model",
               values_to = "auc") %>%
  mutate(model = recode(model,
                        rf_auc  = "Random Forest",
                        xgb_auc = "XGBoost"))

cv_means <- cv_long %>%
  group_by(model) %>%
  summarise(mean_auc = mean(auc), .groups = "drop")

p_cv <- ggplot(cv_long, aes(x = fold, y = auc, colour = model, group = model)) +
  geom_point(size = 3) +
  geom_line(linetype = "dashed") +
  geom_hline(data = cv_means,
             aes(yintercept = mean_auc, colour = model),
             linewidth = 0.8) +
  scale_x_continuous(breaks = 1:K) +
  labs(x = "Fold", y = "AUC", colour = "Model",
       title = "5-Fold Cross-Validation AUC") +
  theme_minimal(base_size = 13)

ggsave("cv_auc_by_fold.pdf", plot = p_cv, width = 7, height = 5)

# ROC curves for all 5 folds with optimal threshold marked
roc_curve_df <- map_dfr(1:K, function(k) {
  ro <- cv_results$xgb_roc[[k]]
  tibble(fold = k, fpr = 1 - ro$specificities, tpr = ro$sensitivities)
})

opt_points_df <- cv_results %>%
  transmute(fold, fpr = 1 - xgb_spec, tpr = xgb_sens)

p_roc_cv <- ggplot(roc_curve_df,
                   aes(x = fpr, y = tpr, colour = factor(fold), group = fold)) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_point(data = opt_points_df,
             aes(x = fpr, y = tpr, colour = factor(fold)),
             shape = 4, size = 3, stroke = 1.5, show.legend = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  scale_colour_brewer(palette = "Set1", name = "Fold") +
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)",
       title = "XGBoost 5-Fold CV ROC Curves",
       subtitle = "× = optimal threshold (min Euclidean distance to (0,1))") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_bw(base_size = 12)

ggsave("cv_xgb_roc_curves.pdf", plot = p_roc_cv, width = 6, height = 6)
cat("Cross-validation complete. Results saved.\n")

set.seed(333)
donor_DCD = unique(df1_DCD_vital$DONOR_ID)
test_ids_DCD <- sample(donor_DCD, size = round(0.3 * length(donor_DCD)))

train_df_DCD <- df1_DCD_vital %>% filter(!(DONOR_ID %in% test_ids_DCD))
test_df_DCD  <- df1_DCD_vital %>% filter(DONOR_ID %in% test_ids_DCD)

rec <- recipe(~ ., data = train_df_DCD) %>%
  step_impute_mean(all_numeric()) %>%   
  step_impute_mode(all_nominal()) %>%   
  prep()

train_imp_DCD <- bake(rec, new_data = NULL)     
test_imp_DCD  <- bake(rec, new_data = test_df_DCD)  

train_DCD = train_imp_DCD %>% select(-DONOR_ID, -CTR_CODE)
test_DCD = test_imp_DCD %>% select(-DONOR_ID, -CTR_CODE)

train_DCD$FIRST_WK_DIAL <- factor(train_DCD$FIRST_WK_DIAL, levels = c("N", "Y"), labels = c("No", "Yes"))
test_DCD$FIRST_WK_DIAL  <- factor(test_DCD$FIRST_WK_DIAL,  levels = c("N", "Y"), labels = c("No", "Yes"))

rf_ranger_DCD <- ranger(
  FIRST_WK_DIAL ~ ., 
  data = train_DCD,
  num.trees = 1000,
  probability = TRUE,
  importance = "permutation"
)
rf_probs_test <- predict(rf_ranger_DCD, data = test_DCD)$predictions[, 2]
rf_labels_test <- test_imp_DCD$FIRST_WK_DIAL
roc_test_DCD <- roc(rf_labels_test, rf_probs_test)
auc_test_DCD <- auc(roc_test_DCD)
auc_test_DCD

write.csv(train_DCD, "train_DCD.csv", row.names = FALSE)
write.csv(test_DCD, "test_DCD.csv", row.names = FALSE)
print("Random forest DCD group modeling with vital signs completed. ")

##########################################################
################### Xgboost ##############################
##########################################################
print("Xgboost Modeling DCD group with vital signs added...")
mm_main <- prepare_matrix_data(train_DCD, test_DCD, outcome = "FIRST_WK_DIAL")

train_DCD <- mm_main$train_df
test_DCD  <- mm_main$test_df
train_mat <- mm_main$train_mat
test_mat  <- mm_main$test_mat

label_train <- as.numeric(train_DCD$FIRST_WK_DIAL == "Yes")
label_test  <- as.numeric(test_DCD$FIRST_WK_DIAL == "Yes")

pos_weight <- sum(label_train == 0) / sum(label_train == 1)

dtrain <- xgb.DMatrix(train_mat, label = label_train)
dtest  <- xgb.DMatrix(test_mat,  label = label_test)

params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  eta = 0.03,
  max_depth = 5,
  subsample = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 5,
  scale_pos_weight = pos_weight,
  lambda = 1,
   alpha = 0.5,
  tree_method = "hist"
)

set.seed(123)
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 500,
  print_every_n = 50
)

pred <- predict(xgb_model, newdata = dtest)
auc_xgb <- pROC::auc(pROC::roc(test_DCD$FIRST_WK_DIAL, pred))
auc_xgb

## NOTE: this 70/30 model (no calibration holdout) is NOT the one the figures
## use. plot_objects.RData is written at the end of the calibration section
## below, from the 70/15/15 model, so that all three figure scripts report the
## same fit:
##   - compute_feature_importance.R / compute_pdp.R load plot_objects.RData
##   - compute_permutation_importance.R independently reproduces the 70/15/15
##     split and retrains with identical params and seed
## The AUC printed just above is therefore a cross-check on the simpler split,
## not the number the figures correspond to.

print("Xgboost DCD group modeling with vital signs completed. ")

#################################################################
################### add calibration #############################
#################################################################

df1_DCD_vital = df1_DCD %>% inner_join(df_long_flat_add, by = "DONOR_ID")

set.seed(333)
donor_DCD = unique(df1_DCD_vital$DONOR_ID)
test_ids_DCD <- sample(donor_DCD, size = round(0.3 * length(donor_DCD)))

train_df_DCD <- df1_DCD_vital %>% filter(!(DONOR_ID %in% test_ids_DCD))
test_df_DCD  <- df1_DCD_vital %>% filter(DONOR_ID %in% test_ids_DCD)

donor_test <- unique(test_df_DCD$DONOR_ID)
calib_prop <- 0.5
set.seed(333)
calib_ids_DCD <- sample(donor_test, size = round(calib_prop * length(donor_test)))

calib_df_DCD <- test_df_DCD %>% filter(DONOR_ID %in% calib_ids_DCD)
test_df_DCD  <- test_df_DCD %>% filter(!(DONOR_ID %in% calib_ids_DCD))

rec <- recipe(~ ., data = train_df_DCD) %>%
  step_impute_mean(all_numeric()) %>%   
  step_impute_mode(all_nominal()) %>%   
  prep()

train_imp_DCD <- bake(rec, new_data = NULL)
calib_imp_DCD <- bake(rec, new_data = calib_df_DCD)
test_imp_DCD  <- bake(rec, new_data = test_df_DCD)

train_DCD = train_imp_DCD %>% select(-DONOR_ID, -CTR_CODE)
calib_DCD = calib_imp_DCD %>% select(-DONOR_ID, -CTR_CODE)
test_DCD  = test_imp_DCD %>% select(-DONOR_ID, -CTR_CODE)

test_ids_df <- test_imp_DCD %>% select(DONOR_ID, CTR_CODE)

train_DCD$FIRST_WK_DIAL <- factor(train_DCD$FIRST_WK_DIAL, levels = c("N", "Y"), labels = c("No", "Yes"))
calib_DCD$FIRST_WK_DIAL <- factor(calib_DCD$FIRST_WK_DIAL, levels = c("N", "Y"), labels = c("No", "Yes"))
test_DCD$FIRST_WK_DIAL  <- factor(test_DCD$FIRST_WK_DIAL,  levels = c("N", "Y"), labels = c("No", "Yes"))

mm_cal <- prepare_matrix_data(train_DCD, test_DCD, calib_DCD, outcome = "FIRST_WK_DIAL")

train_DCD <- mm_cal$train_df
test_DCD  <- mm_cal$test_df
calib_DCD <- mm_cal$calib_df
train_mat <- mm_cal$train_mat
test_mat  <- mm_cal$test_mat
calib_mat <- mm_cal$calib_mat

label_train <- as.numeric(train_DCD$FIRST_WK_DIAL == "Yes")
label_calib <- as.numeric(calib_DCD$FIRST_WK_DIAL == "Yes")
label_test  <- as.numeric(test_DCD$FIRST_WK_DIAL == "Yes")

pos_weight <- sum(label_train == 0) / sum(label_train == 1)

dtrain <- xgb.DMatrix(train_mat, label = label_train)
dcalib <- xgb.DMatrix(calib_mat, label = label_calib)
dtest  <- xgb.DMatrix(test_mat,  label = label_test)

params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  eta = 0.03,
  max_depth = 5,
  subsample = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 5,
  scale_pos_weight = pos_weight,
  lambda = 1,
  alpha = 0.5,
  tree_method = "hist"
)

set.seed(333)
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 500,
  print_every_n = 50
)

pred <- predict(xgb_model, newdata = dtest)
auc_xgb <- pROC::auc(pROC::roc(test_DCD$FIRST_WK_DIAL, pred))
auc_xgb

pred_calib_raw <- predict(xgb_model, newdata = dcalib)
pred_test_raw  <- predict(xgb_model, newdata = dtest)

#### isotonic regression #####
ord <- order(pred_calib_raw)

x_sorted <- pred_calib_raw[ord]
y_sorted <- label_calib[ord]

iso_fit <- isoreg(x_sorted, y_sorted)

calib_fun <- approxfun(
  iso_fit$x,
  iso_fit$yf,
  method = "linear",
  rule = 2
)

p_test_iso <- calib_fun(pred_test_raw)
p_test_iso <- pmin(pmax(p_test_iso, 0), 1)

save(xgb_model, test_mat, test_ids_df, test_DCD, calib_fun,
     file = "center_calibration_inputs.RData")

roc_raw <- roc(label_test, pred_test_raw)
auc_raw <- auc(roc_raw)
roc_iso <- roc(label_test, p_test_iso)
auc_iso <- auc(roc_iso)

auc_raw
auc_iso

# Brier score calculation
brier_raw <- mean((pred_test_raw - label_test)^2)
brier_cal <- mean((p_test_iso   - label_test)^2)
cat(sprintf("Brier score (raw XGBoost):        %.4f\n", brier_raw))
cat(sprintf("Brier score (calibrated XGBoost): %.4f\n", brier_cal))

auprc_raw <- compute_auprc(pred_test_raw, label_test)
auprc_cal <- compute_auprc(p_test_iso,   label_test)
cat(sprintf("AUPRC (raw XGBoost):        %.4f\n", auprc_raw))
cat(sprintf("AUPRC (calibrated XGBoost): %.4f\n", auprc_cal))
cat(sprintf("AUPRC no-skill baseline:    %.4f  (= outcome prevalence)\n", mean(label_test)))

make_calib_plot_qtile <- function(y, p, title, n_bins = 10) {
  df <- data.frame(y = y, p = p) %>%
    mutate(bin = ntile(p, n_bins)) %>%
    group_by(bin) %>%
    summarise(
      mean_pred = mean(p),
      obs_rate  = mean(y),
      n         = n(),
      .groups   = "drop"
    )
  
  ggplot(df, aes(x = mean_pred, y = obs_rate)) +
    geom_point() +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      x = "Mean predicted probability (by quantile bin)",
      y = "Observed event rate",
      title = title
    ) +
    theme_minimal() +
    guides(size = "none")
}

p_raw_q <- make_calib_plot_qtile(
  label_test, pred_test_raw,
  "Calibration (test set, raw XGBoost)", n_bins = 10
)

p_cal_q <- make_calib_plot_qtile(
  label_test, p_test_iso,
  "Calibration (test set, calibrated XGBoost)", n_bins = 10
)

p_raw_q
p_cal_q

combined_plot <- ggarrange(
  p_raw_q,
  p_cal_q,
  ncol = 2,
  nrow = 1
)

ggsave(
  filename = "calibration_xgboost.pdf",
  plot = combined_plot,
  width = 12,
  height = 5,
  units = "in"
)

hl_raw <- hoslem.test(
  x = label_test,
  y = pred_test_raw,
  g = 10
)

hl_iso <- hoslem.test(
  x = label_test,
  y = p_test_iso,
  g = 10
)

hl_raw
hl_iso

# Save test_DCD for future reuse
write.csv(test_DCD, file.path(OUT_DIR, "test_DCD.csv"), row.names = FALSE)

# Subgroup summary
# NOTE: KDPI-based subgroups removed 2026-08 -- KDPI dropped from df_subset2,
# no longer available on test_DCD. AGE_DON subgroups unaffected, kept as-is.
subgroups <- list(
  "AGE_DON <= 60"  = test_DCD %>% filter(AGE_DON <= 60),
  "AGE_DON > 60"   = test_DCD %>% filter(AGE_DON >  60)
)

for (label in names(subgroups)) {
  sg <- subgroups[[label]]
  n_rows   <- nrow(sg)
  n_events <- sum(sg$FIRST_WK_DIAL == "Yes", na.rm = TRUE)
  cat(sprintf("%s: n = %d, events = %d, DGF rate = %.1f%%\n",
              label, n_rows, n_events, n_events / n_rows * 100))
}

# ---- cached inputs for the figure scripts ----------------------------------
# xgb_model/train_DCD/test_DCD/train_mat at this point are the CALIBRATION
# (70/15/15) fit. Saving here -- rather than after the 70/30 fit above -- is
# what keeps the feature-importance, PDP and permutation-importance figures on
# a single model. compute_feature_importance.R and compute_pdp.R read this file.
save(xgb_model, train_DCD, test_DCD, train_mat, file = "plot_objects.RData")
cat("Saved plot_objects.RData (calibration-section model).\n")

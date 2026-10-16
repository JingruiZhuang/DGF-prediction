##############################################################
## Permutation importance (full test set) -- COMPUTE step only.
##
## This is the expensive part (100-repeat permutation importance on the full
## test set) -- run this rarely, only when the underlying model/data actually
## changes. Saves the full results table to
## <PROJECT_DIR>/permutation_importance.csv.
##
## The model is LOADED from plot_objects.RData (written by build_model.R),
## not retrained here. Retraining would not reproduce it bit-for-bit:
## set.seed() does not control XGBoost's internal RNG, and tree_method="hist"
## accumulates histograms in a thread-dependent order. Loading is what keeps
## this figure on the same model as the feature-importance and PDP figures.
##
## For re-plotting (colors, labels, layout, etc.) use
## plot_permutation_importance.R instead -- it only reads that CSV and
## never touches the model, so it's safe to rerun as often as you like.
##############################################################

## ---- paths -----------------------------------------------------------------
## Defaults to the repository root, so running the scripts from the repo
## folder needs no editing. Override it if the data lives elsewhere.
PROJECT_DIR <- "."
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)
setwd(file.path(PROJECT_DIR, "data"))

library(tidyverse)
library(xgboost)
library(pROC)

## ── load the fitted model + data written by build_model.R ───────────────────
load("plot_objects.RData")   # xgb_model, train_DCD, test_DCD, train_mat

cat(sprintf("Loaded plot_objects.RData: train n = %d, test n = %d\n",
            nrow(train_DCD), nrow(test_DCD)))

## prepare_matrix_data() -- verbatim from build_model.R
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

## ── permutation importance helpers -- verbatim from build_model.R ───────────
prep_perm_base <- function(subgroup_df, train_df_ref, model, outcome = "FIRST_WK_DIAL") {
  mm          <- prepare_matrix_data(train_df_ref, subgroup_df, outcome = outcome)
  base_mat    <- mm$test_mat
  base_labels <- as.numeric(mm$test_df[[outcome]] == "Yes")
  base_pred   <- predict(model, newdata = xgb.DMatrix(base_mat))
  baseline_auc <- as.numeric(suppressMessages(
    pROC::auc(pROC::roc(base_labels, base_pred))
  ))
  list(base_mat = base_mat, base_labels = base_labels,
       baseline_auc = baseline_auc, train_df = mm$train_df)
}

auc_from_mat <- function(perm_mat, base_labels, model) {
  pred <- predict(model, newdata = xgb.DMatrix(perm_mat))
  as.numeric(suppressMessages(pROC::auc(pROC::roc(base_labels, pred))))
}

permute_vital_block_auc <- function(subgroup_df, train_df_ref, model,
                                    n_repeats = 100, outcome = "FIRST_WK_DIAL") {
  vital_cols <- grep("^(SYSTOLIC_PRESSURE|DIASTOLIC_PRESSURE|O2_SATURATION)_", names(subgroup_df), value = TRUE)
  base       <- prep_perm_base(subgroup_df, train_df_ref, model, outcome)
  vital_idx  <- which(colnames(base$base_mat) %in% vital_cols)
  drops <- vapply(seq_len(n_repeats), function(i) {
    perm_mat <- base$base_mat
    idx <- sample(nrow(perm_mat))
    perm_mat[, vital_idx] <- perm_mat[idx, vital_idx, drop = FALSE]
    base$baseline_auc - auc_from_mat(perm_mat, base$base_labels, model)
  }, numeric(1))
  list(baseline_auc = base$baseline_auc, mean_drop = mean(drops), sd_drop = sd(drops))
}

permute_individual_features_auc <- function(subgroup_df, train_df_ref, model,
                                             n_repeats = 100, outcome = "FIRST_WK_DIAL") {
  vital_cols   <- grep("^(SYSTOLIC_PRESSURE|DIASTOLIC_PRESSURE|O2_SATURATION)_", names(subgroup_df), value = TRUE)
  feature_cols <- setdiff(names(subgroup_df), c(outcome, vital_cols))
  base         <- prep_perm_base(subgroup_df, train_df_ref, model, outcome)

  results <- lapply(feature_cols, function(feat) {
    col <- base$train_df[[feat]]
    mat_cols <- if (is.factor(col)) {
      intersect(paste0(feat, levels(col)), colnames(base$base_mat))
    } else {
      intersect(feat, colnames(base$base_mat))
    }
    if (length(mat_cols) == 0) return(NULL)
    mat_idx <- which(colnames(base$base_mat) %in% mat_cols)
    drops <- vapply(seq_len(n_repeats), function(i) {
      perm_mat <- base$base_mat
      idx <- sample(nrow(perm_mat))
      perm_mat[, mat_idx] <- perm_mat[idx, mat_idx, drop = FALSE]
      base$baseline_auc - auc_from_mat(perm_mat, base$base_labels, model)
    }, numeric(1))
    data.frame(feature = feat, mean_drop = mean(drops), sd_drop = sd(drops),
               baseline_auc = base$baseline_auc, stringsAsFactors = FALSE)
  })
  do.call(rbind, Filter(Negate(is.null), results))
}

## ── run on the full test set ────────────────────────────────────────────────
subgroup_list <- list("Main (full test set)" = test_DCD)

all_results <- map_dfr(names(subgroup_list), function(sg_label) {
  sg_df <- subgroup_list[[sg_label]]
  cat(sprintf("Running permutation importance for %s (n = %d)...\n",
              sg_label, nrow(sg_df)))

  block_res <- permute_vital_block_auc(sg_df, train_DCD, xgb_model, n_repeats = 100)
  indiv_res <- permute_individual_features_auc(sg_df, train_DCD, xgb_model, n_repeats = 100)

  block_row <- data.frame(
    type = "vital_block", subgroup = sg_label, feature = "vital_block",
    baseline_auc = block_res$baseline_auc,
    mean_drop = block_res$mean_drop, sd_drop = block_res$sd_drop,
    stringsAsFactors = FALSE
  )
  indiv_res$type     <- "individual_feature"
  indiv_res$subgroup <- sg_label

  rbind(block_row,
        indiv_res[, c("type", "subgroup", "feature", "baseline_auc", "mean_drop", "sd_drop")])
})

write.csv(all_results,
          file.path(PROJECT_DIR, "permutation_importance.csv"),
          row.names = FALSE)
cat("Saved permutation_importance.csv -- rerun plot_permutation_importance.R to (re)draw the plots from it.\n")

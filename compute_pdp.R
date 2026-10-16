##############################################################
## Partial dependence data -- COMPUTE step only.
##
## Loads the trained model (plot_objects.RData) and computes raw PDP grids
## (30 points each) for the 13 continuous features, plus the HMP/PUMP_KI
## bar data. Saves everything to data/pdp_13_features_data.csv and
## data/pdp_pump_ki_data.csv. Run this only when the underlying model
## changes (i.e. after build_model.R retrains) -- this is the slow step
## (repeated predict() calls across 30 grid points x 13 features).
##
## For re-plotting (trimming threshold, colors, layout, etc.) use
## plot_pdp.R instead -- it only reads these two CSVs and never touches
## the model.
##############################################################

## ---- paths -----------------------------------------------------------------
## Defaults to the repository root, so running the scripts from the repo
## folder needs no editing. Override it if the data lives elsewhere.
PROJECT_DIR <- "."
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)
setwd(file.path(PROJECT_DIR, "data"))

library(tidyverse)
library(xgboost)
library(pdp)

load("plot_objects.RData")  # loads xgb_model, train_DCD, test_DCD, train_mat

train_mat_df <- as.data.frame(train_mat)

pred_fun_mat <- function(object, newdata) {
  mean(predict(object, xgb.DMatrix(as.matrix(newdata))))
}

features_to_plot <- list(
  list(var = "DIALYSIS_DURATION",       label = "Dialysis duration, days"),
  list(var = "CREAT_TRR",               label = "Recipient serum creatinine at transplant, mg/dL"),
  list(var = "BUN_DON",                 label = "Donor terminal blood urea nitrogen, mg/dL"),
  list(var = "AGE_DON",                 label = "Donor age, years"),
  list(var = "BMI_DON_CALC",            label = "Donor body mass index, kg/m2"),
  list(var = "CREAT_DON",               label = "Donor terminal serum creatinine, mg/dL"),
  list(var = "COLD_ISCH_KI",            label = "Cold ischemia time, hours"),
  list(var = "BMI_CALC",                label = "Recipient body mass index, kg/m2"),
  list(var = "SGPT_DON",                label = "Donor terminal serum alanine aminotransferase, U/L"),
  list(var = "SYSTOLIC_PRESSURE_slope", label = "Systolic blood pressure slope, mmHg/min"),
  list(var = "O2_SATURATION_slope",     label = "Oxygen saturation slope, %/min"),
  list(var = "AGE",                     label = "Recipient age, years")
)

## ── compute raw PDP (30 grid points, untrimmed) for each feature ───────────
pdp_long <- map_dfr(features_to_plot, function(f) {
  cat("Computing PDP for:", f$var, "\n")
  pd <- partial(
    object          = xgb_model,
    pred.var        = f$var,
    train           = train_mat_df,
    pred.fun        = pred_fun_mat,
    grid.resolution = 30,
    progress        = "none"
  )
  tibble(var = f$var, label = f$label, x_value = pd[[f$var]], yhat = pd$yhat)
})

write.csv(pdp_long, "pdp_12_features_data.csv", row.names = FALSE)
cat("Saved pdp_12_features_data.csv\n")

## ── HMP / PUMP_KI bar data ──────────────────────────────────────────────────
## PUMP_KI is one-hot encoded as two dummy columns (PUMP_KIN / PUMP_KIY);
## toggle both together so exactly one level is "on" per prediction.
pump_preds <- sapply(c("N", "Y"), function(lvl) {
  tmp <- train_mat_df
  tmp[["PUMP_KIN"]] <- as.numeric(lvl == "N")
  tmp[["PUMP_KIY"]] <- as.numeric(lvl == "Y")
  mean(predict(xgb_model, xgb.DMatrix(as.matrix(tmp))))
})
org_df <- data.frame(
  group     = factor(c("No", "Yes"), levels = c("No", "Yes")),
  mean_pred = pump_preds
)
write.csv(org_df, "pdp_pump_ki_data.csv", row.names = FALSE)
cat("Saved pdp_pump_ki_data.csv -- rerun plot_pdp.R to (re)draw the plots from these CSVs.\n")

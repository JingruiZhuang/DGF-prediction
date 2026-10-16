##############################################################
## XGBoost feature importance (Gain) -- COMPUTE step only.
##
## Loads the trained model (plot_objects.RData) and computes the top-30
## Gain-based importance table, saving it to
## data/xgb_feature_importance_top30.csv. Run this only when the
## underlying model changes (i.e. after build_model.R retrains).
##
## For re-plotting (colors, labels, layout, etc.) use
## plot_feature_importance.R instead -- it only reads that CSV.
##############################################################

## ---- paths -----------------------------------------------------------------
## Defaults to the repository root, so running the scripts from the repo
## folder needs no editing. Override it if the data lives elsewhere.
PROJECT_DIR <- "."
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)
setwd(file.path(PROJECT_DIR, "data"))

library(xgboost)
library(tidyverse)

load("plot_objects.RData")  # loads xgb_model, train_DCD, test_DCD, train_mat

source(file.path(PROJECT_DIR, "translate_varname.R"))

imp_xgb <- xgb.importance(
  model         = xgb_model,
  feature_names = colnames(train_mat)
)

imp_top30 <- imp_xgb %>%
  arrange(desc(Gain)) %>%
  slice_head(n = 30) %>%
  mutate(label = translate_varname(Feature)) %>%
  select(Feature, label, Gain, Cover, Frequency)

write.csv(imp_top30, "xgb_feature_importance_top30.csv", row.names = FALSE)
cat("Saved xgb_feature_importance_top30.csv -- rerun plot_feature_importance.R to (re)draw the plot from it.\n")

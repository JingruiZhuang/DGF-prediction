##############################################################
## XGBoost feature importance (Gain) -- PLOT step only.
##
## Reads the already-computed data/xgb_feature_importance_top30.csv
## (produced by compute_feature_importance.R) and draws the plot. Does
## NOT load the model -- safe to rerun as often as you like when tweaking
## colors/labels/layout.
##
## If xgb_feature_importance_top30.csv is missing or needs recomputing
## because the underlying model changed, run compute_feature_importance.R
## first (that's the step that touches the model).
##############################################################

## ---- paths -----------------------------------------------------------------
## Defaults to the repository root, so running the scripts from the repo
## folder needs no editing. Override it if the data lives elsewhere.
PROJECT_DIR <- "."
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)
setwd(file.path(PROJECT_DIR, "data"))

library(tidyverse)

imp_top30 <- read.csv("xgb_feature_importance_top30.csv")

p_imp <- ggplot(imp_top30, aes(x = reorder(label, Gain), y = Gain)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 30 Feature Importance (XGBoost, Gain)",
    x     = "Feature",
    y     = "Gain"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("xgb_feature_importance_top30.pdf", plot = p_imp, width = 10, height = 8)

cat("Saved xgb_feature_importance_top30.pdf from xgb_feature_importance_top30.csv\n")

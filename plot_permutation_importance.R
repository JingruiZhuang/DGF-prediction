##############################################################
## Permutation importance (full test set) -- PLOT step only.
##
## Reads the already-computed <PROJECT_DIR>/
## permutation_importance.csv (produced by compute_permutation_importance.R)
## and draws the plot. Does NOT load any model or retrain anything --
## safe to rerun as often as you like when tweaking colors/labels/layout.
##
## If permutation_importance.csv is missing or you need it recomputed
## because the underlying model/data changed, run
## compute_permutation_importance.R first (that's the expensive step).
##############################################################

## ---- paths -----------------------------------------------------------------
## Defaults to the repository root, so running the scripts from the repo
## folder needs no editing. Override it if the data lives elsewhere.
PROJECT_DIR <- "."
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)

library(tidyverse)
library(ggpubr)

all_results <- read.csv(file.path(PROJECT_DIR, "permutation_importance.csv"))

source(file.path(PROJECT_DIR, "translate_varname.R"))

make_perm_df <- function(results_df, subgroup_val) {
  results_df %>%
    filter(subgroup == subgroup_val) %>%
    mutate(
      importance    = mean_drop,
      sd            = sd_drop,
      feature       = ifelse(type == "vital_block", "Hemodynamic data block", feature),
      feature_label = ifelse(feature == "Hemodynamic data block", feature, translate_varname(feature))
    )
}

plot_perm_importance <- function(df, title) {
  top20 <- df %>%
    arrange(desc(importance)) %>%
    slice_head(n = 20) %>%
    mutate(
      feature_label  = factor(feature_label, levels = rev(feature_label)),
      is_vital_block = feature == "Hemodynamic data block"
    )

  ggplot(top20, aes(x = importance, y = feature_label, fill = is_vital_block)) +
    geom_col() +
    scale_fill_manual(
      values = c("FALSE" = "#4472C4", "TRUE" = "#ED7D31"),
      labels = c("FALSE" = "Individual feature", "TRUE" = "Hemodynamic data block"),
      name   = NULL
    ) +
    labs(x = "Mean AUC drop", y = NULL, title = title) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
}

## Main plot
p_main <- plot_perm_importance(make_perm_df(all_results, "Main (full test set)"),
                                NULL)
ggsave(file.path(PROJECT_DIR, "xgb_permutation_importance_top30_with_vital_block.pdf"),
       plot = p_main, width = 11, height = 7)
ggsave(file.path(PROJECT_DIR, "xgb_permutation_importance_top30_with_vital_block.png"),
       plot = p_main, width = 11, height = 7, dpi = 300)

cat("Saved permutation importance plot from permutation_importance.csv.\n")

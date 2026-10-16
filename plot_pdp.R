##############################################################
## Partial dependence plots -- PLOT step only.
##
## Reads the already-computed data/pdp_12_features_data.csv and
## data/pdp_pump_ki_data.csv (produced by compute_pdp.R) and draws both
## PDP figures. Does NOT load the model -- safe to rerun as often as you
## like when tweaking the trimming threshold, colors, or layout.
##
## If those CSVs are missing or need recomputing because the underlying
## model changed, run compute_pdp.R first (that's the step that touches
## the model).
##############################################################

## ---- paths -----------------------------------------------------------------
## Defaults to the repository root, so running the scripts from the repo
## folder needs no editing. Override it if the data lives elsewhere.
PROJECT_DIR <- "."
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)
setwd(file.path(PROJECT_DIR, "data"))

library(tidyverse)
library(ggpubr)

pdp_long <- read.csv("pdp_12_features_data.csv")
org_df   <- read.csv("pdp_pump_ki_data.csv") %>%
  mutate(group = factor(group, levels = c("No", "Yes")))

features_to_plot <- pdp_long %>%
  distinct(var, label) %>%
  { setNames(split(., seq_len(nrow(.))), NULL) } %>%
  lapply(function(r) list(var = r$var, label = r$label))

pd_list <- split(pdp_long, pdp_long$var)[sapply(features_to_plot, `[[`, "var")]

# ── trim leading/trailing flat segments ─────────────────────────────────────
# A step is "flat" if |Δyhat| < 5 % of that feature's own yhat range.
# Keeps one context point at each trimmed edge so the curve reaches the border.
trim_flat_tails <- function(pd, tol_frac = 0.05) {
  yhat <- pd$yhat
  n    <- length(yhat)
  if (n < 3) return(pd)
  feat_range <- diff(range(yhat))
  if (feat_range == 0) return(pd)          # entirely flat — nothing to trim
  thresh <- tol_frac * feat_range
  diffs  <- abs(diff(yhat))               # length n-1

  left <- 1L
  for (i in seq_len(n - 1L)) {
    if (diffs[i] > thresh) break
    left <- i + 1L
  }
  left <- max(1L, left - 1L)

  right <- n
  for (i in rev(seq_len(n - 1L))) {
    if (diffs[i] > thresh) break
    right <- i
  }
  right <- min(n, right + 1L)

  if (left >= right) return(pd)
  pd[left:right, ]
}

pd_trimmed <- lapply(pd_list, trim_flat_tails)

# ── shared y-axis from the 13 continuous features only ──────────────────────
all_yhats <- unlist(lapply(pd_trimmed, `[[`, "yhat"))
y_lo  <- min(all_yhats)
y_hi  <- max(all_yhats)
y_pad <- (y_hi - y_lo) * 0.05
y_lo  <- y_lo - y_pad
y_hi  <- y_hi + y_pad

# ── continuous-feature plots (shared y-axis) ─────────────────────────────────
pdp_plots <- mapply(function(pd, f) {
  ggplot(pd, aes(x = x_value, y = yhat)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    labs(x = f$label, y = "Mean predicted P(DGF)") +
    coord_cartesian(ylim = c(y_lo, y_hi)) +
    theme_minimal(base_size = 10)
}, pd_trimmed, features_to_plot, SIMPLIFY = FALSE)

# ── HMP (PUMP_KI) bar plot (own scale) ───────────────────────────────────────
p_org <- ggplot(org_df, aes(x = group, y = mean_pred)) +
  geom_col(width = 0.5, fill = "steelblue") +
  labs(x = "Hypothermic machine perfusion", y = "Mean predicted P(DGF)") +
  theme_minimal(base_size = 10)

all_plots <- c(pdp_plots, list(p_org))

## 13 panels total (12 continuous features + 1 PUMP_KI bar plot) -- keep the
## 4x4 grid (16 slots) so nothing gets silently dropped; a 4x3 grid (12 slots)
## would not fit all 13.
p_multi <- ggarrange(plotlist = all_plots, ncol = 4, nrow = 4)

ggsave("pdp_12_features.pdf", plot = p_multi, width = 16, height = 14)
cat("Saved: data/pdp_12_features.pdf\n")

# ── individual-scale version (flat tails trimmed, each feature its own y-axis) ─
pdp_plots_indiv <- mapply(function(pd, f) {
  ggplot(pd, aes(x = x_value, y = yhat)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    labs(x = f$label, y = "Mean predicted P(DGF)") +
    theme_minimal(base_size = 10)
}, pd_trimmed, features_to_plot, SIMPLIFY = FALSE)

all_plots_indiv <- c(pdp_plots_indiv, list(p_org))

p_multi_indiv <- ggarrange(plotlist = all_plots_indiv, ncol = 4, nrow = 4)

ggsave("pdp_12_features_indiv_scale.pdf", plot = p_multi_indiv, width = 16, height = 14)
cat("Saved: data/pdp_12_features_indiv_scale.pdf\n")

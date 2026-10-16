###############################################################################
## translate_varname(): raw column name -> plain-English figure label.
##
## Shared by compute_feature_importance.R and plot_permutation_importance.R so
## the feature-importance and permutation-importance figures label the same
## variable identically. Source it with:
##   source(file.path(PROJECT_DIR, "translate_varname.R"))
##
## Handles three cases, in order:
##   1. exact match in static_labels
##   2. a vital-sign summary, <BASE>_<stat>  e.g. SYSTOLIC_PRESSURE_slope
##   3. a static label with a factor level appended, e.g. ABO_DONA
## Anything unmatched is returned unchanged.
###############################################################################

translate_varname <- function(x) {
  static_labels <- c(
    DIALYSIS_DURATION          = "Recipient dialysis duration",
    CREAT_TRR                  = "Recipient serum creatinine at transplant",
    BUN_DON                    = "Donor terminal blood urea nitrogen",
    AGE_DON                    = "Donor age",
    AGE                        = "Recipient age",
    BMI_DON_CALC               = "Donor body mass index",
    BMI_CALC                   = "Recipient body mass index",
    CREAT_DON                  = "Donor terminal serum creatinine",
    COLD_ISCH_KI               = "Cold ischemia time",
    SGPT_DON                   = "Donor terminal alanine aminotransferase",
    KDPI                       = "Kidney donor profile index",
    ETHCAT_DON                 = "Donor race/ethnicity",
    ETHCAT                     = "Recipient race/ethnicity",
    GENDER_DON                 = "Donor sex",
    GENDER                     = "Recipient sex",
    NUM_PREV_TX                = "Recipient number of previous transplants",
    PREV_PREGNANCY             = "Recipient previous pregnancy",
    PREV_TX                    = "Recipient previous transplant",
    DIABETES_DON               = "Donor diabetes",
    DIABETES                   = "Recipient diabetes",
    INIT_CPRA                  = "Recipient CPRA at listing",
    HLAMIS                     = "HLA mismatch",
    ABO_MAT                    = "ABO compatibility",
    ABO_DON                    = "Donor ABO blood type",
    ABO                        = "Recipient ABO blood type",
    COD_CAD_DON                = "Donor cause of death",
    DIAG_KI                    = "Recipient kidney disease diagnosis",
    HIST_HYPERTENS_DON         = "Donor history of hypertension",
    HIST_CAD                   = "Donor history of coronary artery disease",
    HIST_CIG_DON               = "Donor smoking history",
    ALCOHOL_HEAVY_DON          = "Donor heavy alcohol use",
    NON_HRT_DON                = "DCD donor",
    ORG_REC_ON2                = "Hypothermic machine perfusion",
    PUMP_KI                    = "Hypothermic machine perfusion",
    BIOPSY_PERFORMED           = "Biopsy performed",
    GLOMERULI_ADEQUATE_MISSING     = "Glomeruli adequate (missing indicator)",
    GLOMERULOSCLEROSIS_ORD_MISSING = "Glomerulosclerosis (missing indicator)",
    FIBROSIS_ORD_MISSING           = "Interstitial fibrosis (missing indicator)",
    VASCULAR_ORD_MISSING           = "Vascular changes (missing indicator)",
    GLOMERULI_ADEQUATE         = "Glomeruli adequate",
    GLOMERULOSCLEROSIS_ORD     = "Glomerulosclerosis",
    FIBROSIS_ORD               = "Interstitial fibrosis",
    VASCULAR_ORD               = "Vascular changes",
    SYSTOLIC_PRESSURE_slope_missing = "Systolic BP slope (monitoring absent)",
    O2_SATURATION_slope_missing     = "O2 saturation slope (monitoring absent)"
  )
  vital_bases <- c(
    SYSTOLIC_PRESSURE  = "Donor SBP",
    DIASTOLIC_PRESSURE = "Donor DBP",
    O2_SATURATION      = "Oxygen saturation"
  )
  stat_labels <- c(
    max   = "maximum",
    sd    = "SD",
    range = "range",
    last  = "last recorded",
    slope = "slope"
  )
  sorted_static <- names(static_labels)[order(nchar(names(static_labels)), decreasing = TRUE)]
  sorted_vitals <- names(vital_bases)[order(nchar(names(vital_bases)), decreasing = TRUE)]

  sapply(x, function(v) {
    if (v %in% names(static_labels)) return(unname(static_labels[v]))
    for (base in sorted_vitals)
      for (stat in names(stat_labels)) {
        if (v == paste0(base, "_", stat)) {
          # SYSTOLIC_PRESSURE/DIASTOLIC_PRESSURE already carry "Donor" in
          # vital_bases; O2_SATURATION doesn't, so add it here without
          # double-prefixing the other two.
          base_label <- vital_bases[base]
          if (base == "O2_SATURATION") {
            base_label <- paste0("Donor ", tolower(substring(base_label, 1, 1)),
                                  substring(base_label, 2))
          }
          return(paste0(base_label, " (", stat_labels[stat], ")"))
        }
      }
    for (base in sorted_static)
      if (startsWith(v, base) && nchar(v) > nchar(base))
        return(paste0(static_labels[base], " (", substring(v, nchar(base) + 1), ")"))
    v
  }, USE.NAMES = FALSE)
}

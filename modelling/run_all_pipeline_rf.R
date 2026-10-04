# ==============================================================================
# run_all_pipeline_rf.R
#
# Random-forest pipeline: run all the steps, in order
#
# Steps  : each one can be switched off below
#          10) nested CV: tuning, training and all the diagnostics
#          11) exact SHAP values (TreeSHAP)
#          12) SHAP values for dependent covariates (Aas et al. 2019) -- slow
#          13) t-SNE
#          14) daily and monthly prediction maps -- very long (one prediction
#              per day of the grid)
#          Steps 11 to 14 process one frequency and one scheme: target_freq
#          and target_scheme if they are set, otherwise the defaults of
#          00_model_config.R.
#
# Usage  : from the project root (the folder of config.R)
#            source("modelling/run_all_pipeline_rf.R")
# ==============================================================================


# ---- Configuration -----------------------------------------------------------

if (!file.exists("config.R")) {
  stop("Set the working directory to the project root (the folder of ",
       "config.R) before running this script.")
}
source("config.R")   # modelling_dir

run_nested_cv       <- TRUE   # 10
run_shap_treeshap   <- TRUE   # 11
run_shap_aas        <- TRUE   # 12
run_tsne            <- TRUE   # 13
run_prediction_maps <- TRUE   # 14


# ---- Functions ---------------------------------------------------------------

# Source one script of the pipeline and print its duration.
#   enabled : FALSE = skip the step
#   script  : file name of the script, in <modelling_dir>
#   name    : name of the step, for the messages
run_step <- function(enabled, script, name) {
  if (!isTRUE(enabled)) {
    cat("\n[SKIPPED]", name, "\n")
    return(invisible(NULL))
  }
  cat("\n====================\nSTART:", name, "--", format(Sys.time()),
      "\n====================\n")
  t0 <- Sys.time()
  source(file.path(modelling_dir, script))
  cat(sprintf("END: %s (%.1f min)\n", name,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}


# ---- Pipeline ----------------------------------------------------------------

run_step(run_nested_cv,       "10_run_nested_cv_training.R", "Nested CV")
run_step(run_shap_treeshap,   "11_run_shap_treeshap.R",      "SHAP (TreeSHAP)")
run_step(run_shap_aas,        "12_run_shap_dependent_aas.R",
         "SHAP (Aas, dependent covariates)")
run_step(run_tsne,            "13_run_tsne.R",               "t-SNE")
run_step(run_prediction_maps, "14_run_prediction_maps.R",    "Prediction maps")

cat("\nRandom-forest pipeline finished.\n")

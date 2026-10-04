# ==============================================================================
# 03_models.R
#
# Random-forest pipeline: model backend (ranger)
#
# Defines : make_formula(), fit_rf(), predict_rf(), get_importance_rf()
# ==============================================================================

# Model formula: NASC ~ all the covariates.
make_formula <- function() {
  as.formula(paste("NASC ~", paste(covariates_all, collapse = " + ")))
}

# Fit a random forest.
#   train_df : training data
#   params   : list with mtry, min.node.size and num.trees
fit_rf <- function(train_df, params) {
  ranger::ranger(
    formula       = make_formula(),
    data          = train_df,
    mtry          = params$mtry,
    min.node.size = params$min.node.size,
    num.trees     = params$num.trees,
    importance    = "permutation"
  )
}

# Predictions of a model on new data.
predict_rf <- function(model, newdata) {
  predict(model, data = newdata)$predictions
}

# Permutation importance of a model, as a tibble (variable, importance).
get_importance_rf <- function(model) {
  imp <- ranger::importance(model)
  tibble(variable = names(imp), importance = unname(imp))
}

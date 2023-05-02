# This script includes all functions used in the dementia genetic prediction project
# Listed below
# 1. '%!in%': not in
# 2. zeroVar: find variables of zero variance in a dataset
# 3. extract_cv_results: extract cv results from h2o model object

#=================================================
'%!in%' = function(x,y)!('%in%'(x,y))

#=================================================
zeroVar = function(data, useNA = 'ifany') {
  out = apply(data, 2, function(x) {length(table(x, useNA = useNA))})
  which(out == 1)
}

#=================================================
extract_cv_results = function(model_name) {
  result_lst = list()
  cv_metrics = model_name@model$cross_validation_metrics@metrics
  
  MSE = cv_metrics$MSE
  RMSE = cv_metrics$RMSE
  LogLoss = cv_metrics$logloss
  mean_per_class_error = cv_metrics$mean_per_class_error
  AUC = cv_metrics$AUC
  AUCPR = cv_metrics$pr_auc
  Gini = cv_metrics$Gini
  R_2 = cv_metrics$r2
  max_criteria_and_metric_scores = cv_metrics$max_criteria_and_metric_scores
  threshold_f1_index = max_criteria_and_metric_scores %>% filter(metric == "max absolute_mcc") %>% pull(idx)
  max_f1 = max_criteria_and_metric_scores %>% filter(metric == "max absolute_mcc") %>% pull(value)
  # Prediction with threshold F1
  threshold_results = cv_metrics$thresholds_and_metric_scores[(threshold_f1_index+1),]
  
  model_vec = c(MSE, RMSE, LogLoss, mean_per_class_error, AUC, AUCPR, Gini, R_2)
  names(model_vec) = c("MSE", "RMSE", "LogLoss", "mean_per_class_error", "AUC", 
                       "AUCPR", "Gini", "R_2")
  result_lst = c(result_lst, list(model_vec))
  
  perform_vec = list(threshold_results)
  result_lst = c(result_lst, perform_vec)
  return(result_lst)
}

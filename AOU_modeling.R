# Modeling of AOU data
# AD AFR indsig PRS
extract_columns = setdiff(coef_prs_AFR_indsig$variable, "Intercept")
aou_prs_matrix = sample_AFR_AOU[, c(coef_prs_AFR_indsig$variable, "offset_demo")]

coefficients = c(coef_prs_AFR_indsig$mean_coef, 1)
length(coefficients) # 3
# Calculate linear predictor (log-odds) using the coefficients
linear_predictor = as.matrix(aou_prs_matrix) %*% coefficients

# Convert linear predictor to probabilities using the logistic function
predicted_probs = plogis(linear_predictor)
pred_df = cbind(sample_AFR_AOU$dementia, as.data.frame(predicted_probs))
names(pred_df) = c("true_label", "p1")
pr_PRS.AFR.indsig = pr.curve(scores.class0 = pred_df$p1, 
                             weights.class0 = pred_df$true_label)$auc.integral
roc_PRS.AFR.indsig = roc(pred_df$true_label, pred_df$p1) %>% auc()

# AD AFR map PRS
extract_columns = setdiff(coef_prs_AFR_map$variable, "Intercept")
aou_prs_matrix = sample_AFR_AOU[, c(coef_prs_AFR_map$variable, "offset_demo")]

coefficients = c(coef_prs_AFR_map$mean_coef, 1)
length(coefficients) # 3
# Calculate linear predictor (log-odds) using the coefficients
linear_predictor = as.matrix(aou_prs_matrix) %*% coefficients

# Convert linear predictor to probabilities using the logistic function
predicted_probs = plogis(linear_predictor)
pred_df = cbind(sample_AFR_AOU$dementia, as.data.frame(predicted_probs))
names(pred_df) = c("true_label", "p1")
pr_PRS.AFR.map = pr.curve(scores.class0 = pred_df$p1, 
                          weights.class0 = pred_df$true_label)$auc.integral
roc_PRS.AFR.map = roc(pred_df$true_label, pred_df$p1) %>% auc()

# APOE e4 count
extract_columns = setdiff(coef_apoe$variable, "Intercept")
aou_apoe_matrix = sample_AFR_AOU[, c(coef_apoe$variable, "offset_demo")]

coefficients = c(coef_prs_AFR_indsig$mean_coef, 1)
length(coefficients) # 3
# Calculate linear predictor (log-odds) using the coefficients
linear_predictor = as.matrix(aou_apoe_matrix) %*% coefficients

# Convert linear predictor to probabilities using the logistic function
predicted_probs = plogis(linear_predictor)
pred_df = cbind(sample_AFR_AOU$dementia, as.data.frame(predicted_probs))
names(pred_df) = c("true_label", "p1")
pr_apoe = pr.curve(scores.class0 = pred_df$p1, 
                   weights.class0 = pred_df$true_label)$auc.integral
roc_apoe = roc(pred_df$true_label, pred_df$p1) %>% auc()

extract_columns = setdiff(coef_map_summary$variable, "Intercept")
aou_snp_map_matrix = sample_AFR_AOU[, c(coef_map_summary$variable, "offset_demo")]

coefficients = c(coef_map_summary$mean_coef, 1)
length(coefficients) # 3
# Calculate linear predictor (log-odds) using the coefficients
linear_predictor = as.matrix(aou_snp_map_matrix) %*% coefficients

# Convert linear predictor to probabilities using the logistic function
predicted_probs = plogis(linear_predictor)
pred_df = cbind(sample_AFR_AOU$dementia, as.data.frame(predicted_probs))
names(pred_df) = c("true_label", "p1")
pr_SNP.all.map = pr.curve(scores.class0 = pred_df$p1, 
                          weights.class0 = pred_df$true_label)$auc.integral
roc_SNP.all.map = roc(pred_df$true_label, pred_df$p1) %>% auc()

model_list = c("PRS.AFR.indsig",  "PRS.AFR.map", "e4count", "SNP.all.map")
pr_list = c(pr_PRS.AFR.indsig, pr_PRS.AFR.map, pr_apoe, pr_SNP.all.map)
auc_list = c(roc_PRS.AFR.indsig, roc_PRS.AFR.map, roc_apoe, roc_SNP.all.map)
model_summary = cbind(model_list, pr_list, auc_list) %>% 
  as.data.frame() %>% mutate(pr_list = as.numeric(pr_list),
                             auc_list = as.numeric(auc_list))
names(model_summary) = c("model", "pr", "auc")
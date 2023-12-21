# Set response variable
Y = "dementia"

# 1) Model setups
target.train = as.h2o(snp_df_full)
target.train[,Y] = as.factor(target.train[,Y])

# Run two-step linear models
# AD indsig SNPs
X = c(AD.indsig.combine.id)
SNP.lasso.feature = h2o.glm(training_frame = target.train, x = X, y = Y, 
                            balance_classes = T, model_id = "LASSO features", 
                            nfolds = 5, fold_assignment = "Stratified", 
                            offset_column = "offset_demo",
                            seed = seed_num, family = "binomial", 
                            alpha = alpha_feature, lambda_search = TRUE, standardize = T, 
                            keep_cross_validation_models = F)
snps_AD_indsig = SNP.lasso.feature@model$variable_importances %>% 
  as.data.frame() %>% filter(scaled_importance > 0) %>% pull(variable) 
X = c(snps_AD_indsig)
SNP.AD.indsig = h2o.glm(training_frame = target.train, x = X, y = Y, 
                        balance_classes = T, model_id = "LASSO select", 
                        nfolds = 5, fold_assignment = "Stratified", 
                        offset_column = "offset_demo",
                        seed = seed_num, family = "binomial", 
                        alpha = 0, lambda_search = TRUE, standardize = T, 
                        keep_cross_validation_models = F)
# AD lead SNPs
# Step 1: feature selection
X = c(AD.lead.combine.id)
SNP.lasso.feature = h2o.glm(training_frame = target.train, x = X, y = Y, 
                            balance_classes = T, model_id = "LASSO features", 
                            nfolds = 5, fold_assignment = "Stratified", 
                            offset_column = "offset_demo",
                            seed = seed_num, family = "binomial", 
                            alpha = alpha_feature, lambda_search = TRUE, standardize = T, 
                            keep_cross_validation_models = F)
snps_AD_lead = SNP.lasso.feature@model$variable_importances %>% 
  as.data.frame() %>% filter(scaled_importance > 0) %>% pull(variable) 
X = c(snps_AD_lead)
SNP.AD.lead = h2o.glm(training_frame = target.train, x = X, y = Y, 
                      balance_classes = T, model_id = "LASSO select", 
                      nfolds = 5, fold_assignment = "Stratified", 
                      offset_column = "offset_demo",
                      seed = seed_num, family = "binomial", 
                      alpha = 0, lambda_search = TRUE, standardize = T, 
                      keep_cross_validation_models = F)
# AD map SNPs
X = c(AD.map.combine.id)
SNP.lasso.feature = h2o.glm(training_frame = target.train, x = X, y = Y, 
                            balance_classes = T, model_id = "LASSO features", 
                            nfolds = 5, fold_assignment = "Stratified", 
                            offset_column = "offset_demo",
                            seed = seed_num, family = "binomial", 
                            alpha = alpha_feature, lambda_search = TRUE, standardize = T, 
                            keep_cross_validation_models = F)
snps_AD_map = SNP.lasso.feature@model$variable_importances %>% 
  as.data.frame() %>% filter(scaled_importance > 0) %>% pull(variable) 
X = c(snps_AD_map)
SNP.AD.map = h2o.glm(training_frame = target.train, x = X, y = Y, 
                     balance_classes = T, model_id = "LASSO select", 
                     nfolds = 5, fold_assignment = "Stratified", 
                     offset_column = "offset_demo",
                     seed = seed_num, family = "binomial", 
                     alpha = 0, lambda_search = TRUE, standardize = T, 
                     keep_cross_validation_models = F)
# AD map SNPs CADD
X = c(AD.map.CADD.combine.id)
SNP.lasso.feature = h2o.glm(training_frame = target.train, x = X, y = Y, 
                            balance_classes = T, model_id = "LASSO features", 
                            nfolds = 5, fold_assignment = "Stratified", 
                            offset_column = "offset_demo",
                            seed = seed_num, family = "binomial", 
                            alpha = alpha_feature, lambda_search = TRUE, standardize = T, 
                            keep_cross_validation_models = F)
snps_AD_map_CADD = SNP.lasso.feature@model$variable_importances %>% 
  as.data.frame() %>% filter(scaled_importance > 0) %>% pull(variable) 
X = c(snps_AD_map_CADD)
SNP.AD.map.CADD = h2o.glm(training_frame = target.train, x = X, y = Y, 
                          balance_classes = T, model_id = "LASSO select", 
                          nfolds = 5, fold_assignment = "Stratified", 
                          offset_column = "offset_demo",
                          seed = seed_num, family = "binomial", 
                          alpha = 0, lambda_search = TRUE, standardize = T, 
                          keep_cross_validation_models = F)

# All indsig SNPs
# Step 1: feature selection
X = c(all.indsig.combine.id)
SNP.lasso.feature = h2o.glm(training_frame = target.train, x = X, y = Y, 
                            balance_classes = T, model_id = "LASSO features", 
                            nfolds = 5, fold_assignment = "Stratified", 
                            offset_column = "offset_demo",
                            seed = seed_num, family = "binomial", 
                            alpha = alpha_feature, lambda_search = TRUE, standardize = T, 
                            keep_cross_validation_models = F)
snps_all_indsig = SNP.lasso.feature@model$variable_importances %>% 
  as.data.frame() %>% filter(scaled_importance > 0) %>% pull(variable) 
X = c(snps_all_indsig)
SNP.all.indsig = h2o.glm(training_frame = target.train, x = X, y = Y, 
                         balance_classes = T, model_id = "LASSO select", 
                         nfolds = 5, fold_assignment = "Stratified", 
                         offset_column = "offset_demo",
                         seed = seed_num, family = "binomial", 
                         alpha = 0, lambda_search = TRUE, standardize = T, 
                         keep_cross_validation_models = F)
# All lead SNPs
# Step 1: feature selection
X = c(all.lead.combine.id)
SNP.lasso.feature = h2o.glm(training_frame = target.train, x = X, y = Y, 
                            balance_classes = T, model_id = "LASSO features", 
                            nfolds = 5, fold_assignment = "Stratified", 
                            offset_column = "offset_demo",
                            seed = seed_num, family = "binomial", 
                            alpha = alpha_feature, lambda_search = TRUE, standardize = T, 
                            keep_cross_validation_models = F)
snps_all_lead = SNP.lasso.feature@model$variable_importances %>% 
  as.data.frame() %>% filter(scaled_importance > 0) %>% pull(variable) 
X = c(snps_all_lead)
SNP.all.lead = h2o.glm(training_frame = target.train, x = X, y = Y, 
                       balance_classes = T, model_id = "LASSO select", 
                       nfolds = 5, fold_assignment = "Stratified", 
                       offset_column = "offset_demo",
                       seed = seed_num, family = "binomial", 
                       alpha = 0, lambda_search = TRUE, standardize = T, 
                       keep_cross_validation_models = F)
# All mapped SNPs
X = c(all.map.combine.id)
SNP.lasso.feature = h2o.glm(training_frame = target.train, x = X, y = Y, 
                            balance_classes = T, model_id = "LASSO features", 
                            nfolds = 5, fold_assignment = "Stratified", 
                            offset_column = "offset_demo",
                            seed = seed_num, family = "binomial", 
                            alpha = alpha_feature, lambda_search = TRUE, standardize = T, 
                            keep_cross_validation_models = F)
snps_all_map = SNP.lasso.feature@model$variable_importances %>% 
  as.data.frame() %>% filter(scaled_importance > 0) %>% pull(variable) 
X = c(snps_all_map)
SNP.all.map = h2o.glm(training_frame = target.train, x = X, y = Y, 
                      balance_classes = T, model_id = "LASSO select", 
                      nfolds = 5, fold_assignment = "Stratified", 
                      offset_column = "offset_demo",
                      seed = seed_num, family = "binomial", 
                      alpha = 0, lambda_search = TRUE, standardize = T, 
                      keep_cross_validation_models = F)
# All mapped SNPs CADD
X = c(all.map.CADD.combine.id)
SNP.lasso.feature = h2o.glm(training_frame = target.train, x = X, y = Y, 
                            balance_classes = T, model_id = "LASSO features", 
                            nfolds = 5, fold_assignment = "Stratified", 
                            offset_column = "offset_demo",
                            seed = seed_num, family = "binomial", 
                            alpha = alpha_feature, lambda_search = TRUE, standardize = T, 
                            keep_cross_validation_models = F)
snps_all_map_CADD = SNP.lasso.feature@model$variable_importances %>% 
  as.data.frame() %>% filter(scaled_importance > 0) %>% pull(variable) 
X = c(snps_all_map_CADD)
SNP.all.map.CADD = h2o.glm(training_frame = target.train, x = X, y = Y, 
                      balance_classes = T, model_id = "LASSO select", 
                      nfolds = 5, fold_assignment = "Stratified", 
                      offset_column = "offset_demo",
                      seed = seed_num, family = "binomial", 
                      alpha = 0, lambda_search = TRUE, standardize = T, 
                      keep_cross_validation_models = F)

all_model_metric = c()
all_model_performance = c()
# Extract results
model_list = c("SNP.AD.indsig", "SNP.AD.lead", "SNP.AD.map", "SNP.AD.map.CADD", 
               "SNP.all.indsig", "SNP.all.lead", "SNP.all.map", "SNP.all.map.CADD")
for (i in 1:length(model_list)) {
  print(paste0("Model: ", model_list[i]))
  # Run extraction function
  extract_cv = extract_cv_results(get(model_list[i]))
  model_metric = t(extract_cv[[1]]) %>% as.data.frame() %>% 
    mutate(mark = model_list[i])
  model_performance = extract_cv[[2]] %>% as.data.frame() %>% 
    mutate(mark = model_list[i])
  if (i == 1) {
    all_model_metric = model_metric
    all_model_performance = model_performance
  } else {
    all_model_metric = rbind(all_model_metric, model_metric)
    all_model_performance = rbind(all_model_performance, model_performance)
  }
}
# Join results together
snp_model_metric = rbind(all_model_metric)
snp_model_performance = rbind(all_model_performance)
snp_model_summary = snp_model_metric %>% inner_join(snp_model_performance)


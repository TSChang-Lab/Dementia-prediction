rm(list = ls())
# lapply(paste('package:', names(sessionInfo()$otherPkgs), sep = ""),
#        detach, character.only = TRUE, unload = TRUE)
pacman::p_load(tidyverse, h2o, caret, pROC)
raw_data_path = "/Users/Mingzhou/Desktop/Projects/Dementia-prediction/data/"
output_path = "/Users/Mingzhou/Desktop/Projects/Dementia-prediction/output/"
# Source in useful functions
source("/Users/Mingzhou/Desktop/Projects/Dementia-prediction/code/functions.R")
# Load in PRS beta data
load(file = paste0(raw_data_path, "prs/mod/atlas_prs_final.rda"))
load(file = paste0(raw_data_path, "modeling/AFR/sample_AFR_full.rda"))

load(file = paste0(raw_data_path, "modeling/AFR/feature_PRS_full.rda"))
h2o.init(nthreads = -1)
seed_num = 134766
# Set response variable
Y = "dementia"


# get bootstrapping samples
for (k in 966:1000) {
  print(k)
  AFR_case = sample_AFR_full %>% filter(dementia == 1)
  AFR_control = sample_AFR_full %>% filter(dementia == 0) %>% 
    sample_n(nrow(AFR_case)*3)
  sample_AFR_final = rbind(AFR_case, AFR_control) %>% 
    select(UniqueSampleId, age, female, dementia, PC1, PC2, PC3, PC4, PC5, PC6, 
           PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18,
           PC19, PC20) 
  # Join geno data
  sample_AFR_prs = sample_AFR_final %>% 
    mutate(UniqueSampleId = as.character(UniqueSampleId)) %>% 
    left_join(atlas_prs_final) %>% column_to_rownames(var = "UniqueSampleId") %>% 
    select(-c(gen_ancestry, APOE, e2count)) %>% drop_na()
  # fit a logistic regression with demographic sets only
  demo_glm = glm(dementia ~ age + female + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                   PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
                   PC16 + PC17 + PC18 + PC19 + PC20, data = sample_AFR_prs, 
                 family = binomial(link = "logit"))
  demo_glm_pred = predict(demo_glm, type = "link")
  sample_AFR_prs$offset_demo = exp(demo_glm_pred - coef(demo_glm)[1]) 
  # run PRS model
  # 1) Model setups
  target.train = as.h2o(sample_AFR_prs)
  target.train[,Y] = as.factor(target.train[,Y])
  # 2) Run model
  # AD PRS AFR indsig
  X = c("AD_PRS_AFR_indsig")
  PRS.AFR.indsig = h2o.glm(training_frame = target.train, x = X, y = Y, 
                           balance_classes = T, model_id = "AD PRS AFR indsig", 
                           nfolds = 5, fold_assignment = "Stratified", 
                           offset_column = "offset_demo",
                           seed = seed_num, family = "binomial", lambda = 0, 
                           standardize = T, keep_cross_validation_models = F)
  prs_AFR_indsig = PRS.AFR.indsig@model$variable_importances %>%
    as.data.frame() %>% filter(scaled_importance > 0) 
  coef_AFR_indsig = PRS.AFR.indsig@model$coefficients %>% as.data.frame() %>% 
    rownames_to_column() %>% rename("variable" = "rowname", "coef" = ".") %>% 
    filter(abs(coef) > 0) 
  feature_AFR_indsig = coef_AFR_indsig %>% 
    left_join(prs_AFR_indsig, by = "variable") %>% 
    mutate(bootstrap = k, model = "PRS_AFR_indsig")
  # AD PRS AFR map
  X = c("AD_PRS_AFR_map")
  PRS.AFR.map = h2o.glm(training_frame = target.train, x = X, y = Y, 
                           balance_classes = T, model_id = "AD PRS AFR map", 
                           nfolds = 5, fold_assignment = "Stratified", 
                           offset_column = "offset_demo",
                           seed = seed_num, family = "binomial", lambda = 0, 
                           standardize = T, keep_cross_validation_models = F)
  prs_AFR_map = PRS.AFR.map@model$variable_importances %>%
    as.data.frame() %>% filter(scaled_importance > 0) 
  coef_AFR_map = PRS.AFR.map@model$coefficients %>% as.data.frame() %>% 
    rownames_to_column() %>% rename("variable" = "rowname", "coef" = ".") %>% 
    filter(abs(coef) > 0) 
  feature_AFR_map = coef_AFR_map %>% 
    left_join(prs_AFR_map, by = "variable") %>% 
    mutate(bootstrap = k, model = "PRS_AFR_map")
  # e4count
  X = c("e4count")
  e4count = h2o.glm(training_frame = target.train, x = X, y = Y, 
                    balance_classes = T, model_id = "e4count", 
                    nfolds = 5, fold_assignment = "Stratified", 
                    offset_column = "offset_demo",
                    seed = seed_num, family = "binomial", lambda = 0, 
                    standardize = T, keep_cross_validation_models = F)
  var_e4count = e4count@model$variable_importances %>%
    as.data.frame() %>% filter(scaled_importance > 0) 
  coef_e4count = e4count@model$coefficients %>% as.data.frame() %>% 
    rownames_to_column() %>% rename("variable" = "rowname", "coef" = ".") %>% 
    filter(abs(coef) > 0) 
  feature_e4count = coef_e4count %>% 
    left_join(var_e4count, by = "variable") %>% 
    mutate(bootstrap = k, model = "e4count")
  
  feature_sum = rbind(feature_AFR_indsig, feature_AFR_map, feature_e4count)
  # save results
  if (k == 1) {
    feature_PRS_full = feature_sum
  } else {
    feature_PRS_full = rbind(feature_PRS_full, feature_sum)
  }
  # save results
  save(feature_PRS_full,
       file = paste0(raw_data_path, "modeling/AFR/feature_PRS_full.rda") )
}





rm(list = ls())
# lapply(paste('package:', names(sessionInfo()$otherPkgs), sep = ""),
#        detach, character.only = TRUE, unload = TRUE)
pacman::p_load(tidyverse, h2o, caret, pROC)
raw_data_path = "/Users/Mingzhou/Desktop/Projects/Dementia-prediction/data/"
output_path = "/Users/Mingzhou/Desktop/Projects/Dementia-prediction/output/"
# Source in useful functions
source("/Users/Mingzhou/Desktop/Projects/Dementia-prediction/code/functions.R")
# Load in geno beta data
load(file = paste0(raw_data_path, "modeling/freeze60k_geno_freq_full_SNP.rda"))
# Read in SNP.id data
ids_folder = paste0(raw_data_path, "modeling/snp_ids") 
ids_files = list.files(ids_folder, pattern = ".rda$")
for (i in ids_files) {
  load(paste0(raw_data_path, 'modeling/snp_ids/', i))
}

load(file = paste0(raw_data_path, "modeling/AFR/sample_AFR_full.rda"))
load(file = paste0(raw_data_path, "modeling/AFR/feature_map_full.rda"))
h2o.init(nthreads = -1)
seed_num = 134766
# Set response variable
Y = "dementia"
# Set parameters
alpha_feature = 0.005

# get bootstrapping samples
for (k in 901:1000) {
  print(k)
  AFR_case = sample_AFR_full %>% filter(dementia == 1)
  AFR_control = sample_AFR_full %>% filter(dementia == 0) %>% 
    sample_n(nrow(AFR_case)*3)
  sample_AFR_final = rbind(AFR_case, AFR_control) %>% 
    select(UniqueSampleId, age, female, dementia, PC1, PC2, PC3, PC4, PC5, PC6, 
           PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18,
           PC19, PC20) 
  # Join geno data
  snp_df_full = sample_AFR_final %>% 
    mutate(UniqueSampleId = as.character(UniqueSampleId)) %>% 
    left_join(genotype_df_fam, by = c("UniqueSampleId" = "IID")) %>% 
    select(-c(UniqueSampleId)) 
  # fit a logistic regression with demographic sets only
  demo_glm = glm(dementia ~ age + female + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                   PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
                   PC16 + PC17 + PC18 + PC19 + PC20, data = snp_df_full, 
                 family = binomial(link = "logit"))
  demo_glm_pred = predict(demo_glm, type = "link")
  snp_df_full$offset_demo = exp(demo_glm_pred - coef(demo_glm)[1]) 
  # run SNP model
  # 1) Model setups
  target.train = as.h2o(snp_df_full)
  target.train[,Y] = as.factor(target.train[,Y])
  # 2) Run model
  X = c(all.map.combine.id)
  SNP.lasso.feature = h2o.glm(training_frame = target.train, x = X, y = Y, 
                              balance_classes = T, model_id = "LASSO features", 
                              nfolds = 5, fold_assignment = "Stratified", 
                              offset_column = "offset_demo",
                              seed = seed_num, family = "binomial", 
                              alpha = alpha_feature, lambda_search = TRUE, 
                              standardize = T, 
                              keep_cross_validation_models = F)
  snps_all_map = SNP.lasso.feature@model$variable_importances %>%
    as.data.frame() %>% filter(scaled_importance > 0) 
  coef_all_map = SNP.lasso.feature@model$coefficients %>% as.data.frame() %>% 
    rownames_to_column() %>% rename("variable" = "rowname", "coef" = ".") %>% 
    filter(abs(coef) > 0) 
  feature_map = coef_all_map %>% left_join(snps_all_map, by = "variable") %>% 
    mutate(bootstrap = k)
  # save results
  if (k == 1) {
    feature_map_full = feature_map
  } else {
    feature_map_full = rbind(feature_map_full, feature_map)
  }
  # save results
  save(feature_map_full,
       file = paste0(raw_data_path, "modeling/AFR/feature_map_full.rda") )
}





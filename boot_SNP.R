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

load(file = paste0(raw_data_path, "modeling/EAS/sample_EAS_full.rda"))
h2o.init(nthreads = -1)
seed_num = 1347661
# Set parameters
alpha_feature = 0.005

load(file = paste0(raw_data_path, "modeling/EAS/SNP_bootstrap_full_150.rda") )

# get bootstrapping samples
for (k in 1038:1038) {
  print(k)
  EAS_case = sample_EAS_full %>% filter(dementia == 1)
  EAS_control = sample_EAS_full %>% filter(dementia == 0) %>% 
    sample_n(nrow(EAS_case)*3)
  sample_EAS_final = rbind(EAS_case, EAS_control) %>% 
    select(UniqueSampleId, age, female, dementia, PC1, PC2, PC3, PC4, PC5, PC6, 
           PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18,
           PC19, PC20) 
  # Join geno data
  snp_df_full = sample_EAS_final %>% 
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
  # run PRS model
  source("/Users/Mingzhou/Desktop/Projects/Dementia-prediction/code/SNP_modeling.R")
  SNP_bootstrap = snp_model_summary %>% mutate(bootstrap = k)
  # save results
  if (k == 1) {
    SNP_bootstrap_full = SNP_bootstrap
  } else {
    SNP_bootstrap_full = rbind(SNP_bootstrap_full, SNP_bootstrap)
  }
  # save results
  save(SNP_bootstrap_full, file = paste0(raw_data_path, "modeling/EAS/SNP_bootstrap_full_150.rda") )
}




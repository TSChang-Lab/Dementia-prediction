rm(list = ls())
# lapply(paste('package:', names(sessionInfo()$otherPkgs), sep = ""),
#        detach, character.only = TRUE, unload = TRUE)
pacman::p_load(tidyverse, h2o, caret, pROC)
raw_data_path = "/Users/Mingzhou/Desktop/Projects/Dementia-prediction/data/"
output_path = "/Users/Mingzhou/Desktop/Projects/Dementia-prediction/output/"
# Source in useful functions
source("/Users/Mingzhou/Desktop/Projects/Dementia-prediction/code/functions.R")
# Load in PRS + geno beta data
load(file = paste0(raw_data_path, "prs/mod/atlas_prs_final.rda"))
# Read in SNP.id data
ids_folder = paste0(raw_data_path, "modeling/snp_ids") 
ids_files = list.files(ids_folder, pattern = ".rda$")
for (i in ids_files) {
  load(paste0(raw_data_path, 'modeling/snp_ids/', i))
}

load(file = paste0(raw_data_path, "modeling/AMR/sample_AMR_full.rda"))
# load(file = paste0(raw_data_path, "modeling/AMR/PRS_bootstrap_full.rda"))
h2o.init(nthreads = -1)
seed_num = 134766

# get bootstrapping samples
for (k in 1:1000) {
  print(k)
  AMR_case = sample_AMR_full %>% filter(dementia == 1)
  AMR_control = sample_AMR_full %>% filter(dementia == 0) %>% 
    sample_n(nrow(AMR_case)*3)
  sample_AMR_final = rbind(AMR_case, AMR_control) %>% 
    select(UniqueSampleId, age, female, dementia, PC1, PC2, PC3, PC4, PC5, PC6, 
           PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18,
           PC19, PC20) 
  # Join PRSs
  sample_AMR_prs = sample_AMR_final %>% 
    mutate(UniqueSampleId = as.character(UniqueSampleId)) %>% 
    left_join(atlas_prs_final) %>% column_to_rownames(var = "UniqueSampleId") %>% 
    select(-c(gen_ancestry, APOE, e2count)) %>% drop_na()
  # fit a logistic regression with demographic sets only
  demo_glm = glm(dementia ~ age + female + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                   PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
                   PC16 + PC17 + PC18 + PC19 + PC20, data = sample_AMR_prs, 
                 family = binomial(link = "logit"))
  demo_glm_pred = predict(demo_glm, type = "link")
  sample_AMR_prs$offset_demo = exp(demo_glm_pred - coef(demo_glm)[1]) 
  # run PRS model
  source("/Users/Mingzhou/Desktop/Projects/Dementia-prediction/code/PRS_modeling.R")
  PRS_bootstrap = prs_model_summary %>% mutate(bootstrap = k)
  # save results
  if (k == 1) {
    PRS_bootstrap_full = PRS_bootstrap
  } else {
    PRS_bootstrap_full = rbind(PRS_bootstrap_full, PRS_bootstrap)
  }
  # save results
  save(PRS_bootstrap_full, file = paste0(raw_data_path, "modeling/AMR/PRS_bootstrap_full.rda"))
}



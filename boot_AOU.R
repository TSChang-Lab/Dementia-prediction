rm(list = ls())
# lapply(paste('package:', names(sessionInfo()$otherPkgs), sep = ""),
#        detach, character.only = TRUE, unload = TRUE)
pacman::p_load(tidyverse, h2o, caret, pROC, PRROC)
raw_data_path = "/Users/Mingzhou/Desktop/Projects/Dementia-prediction/data/"
output_path = "/Users/Mingzhou/Desktop/Projects/Dementia-prediction/output/"
# Source in useful functions
source("/Users/Mingzhou/Desktop/Projects/Dementia-prediction/code/functions.R")
# Load in AOU data
load(file = paste0(raw_data_path, "modeling/AOU/aou.RData"))

# AFR sample
load(file = paste0(raw_data_path, "AOU/AOU_AFR_PRS.rda"))
AOU_AFR_short = AOU_AFR_PRS %>% 
  select(PatientID, age, female, dementia, record_length, enc_per_yr,
         paste("PC", 1:16, sep = ""), AD_PRS_AFR_indsig, 
         AD_PRS_AFR_map, e4count) %>% drop_na() 
dim(AOU_AFR_short) # (2467, 25) (3957,25)

extract_columns = setdiff(coef_map_summary$variable, "Intercept")
load(file = paste0(raw_data_path, "AOU/aou_geno_freq_full_SNP.rda"))
aou_geno_freq_short = aou_geno_freq_full_SNP %>% 
  select(PatientID, all_of(extract_columns)) %>% drop_na() %>% 
  mutate(Intercept = 1)
dim(aou_geno_freq_short) # dim = (3591, 17) (7836,12)

# combine with SNP columns
sample_AFR_AOU_full = AOU_AFR_short %>% inner_join(aou_geno_freq_short)
dim(sample_AFR_AOU_full) # (3955,36)

load(file = paste0(raw_data_path, "modeling/AOU/AFR_AOU_bootstrap_full.rda"))

# get bootstrapping samples
for (k in 1383:1383) {
  print(k)
  AFR_case = sample_AFR_AOU_full %>% filter(dementia == 1 & age >= 60) %>% 
    sample_n(163, replace = FALSE) # 163 = 181*0.9 of cases
  AFR_control = sample_AFR_AOU_full %>% filter(dementia == 0 & age >= 65) 
  AFR_AOU_df = rbind(AFR_case, AFR_control) %>% as.data.frame()
                           
  sample_AFR_AOU = matchit(dementia ~ age + female, data = AFR_AOU_df, 
                           method = "nearest", ratio = 3) %>% match.data()
  # Regress out demographics
  demo_glm_aou = glm(dementia ~ age + female + PC1 + PC2 + PC3 + PC4 + PC5 +
                       PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 +
                       PC14 + PC15 + PC16, data = sample_AFR_AOU,
                     family = binomial(link = "logit"))
  demo_glm_aou_pred = predict(demo_glm_aou, type = "link")                       
  sample_AFR_AOU$offset_demo = exp(demo_glm_aou_pred - coef(demo_glm_aou)[1])    
  # test for models
  source("/Users/Mingzhou/Desktop/Projects/Dementia-prediction/code/AOU_modeling.R")
  AOU_bootstrap = model_summary %>% mutate(bootstrap = k)
  # save results
  if (k == 1) {
    AOU_bootstrap_full = AOU_bootstrap
  } else {
    AOU_bootstrap_full = rbind(AOU_bootstrap_full, AOU_bootstrap)
  }
  # save results
  save(AOU_bootstrap_full, file = paste0(raw_data_path, "modeling/AOU/AFR_AOU_bootstrap_full.rda"))
}














# Set response variable
Y = "dementia"

# 1) Model setups
target.train = as.h2o(sample_AMR_prs)
target.train[,Y] = as.factor(target.train[,Y])

# 2) Single PRS models
#========= single PRS model =========
# AD PRS EUR
X = c("AD_PRS_EUR_indsig")
PRS.EUR.indsig = h2o.glm(training_frame = target.train, x = X, y = Y, 
                         balance_classes = T, model_id = "AD PRS EUR indsig", 
                         nfolds = 5, fold_assignment = "Stratified", 
                         offset_column = "offset_demo",
                         seed = seed_num, family = "binomial", lambda = 0, 
                         standardize = T, keep_cross_validation_models = F)
X = c("AD_PRS_EUR_lead")
PRS.EUR.lead = h2o.glm(training_frame = target.train, x = X, y = Y, 
                       balance_classes = T, model_id = "AD PRS EUR lead", 
                       nfolds = 5, fold_assignment = "Stratified", 
                       offset_column = "offset_demo",
                       seed = seed_num, family = "binomial", lambda = 0, 
                       standardize = T, keep_cross_validation_models = F)
X = c("AD_PRS_EUR_map")
PRS.EUR.map = h2o.glm(training_frame = target.train, x = X, y = Y, 
                      balance_classes = T, model_id = "AD PRS EUR map", 
                      nfolds = 5, fold_assignment = "Stratified",
                      offset_column = "offset_demo",
                      seed = seed_num, family = "binomial", lambda = 0, 
                      standardize = T, keep_cross_validation_models = F)
# AD PRS AFR
X = c("AD_PRS_AFR_indsig")
PRS.AFR.indsig = h2o.glm(training_frame = target.train, x = X, y = Y, 
                         balance_classes = T, model_id = "AD PRS AFR indsig", 
                         nfolds = 5, fold_assignment = "Stratified", 
                         offset_column = "offset_demo",
                         seed = seed_num, family = "binomial", lambda = 0, 
                         standardize = T, keep_cross_validation_models = F)
X = c("AD_PRS_AFR_lead")
PRS.AFR.lead = h2o.glm(training_frame = target.train, x = X, y = Y, 
                       balance_classes = T, model_id = "AD PRS AFR lead", 
                       nfolds = 5, fold_assignment = "Stratified", 
                       offset_column = "offset_demo",
                       seed = seed_num, family = "binomial", lambda = 0, 
                       standardize = T, keep_cross_validation_models = F)
X = c("AD_PRS_AFR_map")
PRS.AFR.map = h2o.glm(training_frame = target.train, x = X, y = Y, 
                      balance_classes = T, model_id = "AD PRS AFR map", 
                      nfolds = 5, fold_assignment = "Stratified", 
                      offset_column = "offset_demo",
                      seed = seed_num, family = "binomial", lambda = 0, 
                      standardize = T, keep_cross_validation_models = F)
# AD PRS Trans 
X = c("AD_PRS_Trans_indsig")
PRS.Trans.indsig = h2o.glm(training_frame = target.train, x = X, y = Y, 
                           balance_classes = T, model_id = "AD PRS Trans indsig", 
                           nfolds = 5, fold_assignment = "Stratified", 
                           offset_column = "offset_demo",
                           seed = seed_num, family = "binomial", lambda = 0, 
                           standardize = T, keep_cross_validation_models = F)
X = c("AD_PRS_Trans_lead")
PRS.Trans.lead = h2o.glm(training_frame = target.train, x = X, y = Y, 
                         balance_classes = T, model_id = "AD PRS Trans lead", 
                         nfolds = 5, fold_assignment = "Stratified", 
                         offset_column = "offset_demo",
                         seed = seed_num, family = "binomial", lambda = 0, 
                         standardize = T, keep_cross_validation_models = F)
X = c("AD_PRS_Trans_map")
PRS.Trans.map = h2o.glm(training_frame = target.train, x = X, y = Y, 
                        balance_classes = T, model_id = "AD PRS Trans map", 
                        nfolds = 5, fold_assignment = "Stratified", 
                        offset_column = "offset_demo",
                        seed = seed_num, family = "binomial", lambda = 0, 
                        standardize = T, keep_cross_validation_models = F)
# e4count
X = c("e4count")
e4count = h2o.glm(training_frame = target.train, x = X, y = Y, 
                  balance_classes = T, model_id = "e4count", 
                  nfolds = 5, fold_assignment = "Stratified", 
                  offset_column = "offset_demo",
                  seed = seed_num, family = "binomial", lambda = 0, 
                  standardize = T, keep_cross_validation_models = F)

all_model_metric = c()
all_model_performance = c()
# Extract results
model_list = c("e4count", "PRS.EUR.indsig", "PRS.EUR.lead", "PRS.EUR.map", 
               "PRS.AFR.indsig", "PRS.AFR.lead", "PRS.AFR.map", 
               "PRS.Trans.indsig", "PRS.Trans.lead", "PRS.Trans.map")
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
single_prs_model_metric = rbind(all_model_metric)
single_prs_model_performance = rbind(all_model_performance)
single_prs_model_summary = single_prs_model_metric %>% 
  inner_join(single_prs_model_performance)

# 3) Multiple PRS model
#========= multiple PRS model =========
# AD PRS combinations
X = c("AD_PRS_EUR_indsig", "AD_PRS_Trans_indsig", "AD_PRS_AFR_indsig")
PRS.AD.indsig.combine = h2o.glm(training_frame = target.train, x = X, y = Y, 
                                balance_classes = T, 
                                model_id = "AD PRS ind.sig combine", 
                                nfolds = 5, fold_assignment = "Stratified", 
                                offset_column = "offset_demo",
                                seed = seed_num, family = "binomial", lambda = 0, 
                                standardize = T, keep_cross_validation_models = F)
X = c("AD_PRS_EUR_lead", "AD_PRS_Trans_lead", "AD_PRS_AFR_lead")
PRS.AD.lead.combine = h2o.glm(training_frame = target.train, x = X, y = Y, 
                              balance_classes = T, 
                              model_id = "AD PRS lead combine", 
                              nfolds = 5, fold_assignment = "Stratified", 
                              offset_column = "offset_demo",
                              seed = seed_num, family = "binomial", lambda = 0, 
                              standardize = T, keep_cross_validation_models = F)
X = c("AD_PRS_EUR_map", "AD_PRS_Trans_map", "AD_PRS_AFR_map")
PRS.AD.map.combine = h2o.glm(training_frame = target.train, x = X, y = Y, 
                             balance_classes = T, 
                             model_id = "AD PRS map combine", 
                             nfolds = 5, fold_assignment = "Stratified", 
                             offset_column = "offset_demo",
                             seed = seed_num, family = "binomial", lambda = 0, 
                             standardize = T, keep_cross_validation_models = F)
# AD + neuro combinations
X = c("AD_PRS_EUR_indsig", "AD_PRS_Trans_indsig", "AD_PRS_AFR_indsig", 
      "PD_PRS_indsig", "LBD_PRS_indsig", "PSP_PRS_indsig", 
      "Stroke_PRS_indsig")
PRS.neuro.indsig.combine = h2o.glm(training_frame = target.train, x = X, y = Y, 
                                   balance_classes = T, 
                                   model_id = "Neuro PRS ind.sig combine", 
                                   nfolds = 5, fold_assignment = "Stratified", 
                                   offset_column = "offset_demo",
                                   seed = seed_num, family = "binomial", lambda = 0, 
                                   standardize = T, keep_cross_validation_models = F)
X = c("AD_PRS_EUR_lead", "AD_PRS_Trans_lead", "AD_PRS_AFR_lead", 
      "PD_PRS_lead", "LBD_PRS_lead", "PSP_PRS_lead", 
      "Stroke_PRS_lead")
PRS.neuro.lead.combine = h2o.glm(training_frame = target.train, x = X, y = Y, 
                                 balance_classes = T, 
                                 model_id = "Neuro PRS lead combine", 
                                 nfolds = 5, fold_assignment = "Stratified", 
                                 offset_column = "offset_demo",
                                 seed = seed_num, family = "binomial", lambda = 0, 
                                 standardize = T, keep_cross_validation_models = F)
X = c("AD_PRS_EUR_map", "AD_PRS_Trans_map", "AD_PRS_AFR_map", 
      "PD_PRS_map", "LBD_PRS_map", "PSP_PRS_map", "Stroke_PRS_map")
PRS.neuro.map.combine = h2o.glm(training_frame = target.train, x = X, y = Y, 
                                balance_classes = T,
                                model_id = "Neuro PRS map combine", 
                                nfolds = 5, fold_assignment = "Stratified", 
                                offset_column = "offset_demo",
                                seed = seed_num, family = "binomial", lambda = 0, 
                                standardize = T, keep_cross_validation_models = F)

all_model_metric = c()
all_model_performance = c()
# Extract results
model_list = c("PRS.AD.indsig.combine", "PRS.AD.lead.combine", 
               "PRS.AD.map.combine", "PRS.neuro.indsig.combine", 
               "PRS.neuro.lead.combine", "PRS.neuro.map.combine")
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
multi_prs_model_metric = rbind(all_model_metric)
multi_prs_model_performance = rbind(all_model_performance)
multi_prs_model_summary = multi_prs_model_metric %>% 
  inner_join(multi_prs_model_performance)

prs_model_summary = rbind(single_prs_model_summary, multi_prs_model_summary) %>% 
  as.data.frame()


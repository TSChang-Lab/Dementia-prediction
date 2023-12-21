raw_data_path = "/Users/Mingzhou/Desktop/Projects/Dementia-prediction/data/"

eqtl_EUR = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_AD_EUR/eqtl.txt"), 
                      header = T, sep = "\t", fill = T) %>% 
  select(uniqID, symbol, alignedDirection) %>% unique()
eqtl_AFR = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_AD_AFR/eqtl.txt"), 
                      header = T, sep = "\t", fill = T) %>% 
  select(uniqID, symbol, alignedDirection) %>% unique()
eqtl_Trans = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_AD_trans/eqtl.txt"), 
                        header = T, sep = "\t", fill = T) %>% 
  select(uniqID, symbol, alignedDirection) %>% unique()
eqtl_PD = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_PD/eqtl.txt"), 
                     header = T, sep = "\t", fill = T) %>% 
  select(uniqID, symbol, alignedDirection) %>% unique()
eqtl_PSP = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_PSP/eqtl.txt"), 
                      header = T, sep = "\t", fill = T) %>% 
  select(uniqID, symbol, alignedDirection) %>% unique()
eqtl_LBD = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_LBD/eqtl.txt"), 
                      header = T, sep = "\t", fill = T) %>% 
  select(uniqID, symbol, alignedDirection) %>% unique()
eqtl_STROKE = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_STROKE/eqtl.txt"), 
                      header = T, sep = "\t", fill = T) %>% 
  select(uniqID, symbol, alignedDirection) %>% unique()


AD_EUR_ci = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_AD_EUR/ci.txt"), 
                       header = T, sep = "\t", fill = T) %>% 
  select(SNPs, genes) %>% unique()
AD_AFR_ci = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_AD_AFR/ci.txt"), 
                       header = T, sep = "\t", fill = T) %>% 
  select(SNPs, genes) %>% unique()
AD_trans_ci = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_AD_trans/ci.txt"), 
                         header = T, sep = "\t", fill = T) %>% 
  select(SNPs, genes) %>% unique()
LBD_ci = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_LBD/ci.txt"), 
                    header = T, sep = "\t", fill = T) %>% 
  select(SNPs, genes) %>% unique()
PD_ci = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_PD/ci.txt"), 
                   header = T, sep = "\t", fill = T) %>% 
  select(SNPs, genes) %>% unique()
PSP_ci = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_PSP/ci.txt"), 
                    header = T, sep = "\t", fill = T) %>% 
  select(SNPs, genes) %>% unique()
STROKE_ci = read.table(file = paste0(raw_data_path, "FUMA_outputs/FUMA_STROKE/ci.txt"), 
                       header = T, sep = "\t", fill = T) %>% 
  select(SNPs, genes) %>% unique()
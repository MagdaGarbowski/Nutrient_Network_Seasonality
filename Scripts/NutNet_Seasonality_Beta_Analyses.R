
# Analyses code associated with:
# "Nutrient enrichment alters seasonal beta diversity in global grasslands" 
# Garbowski et al. 2023 Journal of Ecology 
# FigShare link for data: 

# libraries
library(brms)

# load data 
beta_data <- read.csv("/Users/MagdaGarbowski 1/Nutrient_Network_Seasonality/Generated_Data/seasonalbeta_observed_null_deviation.csv")
turnover_nestedness <- read.csv("/Users/MagdaGarbowski 1/Nutrient_Network_Seasonality/Generated_Data/turnover_nestedness.csv")
beta_climate <- read.csv("/Users/MagdaGarbowski 1/Nutrient_Network_Seasonality/Generated_Data/seasonalbeta_climate.csv")

#############################  Climate models #####################################################

beta_dev_temp_var_bm_ss_years <- brm(z_score ~ temp_CV, 
                                     data = beta_climate, chains = 4, cores = 4,
                                     warmup = 2500, iter = 3500,
                                     control = list(adapt_delta = 0.9))

summary(beta_dev_temp_var_bm_ss_years)

beta_dev_precip_var_bm_ss_years <- brm(z_score ~ precip_CV, 
                                       data = beta_climate, warmup = 2500, iter = 3500, 
                                       control = list(adapt_delta = 0.9))

summary(beta_dev_precip_var_bm_ss_years)
  
##########################  Dissimilarity models ##################################################

# deviation (z-score)
deviation <- brm(z_score ~ trt + (1|site_code/year/block),
                 data = beta_data, chains = 4, cores = 4,
                 warmup = 2500, iter = 3500, control = list(adapt_delta = 0.99))
summary(deviation)
hypothesis(deviation, "Intercept = (Intercept + trtNPK)", alpha = 0.1)

# nestedness
nestedness <- brm(nestedness ~ trt + (1|site_code/year/block),
                 data = turnover_nestedness, chains = 4, cores = 4,
                 warmup = 2500, iter = 3500, control = list(adapt_delta = 0.99))
summary(nestedness)
hypothesis(nestedness, "Intercept = (Intercept + trtNPK)", alpha = 0.1)

# turnover
turnover <- brm(turnover ~ trt + (1|site_code/year/block),
                  data = turnover_nestedness, chains = 4, cores = 4,
                warmup = 2500, iter = 3500, control = list(adapt_delta = 0.99))
summary(turnover)
hypothesis(turnover, "Intercept = (Intercept + trtNPK)", alpha = 0.1)





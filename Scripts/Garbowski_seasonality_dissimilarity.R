# code associate with "Temperature seasonality and nutrient enrichment drive intra-annual community turnover in global grasslands" 
# obtain seasonal beta observed, null, and deviation (i.e., z-score) estimates

setwd("/Users/MagdaGarbowski 1/Nutrient_Network_Seasonality/")
source("Scripts/Garbowski_seasonality_dissimilarity_functions.R")
library(vegan)

# ------------------------------- data prep  --------------------------------

data <- read.csv("/Users/MagdaGarbowski 1/NutrientNetwork_Seasonality/Data/prepped_data.csv")
seasonal_site_splits <- split(data, data$site_code)

# prep data 
dat_prepped_out <- lapply(seasonal_site_splits, dat_prep)
dat_prepped_out_2 <- unlist(unlist(dat_prepped_out, recursive = FALSE), recursive = FALSE)

# drop communities with less than 4 species - cannot permute these (835 vs. 888)
dat_prepped_out_3 <- dat_prepped_out_2[sapply(dat_prepped_out_2, function(x) length(x) >3)]

# apply early vs. late to communities of 4+ species 
dat_prepped_out_4 <- lapply(dat_prepped_out_3, early_vs_late_function )

# create list of pairwise combinations from list of communities 
combos_obs <- combn(dat_prepped_out_4, 2,simplify = FALSE)

# return list where plot id in list 1 = plot id in list 2 
combos_obs_out <- lapply(combos_obs, combos_function)
combos_obs_out_2 <- combos_obs_out[sapply(combos_obs_out, function(x) length(x) >1)]

# add missing species from one of the two sampling points 
combos_obs_out_3 <- lapply(combos_obs_out_2, missing_sps_function)

# extract "kept" plots  
plots_kept <- data.frame(plot_id = do.call(rbind, lapply(combos_obs_out_3, function(x) plot_id = x[[1]]$plot_id)))
kept_df <- as.data.frame(t(sapply(strsplit(plots_kept$plot_id, "_|\\_"), "[")))
kept_df$site_code <- paste(kept_df$V1, kept_df$V2, sep = ".")
colnames(kept_df)<-c("site_code", "trt", "year", "block", "plot")

# drop year/site combos that don't have matches: chilcas 2017, 2018 and sevi 2015 
combos_obs_out_4 <- lapply(combos_obs_out_3, function (list) 
  if(!grepl("chilcas.ar_Control_2017|chilcas.ar_Control_2018|chilcas.ar_NPK_2017|chilcas.ar_NPK_2018|sevi.us_Control_2015|sevi.us_NPK_2015",list[[1]]$plot_id[1])) {
    return (list)})

combos_obs_out_4 <- combos_obs_out_4[sapply(combos_obs_out_4, function(x) length(x) >1)]

# ----------------------------- bray observed -----------------------------------

bray_obs_out <- do.call(rbind, lapply(combos_obs_out_4, bray_obs_function))

# ------------------------- permutations for null bray  -------------------------

# create 100 permuted communities 
# drop year/site combos that don't have matches as above: chilcas 2017, 2018 and sevi 2015 
permutations_out <-  lapply(dat_prepped_out_3, permutations_function, 100)
permutations_out_1 <- lapply(permutations_out, function (x) 
  if(!grepl("chilcas.ar_Control_2017|chilcas.ar_Control_2018|chilcas.ar_NPK_2017|chilcas.ar_NPK_2018|sevi.us_Control_2015|sevi.us_NPK_2015",x$plot_id[1])) {return (x)})
permutations_out_1 <- permutations_out_1[sapply(permutations_out_1, function(x) length(x) >1)]

# create "paired" datasets to get beta-dissimilarity on permuted communities 
combos <- combn(permutations_out_1, 2, simplify = FALSE)
combos_out <- lapply(combos, combos_function)
combos_out_2 <- combos_out[sapply(combos_out, function(x) length(x) > 1)]
combos_out_3 <- lapply(combos_out_2, missing_sps_function)

# get bray on permuted communities
bray_null_out <- lapply(combos_out_3, bray_null_function)
bray_null_out_df <- do.call(rbind, bray_null_out)

# get null mean and sd from permutations
bray_null_mean_sd <- do.call(rbind, lapply(bray_null_out, null_function))
bray_null_obs <- merge(bray_null_mean_sd, bray_obs_out, by = "plot_id")
bray_null_obs[,2:4] <- apply(bray_null_obs[,2:4], 2, as.numeric)

# get z-scores
bray_null_obs$z_score <- (bray_null_obs$bray - bray_null_obs$bray_mean_null)/(bray_null_obs$bray_sd_null)

# add back column names for analyses
cols = as.data.frame(t(sapply(strsplit(bray_null_obs$plot_id, "_|\\_"), "[")))
colnames(cols)<-c("site_code", "trt", "year", "block", "plot")

bray_null_obs <- cbind(cols, bray_null_obs)

# ------------------------- turnover and nestedness ------------------------------

turnover_nestedness_function <- function(ls){
  df = rbind(ls[[1]], ls[[2]])
  df2 = df[3:length(df)]
  
  beta_part_out <- betapart::beta.pair.abund(df2, index.family = "bray")
  bray = beta_part_out$beta.bray
  turnover = beta_part_out$beta.bray.bal
  nestedness = beta_part_out$beta.bray.gra
  
  df3 = data.frame(index= c("bray", "turnover", "nestedness"),
                   value = rbind(bray, turnover, nestedness))
  
  df3$plot = rep(df$plot_id[1],3)
  return(df3)
}

turnover_nestedness_out <- do.call(rbind, lapply(combos_obs_out_4, turnover_nestedness_function))
cols <- as.data.frame(t(sapply(strsplit(turnover_nestedness_out$plot, "_|\\_"), "[")))
colnames(cols)<-c("site_code", "trt", "year", "block", "plot")
turnover_nestedness_out_2 <- cbind(turnover_nestedness_out, cols)
colnames(turnover_nestedness_out_2)[3] <- "plot_id"

turnover_nestedness_out_wide <- reshape(turnover_nestedness_out_2, 
                                        idvar = c("site_code", "plot_id", "block", "plot", "trt", "year"),
                                        timevar = "index", 
                                        direction = "wide")

colnames(turnover_nestedness_out_wide)[7:9] <- c("bray", "turnover", "nestedness")

# ------------------------- write datasets for analyses ------------------------------

write.csv(bray_null_obs, "/Users/MagdaGarbowski 1/Nutrient_Network_Seasonality/Generated_Data/beta_observed_null_deviation.csv")
write.csv(turnover_nestedness_out_wide, "/Users/MagdaGarbowski 1/Nutrient_Network_Seasonality/Generated_Data/turnover_nestedness.csv")

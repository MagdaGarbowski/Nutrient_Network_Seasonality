
# Garbowski -  "Temperature seasonality and nutrient enrichment drive intra-annual community turnover in global grasslands" 
# October 11, 2022

# functions take community data from two time periods each growing season (early, late) to calculate 
# observed, null, and deviations from expected Bray-Curtis dissimilarity

# general data prep
dat_prep <- function(df){
  df_splits <- split(df, list(df$trt, df$year, df$block, df$plot, df$sampling), drop = TRUE)
  prep_function <- function(df1){
    df2 <- df1[c("Taxon", "site_sampling_trt_yr_blk_plot", "max_cover")]
    df3 <- aggregate(list(max_cover = df2$max_cover), by = list(Taxon = df2$Taxon,
                                                                site_sampling_trt_yr_blk_plot = df2$site_sampling_trt_yr_blk_plot),
                     FUN = "sum")
    df4 <- unique(df3)
    df_wide <- reshape(df4, idvar = "site_sampling_trt_yr_blk_plot", timevar = "Taxon", direction = "wide")
    return(df_wide)
  }
  out = list(lapply(df_splits, prep_function))
  out
}

# set up early to late communities to be compared
early_vs_late_function <- function(df){
  df$plot_id = gsub("_Early_|_Late_", "_", df$site_sampling_trt_yr_blk_plot)
  df$sampling = ifelse(grepl("_Early_", df$site_sampling_trt_yr_blk_plot), "Early", "Late")
  return(df)
}

# combine appropriate plots for comparison (i.e., early vs late from a given year)
combos_function <- function(ls){
  list = ls
  if(list[[1]]$plot_id[1] == list[[2]]$plot_id[1]){
    list_out = list(list[[1]], list[[2]])
    return(list_out)
  }
} 

# add missing species to early vs. late communities and fill with "0" 
missing_sps_function <- function(list){
  mat_1 = list[[1]] # early community
  mat_2 = list[[2]] # late community 
  
  # get species names from both early and late communities 
  sps.names <- unique(c(colnames(mat_1)[2:(length(mat_1)-2)], colnames(mat_2)[2:(length(mat_2)-2)]))
  columns <- c("plot_id", "sampling", sps.names)
  
  # add missing columns (i.e., species) to matrix 1 
  missing_mat_1 <- setdiff(sps.names, names(mat_1))
  mat_1[missing_mat_1] <- 0 
  mat_1 <- mat_1[columns]
  
  # add missing columns (i.e., species) to matrix 2 
  missing_mat_2 <- setdiff(sps.names, names(mat_2))
  mat_2[missing_mat_2] <- 0 
  mat_2 <- mat_2[columns]
  
  list_out <- list(mat_1, mat_2)
  return(list_out)
}

# bray observed function 
bray_obs_function <- function(list){
  m1 = list[[1]]
  m2 = list[[2]]
  
  fun =  function(m_1, m_2){
    rbind(m_1, m_2)
  }
  
  df <-  as.data.frame(mapply(fun, m1, m2))
  
  df_2 <- data.frame(apply(df[3:length(df)], 2, function(x) as.numeric(as.character(x))))
  out = vegdist(df_2, method = "bray")
  out_dat = as.data.frame(t(c(out, df$plot_id[1])))
  colnames(out_dat) <- c("bray", "plot_id")
  return(out_dat)
}

# permutations function 
### permute within communities however many times we decide on
### get beta expected by comparing permuted rows early to permuted rows late

permutations_function <- function(df, n_permutations){
  empty_matrix <- data.frame(matrix(ncol = length(df)-1, nrow = n_permutations)) # create empty matrix of n permutations
  values <- as.numeric(as.vector(df[1,2:length(df)])) # get abundance values 
  
  out_matrix <- as.data.frame(t(apply(empty_matrix, 1, function(x){sample(values)}))) # randomly sample values into n communities 
  colnames(out_matrix) <- names(df)[2:length(df)] # add species names to matrix
  out_matrix = cbind(plot_id_sampling = as.character(df$site_sampling_trt_yr_blk_plot), out_matrix) # add plot names
  out_matrix$plot_id = gsub("_Early_|_Late_", "_", out_matrix$plot_id_sampling) # remove "early" and "late" from plot_id
  out_matrix$sampling = ifelse(grepl("_Early_", df$site_sampling_trt_yr_blk_plot), "Early", "Late") # create sampling column in matrix 
  return(out_matrix)
}

# calculate dissimilarity of rows in matrix 1 to rows in matrix 2 
bray_null_function <- function(list){
  m1 = list[[1]]
  m2 = list[[2]]
  
  fun =  function(m_1, m_2){
    rbind(m_1, m_2)
  }
  
  out <-  as.data.frame(mapply(fun, m1, m2))
  splits <- split(out,rep(1:100,each=2))
  
  veg_dist_function <- function(df){
    df2 = df[3:length(df)]
    df3 <- data.frame(apply(df2, 2, function(x) as.numeric(as.character(x))))
    out = vegdist(df3, method = "bray")
    out_dat = as.data.frame(t(c(out, df$plot_id[1])))
    colnames(out_dat) <- c("bray", "plot_id")
    return(out_dat)
  }
  
  out_all <- do.call(rbind, lapply(splits, veg_dist_function))
  return(out_all)
}

# null function to get z-scores for permutations

null_function <- function(df){
  df$bray <- as.numeric(as.character(df$bray))
  bray_mean_null = mean(df$bray) # gets mean of the expected bray values from randomizations 
  bray_sd_null = sd(df$bray)
  bray_null = as.data.frame(cbind(plot_id = df$plot_id[1], bray_mean_null, bray_sd_null))
  return(bray_null)
}

# get z-scores
bray_null_obs$z_score <- (bray_null_obs$bray - bray_null_obs$bray_mean_null)/(bray_null_obs$bray_sd_null)




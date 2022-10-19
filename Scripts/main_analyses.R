# This script is meant to be used within an HPC environment (i.e. slurm), where an argument indicating the specific
# score choice are provided when the script is ran. The argument should be provided as a numeric,
# ranging between 1:length(score_choices_useful). This argument is read by this script through
# "commandArgs(TRUE)". An example script is provided in template_sbatch_script_for_main_analyses.sh

# Assign the processed data directory path to data_folder_path
procdata_folder_path <- "processed/data/folder/path/"

# Assign the results directory path to res_folder_path
res_folder_path <- "results/folder/path/"

# Dataset names -----------------------------------------------------------

dataset_totint_names <- c("CSF_1", "CSF_2", "CSF_3", "CSF_4", "CSF_5", "CSF_6", "CSF_7",
                          "Dataset1", "Dataset2", "Dataset3", "Dataset4", "Dataset5",
                          "Group 1_Standards", "Group 2_Standards", "Group 3_Standards", "Group 4_Standards", "Group 6_Standards",
                          "Plasma_1", "Plasma_2", "Plasma_Ref_2", "Standard_AminoacidMix", "Standard_Mix_21Aug20", "Standard_Mix_R1", 
                          "Standard_Mix_R2", "Standard_Mix_R3", "STD_Mix2", "UDN_CSF_metab", "Urine_01", "Urine_2")

# Packages ----------------------------------------------------------------

if(!require(dplyr)){
  install.packages('dplyr', repos = "http://cran.r-project.org")
}
library(dplyr)

if(!require(magrittr)){
  install.packages('magrittr', repos = "http://cran.r-project.org")
}
library(magrittr)

if(!require(parallel)){
  install.packages('parallel', repos = "http://cran.r-project.org")
}
library(parallel)

if(!require(rlist)){
  install.packages('rlist', repos = "http://cran.r-project.org")
}

if(!require(gamlss.mx)){
  install.packages('gamlss.mx', repos = "http://cran.r-project.org")
}


# Score Info ---------------------------------------------------------------

score_choices <- c("Similarity Score", "Spectral Similarity Score", 
                   "Weighted Cosine Correlation", "Cosine Correlation",
                   "Stein Scott Similarity", "Stein Scott Similarity Nist", 
                   "Pearson Correlation", "Spearman Correlation",                   
                   "Kendall Tau Correlation", "Euclidean Distance", 
                   "Manhattan Distance", "Jaccard Distance",
                   "DWT Correlation", "DFT Correlation",
                   "Chebyshev Distance", "Squared Euclidean Distance",
                   "Fidelity Similarity", "Matusita Distance",
                   "Squared-chord Distance", "Harmonic mean Distance",
                   "Pearson Chi Squared Distance", "Neyman Chi Squared Distance",
                   "Probabilistic symmetric X2 Distance", "Topsoe Distance",
                   "Chernoff Distance", "Ruzicka Distance",
                   "Roberts Distance", "Motyka Distance",
                   "Canberra Distance", "Canberra Metric",
                   "Kulczynski 1 Distance", "Lorentzian Distance",
                   "Clark Distance", "Hellinger Distance",
                   "Whittaker index of association Distance", "Spectral Contrast Angle",
                   "Wave Hedges Distance", "Dice Similarity",
                   "Inner Product Distance", "Divergence Distance",
                   "Jensen Differences Distance", "Kumar Johnson Distance",
                   "Avg (L1, L8) Distance", "Vicis Wave Hadges Distance",
                   "Vicis-Symmetric X2 1 Distance", "Vicis-Symmetric X2 2 Distance",
                   "Vicis-Symmetric X2 3 Distance", "Max Symmetric Chi Squared Distance",
                   "Min Symmetric Chi Squared Distance", "Additive Symmetric Chi Squared",
                   "Battacharya Distance", "Generalized Ochai Index",
                   "Gower Distance", "Improved Square Root Cosine Similarity",
                   "Intersection Similarity", "J Divergence",
                   "Jensen Shannon Index", "K Divergence",
                   "VW6", "VW5", "VW4", "VW3", "VW2", "VW1",
                   "Taneja Divergence", "Symmetric Chi Squared Distance",
                   "Squared Chi Squared Distance", "Square Root Cosine Correlation",
                   "Sorensen Distance", "Minokowski 3 Distance",
                   "Minokowski 4 Distance", "Kumar Johnson Divergence",
                   "Kumar Hassebrook Similarity", "Kullback Leibler Divergence",
                   "Soergel Distance")
# useful scores idxs
useful_score_idx <- which(score_choices %in% 
                            c("Similarity Score", "Spectral Similarity Score", 
                              "Weighted Cosine Correlation", "Cosine Correlation",
                              "Stein Scott Similarity", "Stein Scott Similarity Nist", 
                              "Pearson Correlation", "Spearman Correlation",                   
                              
                              "DWT Correlation", "DFT Correlation",
                              "Chebyshev Distance", "Squared Euclidean Distance",
                              "Fidelity Similarity", 
                              
                              "Squared-chord Distance", 
                              
                              "Topsoe Distance",
                              
                              "Canberra Metric",
                              
                              "Hellinger Distance",
                              
                              "Dice Similarity",
                              
                              "Jensen Differences Distance", 
                              
                              "Battacharya Distance", 
                              
                              "Gower Distance", "Improved Square Root Cosine Similarity",
                              "Intersection Similarity", 
                              
                              "Squared Chi Squared Distance", "Square Root Cosine Correlation",
                              
                              "Minokowski 4 Distance", 
                              
                              "Kumar Hassebrook Similarity", 
                              
                              "Soergel Distance"))

score_choices_useful <- score_choices[useful_score_idx]

# Assumed all 'distance' metrics should be minimized and all 'similarity' metrics maximized.
score_directions <- c("bigger", "bigger",
                      "bigger", "bigger",
                      "bigger", "bigger",
                      "bigger", "bigger",
                      "bigger", "smaller",
                      "smaller", "smaller",
                      "bigger", "bigger",
                      "smaller", "smaller",
                      "bigger", "smaller",
                      "smaller", "smaller",
                      "smaller", "smaller",
                      "smaller", "smaller",
                      "smaller", "smaller",
                      "smaller", "smaller",
                      "smaller", "smaller",
                      "smaller", "smaller",
                      "smaller", "smaller",
                      "smaller", "smaller",
                      "smaller", "bigger",
                      "smaller", "smaller",
                      "smaller", "smaller",
                      "smaller", "smaller",
                      "smaller", "smaller",
                      "smaller", "smaller",
                      "smaller", "smaller",
                      "smaller", "smaller",
                      "smaller", "bigger",
                      "bigger", "smaller",
                      "smaller", "smaller",
                      "smaller", "smaller", "smaller", "smaller", "smaller", "smaller",
                      "smaller", "smaller",
                      "smaller", "bigger",
                      "smaller", "smaller",
                      "smaller", "smaller",
                      "bigger", "smaller",
                      "smaller")

score_directions_useful <- score_directions[useful_score_idx]

# Functions ---------------------------------------------------------------

# Define function to perform score-variant processing steps on each dataset.
dset_improc_score <- function(gcdat_full, score_choice, score_direction){
  
  gcdat <- gcdat_full
  
  if(score_direction == "bigger"){
    
    # Select the top scoring candidate from each unique spectrum
    gcdat_top <- gcdat %>% 
      dplyr::select(SampID_num, SampID_chr, `Sample name`, `Retention Time`, 
                    `Compound Name`, .data[[score_choice]], Truth.Annotation,
                    SpctCndtsCount, PoMatchCount, Interference, `Peak Height`, refavg_css) %>%
      dplyr::group_by(SampID_num) %>% 
      dplyr::slice_max(order_by = .data[[score_choice]], n = 1) %>%
      dplyr::ungroup()
    
    # Determine whether there are any duplicate entries
    iso_dupes <- gcdat_top %>% 
      dplyr::count(SampID_num) %>% 
      dplyr::filter(n > 1) %>% 
      dplyr::left_join(gcdat_top, by = "SampID_num")
    
    # Duplicate entries are typically 'isomers' of one another. We remove these from the data. 
    gcdat_top <- gcdat_top %>% dplyr::filter(!(SampID_num %in% iso_dupes$SampID_num))
    
    # Create new variable that measures the number of times a reference compound was the top match 
    temp <- gcdat_top %>%
      dplyr::count(`Compound Name`) %>% 
      dplyr::rename(MatchCount = n)
    gcdat_top <- dplyr::left_join(gcdat_top, temp, by = "Compound Name")
    rm(temp)
    
    # Create new variable that measures the difference between top 2 scores among potential matches
    temp <- gcdat %>% 
      dplyr::filter(!(SampID_num %in% iso_dupes$SampID_num), SpctCndtsCount > 1) %>%
      dplyr::select(SampID_num, SampID_chr, `Sample name`, `Retention Time`, 
                    `Compound Name`, .data[[score_choice]], Truth.Annotation,
                    SpctCndtsCount, PoMatchCount, Interference, `Peak Height`, refavg_css) %>%
      dplyr::group_by(SampID_num) %>%
      dplyr::slice_max(order_by = .data[[score_choice]], n = 2) %>%
      dplyr::slice(2) %>%
      dplyr::rename(ScoreNext = .data[[score_choice]]) %>%
      dplyr::select(SampID_num, ScoreNext)
    
    
    gcdat_top <- dplyr::left_join(gcdat_top, temp, by = "SampID_num")
    gcdat_top <- gcdat_top %>% dplyr::mutate(Top2check = ifelse(.data[[score_choice]] >= ScoreNext, 1, 0),
                                             ScoreDiff = ifelse(is.na(ScoreNext), .data[[score_choice]],
                                                                .data[[score_choice]] - ScoreNext)) 
    rm(temp)
    
  } else if(score_direction == "smaller"){
    
    # Select the top scoring candidate from each unique spectrum
    gcdat_top <- gcdat %>% 
      dplyr::select(SampID_num, SampID_chr, `Sample name`, `Retention Time`, 
                    `Compound Name`, .data[[score_choice]], Truth.Annotation,
                    SpctCndtsCount, PoMatchCount, Interference, `Peak Height`, refavg_css) %>%
      dplyr::group_by(SampID_num) %>% 
      dplyr::slice_min(order_by = .data[[score_choice]], n = 1) %>%
      dplyr::ungroup()
    
    # Determine whether there are any duplicate entries
    iso_dupes <- gcdat_top %>% 
      dplyr::count(SampID_num) %>% 
      dplyr::filter(n > 1) %>% 
      dplyr::left_join(gcdat_top)
    
    # Duplicate entries are typically 'isomers' of one another. We remove these from the data. 
    gcdat_top <- gcdat_top %>% dplyr::filter(!(SampID_num %in% iso_dupes$SampID_num))
    
    # Create new variable that measures the number of times a reference compound was the top match 
    temp <- gcdat_top %>%
      dplyr::count(`Compound Name`) %>% 
      dplyr::rename(MatchCount = n)
    gcdat_top <- dplyr::left_join(gcdat_top, temp)
    rm(temp)
    
    # Create new variable that measures the difference between top 2 scores among potential matches
    temp <- gcdat %>% 
      dplyr::filter(!(SampID_num %in% iso_dupes$SampID_num), SpctCndtsCount > 1) %>%
      dplyr::select(SampID_num, SampID_chr, `Sample name`, `Retention Time`, 
                    `Compound Name`, .data[[score_choice]], Truth.Annotation,
                    SpctCndtsCount, PoMatchCount, Interference, `Peak Height`, refavg_css) %>%
      dplyr::group_by(SampID_num) %>%
      dplyr::slice_min(order_by = .data[[score_choice]], n = 2) %>%
      dplyr::slice(2) %>%
      dplyr::rename(ScoreNext = .data[[score_choice]]) %>%
      dplyr::select(SampID_num, ScoreNext)
    
    
    gcdat_top <- dplyr::left_join(gcdat_top, temp)
    largest_score <- max(gcdat[[score_choice]])
    gcdat_top <- gcdat_top %>% dplyr::mutate(Top2check = ifelse(.data[[score_choice]] <= ScoreNext, 1, 0),
                                             ScoreDiff = ifelse(is.na(ScoreNext), .data[[score_choice]] - largest_score,
                                                                .data[[score_choice]] - ScoreNext)) 
    rm(temp)
    
  }
  
  # Toss out all observations with "Unknown"
  gct_nuk <- gcdat_top %>% dplyr::filter(Truth.Annotation != "Unknown")
  gct_nuk <- gct_nuk %>% dplyr::select(-ScoreNext, -Top2check)
  # Standardize all variables.
  gct_nuk <- gct_nuk %>% dplyr::mutate(SpctCndtsCount_z = scale(SpctCndtsCount),
                                       PoMatchCount_z = scale(PoMatchCount),
                                       PeakHeight_z = scale(`Peak Height`),
                                       Interference_z = scale(Interference),
                                       MatchCount_z = scale(MatchCount),
                                       ScoreDiff_z = scale(ScoreDiff),
                                       refavg_css_z = scale(refavg_css))
  return(gct_nuk)
}

# Define function to estimate FDR
estfdr <- function(dataset_of_matches, score_choice, score_direction, thresh, model_pep_name){
  
  if(score_direction == "bigger"){
    tempdat <- dataset_of_matches[dataset_of_matches[[score_choice]] >= thresh,]
    fdr <- mean(tempdat[[model_pep_name]])
    
  } else if(score_direction == "smaller"){
    tempdat <- dataset_of_matches[dataset_of_matches[[score_choice]] <= thresh,]
    fdr <- mean(tempdat[[model_pep_name]])
  }
  
  return(fdr)
}

# Define function to determine actual fdr
actualfdr <- function(dataset_of_matches, score_choice, thresh, score_direction){
  if(score_direction == "bigger"){
    tempdat <- dataset_of_matches[dataset_of_matches[[score_choice]] >= thresh,]
    fdr <- mean(tempdat$Truth.Annotation == "True.Negative")
    
  } else if(score_direction == "smaller"){
    tempdat <- dataset_of_matches[dataset_of_matches[[score_choice]] <= thresh,]
    fdr <- mean(tempdat$Truth.Annotation == "True.Negative")
  }
  return(fdr)
}

# Used to fit GAMLSS models given a specified formula and model type
gmmext_fit <- function(dataset_of_matches, base.formula, formula = NULL, model_type, modseed, modreps){
  
  if(model_type == "baseline"){
    set.seed(modseed)
    result <- tryCatch(gamlssMXfits(n = modreps, base.formula, data = dataset_of_matches,
                                    family = NO, K = 2), error = function(x){return(NULL)})
    
  } else if(model_type == "varypi"){
    set.seed(modseed)
    result <- tryCatch(gamlssMXfits(n = modreps, base.formula, pi.formula = formula,
                                    data = dataset_of_matches,
                                    family = NO, K = 2), error = function(x){return(NULL)})
    
  } else if(model_type == "varymu"){
    set.seed(modseed)
    result <- tryCatch(gamlssMXfits(n = modreps, formula, data = dataset_of_matches,
                                    family = NO, K = 2), error = function(x){return(NULL)})
    
  } else if(model_type == "varyboth"){
    set.seed(modseed)
    result <- tryCatch(gamlssMXfits(n = modreps, formula, pi.formula = formula, data = dataset_of_matches,
                                    family = NO, K = 2), error = function(x){return(NULL)})
    
  }
  
  return(result)
}


# Used to fit GAMLSS models given a specified formula and model type
gmmext_fit_proc <- function(dataset_of_matches, scoreidx_nuk, score_choice, score_direction, base.formula, formula = NULL, model_type, modseed, modreps){
  
  if(model_type == "baseline"){
    set.seed(modseed)
    result <- tryCatch(gamlssMXfits(n = modreps, base.formula, data = dataset_of_matches,
                                    family = NO, K = 2), error = function(x){return(NULL)})
    
  } else if(model_type == "varypi"){
    set.seed(modseed)
    result <- tryCatch(gamlssMXfits(n = modreps, base.formula, pi.formula = formula,
                                    data = dataset_of_matches,
                                    family = NO, K = 2), error = function(x){return(NULL)})
    
  } else if(model_type == "varymu"){
    set.seed(modseed)
    result <- tryCatch(gamlssMXfits(n = modreps, formula, data = dataset_of_matches,
                                    family = NO, K = 2), error = function(x){return(NULL)})
    
  } else if(model_type == "varyboth"){
    set.seed(modseed)
    result <- tryCatch(gamlssMXfits(n = modreps, formula, pi.formula = formula, data = dataset_of_matches,
                                    family = NO, K = 2), error = function(x){return(NULL)})
    
  }
  
  AIC_val <- ifelse(is.null(result), NA, AIC(result))
  BIC_val <- ifelse(is.null(result), NA, GAIC(result, k = log(nrow(dataset_of_matches))))
  
  # To select the col that corresponds to PEP, choose the column with the
  # largest probability at the minimum score index (if score_direction is bigger) or max score index (if score_direction is smaller).
  if(is.null(result)){
    peps_mod <- rep(NA, nrow(dataset_of_matches))
  } else{
    pepidx <- which.max(result$post.prob[scoreidx_nuk,])
    peps_mod <- result$post.prob[,pepidx]
    rm(pepidx)
  }
  
  if(is.null(result)){
    fdr_mod <- rep(NA, nrow(dataset_of_matches))
  } else{
    temp <- dataset_of_matches
    temp$PEP <- peps_mod
    fdr_mod <- sapply(temp[[score_choice]], estfdr,
                      dataset_of_matches = temp, score_choice = score_choice,
                      model_pep_name = "PEP", score_direction = score_direction)
    rm(temp)
  }
  rm(peps_mod)
  
  return(list(AIC  = AIC_val,
              BIC  = BIC_val,
              FDRs = fdr_mod))
}

# Define function to fit all models and output relevant results and summaries
fdrcompare_batch_parallel <- function(dataset_of_matches, score_choice, score_direction, numCores, modseed, modreps = 5){
  
  
  varnames <- c("SpctCndtsCount_z", "PoMatchCount_z", "PeakHeight_z",
                "Interference_z", "MatchCount_z", "ScoreDiff_z", "refavg_css_z")
  
  # Find index of lowest score (if score_direction is bigger)
  # Find index of highest score (if score_direction is smaller)
  scoreidx_nuk <- ifelse(score_direction == "bigger", which.min(dataset_of_matches[[score_choice]]), 
                         which.max(dataset_of_matches[[score_choice]]))
  
  allformula <- NULL
  for(i in 1:length(varnames)){
    allformula <- c(allformula,
                    lapply(combn(varnames, i, simplify = FALSE), 
                           function(x){as.formula(paste0("`", score_choice, "`", " ~ ", paste(x, collapse = " + ")))}))
  }
  modnames <- gsub(paste0("`", score_choice, "`", " ~ "), "", as.character(allformula))
  
  # baseline model
  mod_baseline <- gmmext_fit(dataset_of_matches = dataset_of_matches, base.formula = as.formula(paste0("`", score_choice, "`", " ~ 1")),
                             formula = NULL, model_type = "baseline", modseed = modseed, modreps = modreps)
  
  AIC_baseline <- ifelse(is.null(mod_baseline), NA, AIC(mod_baseline))
  BIC_baseline <- ifelse(is.null(mod_baseline), NA, GAIC(mod_baseline, k = log(nrow(dataset_of_matches))))
  
  # To select the col that corresponds to PEP, choose the column with the
  # largest probability at the minimum score index (if score_direction is bigger) or max score index (if score_direction is smaller).
  if(is.null(mod_baseline)){
    peps_baseline <- rep(NA, nrow(dataset_of_matches))
  } else{
    pepidx <- which.max(mod_baseline$post.prob[scoreidx_nuk,])
    peps_baseline <- mod_baseline$post.prob[,pepidx]
    rm(pepidx)
  }
  
  if(is.null(mod_baseline)){
    fdrs_baseline <- rep(NA, nrow(dataset_of_matches))
  } else{
    temp <- dataset_of_matches
    temp$PEP <- peps_baseline
    fdrs_baseline <- sapply(temp[[score_choice]], estfdr,
                            dataset_of_matches = temp, score_choice = score_choice,
                            model_pep_name = "PEP", score_direction = score_direction)
    rm(temp)
  }
  rm(mod_baseline, peps_baseline)
  gc()
  
  
  clust <- makeCluster(numCores, type = "FORK")
  
  # extension 1: vary pi
  mods_varypi <- parLapply(cl = clust, allformula, gmmext_fit_proc, dataset_of_matches = dataset_of_matches, 
                           scoreidx_nuk = scoreidx_nuk, score_choice = score_choice, score_direction = score_direction,
                           base.formula = as.formula(paste0("`", score_choice, "`", " ~ 1")),
                           model_type = "varypi", modseed = modseed, modreps = modreps)
  
  stopCluster(clust)
  
  clust <- makeCluster(numCores, type = "FORK")
  
  # extension 2: vary mu
  mods_varymu <- parLapply(cl = clust, allformula, gmmext_fit_proc, dataset_of_matches = dataset_of_matches, 
                           scoreidx_nuk = scoreidx_nuk, score_choice = score_choice, score_direction = score_direction,
                           base.formula = as.formula(paste0("`", score_choice, "`", " ~ 1")),
                           model_type = "varymu", modseed = modseed, modreps = modreps)
  stopCluster(clust)
  
  
  clust <- makeCluster(numCores, type = "FORK")
  
  # extension 3: vary both
  mods_varyboth <- parLapply(cl = clust, allformula, gmmext_fit_proc, dataset_of_matches = dataset_of_matches, 
                             scoreidx_nuk = scoreidx_nuk, score_choice = score_choice, score_direction = score_direction,
                             base.formula = as.formula(paste0("`", score_choice, "`", " ~ 1")),
                             model_type = "varyboth", modseed = modseed, modreps = modreps)
  
  stopCluster(clust)
  rm(clust)
  
  modnames <- gsub(paste0("`", score_choice, "`", " ~ "), "", as.character(allformula))
  
  # Compile Model Selection Results
  AIC_varypi   <- Reduce("c", rlist::list.ungroup(rlist::list.select(mods_varypi, AIC)))
  AIC_varymu   <- Reduce("c", rlist::list.ungroup(rlist::list.select(mods_varymu, AIC)))
  AIC_varyboth <- Reduce("c", rlist::list.ungroup(rlist::list.select(mods_varyboth, AIC)))
  
  AIC_results <- data.frame(AIC = c(AIC_baseline, AIC_varypi, AIC_varymu, AIC_varyboth),
                            Model = c("baseline", paste0("Ext #1: ", modnames), 
                                      paste0("Ext #2: ", modnames), paste0("Ext #3: ", modnames))) %>%
    dplyr::arrange(AIC) %>%
    dplyr::mutate(AICdiff = AIC - min(AIC, na.rm = TRUE),
                  AICrank = dplyr::row_number()) %>% 
    dplyr::select(Model, AICdiff, AICrank)
  
  
  BIC_varypi   <- Reduce("c", rlist::list.ungroup(rlist::list.select(mods_varypi, BIC)))
  BIC_varymu   <- Reduce("c", rlist::list.ungroup(rlist::list.select(mods_varymu, BIC)))
  BIC_varyboth <- Reduce("c", rlist::list.ungroup(rlist::list.select(mods_varyboth, BIC)))
  
  BIC_results <- data.frame(BIC = c(BIC_baseline, BIC_varypi, BIC_varymu, BIC_varyboth),
                            Model = c("baseline", paste0("Ext #1: ", modnames), 
                                      paste0("Ext #2: ", modnames), paste0("Ext #3: ", modnames))) %>%
    dplyr::arrange(BIC) %>%
    dplyr::mutate(BICdiff = BIC - min(BIC, na.rm = TRUE),
                  BICrank = dplyr::row_number()) %>% 
    dplyr::select(Model, BICdiff, BICrank)
  
  MS_results <- dplyr::left_join(AIC_results, BIC_results) %>% 
    dplyr::mutate(ScoreChoice = score_choice)
  
  rm(AIC_results, BIC_results, AIC_baseline, AIC_varyboth, AIC_varypi, BIC_baseline, BIC_varyboth, BIC_varypi)
  
  # Compile FDRs results
  fdrs_varypi   <- Reduce("cbind", rlist::list.ungroup(rlist::list.select(mods_varypi, FDRs)))
  fdrs_varymu   <- Reduce("cbind", rlist::list.ungroup(rlist::list.select(mods_varymu, FDRs)))
  fdrs_varyboth <- Reduce("cbind", rlist::list.ungroup(rlist::list.select(mods_varyboth, FDRs)))
  
  fdrs_df <- data.frame(data.frame(fdrs_baseline), fdrs_varypi, fdrs_varymu, fdrs_varyboth)
  colnames(fdrs_df) <- c("FDR_baseline", paste0("FDR_ext1_", modnames), 
                         paste0("FDR_ext2_", modnames), paste0("FDR_ext3_", modnames))
  fdr_names <- colnames(fdrs_df)
  
  dataset_of_matches <- data.frame(dataset_of_matches, fdrs_df, check.names = FALSE)
  
  rm(fdrs_df, fdrs_varypi, fdrs_varymu, fdrs_varyboth)
  gc()
  
  
  # Truth
  dataset_of_matches$truFDR <- sapply(dataset_of_matches[[score_choice]], actualfdr, dataset_of_matches = dataset_of_matches, 
                                      score_choice = score_choice, score_direction = score_direction) 
  
  dataset_of_matches$ScoreChoice <- score_choice
  
  fdr_err <- lapply(as.list(fdr_names), function(x, dataset_of_matches){
    return(
      dataset_of_matches[[x]] - dataset_of_matches$truFDR
    )
  }, dataset_of_matches = dataset_of_matches)
  fdr_err <- Reduce("c", fdr_err)
  
  fdr_abserr <- lapply(as.list(fdr_names), function(x, dataset_of_matches){
    return(
      abs(dataset_of_matches[[x]] - dataset_of_matches$truFDR)
    )
  }, dataset_of_matches = dataset_of_matches)
  fdr_abserr <- Reduce("c", fdr_abserr)
  
  fdr_relabserr <- lapply(as.list(fdr_names), function(x, dataset_of_matches, fdrname_bsln){
    return(
      abs(dataset_of_matches[[x]] - dataset_of_matches$truFDR) - abs(dataset_of_matches[[fdrname_bsln]] - dataset_of_matches$truFDR)
    )
  }, dataset_of_matches = dataset_of_matches, fdrname_bsln = fdr_names[1])
  fdr_relabserr <- Reduce("c", fdr_relabserr)
  
  fdr_sqerr <- lapply(as.list(fdr_names), function(x, dataset_of_matches){
    return(
      (dataset_of_matches[[x]] - dataset_of_matches$truFDR)^2
    )
  }, dataset_of_matches = dataset_of_matches)
  fdr_sqerr <- Reduce("c", fdr_sqerr)
  
  FDRerr_df <- data.frame(Err = fdr_err,
                          AbsErr = fdr_abserr,
                          RelAbsErr = fdr_relabserr,
                          SqErr = fdr_sqerr,
                          ModelDetailed = factor(rep(c("baseline", paste0("Ext #1: ", modnames), 
                                                       paste0("Ext #2: ", modnames), paste0("Ext #3: ", modnames)),
                                                     each = nrow(dataset_of_matches)),
                                                 levels = c("baseline", paste0("Ext #1: ", modnames), 
                                                            paste0("Ext #2: ", modnames), paste0("Ext #3: ", modnames))),
                          Score       = dataset_of_matches[[score_choice]],
                          ScoreChoice = score_choice) 
  rm(fdr_err, fdr_abserr, fdr_sqerr)
  gc()
  
  FDRerr_dfsummary <- FDRerr_df %>% dplyr::group_by(ModelDetailed) %>% dplyr::summarize(MeanErr = mean(Err, na.rm = TRUE),
                                                                                        MedianErr = median(Err, na.rm = TRUE),
                                                                                        MeanAbsErr = mean(AbsErr, na.rm = TRUE),
                                                                                        MedianAbsErr = median(AbsErr, na.rm = TRUE),
                                                                                        MeanRelAbsErr = mean(RelAbsErr, na.rm = TRUE),
                                                                                        MedianRelAbsErr = median(RelAbsErr, na.rm = TRUE),
                                                                                        MeanSqErr = mean(SqErr, na.rm = TRUE),
                                                                                        MedianSqErr = median(SqErr, na.rm = TRUE)) %>% 
    dplyr::mutate(ScoreChoice = score_choice)
  
  output_res <- list(Match_data = dataset_of_matches,
                     MSres = MS_results,
                     FDRerr_df = FDRerr_df,
                     FDRerr_dfsummary = FDRerr_dfsummary,
                     modnames = modnames,
                     scorechoice = score_choice)
  
  return(output_res)
}


# Script ------------------------------------------------------------------

for(i in 1:29){
  # Determine number of available cores, use that number less two as total cores for script
  ncores <- parallel::detectCores() - 2
  print(ncores)
  
  # read in the dataset
  fulldat <- readRDS(file = paste0(procdata_folder_path, dataset_totint_names[i], ".rds"))
  
  # Extract arg that indexes score to use, and define score and related args
  myargs <- as.integer(commandArgs(TRUE))
  score_choice    <- score_choices_useful[myargs]
  score_direction <- score_directions_useful[myargs]
  
  # Generate dataset of matches
  matchdat <- dset_improc_score(gcdat_full = fulldat,
                                score_choice = score_choice,
                                score_direction = score_direction)
  
  # Run main script to obtain fdr estimates and main result set
  result <- fdrcompare_batch_parallel(dataset_of_matches = matchdat,
                                      score_choice = score_choice,
                                      score_direction = score_direction,
                                      numCores = ncores, 
                                      modseed = 1123,
                                      modreps = 5)
  saveRDS(result, file = paste0(res_folder_path,"mqfdr_", dataset_totint_names[i], "_", myargs, "_useful.rds"))
  
}



# This script is meant to be used within an HPC environment (i.e. slurm), where an argument indicating the specific
# score choice are provided when the script is ran. The argument should be provided as a numeric,
# ranging between 1:length(score_choices_useful). This argument is read by this script through
# "commandArgs(TRUE)"

# Assign the data directory path to data_folder_path
data_folder_path <- "/data/folder/path/"

# Assign the path to the directory containing library size result files to res_folder_path
libres_folder_path <- "/libres/folder/path/"

# Assign the misc directory path to misc_folder_path
misc_folder_path <- "/misc/folder/path/"

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

# Dataset Info -----------------------------------------------------------

all_dataset_paths <- list.files(path = data_folder_path)
all_dataset_paths <- sapply(all_dataset_paths, function(x){paste0(data_folder_path,x)},
                            USE.NAMES = FALSE)
print(all_dataset_paths)

dataset_totint_names <- c("CSF_1", "CSF_2", "CSF_3", "CSF_4", "CSF_5", "CSF_6", "CSF_7",
                          "Dataset1", "Dataset2", "Dataset3", "Dataset4", "Dataset5",
                          "Group 1_Standards", "Group 2_Standards", "Group 3_Standards", "Group 4_Standards", "Group 6_Standards",
                          "Plasma_1", "Plasma_2", "Plasma_Ref_2", "Standard_AminoacidMix", "Standard_Mix_21Aug20", "Standard_Mix_R1", 
                          "Standard_Mix_R2", "Standard_Mix_R3", "STD_Mix2", "UDN_CSF_metab", "Urine_01", "Urine_2")


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

# Library Info ------------------------------------------------------------

ll_scoremat <- readRDS(file = paste0(misc_folder_path, "ll_scoremat.rds"))

# Functions ---------------------------------------------------------------

# Define function for loading and score-invariant processing steps on each dataset - adapted to a modified library
dset_improc_batch_slib <- function(dset_paths, dset_totint_names, ll_scoremat, slibcmpds){
  
  alldat <- vector("list", length = length(dset_paths))
  for(i in 1:length(dset_paths)){
    if(dset_totint_names[i] == "Group 6_Standards"){
      
      # This special case is created because using readr::read_delim has trouble with this file in particular.
      
      alldat[[i]] <- data.table::fread(file = dset_paths[i], na.strings = c("", NA), 
                                       colClasses = list(character = c("Sample name", "Compound Name", "Reason", "Truth.Annotation")))
      # For whatever reason, fread interprets some numeric collumns as character. Need to add step below to make sure numeric columns
      # are numeric.
      for(numcols in colnames(alldat[[i]])[-c(1,12, 86, 87)]){
        alldat[[i]][[numcols]] <- as.numeric(alldat[[i]][[numcols]])
      }
      rm(numcols)
    } else{
      
      alldat[[i]] <- readr::read_delim(file = dset_paths[i], delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                                       show_col_types = FALSE)
    }
  }
  gcdat <- Reduce("rbind", alldat)
  rm(alldat, i)
  
  # Retain only the rows where the candidate compound is among the compounds in slibcmpds
  gcdat <- gcdat %>% filter(`Compound Name` %in% slibcmpds)
  
  nrow_init <- nrow(gcdat)
  
  # Create a character-valued and numeric-valued ID specific to each spectrum.
  gcdat <- gcdat %>% 
    filter(!is.na(`Compound Name`)) %>% 
    mutate(SampID_chr = paste0(`Sample name`, "_", round(`Retention Time`, 6))) %>%
    arrange(SampID_chr) %>%
    group_by(SampID_chr) %>%
    mutate(SampID_num = cur_group_id()) %>%
    arrange(SampID_num) %>%
    ungroup()
  
  nrow_next <- nrow(gcdat)
  
  # Convert all Truth.Annotation NA values to "Unknown"
  gcdat <- gcdat %>% 
    mutate(Truth.Annotation = ifelse(is.na(Truth.Annotation), "Unknown", Truth.Annotation))
  
  # Create a variable that measures the number of candidates per spectrum
  temp <- gcdat %>% 
    count(SampID_num) %>% 
    rename(SpctCndtsCount = n)
  gcdat <- left_join(gcdat, temp, by = "SampID_num")
  rm(temp)
  
  # Create a variable that measures the number of times a reference compound is identified as a candidate across spectra
  temp <- gcdat %>% 
    count(`Compound Name`) %>% 
    rename(PoMatchCount = n)
  gcdat <- left_join(gcdat, temp, by = "Compound Name")
  rm(temp)
  
  # Create a variable that measures the spectral? Interference
  totintensity <- readRDS(file = paste0(misc_folder_path,"totintensity.rds"))
  totintensity <- totintensity %>%
    filter(`File` %in% dset_totint_names) %>%
    mutate(SampID_chr = paste0(`Sample name`, "_", round(`Retention Time`, 6))) %>%
    dplyr::select(SampID_chr, `TotalIntensity`)
  gcdat <- left_join(gcdat, totintensity, by = "SampID_chr")
  gcdat <- gcdat %>% 
    mutate(Interference =  `TotalIntensity` - `Peak Height`)
  rm(totintensity)
  
  # Compute the average cosine similarity score for each reference compound relative to all other reference compounds
  sll_scoremat <- ll_scoremat[rownames(ll_scoremat) %in% slibcmpds, colnames(ll_scoremat) %in% slibcmpds]
  avgscoredf <- data.frame(`Compound Name` = rownames(sll_scoremat),
                           refavg_css      = NA,
                           check.names = FALSE)
  for(i in 1:nrow(avgscoredf)){
    avgscoredf$refavg_css[i] <- mean(as.numeric(ll_scoremat[i,-i]))
  }
  gcdat <- left_join(gcdat, avgscoredf, by = "Compound Name")
  rm(avgscoredf, sll_scoremat)
  
  nrow_last <- nrow(gcdat)
  
  attr(gcdat, "nrows") <- c(nrow_init, nrow_next, nrow_last)
  
  return(gcdat)
}

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
    
    nrow_init <- nrow(gcdat_top)
    
    # Determine whether there are any duplicate entries
    iso_dupes <- gcdat_top %>% 
      dplyr::count(SampID_num) %>% 
      dplyr::filter(n > 1) %>% 
      dplyr::left_join(gcdat_top, by = "SampID_num")
    
    # Duplicate entries are typically 'isomers' of one another. we remove these from the data. 
    gcdat_top <- gcdat_top %>% dplyr::filter(!(SampID_num %in% iso_dupes$SampID_num))
    
    nrow_next <- nrow(gcdat_top)
    
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
    
    nrow_next2 <- nrow(gcdat_top)
    
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
    
    nrow_init <- nrow(gcdat_top)
    
    # Determine whether there are any duplicate entries
    iso_dupes <- gcdat_top %>% 
      dplyr::count(SampID_num) %>% 
      dplyr::filter(n > 1) %>% 
      dplyr::left_join(gcdat_top)
    
    # Duplicate entries are typically 'isomers' of one another. we remove these from the data. 
    gcdat_top <- gcdat_top %>% dplyr::filter(!(SampID_num %in% iso_dupes$SampID_num))
    
    nrow_next <- nrow(gcdat_top)
    
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
    nrow_next2 <- nrow(gcdat_top)
    
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
  nrow_last <- nrow(gct_nuk)
  
  attr(gct_nuk, "nrows") <- c(nrow_init, nrow_next, nrow_next2, nrow_last)
  
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
fdrcompare_batch_parallel_slim <- function(dataset_of_matches, score_choice, score_direction, numCores, modseed, modreps = 5){
  
  
  
  # Find index of lowest score (if score_direction is bigger)
  # Find index of highest score (if score_direction is smaller)
  scoreidx_nuk <- ifelse(score_direction == "bigger", which.min(dataset_of_matches[[score_choice]]), 
                         which.max(dataset_of_matches[[score_choice]]))
  
 allformula <- list(as.formula(paste0("`", score_choice, "`", " ~ ", "ScoreDiff_z")),
                     as.formula(paste0("`", score_choice, "`", " ~ ", "SpctCndtsCount_z + ScoreDiff_z")),
                     as.formula(paste0("`", score_choice, "`", " ~ ", "PoMatchCount_z + ScoreDiff_z")))
  modnames <- gsub(paste0("`", score_choice, "`", " ~ "), "", as.character(allformula))
  
  # baseline model
  mod_baseline <- gmmext_fit(dataset_of_matches = dataset_of_matches, base.formula = as.formula(paste0("`", score_choice, "`", " ~ 1")),
                             formula = NULL, model_type = "baseline", modseed = modseed, modreps = modreps)
  
  
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
  
  # Compile FDRs results
  fdrs_varypi   <- Reduce("cbind", rlist::list.ungroup(rlist::list.select(mods_varypi, FDRs)))
 
  fdrs_df <- data.frame(data.frame(fdrs_baseline), fdrs_varypi) #, fdrs_varymu, fdrs_varyboth)
  colnames(fdrs_df) <- c("FDR_baseline", paste0("FDR_ext1_", modnames))#, 
  fdr_names <- colnames(fdrs_df)
  
  dataset_of_matches <- data.frame(dataset_of_matches, fdrs_df, check.names = FALSE)
  
  rm(fdrs_df, fdrs_varypi)
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
                          ModelDetailed = factor(rep(c("baseline", paste0("Ext #1: ", modnames)),#, 
                                                     each = nrow(dataset_of_matches)),
                                                 levels = c("baseline", paste0("Ext #1: ", modnames))),
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
                     FDRerr_df = FDRerr_df,
                     FDRerr_dfsummary = FDRerr_dfsummary,
                     modnames = modnames,
                     scorechoice = score_choice)
  
  return(output_res)
}


# Script ------------------------------------------------------------------

# Design of sim study:
# Consider three percentages - 25%, 50%, and 75%
# Consider libraries that are 25%, 50%, and 75% the size of our original reference library.
# Given that our library contains 1284 compounds, this translates to libraries of size 321, 642, and 963 compounds, respectively.
# For each case, to generate a library, randomly sample compound names from the 1284 reference compounds. 
# Do this nsims times for each case, and for each of those nsims datasets, repeat the same analyses used to determine a 
# best model, but only consider the baseline model, and two models that were identified as best based on the analyses 
# using the full reference library.


libpercs <- c(0.25, 0.50, 0.75)
nsims <- 30
seed <- 1123

for(i in 1:length(libpercs)){
  for(j in 1:29){
    set.seed(seed)
    tempdat <- lapply(1:nsims, function(x, dset_path, dset_totint_name, ll_scoremat, libperc){
      slibcmpds <- sample(x = rownames(ll_scoremat), size = libperc*nrow(ll_scoremat), replace = FALSE)
      simdat <- dset_improc_batch_slib(dset_paths        = dset_path, 
                                       dset_totint_names = dset_totint_name, 
                                       ll_scoremat       = ll_scoremat, 
                                       slibcmpds         = slibcmpds)
      return(list(dset = simdat,
                  slib = slibcmpds))
    }, 
    dset_path        = all_dataset_paths[j], 
    dset_totint_name = dataset_totint_names[j],
    ll_scoremat      = ll_scoremat,
    libperc          = libpercs[i])
    
    # Determine number of available cores, use that number less one as total cores for script
    ncores <- parallel::detectCores() - 1
    
    # Extract seed that indexes score to use, and define score and related args
    myargs <- as.integer(commandArgs(TRUE))
    score_choice    <- score_choices_useful[myargs]
    score_direction <- score_directions_useful[myargs]
    
    tempres <- lapply(tempdat, function(x, score_choice, score_direction, ncores){
      
      # Generate dataset of matches
      matchdat <- dset_improc_score(gcdat_full = x$dset,
                                    score_choice = score_choice,
                                    score_direction = score_direction)
      
      # Run main script to obtain fdr estimates and main result set
      result <- fdrcompare_batch_parallel_slim(dataset_of_matches = matchdat,
                                               score_choice = score_choice,
                                               score_direction = score_direction,
                                               numCores = ncores, 
                                               modseed = 1123,
                                               modreps = 5)
      return(list(results = result,
                  rmhist  = c(attr(x$dset, "nrows"), attr(matchdat, "nrows")),
                  slib    = x$slib))
    },
    score_choice    = score_choice,
    score_direction = score_direction,
    ncores          = ncores)
    
    saveRDS(tempres, file = paste0(libres_folder_path, "lbsz_", 
                                   libpercs[i]*100, "_", dataset_totint_names[j], "_", myargs, "_useful.rds"))
    
    rm(tempdat, tempres)
    gc()
  }
}



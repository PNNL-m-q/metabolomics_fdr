# Required packages: 
# require(dplyr)
# require(magrittr)
# require(rlist)
# require(readr)
# require(data.table)

# Assign the data directory path to data_folder_path
data_folder_path <- "/data/folder/path/"
# Assign the misc directory path to misc_folder_path
misc_folder_path <- "/misc/folder/path/"

## Data and functions -------------------------------------------------------------------------

# Retrieves all data file names in specified directory
all_dataset_paths <- list.files(path = data_folder_path)
all_dataset_paths <- sapply(all_dataset_paths, function(x){paste0(data_folder_path, x)},
                            USE.NAMES = FALSE)

# Shorthand names for datafiles
dataset_totint_names <- c("CSF_1", "CSF_2", "CSF_3", "CSF_4", "CSF_5", "CSF_6", "CSF_7",
                          "Dataset1", "Dataset2", "Dataset3", "Dataset4", "Dataset5",
                          "Group 1_Standards", "Group 2_Standards", "Group 3_Standards", "Group 4_Standards", "Group 6_Standards",
                          "Plasma_1", "Plasma_2", "Plasma_Ref_2", "Standard_AminoacidMix", "Standard_Mix_21Aug20", "Standard_Mix_R1", 
                          "Standard_Mix_R2", "Standard_Mix_R3", "STD_Mix2", "UDN_CSF_metab", "Urine_01", "Urine_2")

# All similarity metrics computed in data
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

# useful score indexes
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
# useful score names
score_choices_useful <- score_choices[useful_score_idx]

# Assumes all 'distance' metrics should be minimized (smaller) and all 'similarity' metrics maximized (bigger).
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

# directions corresponding to 'useful' scores
score_directions_useful <- score_directions[useful_score_idx]

# Define function for loading and score-invariant processing steps on each dataset.
# arguments are the path(s) to the dataset(s) (dset_paths) and the corresponding shorthand name(s) (dset_totint_names)
dset_improc_batch <- function(dset_paths, dset_totint_names){
  
  # Read in all specified datasets
  alldat <- vector("list", length = length(dset_paths))
  for(i in 1:length(dset_paths)){
    if(dset_totint_names[i] == "Group 6_Standards"){
      
      # readr::read_delim has difficulty reading in this file. Use data.table::fread instead
      
      alldat[[i]] <- data.table::fread(file = dset_paths[i], na.strings = c("", NA), 
                                       colClasses = list(character = c("Sample name", "Compound Name", "Reason", "Truth.Annotation")))
      # Step below to makes sure known numeric columns are read in as numeric.
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
  
  # Create a character-valued and numeric-valued ID specific to each spectrum.
  gcdat <- gcdat %>% 
    filter(!is.na(`Compound Name`)) %>% 
    mutate(SampID_chr = paste0(`Sample name`, "_", round(`Retention Time`, 6))) %>%
    arrange(SampID_chr) %>%
    group_by(SampID_chr) %>%
    mutate(SampID_num = cur_group_id()) %>%
    arrange(SampID_num) %>%
    ungroup()
  
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
  
  # Create a variable that measures the spectral Interference
  totintensity <- readr::read_csv(file = paste0(misc_folder_path, "Total_Intensity_for_Interference.csv"),
                                  show_col_types = FALSE)
  totintensity <- totintensity %>%
    filter(`File` %in% dset_totint_names) %>%
    mutate(SampID_chr = paste0(`Sample name`, "_", round(`Retention Time`, 6))) %>%
    dplyr::select(SampID_chr, `TotalIntensity`)
  gcdat <- left_join(gcdat, totintensity, by = "SampID_chr")
  gcdat <- gcdat %>% 
    mutate(Interference =  `TotalIntensity` - `Peak Height`)
  rm(totintensity)
  
  # Compute the average cosine similarity score for each reference compound relative to all other reference compounds
  ll_scoremat <- read.csv(paste0(misc_folder_path, "CosineCorrelationLowResGCMSDataBase.csv"), row.names = 1, header = TRUE)
  rownames(ll_scoremat) <- trimws(rownames(ll_scoremat), which = "both") # Important! Trims leading and trailing spaces
  ll_scoremat[lower.tri(ll_scoremat)] <- t(ll_scoremat)[lower.tri(t(ll_scoremat))]
  colnames(ll_scoremat) <- rownames(ll_scoremat)
  avgscoredf <- data.frame(`Compound Name` = rownames(ll_scoremat),
                           refavg_css      = NA,
                           check.names = FALSE)
  for(i in 1:nrow(avgscoredf)){
    avgscoredf$refavg_css[i] <- mean(as.numeric(ll_scoremat[i,-i]))
  }
  gcdat <- left_join(gcdat, avgscoredf, by = "Compound Name")
  rm(avgscoredf, ll_scoremat)
  
  return(gcdat)
}

# Define function to perform score-variant processing steps on each dataset.
# arguments are the dataset yielded by dset_improc_batch (gcdat_full), and the 
# score choice and its corresponding direction.
dset_improc_score_keepuk <- function(gcdat_full, score_choice, score_direction){
  
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
    
    # Duplicate entries are typically 'isomers' of one another. 
    # These are removed from the data. 
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
    
    # Duplicate entries are typically 'isomers' of one another.
    # These are removed from the data.  
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
  
  
  gct_nuk <- gcdat_top 
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

## Table 1 ---------------------------------------------------

# Data import / processing
reslist <- vector("list", length = length(all_dataset_paths))
for(i in 1:length(all_dataset_paths)){
  tdat <- dset_improc_batch(dset_paths = all_dataset_paths[i],
                            dset_totint_names = dataset_totint_names[i])
  reslist[[i]] <- list(nsamps = length(unique(tdat$`Sample name`)),
                       nspectra = length(unique(tdat$`SampID_chr`)),
                       ntp = sum(tdat$Truth.Annotation == "True.Positive"),
                       ntn = sum(tdat$Truth.Annotation == "True.Negative"),
                       nuk = sum(tdat$Truth.Annotation == "Unknown"))
  rm(tdat)
  gc()
}

# Standards
stds_reslist <- reslist[c(13:17, 21:26)]
sum(Reduce("c", list.ungroup(list.select(stds_reslist, nsamps)))) # number of samples
sum(Reduce("c", list.ungroup(list.select(stds_reslist, nspectra)))) # number of spectra
sum(Reduce("c", list.ungroup(list.select(stds_reslist, ntp)))) # number of TP
sum(Reduce("c", list.ungroup(list.select(stds_reslist, ntn)))) # number of TN
sum(Reduce("c", list.ungroup(list.select(stds_reslist, nuk)))) # number of Unknown

# Human CSF
csf_reslist <- reslist[c(1:7, 27)]
sum(Reduce("c", list.ungroup(list.select(csf_reslist, nsamps)))) # number of samples
sum(Reduce("c", list.ungroup(list.select(csf_reslist, nspectra)))) # number of spectra
sum(Reduce("c", list.ungroup(list.select(csf_reslist, ntp)))) # number of TP
sum(Reduce("c", list.ungroup(list.select(csf_reslist, ntn)))) # number of TN
sum(Reduce("c", list.ungroup(list.select(csf_reslist, nuk)))) # number of Unknown

# Human Blood Plasma
blood_reslist <- reslist[c(18, 19, 20)]
sum(Reduce("c", list.ungroup(list.select(blood_reslist, nsamps)))) # number of samples
sum(Reduce("c", list.ungroup(list.select(blood_reslist, nspectra)))) # number of spectra
sum(Reduce("c", list.ungroup(list.select(blood_reslist, ntp)))) # number of TP
sum(Reduce("c", list.ungroup(list.select(blood_reslist, ntn)))) # number of TN
sum(Reduce("c", list.ungroup(list.select(blood_reslist, nuk)))) # number of Unknown

# Human Urine
urine_reslist <- reslist[c(28, 29)]
sum(Reduce("c", list.ungroup(list.select(urine_reslist, nsamps)))) # number of samples
sum(Reduce("c", list.ungroup(list.select(urine_reslist, nspectra)))) # number of spectra
sum(Reduce("c", list.ungroup(list.select(urine_reslist, ntp)))) # number of TP
sum(Reduce("c", list.ungroup(list.select(urine_reslist, ntn)))) # number of TN
sum(Reduce("c", list.ungroup(list.select(urine_reslist, nuk)))) # number of Unknown

# Fungi
fungi_reslist <- reslist[c(10, 11, 12)]
sum(Reduce("c", list.ungroup(list.select(fungi_reslist, nsamps)))) # number of samples
sum(Reduce("c", list.ungroup(list.select(fungi_reslist, nspectra)))) # number of spectra
sum(Reduce("c", list.ungroup(list.select(fungi_reslist, ntp)))) # number of TP
sum(Reduce("c", list.ungroup(list.select(fungi_reslist, ntn)))) # number of TN
sum(Reduce("c", list.ungroup(list.select(fungi_reslist, nuk)))) # number of Unknown

# Soil
soil_reslist <- reslist[c(8, 9)]
sum(Reduce("c", list.ungroup(list.select(soil_reslist, nsamps)))) # number of samples
sum(Reduce("c", list.ungroup(list.select(soil_reslist, nspectra)))) # number of spectra
sum(Reduce("c", list.ungroup(list.select(soil_reslist, ntp)))) # number of TP
sum(Reduce("c", list.ungroup(list.select(soil_reslist, ntn)))) # number of TN
sum(Reduce("c", list.ungroup(list.select(soil_reslist, nuk)))) # number of Unknown

rm(list = ls())

## Table 2 -------------------------------------------------------------------------

reslist <- vector("list", length = length(all_dataset_paths)*length(score_directions_useful))
pb <- txtProgressBar(min = 0, max = length(all_dataset_paths), style = 3)
for(i in 1:length(all_dataset_paths)){
  setTxtProgressBar(pb, i)
  tdat <- dset_improc_batch(dset_paths = all_dataset_paths[i],
                            dset_totint_names = dataset_totint_names[i])
  
  if(i %in% c(13:17, 21:26)){
    dname <- "Standards"
  } else if(i %in% c(1:7, 27)){
    dname <- "CSF"
  } else if(i %in% c(18, 19, 20)){
    dname <- "Blood"
  } else if(i %in% c(28, 29)){
    dname <- "Urine"
  } else if(i %in% c(10, 11, 12)){
    dname <- "Fungi"
  } else if(i %in% c(8, 9)){
    dname <- "Soil"
  }
  
  for(j in 1:length(score_choices_useful)){
    score_choice    <- score_choices_useful[j]
    score_direction <- score_directions_useful[j]
    
    # Generate dataset of matches (identification list)
    matchdat <- dset_improc_score_keepuk(gcdat_full = tdat,
                                         score_choice = score_choice,
                                         score_direction = score_direction)
    
    reslist[[28*(i-1)+j]] <- list(ntp = sum(matchdat$Truth.Annotation == "True.Positive"),
                                  ntn = sum(matchdat$Truth.Annotation == "True.Negative"),
                                  nuk = sum(matchdat$Truth.Annotation == "Unknown"),
                                  perctp = sum(matchdat$Truth.Annotation == "True.Positive")/nrow(matchdat)*100,
                                  perctn = sum(matchdat$Truth.Annotation == "True.Negative")/nrow(matchdat)*100,
                                  percuk = sum(matchdat$Truth.Annotation == "Unknown")/nrow(matchdat)*100,
                                  dset_name = dname,
                                  score_choice = score_choices_useful[j])
    rm(matchdat, score_choice, score_direction)
    gc()
  }
  
  rm(tdat, dname)
  gc()
}
close(pb)


dsetnames <- Reduce("c", list.ungroup(list.select(reslist, dset_name)))

stds_reslist <- reslist[which(dsetnames == "Standards")]
mean(Reduce("c", list.ungroup(list.select(stds_reslist, ntp)))) # number of TP
mean(Reduce("c", list.ungroup(list.select(stds_reslist, ntn)))) # number of TN
mean(Reduce("c", list.ungroup(list.select(stds_reslist, nuk)))) # number of Unknown
mean(Reduce("c", list.ungroup(list.select(stds_reslist, perctp)))) # percentage of TP
mean(Reduce("c", list.ungroup(list.select(stds_reslist, perctn)))) # percentage of TN
mean(Reduce("c", list.ungroup(list.select(stds_reslist, percuk)))) # percentage of Unknown


csf_reslist <- reslist[which(dsetnames == "CSF")]
mean(Reduce("c", list.ungroup(list.select(csf_reslist, ntp)))) # number of TP
mean(Reduce("c", list.ungroup(list.select(csf_reslist, ntn)))) # number of TN
mean(Reduce("c", list.ungroup(list.select(csf_reslist, nuk)))) # number of Unknown
mean(Reduce("c", list.ungroup(list.select(csf_reslist, perctp)))) # percentage of TP
mean(Reduce("c", list.ungroup(list.select(csf_reslist, perctn)))) # percentage of TN
mean(Reduce("c", list.ungroup(list.select(csf_reslist, percuk)))) # percentage of Unknown

blood_reslist <- reslist[which(dsetnames == "Blood")]
mean(Reduce("c", list.ungroup(list.select(blood_reslist, ntp)))) # number of TP
mean(Reduce("c", list.ungroup(list.select(blood_reslist, ntn)))) # number of TN
mean(Reduce("c", list.ungroup(list.select(blood_reslist, nuk)))) # number of Unknown
mean(Reduce("c", list.ungroup(list.select(blood_reslist, perctp)))) # percentage of TP
mean(Reduce("c", list.ungroup(list.select(blood_reslist, perctn)))) # percentage of TN
mean(Reduce("c", list.ungroup(list.select(blood_reslist, percuk)))) # percentage of Unknown

urine_reslist <- reslist[which(dsetnames == "Urine")]
mean(Reduce("c", list.ungroup(list.select(urine_reslist, ntp)))) # number of TP
mean(Reduce("c", list.ungroup(list.select(urine_reslist, ntn)))) # number of TN
mean(Reduce("c", list.ungroup(list.select(urine_reslist, nuk)))) # number of Unknown
mean(Reduce("c", list.ungroup(list.select(urine_reslist, perctp)))) # percentage of TP
mean(Reduce("c", list.ungroup(list.select(urine_reslist, perctn)))) # percentage of TN
mean(Reduce("c", list.ungroup(list.select(urine_reslist, percuk)))) # percentage of Unknown

fungi_reslist <- reslist[which(dsetnames == "Fungi")]
mean(Reduce("c", list.ungroup(list.select(fungi_reslist, ntp)))) # number of TP
mean(Reduce("c", list.ungroup(list.select(fungi_reslist, ntn)))) # number of TN
mean(Reduce("c", list.ungroup(list.select(fungi_reslist, nuk)))) # number of Unknown
mean(Reduce("c", list.ungroup(list.select(fungi_reslist, perctp)))) # percentage of TP
mean(Reduce("c", list.ungroup(list.select(fungi_reslist, perctn)))) # percentage of TN
mean(Reduce("c", list.ungroup(list.select(fungi_reslist, percuk)))) # percentage of Unknown

soil_reslist <- reslist[which(dsetnames == "Soil")]
mean(Reduce("c", list.ungroup(list.select(soil_reslist, ntp)))) # number of TP
mean(Reduce("c", list.ungroup(list.select(soil_reslist, ntn)))) # number of TN
mean(Reduce("c", list.ungroup(list.select(soil_reslist, nuk)))) # number of Unknown
mean(Reduce("c", list.ungroup(list.select(soil_reslist, perctp)))) # percentage of TP
mean(Reduce("c", list.ungroup(list.select(soil_reslist, perctn)))) # percentage of TN
mean(Reduce("c", list.ungroup(list.select(soil_reslist, percuk)))) # percentage of Unknown

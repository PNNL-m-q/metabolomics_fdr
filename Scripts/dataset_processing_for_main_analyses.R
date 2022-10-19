# Required packages:
# require(tidyverse)
# require(parallel)
# require(readr)
# require(data.table)

# Assign the data directory path to data_folder_path
data_folder_path <- "/data/folder/path/"
# Assign the misc directory path to misc_folder_path
misc_folder_path <- "/misc/folder/path/"
# Assign the processed data directory path to data_folder_path
procdata_folder_path <- "processed/data/folder/path/"

# Data and functions -------------------------------------------------------------------------

all_dataset_paths <- list.files(path = data_folder_path)
all_dataset_paths <- sapply(all_dataset_paths, function(x){paste0(data_folder_path,x)},
                            USE.NAMES = FALSE)
dataset_totint_names <- c("CSF_1", "CSF_2", "CSF_3", "CSF_4", "CSF_5", "CSF_6", "CSF_7",
                          "Dataset1", "Dataset2", "Dataset3", "Dataset4", "Dataset5",
                          "Group 1_Standards", "Group 2_Standards", "Group 3_Standards", "Group 4_Standards", "Group 6_Standards",
                          "Plasma_1", "Plasma_2", "Plasma_Ref_2", "Standard_AminoacidMix", "Standard_Mix_21Aug20", "Standard_Mix_R1", 
                          "Standard_Mix_R2", "Standard_Mix_R3", "STD_Mix2", "UDN_CSF_metab", "Urine_01", "Urine_2")
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

# Define function for loading and score-invariant processing steps on each dataset.
dset_improc_batch <- function(dset_paths, dset_totint_names){
  
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
  
  # Create a variable that measures the spectral? Interference
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

# Single dataset import ---------------------------------------------------

for(i in 1:length(all_dataset_paths)){
  tdat <- dset_improc_batch(dset_paths = all_dataset_paths[i],
                            dset_totint_names = dataset_totint_names[i])
  saveRDS(tdat,
          file = paste0(procdata_folder_path, dataset_totint_names[i], ".rds" ))
  rm(tdat)
  gc()
}
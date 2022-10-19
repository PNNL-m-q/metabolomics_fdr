# Required packages: 
# require(dplyr)
# require(magrittr)
# require(rlist)
# require(ggplot2)
# require(ggpubr)

# Assign the path to the directory containing result files to res_folder_path
res_folder_path <- "/data/folder/path/"

## Score Info --------------------------------------------------------------

scr_info <- data.frame(ScoreChoice = c("Similarity Score", "Spectral Similarity Score", 
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
                                       "Soergel Distance"),
                       Direction = c("bigger", "bigger",
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
                                     "smaller"))

# -------------------------------------------------------------------------

## Dataset Info ------------------------------------------------------------

all_resfile_paths <- list.files(path = res_folder_path)

dataset_totint_names <- c("CSF_1", "CSF_2", "CSF_3", "CSF_4", "CSF_5", "CSF_6", "CSF_7",
                          "Dataset1", "Dataset2", "Dataset3", "Dataset4", "Dataset5",
                          "Group 1_Standards", "Group 2_Standards", "Group 3_Standards", "Group 4_Standards", "Group 6_Standards",
                          "Plasma_1", "Plasma_2", "Plasma_Ref_2", "Standard_AminoacidMix", "Standard_Mix_21Aug20", "Standard_Mix_R1", 
                          "Standard_Mix_R2", "Standard_Mix_R3", "STD_Mix2", "UDN_CSF_metab", "Urine_01", "Urine_2")
# -------------------------------------------------------------------------

## Soil ---------------------------------------------------------------------

# Result Processing ---------------------------------------------------------

# "Dataset1", "Dataset2"
# dataset_totint_names[c(8, 9)]

soil_results <- vector("list", length = length(dataset_totint_names[c(8, 9)]))
names(soil_results) <- dataset_totint_names[c(8, 9)]
for(dset in dataset_totint_names[c(8, 9)]){
  resfile_paths <- all_resfile_paths[which(grepl(dset, all_resfile_paths))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(res_folder_path,
                                   resfile_paths[i]))
  }
  soil_results[[dset]] <- tempres
  for(i in 1:length(soil_results[[dset]])){
    soil_results[[dset]][[i]] <- c(soil_results[[dset]][[i]], dataset = dset)
    soil_results[[dset]][[i]]$Match_data$Dataset       <- dset
    soil_results[[dset]][[i]]$MSres$Dataset            <- dset
    soil_results[[dset]][[i]]$FDRerr_df$Dataset        <- dset
    soil_results[[dset]][[i]]$FDRerr_dfsummary$Dataset <- dset
  }
}
rm(tempres, dset, i)



soil_matchdat <- list.ungroup(list.select(list.ungroup(soil_results), Match_data))

soil_tfdr_dat <- lapply(soil_matchdat, function(x){
  score_choice <- unique(x$ScoreChoice)
  tfdr_df <- x %>% dplyr::select(score_choice, truFDR, ScoreChoice, Dataset)
  colnames(tfdr_df) <- c("ScoreVals", "TrueFDR", "ScoreChoice", "Dataset")
  tfdr_df <- tfdr_df %>% 
    dplyr::mutate(ScoreVals = (ScoreVals-min(ScoreVals))/(max(ScoreVals)-min(ScoreVals))) %>%
    dplyr::arrange(ScoreVals)
  return(tfdr_df)
})

soil_tfdr_dat <- Reduce("rbind", soil_tfdr_dat)
soil_tfdr_dat <- dplyr::left_join(soil_tfdr_dat, scr_info, by = "ScoreChoice") %>%
  dplyr::mutate(Group = paste0(Dataset, ":", ScoreChoice),
                Direction = ifelse(Direction == "bigger", "Larger is More Similar", 
                                   "Smaller is More Similar"))

rm(soil_results, soil_matchdat)
gc()

# -------------------------------------------------------------------------

## Fungi ---------------------------------------------------------------------

# Result Processing ---------------------------------------------------------

# "Dataset3", "Dataset4", "Dataset5"
# dataset_totint_names[c(10, 11, 12)]

fungi_results <- vector("list", length = length(dataset_totint_names[c(10, 11, 12)]))
names(fungi_results) <- dataset_totint_names[c(10, 11, 12)]
for(dset in dataset_totint_names[c(10, 11, 12)]){
  resfile_paths <- all_resfile_paths[which(grepl(dset, all_resfile_paths))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(res_folder_path,
                                   resfile_paths[i]))
  }
  fungi_results[[dset]] <- tempres
  for(i in 1:length(fungi_results[[dset]])){
    fungi_results[[dset]][[i]] <- c(fungi_results[[dset]][[i]], dataset = dset)
    fungi_results[[dset]][[i]]$Match_data$Dataset       <- dset
    fungi_results[[dset]][[i]]$MSres$Dataset            <- dset
    fungi_results[[dset]][[i]]$FDRerr_df$Dataset        <- dset
    fungi_results[[dset]][[i]]$FDRerr_dfsummary$Dataset <- dset
  }
}
rm(tempres, dset, i)


fungi_matchdat <- list.ungroup(list.select(list.ungroup(fungi_results), Match_data))

fungi_tfdr_dat <- lapply(fungi_matchdat, function(x){
  score_choice <- unique(x$ScoreChoice)
  tfdr_df <- x %>% dplyr::select(score_choice, truFDR, ScoreChoice, Dataset)
  colnames(tfdr_df) <- c("ScoreVals", "TrueFDR", "ScoreChoice", "Dataset")
  tfdr_df <- tfdr_df %>% 
    dplyr::mutate(ScoreVals = (ScoreVals-min(ScoreVals))/(max(ScoreVals)-min(ScoreVals))) %>%
    dplyr::arrange(ScoreVals)
  return(tfdr_df)
})

fungi_tfdr_dat <- Reduce("rbind", fungi_tfdr_dat)
fungi_tfdr_dat <- dplyr::left_join(fungi_tfdr_dat, scr_info, by = "ScoreChoice") %>%
  dplyr::mutate(Group = paste0(Dataset, ":", ScoreChoice),
                Direction = ifelse(Direction == "bigger", "Larger is More Similar", 
                                   "Smaller is More Similar"))

rm(fungi_results, fungi_matchdat)
gc()

# -------------------------------------------------------------------------

## CSF ---------------------------------------------------------------------

# Result Processing ---------------------------------------------------------

# "CSF_1", "CSF_2", "CSF_3", "CSF_4", "CSF_5", "CSF_6", "CSF_7", "UDN_CSF_metab"
# dataset_totint_names[c(1:7, 27)]

csf_results <- vector("list", length = length(dataset_totint_names[c(1:7, 27)]))
names(csf_results) <- dataset_totint_names[c(1:7, 27)]
for(dset in dataset_totint_names[c(1:7, 27)]){
  resfile_paths <- all_resfile_paths[which(grepl(dset, all_resfile_paths))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(res_folder_path,
                                   resfile_paths[i]))
  }
  csf_results[[dset]] <- tempres
  for(i in 1:length(csf_results[[dset]])){
    csf_results[[dset]][[i]] <- c(csf_results[[dset]][[i]], dataset = dset)
    csf_results[[dset]][[i]]$Match_data$Dataset       <- dset
    csf_results[[dset]][[i]]$MSres$Dataset            <- dset
    csf_results[[dset]][[i]]$FDRerr_df$Dataset        <- dset
    csf_results[[dset]][[i]]$FDRerr_dfsummary$Dataset <- dset
  }
}
rm(tempres, dset, i)


csf_matchdat <- list.ungroup(list.select(list.ungroup(csf_results), Match_data))

csf_tfdr_dat <- lapply(csf_matchdat, function(x){
  score_choice <- unique(x$ScoreChoice)
  tfdr_df <- x %>% dplyr::select(score_choice, truFDR, ScoreChoice, Dataset)
  colnames(tfdr_df) <- c("ScoreVals", "TrueFDR", "ScoreChoice", "Dataset")
  tfdr_df <- tfdr_df %>% 
    dplyr::mutate(ScoreVals = (ScoreVals-min(ScoreVals))/(max(ScoreVals)-min(ScoreVals))) %>%
    dplyr::arrange(ScoreVals)
  return(tfdr_df)
})

csf_tfdr_dat <- Reduce("rbind", csf_tfdr_dat)
csf_tfdr_dat <- dplyr::left_join(csf_tfdr_dat, scr_info, by = "ScoreChoice") %>%
  dplyr::mutate(Group = paste0(Dataset, ":", ScoreChoice),
                Direction = ifelse(Direction == "bigger", "Larger is More Similar", 
                                   "Smaller is More Similar"))

rm(csf_results, csf_matchdat)
gc()

# -------------------------------------------------------------------------

## Blood ---------------------------------------------------------------------

# Result Processing ---------------------------------------------------------

# "Plasma_1", "Plasma_2", "Plasma_Ref_2"
# dataset_totint_names[c(18, 19, 20)]

blood_results <- vector("list", length = length(dataset_totint_names[c(18, 19, 20)]))
names(blood_results) <- dataset_totint_names[c(18, 19, 20)]
for(dset in dataset_totint_names[c(18, 19, 20)]){
  resfile_paths <- all_resfile_paths[which(grepl(dset, all_resfile_paths))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(res_folder_path,
                                   resfile_paths[i]))
  }
  blood_results[[dset]] <- tempres
  for(i in 1:length(blood_results[[dset]])){
    blood_results[[dset]][[i]] <- c(blood_results[[dset]][[i]], dataset = dset)
    blood_results[[dset]][[i]]$Match_data$Dataset       <- dset
    blood_results[[dset]][[i]]$MSres$Dataset            <- dset
    blood_results[[dset]][[i]]$FDRerr_df$Dataset        <- dset
    blood_results[[dset]][[i]]$FDRerr_dfsummary$Dataset <- dset
  }
}
rm(tempres, dset, i)


blood_matchdat <- list.ungroup(list.select(list.ungroup(blood_results), Match_data))

blood_tfdr_dat <- lapply(blood_matchdat, function(x){
  score_choice <- unique(x$ScoreChoice)
  tfdr_df <- x %>% dplyr::select(score_choice, truFDR, ScoreChoice, Dataset)
  colnames(tfdr_df) <- c("ScoreVals", "TrueFDR", "ScoreChoice", "Dataset")
  tfdr_df <- tfdr_df %>% 
    dplyr::mutate(ScoreVals = (ScoreVals-min(ScoreVals))/(max(ScoreVals)-min(ScoreVals))) %>%
    dplyr::arrange(ScoreVals)
  return(tfdr_df)
})

blood_tfdr_dat <- Reduce("rbind", blood_tfdr_dat)
blood_tfdr_dat <- dplyr::left_join(blood_tfdr_dat, scr_info, by = "ScoreChoice") %>%
  dplyr::mutate(Group = paste0(Dataset, ":", ScoreChoice),
                Direction = ifelse(Direction == "bigger", "Larger is More Similar", 
                                   "Smaller is More Similar"))

rm(blood_results, blood_matchdat)
gc()

# -------------------------------------------------------------------------

## Urine ---------------------------------------------------------------------

# Result Processing ---------------------------------------------------------

# "Urine_01", "Urine_2"
# dataset_totint_names[c(28, 29)]

urine_results <- vector("list", length = length(dataset_totint_names[c(28, 29)]))
names(urine_results) <- dataset_totint_names[c(28, 29)]
for(dset in dataset_totint_names[c(28, 29)]){
  resfile_paths <- all_resfile_paths[which(grepl(dset, all_resfile_paths))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(res_folder_path,
                                   resfile_paths[i]))
  }
  urine_results[[dset]] <- tempres
  for(i in 1:length(urine_results[[dset]])){
    urine_results[[dset]][[i]] <- c(urine_results[[dset]][[i]], dataset = dset)
    urine_results[[dset]][[i]]$Match_data$Dataset       <- dset
    urine_results[[dset]][[i]]$MSres$Dataset            <- dset
    urine_results[[dset]][[i]]$FDRerr_df$Dataset        <- dset
    urine_results[[dset]][[i]]$FDRerr_dfsummary$Dataset <- dset
  }
}
rm(tempres, dset, i)


urine_matchdat <- list.ungroup(list.select(list.ungroup(urine_results), Match_data))

urine_tfdr_dat <- lapply(urine_matchdat, function(x){
  score_choice <- unique(x$ScoreChoice)
  tfdr_df <- x %>% dplyr::select(score_choice, truFDR, ScoreChoice, Dataset)
  colnames(tfdr_df) <- c("ScoreVals", "TrueFDR", "ScoreChoice", "Dataset")
  tfdr_df <- tfdr_df %>% 
    dplyr::mutate(ScoreVals = (ScoreVals-min(ScoreVals))/(max(ScoreVals)-min(ScoreVals))) %>%
    dplyr::arrange(ScoreVals)
  return(tfdr_df)
})

urine_tfdr_dat <- Reduce("rbind", urine_tfdr_dat)
urine_tfdr_dat <- dplyr::left_join(urine_tfdr_dat, scr_info, by = "ScoreChoice") %>%
  dplyr::mutate(Group = paste0(Dataset, ":", ScoreChoice),
                Direction = ifelse(Direction == "bigger", "Larger is More Similar", 
                                   "Smaller is More Similar"))

rm(urine_results, urine_matchdat)
gc()

# -------------------------------------------------------------------------

## Standards ---------------------------------------------------------------------

# Result Processing ---------------------------------------------------------

# "Group 1_Standards", "Group 2_Standards", "Group 3_Standards", "Group 4_Standards", "Group 6_Standards"
# "Standard_AminoacidMix", "Standard_Mix_21Aug20", "Standard_Mix_R1", "Standard_Mix_R2", "Standard_Mix_R3", "STD_Mix2"
# dataset_totint_names[c(13:17, 21:26)]

standards_results <- vector("list", length = length(dataset_totint_names[c(13:17, 21:26)]))
names(standards_results) <- dataset_totint_names[c(13:17, 21:26)]
for(dset in dataset_totint_names[c(13:17, 21:26)]){
  resfile_paths <- all_resfile_paths[which(grepl(dset, all_resfile_paths))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(res_folder_path,
                                   resfile_paths[i]))
  }
  standards_results[[dset]] <- tempres
  for(i in 1:length(standards_results[[dset]])){
    standards_results[[dset]][[i]] <- c(standards_results[[dset]][[i]], dataset = dset)
    standards_results[[dset]][[i]]$Match_data$Dataset       <- dset
    standards_results[[dset]][[i]]$MSres$Dataset            <- dset
    standards_results[[dset]][[i]]$FDRerr_df$Dataset        <- dset
    standards_results[[dset]][[i]]$FDRerr_dfsummary$Dataset <- dset
  }
}
rm(tempres, dset, i)


standards_matchdat <- list.ungroup(list.select(list.ungroup(standards_results), Match_data))

standards_tfdr_dat <- lapply(standards_matchdat, function(x){
  score_choice <- unique(x$ScoreChoice)
  tfdr_df <- x %>% dplyr::select(score_choice, truFDR, ScoreChoice, Dataset)
  colnames(tfdr_df) <- c("ScoreVals", "TrueFDR", "ScoreChoice", "Dataset")
  tfdr_df <- tfdr_df %>% 
    dplyr::mutate(ScoreVals = (ScoreVals-min(ScoreVals))/(max(ScoreVals)-min(ScoreVals))) %>%
    dplyr::arrange(ScoreVals)
  return(tfdr_df)
})

standards_tfdr_dat <- Reduce("rbind", standards_tfdr_dat)
standards_tfdr_dat <- dplyr::left_join(standards_tfdr_dat, scr_info, by = "ScoreChoice") %>%
  dplyr::mutate(Group = paste0(Dataset, ":", ScoreChoice),
                Direction = ifelse(Direction == "bigger", "Larger is More Similar", 
                                   "Smaller is More Similar"))


rm(standards_results, standards_matchdat)
gc()

# -------------------------------------------------------------------------

## Figure S1 ---------------------------------------------------

# Soil
soil_p <- ggplot(data = soil_tfdr_dat, 
                 aes(x = ScoreVals, y = TrueFDR, group = Group)) + xlab("Normalized Score Value") +
  ylab("False Discovery Rate") +
  geom_line(alpha = 0.2, color = "red") + theme_bw() + facet_wrap(~Direction) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
  ggtitle("Soil")

# Fungi
fungi_p <- ggplot(data = fungi_tfdr_dat, 
                  aes(x = ScoreVals, y = TrueFDR, group = Group)) + xlab("Normalized Score Value") +
  ylab("False Discovery Rate") +
  geom_line(alpha = 0.2, color = "red") + theme_bw() + facet_wrap(~Direction) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()) +
  ggtitle("Fungi")

# CSF
csf_p <- ggplot(data = csf_tfdr_dat, 
                aes(x = ScoreVals, y = TrueFDR, group = Group)) + xlab("Normalized Score Value") +
  ylab("False Discovery Rate") +
  geom_line(alpha = 0.2, color = "red") + theme_bw() + facet_wrap(~Direction) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()) +
  ggtitle("Human CSF")

# Blood
blood_p <- ggplot(data = blood_tfdr_dat, 
                  aes(x = ScoreVals, y = TrueFDR, group = Group)) + xlab("Normalized Score Value") +
  ylab("False Discovery Rate") +
  geom_line(alpha = 0.2, color = "red") + theme_bw() + facet_wrap(~Direction) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle("Human Blood Plasma")

# Urine
urine_p <- ggplot(data = urine_tfdr_dat, 
                  aes(x = ScoreVals, y = TrueFDR, group = Group)) + xlab("Normalized Score Value") +
  ylab("False Discovery Rate") +
  geom_line(alpha = 0.2, color = "red") + theme_bw() + facet_wrap(~Direction) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()) +
  ggtitle("Human Urine")

# Standards
standards_p <- ggplot(data = standards_tfdr_dat, 
                      aes(x = ScoreVals, y = TrueFDR, group = Group)) + xlab("Normalized Score Value") +
  ylab("False Discovery Rate") +
  geom_line(alpha = 0.2, color = "red") + theme_bw() + facet_wrap(~Direction) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()) +
  ggtitle("Standards")


plot <- ggpubr::ggarrange(soil_p, fungi_p, csf_p, blood_p, urine_p, standards_p, ncol = 3, nrow = 2, widths = c(1,1,1))
ggpubr::annotate_figure(plot, bottom = ggpubr::text_grob("Normalized Score Value", size = 20, hjust = 0.30),
                        left = ggpubr::text_grob("FDR", size = 20, rot = 90))

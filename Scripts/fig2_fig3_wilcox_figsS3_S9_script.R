# Required packages: 
# require(dplyr)
# require(magrittr)
# require(rlist)
# require(ggplot2)

# Assign the path to the directory containing result files to res_folder_path
res_folder_path <- "/data/folder/path/"

## Score Info --------------------------------------------------------------

scr_info <- data.frame(Score = c("Similarity Score", "Spectral Similarity Score", 
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

# Result Import/Processing ---------------------------------------------------------

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

# Aggregate result sets
soil_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(soil_results), FDRerr_dfsummary)))

# Summarize aggregated results
soil_medMAE_df <- soil_combined_fdrerr_df %>% 
  dplyr::group_by(ModelDetailed) %>% 
  summarise(MedMAE = median(MedianAbsErr, na.rm = TRUE)) %>%
  mutate(`Sample Type` = "Soil") %>%
  dplyr::ungroup()

rm(soil_results)
gc()

# -------------------------------------------------------------------------

## Fungi ---------------------------------------------------------------------

# Result Import/Processing ---------------------------------------------------------

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

# Aggregate result sets
fungi_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(fungi_results), FDRerr_dfsummary)))

# Summarize aggregated results
fungi_medMAE_df <- fungi_combined_fdrerr_df %>% 
  dplyr::group_by(ModelDetailed) %>% 
  summarise(MedMAE = median(MedianAbsErr, na.rm = TRUE)) %>%
  mutate(`Sample Type` = "Fungi") %>%
  dplyr::ungroup()

rm(fungi_results)
gc()

# -------------------------------------------------------------------------

## CSF ---------------------------------------------------------------------

# Result Import/Processing ---------------------------------------------------------

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

# Aggregate result sets
csf_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(csf_results), FDRerr_dfsummary)))

# Summarize aggregated results
csf_medMAE_df <- csf_combined_fdrerr_df %>% 
  dplyr::group_by(ModelDetailed) %>% 
  summarise(MedMAE = median(MedianAbsErr, na.rm = TRUE)) %>%
  mutate(`Sample Type` = "Human CSF") %>%
  dplyr::ungroup()

rm(csf_results)
gc()

# -------------------------------------------------------------------------

## Blood ---------------------------------------------------------------------

# Result Import/Processing ---------------------------------------------------------

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

# Aggregate result sets
blood_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(blood_results), FDRerr_dfsummary)))

# Summarize aggregated results
blood_medMAE_df <- blood_combined_fdrerr_df %>% 
  dplyr::group_by(ModelDetailed) %>% 
  summarise(MedMAE = median(MedianAbsErr, na.rm = TRUE)) %>%
  mutate(`Sample Type` = "Human Blood Plasma") %>%
  dplyr::ungroup()

rm(blood_results)
gc()

# -------------------------------------------------------------------------

## Urine ---------------------------------------------------------------------

# Result Import/Processing ---------------------------------------------------------

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

# Aggregate result sets
urine_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(urine_results), FDRerr_dfsummary)))

# Summarize aggregated results
urine_medMAE_df <- urine_combined_fdrerr_df %>% 
  dplyr::group_by(ModelDetailed) %>% 
  summarise(MedMAE = median(MedianAbsErr, na.rm = TRUE)) %>%
  mutate(`Sample Type` = "Human Urine") %>%
  dplyr::ungroup()

rm(urine_results)
gc()

# -------------------------------------------------------------------------

## Standards ---------------------------------------------------------------------

# Result Import/Processing ---------------------------------------------------------

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

# Aggregate result sets
standards_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(standards_results), FDRerr_dfsummary)))

# Summarize aggregated results
standards_medMAE_df <- standards_combined_fdrerr_df %>% 
  dplyr::group_by(ModelDetailed) %>% 
  summarise(MedMAE = median(MedianAbsErr, na.rm = TRUE)) %>%
  mutate(`Sample Type` = "Standards") %>%
  dplyr::ungroup()

rm(standards_results)
gc()

# -------------------------------------------------------------------------



# (Figure 2) Median MAE distributions by Sample Type ---------------------------------

all_medMAE_df <- rbind.data.frame(standards_medMAE_df,
                                  csf_medMAE_df,
                                  blood_medMAE_df,
                                  urine_medMAE_df,
                                  fungi_medMAE_df,
                                  soil_medMAE_df)

# Standards, Human CSF, Human Blood Plasma, Human Urine, Fungi, Soil
temp_bsln <- all_medMAE_df %>% dplyr::filter(ModelDetailed == "baseline") %>% .$MedMAE

all_medMAE_df2 <- all_medMAE_df %>% 
  dplyr::filter(ModelDetailed != "baseline")
all_medMAE_df2 <- all_medMAE_df2 %>% 
  dplyr::mutate(Extension = factor(ifelse(grepl("Ext #1", all_medMAE_df2$ModelDetailed), "Extension #1",
                                          ifelse(grepl("Ext #2", all_medMAE_df2$ModelDetailed), "Extension #2", 
                                                 "Extension #3")), 
                                   levels = c("Extension #1", "Extension #2", "Extension #3")))


temp_standards <- all_medMAE_df2 %>% dplyr::filter(`Sample Type` == "Standards")
sum(temp_standards$MedMAE[temp_standards$Extension == "Extension #1"] < temp_bsln[1])/127*100
sum(temp_standards$MedMAE[temp_standards$Extension == "Extension #2"] < temp_bsln[1])/127*100
sum(temp_standards$MedMAE[temp_standards$Extension == "Extension #3"] < temp_bsln[1])/127*100

temp_csf <- all_medMAE_df2 %>% dplyr::filter(`Sample Type` == "Human CSF")
sum(temp_csf$MedMAE[temp_csf$Extension == "Extension #1"] < temp_bsln[2])/127*100
sum(temp_csf$MedMAE[temp_csf$Extension == "Extension #2"] < temp_bsln[2])/127*100
sum(temp_csf$MedMAE[temp_csf$Extension == "Extension #3"] < temp_bsln[2])/127*100

temp_blood <- all_medMAE_df2 %>% dplyr::filter(`Sample Type` == "Human Blood Plasma")
sum(temp_blood$MedMAE[temp_blood$Extension == "Extension #1"] < temp_bsln[3])/127*100
sum(temp_blood$MedMAE[temp_blood$Extension == "Extension #2"] < temp_bsln[3])/127*100
sum(temp_blood$MedMAE[temp_blood$Extension == "Extension #3"] < temp_bsln[3])/127*100

temp_urine <- all_medMAE_df2 %>% dplyr::filter(`Sample Type` == "Human Urine")
sum(temp_urine$MedMAE[temp_urine$Extension == "Extension #1"] < temp_bsln[4])/127*100
sum(temp_urine$MedMAE[temp_urine$Extension == "Extension #2"] < temp_bsln[4])/127*100
sum(temp_urine$MedMAE[temp_urine$Extension == "Extension #3"] < temp_bsln[4])/127*100

temp_fungi <- all_medMAE_df2 %>% dplyr::filter(`Sample Type` == "Fungi")
sum(temp_fungi$MedMAE[temp_fungi$Extension == "Extension #1"] < temp_bsln[5])/127*100
sum(temp_fungi$MedMAE[temp_fungi$Extension == "Extension #2"] < temp_bsln[5])/127*100
sum(temp_fungi$MedMAE[temp_fungi$Extension == "Extension #3"] < temp_bsln[5])/127*100

temp_soil <- all_medMAE_df2 %>% dplyr::filter(`Sample Type` == "Soil")
sum(temp_soil$MedMAE[temp_soil$Extension == "Extension #1"] < temp_bsln[6])/127*100
sum(temp_soil$MedMAE[temp_soil$Extension == "Extension #2"] < temp_bsln[6])/127*100
sum(temp_soil$MedMAE[temp_soil$Extension == "Extension #3"] < temp_bsln[6])/127*100

annotation <- data.frame(x = c(rep("Standards", 3), rep("Human CSF", 3), rep("Human Blood Plasma", 3), 
                               rep("Human Urine", 3), rep("Fungi", 3), rep("Soil", 3)),
                         y = rep(0.7, 6*3),
                         y2 = c(rep(temp_bsln[1], 3),
                                rep(temp_bsln[2], 3),
                                rep(temp_bsln[3], 3),
                                rep(temp_bsln[4], 3),
                                rep(temp_bsln[5], 3),
                                rep(temp_bsln[6], 3)),
                         Extension = c(rep(c("Extension #1", "Extension #2", "Extension #3"), times = 6)),
                         label = paste0(round(c(sum(temp_standards$MedMAE[temp_standards$Extension == "Extension #1"] < temp_bsln[1])/127*100,
                                                sum(temp_standards$MedMAE[temp_standards$Extension == "Extension #2"] < temp_bsln[1])/127*100,
                                                sum(temp_standards$MedMAE[temp_standards$Extension == "Extension #3"] < temp_bsln[1])/127*100,
                                                sum(temp_csf$MedMAE[temp_csf$Extension == "Extension #1"] < temp_bsln[2])/127*100,
                                                sum(temp_csf$MedMAE[temp_csf$Extension == "Extension #2"] < temp_bsln[2])/127*100,
                                                sum(temp_csf$MedMAE[temp_csf$Extension == "Extension #3"] < temp_bsln[2])/127*100,
                                                sum(temp_blood$MedMAE[temp_blood$Extension == "Extension #1"] < temp_bsln[3])/127*100,
                                                sum(temp_blood$MedMAE[temp_blood$Extension == "Extension #2"] < temp_bsln[3])/127*100,
                                                sum(temp_blood$MedMAE[temp_blood$Extension == "Extension #3"] < temp_bsln[3])/127*100,
                                                sum(temp_urine$MedMAE[temp_urine$Extension == "Extension #1"] < temp_bsln[4])/127*100,
                                                sum(temp_urine$MedMAE[temp_urine$Extension == "Extension #2"] < temp_bsln[4])/127*100,
                                                sum(temp_urine$MedMAE[temp_urine$Extension == "Extension #3"] < temp_bsln[4])/127*100,
                                                sum(temp_fungi$MedMAE[temp_fungi$Extension == "Extension #1"] < temp_bsln[5])/127*100,
                                                sum(temp_fungi$MedMAE[temp_fungi$Extension == "Extension #2"] < temp_bsln[5])/127*100,
                                                sum(temp_fungi$MedMAE[temp_fungi$Extension == "Extension #3"] < temp_bsln[5])/127*100,
                                                sum(temp_soil$MedMAE[temp_soil$Extension == "Extension #1"] < temp_bsln[6])/127*100,
                                                sum(temp_soil$MedMAE[temp_soil$Extension == "Extension #2"] < temp_bsln[6])/127*100,
                                                sum(temp_soil$MedMAE[temp_soil$Extension == "Extension #3"] < temp_bsln[6])/127*100), 2), 
                                        "%"))

ggplot(all_medMAE_df2, aes(x = `Sample Type`, y = MedMAE)) + geom_violin() + 
  geom_point(data = annotation, aes(x = x, y = y2), size = 2, color = "white", shape = 21, fill = "black") + 
  geom_text(data = annotation, aes(x = x, y = y, label = label), size = 5) + coord_flip() +
  ylab("Median MAE") + 
  theme_bw() + theme(strip.text.x = element_text(size = 20),
                     text = element_text(size = 20)) +
  facet_wrap(~Extension) 


# Wilcoxon tests ----------------------------------------------------------

# Soil
soil_bm1 <- soil_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: ScoreDiff_z")) %>% .$MedianAbsErr
soil_bm2 <- soil_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: SpctCndtsCount_z + ScoreDiff_z")) %>% .$MedianAbsErr
soil_bm3 <- soil_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: PoMatchCount_z + ScoreDiff_z")) %>% .$MedianAbsErr
soil_bs <- soil_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("baseline")) %>% .$MedianAbsErr

(soil_bsltest <- c(wilcox.test(soil_bm1, soil_bs, alternative = "less")$p.value,
                   wilcox.test(soil_bm2, soil_bs, alternative = "less")$p.value,
                   wilcox.test(soil_bm3, soil_bs, alternative = "less")$p.value))

(soil_bmeqtest <- c(wilcox.test(soil_bm1, soil_bm2)$p.value,
                    wilcox.test(soil_bm1, soil_bm3)$p.value,
                    wilcox.test(soil_bm2, soil_bm3)$p.value))

# Fungi
fungi_bm1 <- fungi_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: ScoreDiff_z")) %>% .$MedianAbsErr
fungi_bm2 <- fungi_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: SpctCndtsCount_z + ScoreDiff_z")) %>% .$MedianAbsErr
fungi_bm3 <- fungi_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: PoMatchCount_z + ScoreDiff_z")) %>% .$MedianAbsErr
fungi_bs <- fungi_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("baseline")) %>% .$MedianAbsErr

(fungi_bsltest <- c(wilcox.test(fungi_bm1, fungi_bs, alternative = "less")$p.value,
                    wilcox.test(fungi_bm2, fungi_bs, alternative = "less")$p.value,
                    wilcox.test(fungi_bm3, fungi_bs, alternative = "less")$p.value))

(fungi_bmeqtest <- c(wilcox.test(fungi_bm1, fungi_bm2)$p.value,
                     wilcox.test(fungi_bm1, fungi_bm3)$p.value,
                     wilcox.test(fungi_bm2, fungi_bm3)$p.value))

# Human CSF
csf_bm1 <- csf_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: ScoreDiff_z")) %>% .$MedianAbsErr
csf_bm2 <- csf_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: SpctCndtsCount_z + ScoreDiff_z")) %>% .$MedianAbsErr
csf_bm3 <- csf_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: PoMatchCount_z + ScoreDiff_z")) %>% .$MedianAbsErr
csf_bs <- csf_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("baseline")) %>% .$MedianAbsErr

(csf_bsltest <- c(wilcox.test(csf_bm1, csf_bs, alternative = "less")$p.value,
                  wilcox.test(csf_bm2, csf_bs, alternative = "less")$p.value,
                  wilcox.test(csf_bm3, csf_bs, alternative = "less")$p.value))

(csf_bmeqtest <- c(wilcox.test(csf_bm1, csf_bm2)$p.value,
                   wilcox.test(csf_bm1, csf_bm3)$p.value,
                   wilcox.test(csf_bm2, csf_bm3)$p.value))

# Human Blood Plasma
blood_bm1 <- blood_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: ScoreDiff_z")) %>% .$MedianAbsErr
blood_bm2 <- blood_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: SpctCndtsCount_z + ScoreDiff_z")) %>% .$MedianAbsErr
blood_bm3 <- blood_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: PoMatchCount_z + ScoreDiff_z")) %>% .$MedianAbsErr
blood_bs <- blood_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("baseline")) %>% .$MedianAbsErr

(blood_bsltest <- c(wilcox.test(blood_bm1, blood_bs, alternative = "less")$p.value,
                    wilcox.test(blood_bm2, blood_bs, alternative = "less")$p.value,
                    wilcox.test(blood_bm3, blood_bs, alternative = "less")$p.value))

(blood_bmeqtest <- c(wilcox.test(blood_bm1, blood_bm2)$p.value,
                     wilcox.test(blood_bm1, blood_bm3)$p.value,
                     wilcox.test(blood_bm2, blood_bm3)$p.value))

# Human Urine
urine_bm1 <- urine_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: ScoreDiff_z")) %>% .$MedianAbsErr
urine_bm2 <- urine_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: SpctCndtsCount_z + ScoreDiff_z")) %>% .$MedianAbsErr
urine_bm3 <- urine_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: PoMatchCount_z + ScoreDiff_z")) %>% .$MedianAbsErr
urine_bs <- urine_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("baseline")) %>% .$MedianAbsErr

(urine_bsltest <- c(wilcox.test(urine_bm1, urine_bs, alternative = "less")$p.value,
                    wilcox.test(urine_bm2, urine_bs, alternative = "less")$p.value,
                    wilcox.test(urine_bm3, urine_bs, alternative = "less")$p.value))

(urine_bmeqtest <- c(wilcox.test(urine_bm1, urine_bm2)$p.value,
                     wilcox.test(urine_bm1, urine_bm3)$p.value,
                     wilcox.test(urine_bm2, urine_bm3)$p.value))

# Standards
standards_bm1 <- standards_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: ScoreDiff_z")) %>% .$MedianAbsErr
standards_bm2 <- standards_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: SpctCndtsCount_z + ScoreDiff_z")) %>% .$MedianAbsErr
standards_bm3 <- standards_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("Ext #1: PoMatchCount_z + ScoreDiff_z")) %>% .$MedianAbsErr
standards_bs <- standards_combined_fdrerr_df %>% 
  filter(ModelDetailed %in% c("baseline")) %>% .$MedianAbsErr

(standards_bsltest <- c(wilcox.test(standards_bm1, standards_bs, alternative = "less")$p.value,
                        wilcox.test(standards_bm2, standards_bs, alternative = "less")$p.value,
                        wilcox.test(standards_bm3, standards_bs, alternative = "less")$p.value))

(standards_bsbmeqtest <- c(wilcox.test(standards_bm1, standards_bs)$p.value,
                           wilcox.test(standards_bm2, standards_bs)$p.value,
                           wilcox.test(standards_bm3, standards_bs)$p.value))

(standards_bmeqtest <- c(wilcox.test(standards_bm1, standards_bm2)$p.value,
                         wilcox.test(standards_bm1, standards_bm3)$p.value,
                         wilcox.test(standards_bm2, standards_bm3)$p.value))


# (Figure 3) ------------------------------------------------------------------


soil_msp1 <- ggplot(data = soil_combined_fdrerr_df %>% filter(ModelDetailed %in% c("baseline", 
                                                                                       "Ext #1: ScoreDiff_z",
                                                                                       "Ext #1: SpctCndtsCount_z + ScoreDiff_z",
                                                                                       "Ext #1: PoMatchCount_z + ScoreDiff_z")), 
                    aes(x = ModelDetailed, y = MedianAbsErr)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() + geom_jitter(size = 0.5, width = 0.1) +
  ylab("Median Absolute Estimation Error (MAE)") + theme_bw() + 
  stat_summary(geom = "text", fun = median, aes(label = paste0("", sprintf("%1.2f", ..y..))), color = "black", 
               position = position_nudge(x = 0.5)) +
  scale_x_discrete(labels = rev(c("Ext #1: PoMatchCount + ScoreDiff", "Ext #1: SpctCndtsCount + ScoreDiff", 
                                  "Ext #1: ScoreDiff", "Baseline"))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle("Soil")

fungi_msp1 <- ggplot(data = fungi_combined_fdrerr_df %>% filter(ModelDetailed %in% c("baseline", 
                                                                                         "Ext #1: ScoreDiff_z",
                                                                                         "Ext #1: SpctCndtsCount_z + ScoreDiff_z",
                                                                                         "Ext #1: PoMatchCount_z + ScoreDiff_z")), 
                     aes(x = ModelDetailed, y = MedianAbsErr)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() + geom_jitter(size = 0.5, width = 0.1) +
  ylab("Median Absolute Estimation Error") + theme_bw() + 
  stat_summary(geom = "text", fun = median, aes(label = paste0("", sprintf("%1.2f", ..y..))), color = "black", 
               position = position_nudge(x = 0.5)) +
  scale_x_discrete(labels = rev(c("Ext #1: PoMatchCount + ScoreDiff", "Ext #1: SpctCndtsCount + ScoreDiff", 
                                  "Ext #1: ScoreDiff", "Baseline"))) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()) +
  ggtitle("Fungi")

csf_msp1 <- ggplot(data = csf_combined_fdrerr_df %>% filter(ModelDetailed %in% c("baseline", 
                                                                                     "Ext #1: ScoreDiff_z",
                                                                                     "Ext #1: SpctCndtsCount_z + ScoreDiff_z",
                                                                                     "Ext #1: PoMatchCount_z + ScoreDiff_z")), 
                   aes(x = ModelDetailed, y = MedianAbsErr)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() + geom_jitter(size = 0.5, width = 0.1) +
  ylab("Median Absolute Estimation Error") + theme_bw() + 
  stat_summary(geom = "text", fun = median, aes(label = paste0("", sprintf("%1.2f", ..y..))), color = "black", 
               position = position_nudge(x = 0.5)) +
  scale_x_discrete(labels = rev(c("Ext #1: PoMatchCount + ScoreDiff", "Ext #1: SpctCndtsCount + ScoreDiff", 
                                  "Ext #1: ScoreDiff", "Baseline"))) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()) +
  ggtitle("Human CSF")

blood_msp1 <- ggplot(data = blood_combined_fdrerr_df %>% filter(ModelDetailed %in% c("baseline", 
                                                                                         "Ext #1: ScoreDiff_z",
                                                                                         "Ext #1: SpctCndtsCount_z + ScoreDiff_z",
                                                                                         "Ext #1: PoMatchCount_z + ScoreDiff_z")), 
                     aes(x = ModelDetailed, y = MedianAbsErr)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() + geom_jitter(size = 0.5, width = 0.1) +
  ylab("Median Absolute Estimation Error") + theme_bw() + 
  stat_summary(geom = "text", fun = median, aes(label = paste0("", sprintf("%1.2f", ..y..))), color = "black", 
               position = position_nudge(x = 0.5)) +
  scale_x_discrete(labels = rev(c("Ext #1: PoMatchCount + ScoreDiff", "Ext #1: SpctCndtsCount + ScoreDiff", 
                                  "Ext #1: ScoreDiff", "Baseline"))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle("Human Blood Plasma")

urine_msp1 <- ggplot(data = urine_combined_fdrerr_df %>% filter(ModelDetailed %in% c("baseline", 
                                                                                         "Ext #1: ScoreDiff_z",
                                                                                         "Ext #1: SpctCndtsCount_z + ScoreDiff_z",
                                                                                         "Ext #1: PoMatchCount_z + ScoreDiff_z")), 
                     aes(x = ModelDetailed, y = MedianAbsErr)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() + geom_jitter(size = 0.5, width = 0.1) +
  ylab("Median Absolute Estimation Error") + theme_bw() + 
  stat_summary(geom = "text", fun = median, aes(label = paste0("", sprintf("%1.2f", ..y..))), color = "black",
               position = position_nudge(x = 0.5)) +
  scale_x_discrete(labels = rev(c("Ext #1: PoMatchCount + ScoreDiff", "Ext #1: SpctCndtsCount + ScoreDiff", 
                                  "Ext #1: ScoreDiff", "Baseline"))) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()) +
  ggtitle("Human Urine")

standards_msp1 <- ggplot(data = standards_combined_fdrerr_df %>% filter(ModelDetailed %in% c("baseline", 
                                                                                                 "Ext #1: ScoreDiff_z",
                                                                                                 "Ext #1: SpctCndtsCount_z + ScoreDiff_z",
                                                                                                 "Ext #1: PoMatchCount_z + ScoreDiff_z")), 
                         aes(x = ModelDetailed, y = MedianAbsErr)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() + geom_jitter(size = 0.5, width = 0.1) +
  ylab("Median Absolute Estimation Error") + theme_bw() + 
  stat_summary(geom = "text", fun = median, aes(label = paste0("", sprintf("%1.2f", ..y..))), color = "black", 
               position = position_nudge(x = 0.5)) +
  scale_x_discrete(labels = rev(c("Ext #1: PoMatchCount + ScoreDiff", "Ext #1: SpctCndtsCount + ScoreDiff", 
                                  "Ext #1: ScoreDiff", "Baseline"))) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()) +
  ggtitle("Standards")


plot <- ggpubr::ggarrange(soil_msp1, fungi_msp1, csf_msp1, blood_msp1, urine_msp1, standards_msp1, ncol = 3, nrow = 2, widths = c(2.1,1,1))
ggpubr::annotate_figure(plot, bottom = ggpubr::text_grob("Median Absolute Estimation Error (MAE)", hjust = 0.15),
                        left = ggpubr::text_grob("Model", rot = 90))

# (Figures S3-S9) ---------------------------------------------------------

dataset_totint_names <- c("CSF_1", "CSF_2", "CSF_3", "CSF_4", "CSF_5", "CSF_6", "CSF_7",
                          "Dataset1", "Dataset2", "Dataset3", "Dataset4", "Dataset5",
                          "Group 1_Standards", "Group 2_Standards", "Group 3_Standards", "Group 4_Standards", "Group 6_Standards",
                          "Plasma_1", "Plasma_2", "Plasma_Ref_2", "Standard_AminoacidMix", "Standard_Mix_21Aug20", "Standard_Mix_R1", 
                          "Standard_Mix_R2", "Standard_Mix_R3", "STD_Mix2", "UDN_CSF_metab", "Urine_01", "Urine_2")

alldat <- vector("list", length = length(dataset_totint_names))
for(i in 1:length(dataset_totint_names)){
  # read in the dataset
  alldat[[i]] <- readRDS(file = paste0(dataset_totint_names[i], ".rds")) # These .rds files may be generated through dataset_processing_for_main_analyses.R
  alldat[[i]]$Dataset <- dataset_totint_names[i]
}
alldat <- Reduce("rbind", alldat)


alldat <- alldat %>%
  group_by(Dataset) %>%
  summarise(avg_PoMatchCount = mean(PoMatchCount),
            avg_Interference = mean(Interference),
            avg_PeakHgt = mean(`Peak Height`),
            avg_RT = mean(`Retention Time`),
            TP_prop = mean(Truth.Annotation == "True.Positive"),
            TN_prop = mean(Truth.Annotation == "True.Negative"),
            UK_prop = mean(Truth.Annotation == "Unknown"))

blood_combined_fdrerr_df$`Sample Type` <- "Human Blood Plasma"  
csf_combined_fdrerr_df$`Sample Type` <- "Human CSF"
fungi_combined_fdrerr_df$`Sample Type` <- "Fungi"
soil_combined_fdrerr_df$`Sample Type` <- "Soil"
standards_combined_fdrerr_df$`Sample Type` <- "Standards"
urine_combined_fdrerr_df$`Sample Type` <- "Human Urine"

blood_msd <- blood_combined_fdrerr_df %>%
  group_by(Dataset, ModelDetailed) %>%
  summarise(avgMAE = mean(MedianAbsErr),
            sdMAE = sd(MedianAbsErr)) %>%
  mutate(`Sample Type` = "Human Blood Plasma")

csf_msd <- csf_combined_fdrerr_df %>%
  group_by(Dataset, ModelDetailed) %>%
  summarise(avgMAE = mean(MedianAbsErr),
            sdMAE = sd(MedianAbsErr)) %>%
  mutate(`Sample Type` = "Human CSF")

fungi_msd <- fungi_combined_fdrerr_df %>%
  group_by(Dataset, ModelDetailed) %>%
  summarise(avgMAE = mean(MedianAbsErr),
            sdMAE = sd(MedianAbsErr)) %>%
  mutate(`Sample Type` = "Fungi")

soil_msd <- soil_combined_fdrerr_df %>%
  group_by(Dataset, ModelDetailed) %>%
  summarise(avgMAE = mean(MedianAbsErr),
            sdMAE = sd(MedianAbsErr)) %>%
  mutate(`Sample Type` = "Soil")

standards_msd <- standards_combined_fdrerr_df %>%
  group_by(Dataset, ModelDetailed) %>%
  summarise(avgMAE = mean(MedianAbsErr),
            sdMAE = sd(MedianAbsErr)) %>%
  mutate(`Sample Type` = "Standard")

urine_msd <- urine_combined_fdrerr_df %>%
  group_by(Dataset, ModelDetailed) %>%
  summarise(avgMAE = mean(MedianAbsErr),
            sdMAE = sd(MedianAbsErr)) %>%
  mutate(`Sample Type` = "Human Urine")

allres_msd <- rbind.data.frame(blood_msd,
                               csf_msd,
                               fungi_msd,
                               soil_msd,
                               standards_msd,
                               urine_msd)

comb_dat_msd <- left_join(allres_msd, alldat, by = "Dataset")
comb_dat_msd_filt <- comb_dat_msd %>% filter(ModelDetailed %in% c("baseline", "Ext #1: ScoreDiff_z",
                                                                  "Ext #1: SpctCndtsCount_z + ScoreDiff_z", 
                                                                  "Ext #1: PoMatchCount_z + ScoreDiff_z")) %>%
  mutate(ModelDetailed = as.character(ModelDetailed))
comb_dat_msd_filt <- comb_dat_msd_filt %>% 
  dplyr::mutate(ModelDetailed = case_when(
    ModelDetailed == "baseline" ~ "Baseline",
    ModelDetailed == "Ext #1: ScoreDiff_z" ~ "Ext #1: ScoreDiff",
    ModelDetailed == "Ext #1: SpctCndtsCount_z + ScoreDiff_z" ~ "Ext #1: SpctCndtsCount + ScoreDiff",
    ModelDetailed == "Ext #1: PoMatchCount_z + ScoreDiff_z" ~ "Ext #1: PoMatchCount + ScoreDiff"
  )) 

png(filename = "C:\\Users\\flor829\\OneDrive - PNNL\\Projects\\mq\\Papers\\GCMS_FDR\\Reviewer Comments\\Plots\\samptype_trends_pomatch.png", units="in", width=9, height=6, res=600)
ggplot(data = comb_dat_msd_filt, aes(x = avg_PoMatchCount, y = avgMAE)) + 
  geom_smooth(method = "loess") +
  geom_point(aes(color = `Sample Type`)) + theme_bw() + geom_errorbar(aes(ymin = avgMAE - sdMAE, ymax = avgMAE + sdMAE, color = `Sample Type`)) + 
  facet_wrap(~ModelDetailed) + ylab("Median Absolute Estimation Error (MAE)") + 
  xlab("Average Number of Candidate Matches")
dev.off()

png(filename = "C:\\Users\\flor829\\OneDrive - PNNL\\Projects\\mq\\Papers\\GCMS_FDR\\Reviewer Comments\\Plots\\samptype_trends_interference.png", units="in", width=9, height=6, res=600)
ggplot(data = comb_dat_msd_filt, aes(x = avg_Interference, y = avgMAE)) + 
  geom_smooth(method = "loess") +
  geom_point(aes(color = `Sample Type`)) + theme_bw() + geom_errorbar(aes(ymin = avgMAE - sdMAE, ymax = avgMAE + sdMAE, color = `Sample Type`)) + 
  facet_wrap(~ModelDetailed) + ylab("Median Absolute Estimation Error (MAE)") + 
  xlab("Average Spectral Interference")
dev.off()

png(filename = "C:\\Users\\flor829\\OneDrive - PNNL\\Projects\\mq\\Papers\\GCMS_FDR\\Reviewer Comments\\Plots\\samptype_trends_peakheight.png", units="in", width=9, height=6, res=600)
ggplot(data = comb_dat_msd_filt, aes(x = avg_PeakHgt, y = avgMAE)) + 
  geom_smooth(method = "loess") +
  geom_point(aes(color = `Sample Type`)) + theme_bw() + geom_errorbar(aes(ymin = avgMAE - sdMAE, ymax = avgMAE + sdMAE, color = `Sample Type`)) + 
  facet_wrap(~ModelDetailed) + ylab("Median Absolute Estimation Error (MAE)") + 
  xlab("Average Spectral Peak Height")
dev.off()

png(filename = "C:\\Users\\flor829\\OneDrive - PNNL\\Projects\\mq\\Papers\\GCMS_FDR\\Reviewer Comments\\Plots\\samptype_trends_rt.png", units="in", width=9, height=6, res=600)
ggplot(data = comb_dat_msd_filt, aes(x = avg_RT, y = avgMAE)) + 
  geom_smooth(method = "loess") +
  geom_point(aes(color = `Sample Type`)) + theme_bw() + geom_errorbar(aes(ymin = avgMAE - sdMAE, ymax = avgMAE + sdMAE, color = `Sample Type`)) + 
  facet_wrap(~ModelDetailed) + ylab("Median Absolute Estimation Error (MAE)") + 
  xlab("Average Measured Retention Time")
dev.off()

png(filename = "C:\\Users\\flor829\\OneDrive - PNNL\\Projects\\mq\\Papers\\GCMS_FDR\\Reviewer Comments\\Plots\\samptype_trends_tpprop.png", units="in", width=9, height=6, res=600)
ggplot(data = comb_dat_msd_filt, aes(x = TP_prop, y = avgMAE)) + 
  geom_smooth(method = "loess") +
  geom_point(aes(color = `Sample Type`)) + theme_bw() + geom_errorbar(aes(ymin = avgMAE - sdMAE, ymax = avgMAE + sdMAE, color = `Sample Type`)) + 
  facet_wrap(~ModelDetailed) + ylab("Median Absolute Estimation Error (MAE)") + 
  xlab("Proportion of True Positive Annotations")
dev.off()

png(filename = "C:\\Users\\flor829\\OneDrive - PNNL\\Projects\\mq\\Papers\\GCMS_FDR\\Reviewer Comments\\Plots\\samptype_trends_tnprop.png", units="in", width=9, height=6, res=600)
ggplot(data = comb_dat_msd_filt, aes(x = TN_prop, y = avgMAE)) + 
  geom_smooth(method = "loess") +
  geom_point(aes(color = `Sample Type`)) + theme_bw() + geom_errorbar(aes(ymin = avgMAE - sdMAE, ymax = avgMAE + sdMAE, color = `Sample Type`)) + 
  facet_wrap(~ModelDetailed) + ylab("Median Absolute Estimation Error (MAE)") + 
  xlab("Proportion of True Negative Annotations")
dev.off()

png(filename = "C:\\Users\\flor829\\OneDrive - PNNL\\Projects\\mq\\Papers\\GCMS_FDR\\Reviewer Comments\\Plots\\samptype_trends_ukprop.png", units="in", width=9, height=6, res=600)
ggplot(data = comb_dat_msd_filt, aes(x = UK_prop, y = avgMAE)) + 
  geom_smooth(method = "loess") +
  geom_point(aes(color = `Sample Type`)) + theme_bw() + geom_errorbar(aes(ymin = avgMAE - sdMAE, ymax = avgMAE + sdMAE, color = `Sample Type`)) + 
  facet_wrap(~ModelDetailed) + ylab("Median Absolute Estimation Error (MAE)") + 
  xlab("Proportion of Unknown Annotations")
dev.off()

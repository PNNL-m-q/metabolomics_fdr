# Required packages: 
# require(dplyr)
# require(magrittr)
# require(rlist)
# require(ggplot2)

# Assign the path to the directory containing library size result files to res_folder_path
libres_folder_path <- "/libres/folder/path/"

## Dataset Info ------------------------------------------------------------

all_resfile_paths <- list.files(path = libres_folder_path)
resfile_paths_lb25 <- all_resfile_paths[which(grepl("lbsz_25", all_resfile_paths))]
resfile_paths_lb50 <- all_resfile_paths[which(grepl("lbsz_50", all_resfile_paths))]
resfile_paths_lb75 <- all_resfile_paths[which(grepl("lbsz_75", all_resfile_paths))]

dataset_totint_names <- c("CSF_1", "CSF_2", "CSF_3", "CSF_4", "CSF_5", "CSF_6", "CSF_7",
                          "Dataset1", "Dataset2", "Dataset3", "Dataset4", "Dataset5",
                          "Group 1_Standards", "Group 2_Standards", "Group 3_Standards", "Group 4_Standards", "Group 6_Standards",
                          "Plasma_1", "Plasma_2", "Plasma_Ref_2", "Standard_AminoacidMix", "Standard_Mix_21Aug20", "Standard_Mix_R1", 
                          "Standard_Mix_R2", "Standard_Mix_R3", "STD_Mix2", "UDN_CSF_metab", "Urine_01", "Urine_2")

# -------------------------------------------------------------------------

## Soil --------------------------------------------------------------------

# Results Processing (25% set) ------------------------------------------------------

# "Dataset1", "Dataset2"
# dataset_totint_names[c(8, 9)]

soil25_results <- vector("list", length = length(dataset_totint_names[c(8, 9)]))
names(soil25_results) <- dataset_totint_names[c(8, 9)]
for(dset in dataset_totint_names[c(8, 9)]){
  resfile_paths <- resfile_paths_lb25[which(grepl(dset, resfile_paths_lb25))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  soil25_results[[dset]] <- tempres
  for(i in 1:length(soil25_results[[dset]])){
    for(j in 1:length(soil25_results[[dset]][[i]])){
      soil25_results[[dset]][[i]][[j]] <- c(soil25_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      soil25_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      soil25_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      soil25_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      soil25_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      soil25_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      soil25_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

soil25_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(soil25_results)), results)), FDRerr_dfsummary)))
soil25_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(soil25_results)), rmhist)))
soil25_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(soil25_results)), dataset)))
soil25_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(soil25_results)), sim)))
soil25_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(soil25_results)), results)), scorechoice)))

soil25_combined_rmhist_summary <- data.frame(Dataset     = soil25_combined_dset,
                                             Score       = soil25_combined_score,
                                             Sim         = soil25_combined_sim,
                                             Datsize     = soil25_combined_rmhist[,1],
                                             MDatsize    = soil25_combined_rmhist[,4],
                                             Dupsrm      = soil25_combined_rmhist[,4] - soil25_combined_rmhist[,5],
                                             Dupsrm_perc = (soil25_combined_rmhist[,4] - soil25_combined_rmhist[,5])/soil25_combined_rmhist[,4]*100,
                                             Unksrm      = soil25_combined_rmhist[,6] - soil25_combined_rmhist[,7],
                                             Unksrm_perc = (soil25_combined_rmhist[,6] - soil25_combined_rmhist[,7])/soil25_combined_rmhist[,6]*100)
rm(soil25_results, soil25_combined_rmhist, soil25_combined_dset, soil25_combined_sim, soil25_combined_score)
gc()

# Results Processing (50% set) ------------------------------------------------------

soil50_results <- vector("list", length = length(dataset_totint_names[c(8, 9)]))
names(soil50_results) <- dataset_totint_names[c(8, 9)]
for(dset in dataset_totint_names[c(8, 9)]){
  resfile_paths <- resfile_paths_lb50[which(grepl(dset, resfile_paths_lb50))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  soil50_results[[dset]] <- tempres
  for(i in 1:length(soil50_results[[dset]])){
    for(j in 1:length(soil50_results[[dset]][[i]])){
      soil50_results[[dset]][[i]][[j]] <- c(soil50_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      soil50_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      soil50_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      soil50_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      soil50_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      soil50_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      soil50_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

soil50_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(soil50_results)), results)), FDRerr_dfsummary)))
soil50_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(soil50_results)), rmhist)))
soil50_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(soil50_results)), dataset)))
soil50_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(soil50_results)), sim)))
soil50_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(soil50_results)), results)), scorechoice)))

soil50_combined_rmhist_summary <- data.frame(Dataset     = soil50_combined_dset,
                                             Score       = soil50_combined_score,
                                             Sim         = soil50_combined_sim,
                                             Datsize     = soil50_combined_rmhist[,1],
                                             MDatsize    = soil50_combined_rmhist[,4],
                                             Dupsrm      = soil50_combined_rmhist[,4] - soil50_combined_rmhist[,5],
                                             Dupsrm_perc = (soil50_combined_rmhist[,4] - soil50_combined_rmhist[,5])/soil50_combined_rmhist[,4]*100,
                                             Unksrm      = soil50_combined_rmhist[,6] - soil50_combined_rmhist[,7],
                                             Unksrm_perc = (soil50_combined_rmhist[,6] - soil50_combined_rmhist[,7])/soil50_combined_rmhist[,6]*100)
rm(soil50_results, soil50_combined_rmhist, soil50_combined_dset, soil50_combined_sim, soil50_combined_score)
gc()

# Results Processing (75% set) ------------------------------------------------------

soil75_results <- vector("list", length = length(dataset_totint_names[c(8, 9)]))
names(soil75_results) <- dataset_totint_names[c(8, 9)]
for(dset in dataset_totint_names[c(8, 9)]){
  resfile_paths <- resfile_paths_lb75[which(grepl(dset, resfile_paths_lb75))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  soil75_results[[dset]] <- tempres
  for(i in 1:length(soil75_results[[dset]])){
    for(j in 1:length(soil75_results[[dset]][[i]])){
      soil75_results[[dset]][[i]][[j]] <- c(soil75_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      soil75_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      soil75_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      soil75_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      soil75_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      soil75_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      soil75_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

soil75_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(soil75_results)), results)), FDRerr_dfsummary)))
soil75_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(soil75_results)), rmhist)))
soil75_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(soil75_results)), dataset)))
soil75_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(soil75_results)), sim)))
soil75_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(soil75_results)), results)), scorechoice)))

soil75_combined_rmhist_summary <- data.frame(Dataset     = soil75_combined_dset,
                                             Score       = soil75_combined_score,
                                             Sim         = soil75_combined_sim,
                                             Datsize     = soil75_combined_rmhist[,1],
                                             MDatsize    = soil75_combined_rmhist[,4],
                                             Dupsrm      = soil75_combined_rmhist[,4] - soil75_combined_rmhist[,5],
                                             Dupsrm_perc = (soil75_combined_rmhist[,4] - soil75_combined_rmhist[,5])/soil75_combined_rmhist[,4]*100,
                                             Unksrm      = soil75_combined_rmhist[,6] - soil75_combined_rmhist[,7],
                                             Unksrm_perc = (soil75_combined_rmhist[,6] - soil75_combined_rmhist[,7])/soil75_combined_rmhist[,6]*100)
rm(soil75_results, soil75_combined_rmhist, soil75_combined_dset, soil75_combined_sim, soil75_combined_score)
gc()


# Combined ------------------------------------------------------

soil25_combined_fdrerr_df <- soil25_combined_fdrerr_df %>% mutate(`Library Size` = "25%")
soil50_combined_fdrerr_df <- soil50_combined_fdrerr_df %>% mutate(`Library Size` = "50%")
soil75_combined_fdrerr_df <- soil75_combined_fdrerr_df %>% mutate(`Library Size` = "75%")

soil_combined_fdrerr_df <- rbind.data.frame(soil25_combined_fdrerr_df,
                                            soil50_combined_fdrerr_df,
                                            soil75_combined_fdrerr_df)
soil_combined_fdrerr_df <- soil_combined_fdrerr_df %>% 
  mutate(`Library Size` = factor(`Library Size`, levels = c("75%", "50%", "25%")))

# -------------------------------------------------------------------------


## Fungi -------------------------------------------------------------------

# Results Processing (25% set) ------------------------------------------------------

# "Dataset3", "Dataset4", "Dataset5"
# dataset_totint_names[c(10, 11, 12)]

fungi25_results <- vector("list", length = length(dataset_totint_names[c(10, 11, 12)]))
names(fungi25_results) <- dataset_totint_names[c(10, 11, 12)]
for(dset in dataset_totint_names[c(10, 11, 12)]){
  resfile_paths <- resfile_paths_lb25[which(grepl(dset, resfile_paths_lb25))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  fungi25_results[[dset]] <- tempres
  for(i in 1:length(fungi25_results[[dset]])){
    for(j in 1:length(fungi25_results[[dset]][[i]])){
      fungi25_results[[dset]][[i]][[j]] <- c(fungi25_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      fungi25_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      fungi25_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      fungi25_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      fungi25_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      fungi25_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      fungi25_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

fungi25_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(fungi25_results)), results)), FDRerr_dfsummary)))
fungi25_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(fungi25_results)), rmhist)))
fungi25_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(fungi25_results)), dataset)))
fungi25_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(fungi25_results)), sim)))
fungi25_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(fungi25_results)), results)), scorechoice)))

fungi25_combined_rmhist_summary <- data.frame(Dataset     = fungi25_combined_dset,
                                              Score       = fungi25_combined_score,
                                              Sim         = fungi25_combined_sim,
                                              Datsize     = fungi25_combined_rmhist[,1],
                                              MDatsize    = fungi25_combined_rmhist[,4],
                                              Dupsrm      = fungi25_combined_rmhist[,4] - fungi25_combined_rmhist[,5],
                                              Dupsrm_perc = (fungi25_combined_rmhist[,4] - fungi25_combined_rmhist[,5])/fungi25_combined_rmhist[,4]*100,
                                              Unksrm      = fungi25_combined_rmhist[,6] - fungi25_combined_rmhist[,7],
                                              Unksrm_perc = (fungi25_combined_rmhist[,6] - fungi25_combined_rmhist[,7])/fungi25_combined_rmhist[,6]*100)
rm(fungi25_results, fungi25_combined_rmhist, fungi25_combined_dset, fungi25_combined_sim, fungi25_combined_score)
gc()

# Results Processing (50% set) ------------------------------------------------------

fungi50_results <- vector("list", length = length(dataset_totint_names[c(10, 11, 12)]))
names(fungi50_results) <- dataset_totint_names[c(10, 11, 12)]
for(dset in dataset_totint_names[c(10, 11, 12)]){
  resfile_paths <- resfile_paths_lb50[which(grepl(dset, resfile_paths_lb50))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  fungi50_results[[dset]] <- tempres
  for(i in 1:length(fungi50_results[[dset]])){
    for(j in 1:length(fungi50_results[[dset]][[i]])){
      fungi50_results[[dset]][[i]][[j]] <- c(fungi50_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      fungi50_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      fungi50_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      fungi50_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      fungi50_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      fungi50_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      fungi50_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

fungi50_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(fungi50_results)), results)), FDRerr_dfsummary)))
fungi50_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(fungi50_results)), rmhist)))
fungi50_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(fungi50_results)), dataset)))
fungi50_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(fungi50_results)), sim)))
fungi50_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(fungi50_results)), results)), scorechoice)))

fungi50_combined_rmhist_summary <- data.frame(Dataset     = fungi50_combined_dset,
                                              Score       = fungi50_combined_score,
                                              Sim         = fungi50_combined_sim,
                                              Datsize     = fungi50_combined_rmhist[,1],
                                              MDatsize    = fungi50_combined_rmhist[,4],
                                              Dupsrm      = fungi50_combined_rmhist[,4] - fungi50_combined_rmhist[,5],
                                              Dupsrm_perc = (fungi50_combined_rmhist[,4] - fungi50_combined_rmhist[,5])/fungi50_combined_rmhist[,4]*100,
                                              Unksrm      = fungi50_combined_rmhist[,6] - fungi50_combined_rmhist[,7],
                                              Unksrm_perc = (fungi50_combined_rmhist[,6] - fungi50_combined_rmhist[,7])/fungi50_combined_rmhist[,6]*100)
rm(fungi50_results, fungi50_combined_rmhist, fungi50_combined_dset, fungi50_combined_sim, fungi50_combined_score)
gc()


# Results Processing (75% set) ------------------------------------------------------

fungi75_results <- vector("list", length = length(dataset_totint_names[c(10, 11, 12)]))
names(fungi75_results) <- dataset_totint_names[c(10, 11, 12)]
for(dset in dataset_totint_names[c(10, 11, 12)]){
  resfile_paths <- resfile_paths_lb75[which(grepl(dset, resfile_paths_lb75))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  fungi75_results[[dset]] <- tempres
  for(i in 1:length(fungi75_results[[dset]])){
    for(j in 1:length(fungi75_results[[dset]][[i]])){
      fungi75_results[[dset]][[i]][[j]] <- c(fungi75_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      fungi75_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      fungi75_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      fungi75_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      fungi75_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      fungi75_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      fungi75_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

fungi75_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(fungi75_results)), results)), FDRerr_dfsummary)))
fungi75_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(fungi75_results)), rmhist)))
fungi75_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(fungi75_results)), dataset)))
fungi75_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(fungi75_results)), sim)))
fungi75_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(fungi75_results)), results)), scorechoice)))

fungi75_combined_rmhist_summary <- data.frame(Dataset     = fungi75_combined_dset,
                                              Score       = fungi75_combined_score,
                                              Sim         = fungi75_combined_sim,
                                              Datsize     = fungi75_combined_rmhist[,1],
                                              MDatsize    = fungi75_combined_rmhist[,4],
                                              Dupsrm      = fungi75_combined_rmhist[,4] - fungi75_combined_rmhist[,5],
                                              Dupsrm_perc = (fungi75_combined_rmhist[,4] - fungi75_combined_rmhist[,5])/fungi75_combined_rmhist[,4]*100,
                                              Unksrm      = fungi75_combined_rmhist[,6] - fungi75_combined_rmhist[,7],
                                              Unksrm_perc = (fungi75_combined_rmhist[,6] - fungi75_combined_rmhist[,7])/fungi75_combined_rmhist[,6]*100)
rm(fungi75_results, fungi75_combined_rmhist, fungi75_combined_dset, fungi75_combined_sim, fungi75_combined_score)
gc()



# Combined ------------------------------------------------------

fungi25_combined_fdrerr_df <- fungi25_combined_fdrerr_df %>% mutate(`Library Size` = "25%")
fungi50_combined_fdrerr_df <- fungi50_combined_fdrerr_df %>% mutate(`Library Size` = "50%")
fungi75_combined_fdrerr_df <- fungi75_combined_fdrerr_df %>% mutate(`Library Size` = "75%")

fungi_combined_fdrerr_df <- rbind.data.frame(fungi25_combined_fdrerr_df,
                                             fungi50_combined_fdrerr_df,
                                             fungi75_combined_fdrerr_df)
fungi_combined_fdrerr_df <- fungi_combined_fdrerr_df %>% 
  mutate(`Library Size` = factor(`Library Size`, levels = c("75%", "50%", "25%")))


# -------------------------------------------------------------------------


## CSF ---------------------------------------------------------------------

# Results Processing (25% set) ------------------------------------------------------

# "CSF_1", "CSF_2", "CSF_3", "CSF_4", "CSF_5", "CSF_6", "CSF_7", "UDN_CSF_metab"
# dataset_totint_names[c(1:7, 27)]

csf25_results <- vector("list", length = length(dataset_totint_names[c(1:7, 27)]))
names(csf25_results) <- dataset_totint_names[c(1:7, 27)]
for(dset in dataset_totint_names[c(1:7, 27)]){
  resfile_paths <- resfile_paths_lb25[which(grepl(dset, resfile_paths_lb25))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  csf25_results[[dset]] <- tempres
  for(i in 1:length(csf25_results[[dset]])){
    for(j in 1:length(csf25_results[[dset]][[i]])){
      csf25_results[[dset]][[i]][[j]] <- c(csf25_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      csf25_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      csf25_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      csf25_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      csf25_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      csf25_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      csf25_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

csf25_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(csf25_results)), results)), FDRerr_dfsummary)))
csf25_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(csf25_results)), rmhist)))
csf25_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(csf25_results)), dataset)))
csf25_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(csf25_results)), sim)))
csf25_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(csf25_results)), results)), scorechoice)))


csf25_combined_rmhist_summary <- data.frame(Dataset     = csf25_combined_dset,
                                            Score       = csf25_combined_score,
                                            Sim         = csf25_combined_sim,
                                            Datsize     = csf25_combined_rmhist[,1],
                                            MDatsize    = csf25_combined_rmhist[,4],
                                            Dupsrm      = csf25_combined_rmhist[,4] - csf25_combined_rmhist[,5],
                                            Dupsrm_perc = (csf25_combined_rmhist[,4] - csf25_combined_rmhist[,5])/csf25_combined_rmhist[,4]*100,
                                            Unksrm      = csf25_combined_rmhist[,6] - csf25_combined_rmhist[,7],
                                            Unksrm_perc = (csf25_combined_rmhist[,6] - csf25_combined_rmhist[,7])/csf25_combined_rmhist[,6]*100)
rm(csf25_results, csf25_combined_rmhist, csf25_combined_dset, csf25_combined_sim, csf25_combined_score)
gc()

# Results Processing (50% set) ------------------------------------------------------

csf50_results <- vector("list", length = length(dataset_totint_names[c(1:7, 27)]))
names(csf50_results) <- dataset_totint_names[c(1:7, 27)]
for(dset in dataset_totint_names[c(1:7, 27)]){
  resfile_paths <- resfile_paths_lb50[which(grepl(dset, resfile_paths_lb50))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  csf50_results[[dset]] <- tempres
  for(i in 1:length(csf50_results[[dset]])){
    for(j in 1:length(csf50_results[[dset]][[i]])){
      csf50_results[[dset]][[i]][[j]] <- c(csf50_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      csf50_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      csf50_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      csf50_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      csf50_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      csf50_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      csf50_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

csf50_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(csf50_results)), results)), FDRerr_dfsummary)))
csf50_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(csf50_results)), rmhist)))
csf50_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(csf50_results)), dataset)))
csf50_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(csf50_results)), sim)))
csf50_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(csf50_results)), results)), scorechoice)))


csf50_combined_rmhist_summary <- data.frame(Dataset     = csf50_combined_dset,
                                            Score       = csf50_combined_score,
                                            Sim         = csf50_combined_sim,
                                            Datsize     = csf50_combined_rmhist[,1],
                                            MDatsize    = csf50_combined_rmhist[,4],
                                            Dupsrm      = csf50_combined_rmhist[,4] - csf50_combined_rmhist[,5],
                                            Dupsrm_perc = (csf50_combined_rmhist[,4] - csf50_combined_rmhist[,5])/csf50_combined_rmhist[,4]*100,
                                            Unksrm      = csf50_combined_rmhist[,6] - csf50_combined_rmhist[,7],
                                            Unksrm_perc = (csf50_combined_rmhist[,6] - csf50_combined_rmhist[,7])/csf50_combined_rmhist[,6]*100)
rm(csf50_results, csf50_combined_rmhist, csf50_combined_dset, csf50_combined_sim, csf50_combined_score)
gc()

# Results Processing (75% set) ------------------------------------------------------

csf75_results <- vector("list", length = length(dataset_totint_names[c(1:7, 27)]))
names(csf75_results) <- dataset_totint_names[c(1:7, 27)]
for(dset in dataset_totint_names[c(1:7, 27)]){
  resfile_paths <- resfile_paths_lb75[which(grepl(dset, resfile_paths_lb75))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  csf75_results[[dset]] <- tempres
  for(i in 1:length(csf75_results[[dset]])){
    for(j in 1:length(csf75_results[[dset]][[i]])){
      csf75_results[[dset]][[i]][[j]] <- c(csf75_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      csf75_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      csf75_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      csf75_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      csf75_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      csf75_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      csf75_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

csf75_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(csf75_results)), results)), FDRerr_dfsummary)))
csf75_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(csf75_results)), rmhist)))
csf75_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(csf75_results)), dataset)))
csf75_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(csf75_results)), sim)))
csf75_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(csf75_results)), results)), scorechoice)))


csf75_combined_rmhist_summary <- data.frame(Dataset     = csf75_combined_dset,
                                            Score       = csf75_combined_score,
                                            Sim         = csf75_combined_sim,
                                            Datsize     = csf75_combined_rmhist[,1],
                                            MDatsize    = csf75_combined_rmhist[,4],
                                            Dupsrm      = csf75_combined_rmhist[,4] - csf75_combined_rmhist[,5],
                                            Dupsrm_perc = (csf75_combined_rmhist[,4] - csf75_combined_rmhist[,5])/csf75_combined_rmhist[,4]*100,
                                            Unksrm      = csf75_combined_rmhist[,6] - csf75_combined_rmhist[,7],
                                            Unksrm_perc = (csf75_combined_rmhist[,6] - csf75_combined_rmhist[,7])/csf75_combined_rmhist[,6]*100)
rm(csf75_results, csf75_combined_rmhist, csf75_combined_dset, csf75_combined_sim, csf75_combined_score)
gc()


# Combined ------------------------------------------------------

csf25_combined_fdrerr_df <- csf25_combined_fdrerr_df %>% mutate(`Library Size` = "25%")
csf50_combined_fdrerr_df <- csf50_combined_fdrerr_df %>% mutate(`Library Size` = "50%")
csf75_combined_fdrerr_df <- csf75_combined_fdrerr_df %>% mutate(`Library Size` = "75%")

csf_combined_fdrerr_df <- rbind.data.frame(csf25_combined_fdrerr_df,
                                           csf50_combined_fdrerr_df,
                                           csf75_combined_fdrerr_df)
csf_combined_fdrerr_df <- csf_combined_fdrerr_df %>% 
  mutate(`Library Size` = factor(`Library Size`, levels = c("75%", "50%", "25%")))

# -------------------------------------------------------------------------

## Blood -------------------------------------------------------------------

# Results Processing (25% set) ------------------------------------------------------

# "Plasma_1", "Plasma_2", "Plasma_Ref_2"
# dataset_totint_names[c(18, 19, 20)]

blood25_results <- vector("list", length = length(dataset_totint_names[c(18, 19, 20)]))
names(blood25_results) <- dataset_totint_names[c(18, 19, 20)]
for(dset in dataset_totint_names[c(18, 19, 20)]){
  resfile_paths <- resfile_paths_lb25[which(grepl(dset, resfile_paths_lb25))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  blood25_results[[dset]] <- tempres
  for(i in 1:length(blood25_results[[dset]])){
    for(j in 1:length(blood25_results[[dset]][[i]])){
      blood25_results[[dset]][[i]][[j]] <- c(blood25_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      blood25_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      blood25_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      blood25_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      blood25_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      blood25_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      blood25_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

blood25_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(blood25_results)), results)), FDRerr_dfsummary)))
blood25_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(blood25_results)), rmhist)))
blood25_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(blood25_results)), dataset)))
blood25_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(blood25_results)), sim)))
blood25_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(blood25_results)), results)), scorechoice)))

blood25_combined_rmhist_summary <- data.frame(Dataset     = blood25_combined_dset,
                                              Score       = blood25_combined_score,
                                              Sim         = blood25_combined_sim,
                                              Datsize     = blood25_combined_rmhist[,1],
                                              MDatsize    = blood25_combined_rmhist[,4],
                                              Dupsrm      = blood25_combined_rmhist[,4] - blood25_combined_rmhist[,5],
                                              Dupsrm_perc = (blood25_combined_rmhist[,4] - blood25_combined_rmhist[,5])/blood25_combined_rmhist[,4]*100,
                                              Unksrm      = blood25_combined_rmhist[,6] - blood25_combined_rmhist[,7],
                                              Unksrm_perc = (blood25_combined_rmhist[,6] - blood25_combined_rmhist[,7])/blood25_combined_rmhist[,6]*100)
rm(blood25_results, blood25_combined_rmhist, blood25_combined_dset, blood25_combined_sim, blood25_combined_score)
gc()

# Results Processing (50% set) ------------------------------------------------------

blood50_results <- vector("list", length = length(dataset_totint_names[c(18, 19, 20)]))
names(blood50_results) <- dataset_totint_names[c(18, 19, 20)]
for(dset in dataset_totint_names[c(18, 19, 20)]){
  resfile_paths <- resfile_paths_lb50[which(grepl(dset, resfile_paths_lb50))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  blood50_results[[dset]] <- tempres
  for(i in 1:length(blood50_results[[dset]])){
    for(j in 1:length(blood50_results[[dset]][[i]])){
      blood50_results[[dset]][[i]][[j]] <- c(blood50_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      blood50_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      blood50_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      blood50_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      blood50_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      blood50_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      blood50_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

blood50_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(blood50_results)), results)), FDRerr_dfsummary)))
blood50_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(blood50_results)), rmhist)))
blood50_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(blood50_results)), dataset)))
blood50_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(blood50_results)), sim)))
blood50_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(blood50_results)), results)), scorechoice)))


blood50_combined_rmhist_summary <- data.frame(Dataset     = blood50_combined_dset,
                                              Score       = blood50_combined_score,
                                              Sim         = blood50_combined_sim,
                                              Datsize     = blood50_combined_rmhist[,1],
                                              MDatsize    = blood50_combined_rmhist[,4],
                                              Dupsrm      = blood50_combined_rmhist[,4] - blood50_combined_rmhist[,5],
                                              Dupsrm_perc = (blood50_combined_rmhist[,4] - blood50_combined_rmhist[,5])/blood50_combined_rmhist[,4]*100,
                                              Unksrm      = blood50_combined_rmhist[,6] - blood50_combined_rmhist[,7],
                                              Unksrm_perc = (blood50_combined_rmhist[,6] - blood50_combined_rmhist[,7])/blood50_combined_rmhist[,6]*100)
rm(blood50_results, blood50_combined_rmhist, blood50_combined_dset, blood50_combined_sim, blood50_combined_score)
gc()

# Results Processing (75% set) ------------------------------------------------------

blood75_results <- vector("list", length = length(dataset_totint_names[c(18, 19, 20)]))
names(blood75_results) <- dataset_totint_names[c(18, 19, 20)]
for(dset in dataset_totint_names[c(18, 19, 20)]){
  resfile_paths <- resfile_paths_lb75[which(grepl(dset, resfile_paths_lb75))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  blood75_results[[dset]] <- tempres
  for(i in 1:length(blood75_results[[dset]])){
    for(j in 1:length(blood75_results[[dset]][[i]])){
      blood75_results[[dset]][[i]][[j]] <- c(blood75_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      blood75_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      blood75_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      blood75_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      blood75_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      blood75_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      blood75_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

blood75_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(blood75_results)), results)), FDRerr_dfsummary)))
blood75_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(blood75_results)), rmhist)))
blood75_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(blood75_results)), dataset)))
blood75_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(blood75_results)), sim)))
blood75_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(blood75_results)), results)), scorechoice)))

blood75_combined_rmhist_summary <- data.frame(Dataset     = blood75_combined_dset,
                                              Score       = blood75_combined_score,
                                              Sim         = blood75_combined_sim,
                                              Datsize     = blood75_combined_rmhist[,1],
                                              MDatsize    = blood75_combined_rmhist[,4],
                                              Dupsrm      = blood75_combined_rmhist[,4] - blood75_combined_rmhist[,5],
                                              Dupsrm_perc = (blood75_combined_rmhist[,4] - blood75_combined_rmhist[,5])/blood75_combined_rmhist[,4]*100,
                                              Unksrm      = blood75_combined_rmhist[,6] - blood75_combined_rmhist[,7],
                                              Unksrm_perc = (blood75_combined_rmhist[,6] - blood75_combined_rmhist[,7])/blood75_combined_rmhist[,6]*100)
rm(blood75_results, blood75_combined_rmhist, blood75_combined_dset, blood75_combined_sim, blood75_combined_score)
gc()



# Combined ------------------------------------------------------

blood25_combined_fdrerr_df <- blood25_combined_fdrerr_df %>% mutate(`Library Size` = "25%")
blood50_combined_fdrerr_df <- blood50_combined_fdrerr_df %>% mutate(`Library Size` = "50%")
blood75_combined_fdrerr_df <- blood75_combined_fdrerr_df %>% mutate(`Library Size` = "75%")

blood_combined_fdrerr_df <- rbind.data.frame(blood25_combined_fdrerr_df,
                                             blood50_combined_fdrerr_df,
                                             blood75_combined_fdrerr_df)
blood_combined_fdrerr_df <- blood_combined_fdrerr_df %>% 
  mutate(`Library Size` = factor(`Library Size`, levels = c("75%", "50%", "25%")))

# -------------------------------------------------------------------------

## Urine -------------------------------------------------------------------

# Results Processing (25% set) ------------------------------------------------------

# "Urine_01", "Urine_2"
# dataset_totint_names[c(28, 29)]

urine25_results <- vector("list", length = length(dataset_totint_names[c(28, 29)]))
names(urine25_results) <- dataset_totint_names[c(28, 29)]
for(dset in dataset_totint_names[c(28, 29)]){
  resfile_paths <- resfile_paths_lb25[which(grepl(dset, resfile_paths_lb25))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  urine25_results[[dset]] <- tempres
  for(i in 1:length(urine25_results[[dset]])){
    for(j in 1:length(urine25_results[[dset]][[i]])){
      urine25_results[[dset]][[i]][[j]] <- c(urine25_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      urine25_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      urine25_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      urine25_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      urine25_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      urine25_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      urine25_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

urine25_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(urine25_results)), results)), FDRerr_dfsummary)))
urine25_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(urine25_results)), rmhist)))
urine25_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(urine25_results)), dataset)))
urine25_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(urine25_results)), sim)))
urine25_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(urine25_results)), results)), scorechoice)))

urine25_combined_rmhist_summary <- data.frame(Dataset     = urine25_combined_dset,
                                              Score       = urine25_combined_score,
                                              Sim         = urine25_combined_sim,
                                              Datsize     = urine25_combined_rmhist[,1],
                                              MDatsize    = urine25_combined_rmhist[,4],
                                              Dupsrm      = urine25_combined_rmhist[,4] - urine25_combined_rmhist[,5],
                                              Dupsrm_perc = (urine25_combined_rmhist[,4] - urine25_combined_rmhist[,5])/urine25_combined_rmhist[,4]*100,
                                              Unksrm      = urine25_combined_rmhist[,6] - urine25_combined_rmhist[,7],
                                              Unksrm_perc = (urine25_combined_rmhist[,6] - urine25_combined_rmhist[,7])/urine25_combined_rmhist[,6]*100)
rm(urine25_results, urine25_combined_rmhist, urine25_combined_dset, urine25_combined_sim, urine25_combined_score)
gc()

# Results Processing (50% set) ------------------------------------------------------

urine50_results <- vector("list", length = length(dataset_totint_names[c(28, 29)]))
names(urine50_results) <- dataset_totint_names[c(28, 29)]
for(dset in dataset_totint_names[c(28, 29)]){
  resfile_paths <- resfile_paths_lb50[which(grepl(dset, resfile_paths_lb50))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  urine50_results[[dset]] <- tempres
  for(i in 1:length(urine50_results[[dset]])){
    for(j in 1:length(urine50_results[[dset]][[i]])){
      urine50_results[[dset]][[i]][[j]] <- c(urine50_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      urine50_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      urine50_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      urine50_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      urine50_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      urine50_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      urine50_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

urine50_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(urine50_results)), results)), FDRerr_dfsummary)))
urine50_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(urine50_results)), rmhist)))
urine50_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(urine50_results)), dataset)))
urine50_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(urine50_results)), sim)))
urine50_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(urine50_results)), results)), scorechoice)))

urine50_combined_rmhist_summary <- data.frame(Dataset     = urine50_combined_dset,
                                              Score       = urine50_combined_score,
                                              Sim         = urine50_combined_sim,
                                              Datsize     = urine50_combined_rmhist[,1],
                                              MDatsize    = urine50_combined_rmhist[,4],
                                              Dupsrm      = urine50_combined_rmhist[,4] - urine50_combined_rmhist[,5],
                                              Dupsrm_perc = (urine50_combined_rmhist[,4] - urine50_combined_rmhist[,5])/urine50_combined_rmhist[,4]*100,
                                              Unksrm      = urine50_combined_rmhist[,6] - urine50_combined_rmhist[,7],
                                              Unksrm_perc = (urine50_combined_rmhist[,6] - urine50_combined_rmhist[,7])/urine50_combined_rmhist[,6]*100)
rm(urine50_results, urine50_combined_rmhist, urine50_combined_dset, urine50_combined_sim, urine50_combined_score)
gc()

# Results Processing (75% set) ------------------------------------------------------

urine75_results <- vector("list", length = length(dataset_totint_names[c(28, 29)]))
names(urine75_results) <- dataset_totint_names[c(28, 29)]
for(dset in dataset_totint_names[c(28, 29)]){
  resfile_paths <- resfile_paths_lb75[which(grepl(dset, resfile_paths_lb75))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  urine75_results[[dset]] <- tempres
  for(i in 1:length(urine75_results[[dset]])){
    for(j in 1:length(urine75_results[[dset]][[i]])){
      urine75_results[[dset]][[i]][[j]] <- c(urine75_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      urine75_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      urine75_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      urine75_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      urine75_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      urine75_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      urine75_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

urine75_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(urine75_results)), results)), FDRerr_dfsummary)))
urine75_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(urine75_results)), rmhist)))
urine75_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(urine75_results)), dataset)))
urine75_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(urine75_results)), sim)))
urine75_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(urine75_results)), results)), scorechoice)))

urine75_combined_rmhist_summary <- data.frame(Dataset     = urine75_combined_dset,
                                              Score       = urine75_combined_score,
                                              Sim         = urine75_combined_sim,
                                              Datsize     = urine75_combined_rmhist[,1],
                                              MDatsize    = urine75_combined_rmhist[,4],
                                              Dupsrm      = urine75_combined_rmhist[,4] - urine75_combined_rmhist[,5],
                                              Dupsrm_perc = (urine75_combined_rmhist[,4] - urine75_combined_rmhist[,5])/urine75_combined_rmhist[,4]*100,
                                              Unksrm      = urine75_combined_rmhist[,6] - urine75_combined_rmhist[,7],
                                              Unksrm_perc = (urine75_combined_rmhist[,6] - urine75_combined_rmhist[,7])/urine75_combined_rmhist[,6]*100)
rm(urine75_results, urine75_combined_rmhist, urine75_combined_dset, urine75_combined_sim, urine75_combined_score)
gc()


# Combined ------------------------------------------------------

urine25_combined_fdrerr_df <- urine25_combined_fdrerr_df %>% mutate(`Library Size` = "25%")
urine50_combined_fdrerr_df <- urine50_combined_fdrerr_df %>% mutate(`Library Size` = "50%")
urine75_combined_fdrerr_df <- urine75_combined_fdrerr_df %>% mutate(`Library Size` = "75%")

urine_combined_fdrerr_df <- rbind.data.frame(urine25_combined_fdrerr_df,
                                             urine50_combined_fdrerr_df,
                                             urine75_combined_fdrerr_df)
urine_combined_fdrerr_df <- urine_combined_fdrerr_df %>% 
  mutate(`Library Size` = factor(`Library Size`, levels = c("75%", "50%", "25%")))


# -------------------------------------------------------------------------


## Standards ---------------------------------------------------------------

# Results Processing (25% set) ------------------------------------------------------

# "Group 1_Standards", "Group 2_Standards", "Group 3_Standards", "Group 4_Standards", "Group 6_Standards"
# "Standard_AminoacidMix", "Standard_Mix_21Aug20", "Standard_Mix_R1", "Standard_Mix_R2", "Standard_Mix_R3", "STD_Mix2"
# dataset_totint_names[c(13:17, 21:26)]

standards25_results <- vector("list", length = length(dataset_totint_names[c(13:17, 21:26)]))
names(standards25_results) <- dataset_totint_names[c(13:17, 21:26)]
for(dset in dataset_totint_names[c(13:17, 21:26)]){
  resfile_paths <- resfile_paths_lb25[which(grepl(dset, resfile_paths_lb25))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  standards25_results[[dset]] <- tempres
  for(i in 1:length(standards25_results[[dset]])){
    for(j in 1:length(standards25_results[[dset]][[i]])){
      standards25_results[[dset]][[i]][[j]] <- c(standards25_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      standards25_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      standards25_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      standards25_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      standards25_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      standards25_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      standards25_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

standards25_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(standards25_results)), results)), FDRerr_dfsummary)))
standards25_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(standards25_results)), rmhist)))
standards25_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(standards25_results)), dataset)))
standards25_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(standards25_results)), sim)))
standards25_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(standards25_results)), results)), scorechoice)))

standards25_combined_rmhist_summary <- data.frame(Dataset     = standards25_combined_dset,
                                                  Score       = standards25_combined_score,
                                                  Sim         = standards25_combined_sim,
                                                  Datsize     = standards25_combined_rmhist[,1],
                                                  MDatsize    = standards25_combined_rmhist[,4],
                                                  Dupsrm      = standards25_combined_rmhist[,4] - standards25_combined_rmhist[,5],
                                                  Dupsrm_perc = (standards25_combined_rmhist[,4] - standards25_combined_rmhist[,5])/standards25_combined_rmhist[,4]*100,
                                                  Unksrm      = standards25_combined_rmhist[,6] - standards25_combined_rmhist[,7],
                                                  Unksrm_perc = (standards25_combined_rmhist[,6] - standards25_combined_rmhist[,7])/standards25_combined_rmhist[,6]*100)
rm(standards25_results, standards25_combined_rmhist, standards25_combined_dset, standards25_combined_sim, standards25_combined_score)
gc()

# Results Processing (50% set) ------------------------------------------------------

standards50_results <- vector("list", length = length(dataset_totint_names[c(13:17, 21:26)]))
names(standards50_results) <- dataset_totint_names[c(13:17, 21:26)]
for(dset in dataset_totint_names[c(13:17, 21:26)]){
  resfile_paths <- resfile_paths_lb50[which(grepl(dset, resfile_paths_lb50))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  standards50_results[[dset]] <- tempres
  for(i in 1:length(standards50_results[[dset]])){
    for(j in 1:length(standards50_results[[dset]][[i]])){
      standards50_results[[dset]][[i]][[j]] <- c(standards50_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      standards50_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      standards50_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      standards50_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      standards50_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      standards50_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      standards50_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

standards50_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(standards50_results)), results)), FDRerr_dfsummary)))
standards50_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(standards50_results)), rmhist)))
standards50_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(standards50_results)), dataset)))
standards50_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(standards50_results)), sim)))
standards50_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(standards50_results)), results)), scorechoice)))

standards50_combined_rmhist_summary <- data.frame(Dataset     = standards50_combined_dset,
                                                  Score       = standards50_combined_score,
                                                  Sim         = standards50_combined_sim,
                                                  Datsize     = standards50_combined_rmhist[,1],
                                                  MDatsize    = standards50_combined_rmhist[,4],
                                                  Dupsrm      = standards50_combined_rmhist[,4] - standards50_combined_rmhist[,5],
                                                  Dupsrm_perc = (standards50_combined_rmhist[,4] - standards50_combined_rmhist[,5])/standards50_combined_rmhist[,4]*100,
                                                  Unksrm      = standards50_combined_rmhist[,6] - standards50_combined_rmhist[,7],
                                                  Unksrm_perc = (standards50_combined_rmhist[,6] - standards50_combined_rmhist[,7])/standards50_combined_rmhist[,6]*100)
rm(standards50_results, standards50_combined_rmhist, standards50_combined_dset, standards50_combined_sim, standards50_combined_score)
gc()

# Results Processing (75% set) ------------------------------------------------------

standards75_results <- vector("list", length = length(dataset_totint_names[c(13:17, 21:26)]))
names(standards75_results) <- dataset_totint_names[c(13:17, 21:26)]
for(dset in dataset_totint_names[c(13:17, 21:26)]){
  resfile_paths <- resfile_paths_lb75[which(grepl(dset, resfile_paths_lb75))]
  
  tempres <- vector("list", length = length(resfile_paths))
  for(i in 1:length(tempres)){
    tempres[[i]] <- readRDS(paste0(libres_folder_path,
                                   resfile_paths[i]))
  }
  
  standards75_results[[dset]] <- tempres
  for(i in 1:length(standards75_results[[dset]])){
    for(j in 1:length(standards75_results[[dset]][[i]])){
      standards75_results[[dset]][[i]][[j]] <- c(standards75_results[[dset]][[i]][[j]], dataset = dset, sim = j)
      standards75_results[[dset]][[i]][[j]]$results$Match_data$Dataset       <- dset
      standards75_results[[dset]][[i]][[j]]$results$FDRerr_df$Dataset        <- dset
      standards75_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Dataset <- dset
      standards75_results[[dset]][[i]][[j]]$results$Match_data$Sim       <- j
      standards75_results[[dset]][[i]][[j]]$results$FDRerr_df$Sim        <- j
      standards75_results[[dset]][[i]][[j]]$results$FDRerr_dfsummary$Sim <- j
    }
  }
}
rm(resfile_paths, tempres, dset, i, j)

standards75_combined_fdrerr_df <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(standards75_results)), results)), FDRerr_dfsummary)))
standards75_combined_rmhist <- Reduce("rbind", list.ungroup(list.select(list.ungroup(list.ungroup(standards75_results)), rmhist)))
standards75_combined_dset <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(standards75_results)), dataset)))
standards75_combined_sim <- Reduce("c", list.ungroup(list.select(list.ungroup(list.ungroup(standards75_results)), sim)))
standards75_combined_score <- Reduce("c", list.ungroup(list.select(list.ungroup(list.select(list.ungroup(list.ungroup(standards75_results)), results)), scorechoice)))

standards75_combined_rmhist_summary <- data.frame(Dataset     = standards75_combined_dset,
                                                  Score       = standards75_combined_score,
                                                  Sim         = standards75_combined_sim,
                                                  Datsize     = standards75_combined_rmhist[,1],
                                                  MDatsize    = standards75_combined_rmhist[,4],
                                                  Dupsrm      = standards75_combined_rmhist[,4] - standards75_combined_rmhist[,5],
                                                  Dupsrm_perc = (standards75_combined_rmhist[,4] - standards75_combined_rmhist[,5])/standards75_combined_rmhist[,4]*100,
                                                  Unksrm      = standards75_combined_rmhist[,6] - standards75_combined_rmhist[,7],
                                                  Unksrm_perc = (standards75_combined_rmhist[,6] - standards75_combined_rmhist[,7])/standards75_combined_rmhist[,6]*100)
rm(standards75_results, standards75_combined_rmhist, standards75_combined_dset, standards75_combined_sim, standards75_combined_score)
gc()


# Combined -----------------------------------------------------------

standards25_combined_fdrerr_df <- standards25_combined_fdrerr_df %>% mutate(`Library Size` = "25%")
standards50_combined_fdrerr_df <- standards50_combined_fdrerr_df %>% mutate(`Library Size` = "50%")
standards75_combined_fdrerr_df <- standards75_combined_fdrerr_df %>% mutate(`Library Size` = "75%")

standards_combined_fdrerr_df <- rbind.data.frame(standards25_combined_fdrerr_df,
                                                 standards50_combined_fdrerr_df,
                                                 standards75_combined_fdrerr_df)
standards_combined_fdrerr_df <- standards_combined_fdrerr_df %>% 
  mutate(`Library Size` = factor(`Library Size`, levels = c("75%", "50%", "25%")))

# -------------------------------------------------------------------------


## Figure 4 ----------------------------------------------------------------

soilp <- ggplot(data = soil_combined_fdrerr_df, 
                aes(x = ModelDetailed, y = MedianAbsErr, fill = `Library Size`)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() +
  stat_summary(fun.data = function(x){return(c(y = median(x)*1.055, label = length(x)))}, 
               geom = "text", fun = median, aes(label = paste0("", sprintf("%1.2f", ..y../1.055))), color = "black",  size = 4, 
               position = position_dodge(width = 0.75)) + xlab("Model") +
  ylab("Median Absolute Estimation Error") + theme_bw() + 
  scale_x_discrete(labels = rev(c("Ext #1: PoMatchCount + ScoreDiff", "Ext #1: SpctCndtsCount + ScoreDiff", 
                                  "Ext #1: ScoreDiff", "Baseline"))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  scale_fill_brewer(palette = "Oranges",
                    guide = guide_legend(reverse=T)) +
  ggtitle("Soil")

fungip <- ggplot(data = fungi_combined_fdrerr_df, 
                 aes(x = ModelDetailed, y = MedianAbsErr, fill = `Library Size`)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() +
  stat_summary(fun.data = function(x){return(c(y = median(x)*1.075, label = length(x)))}, 
               geom = "text", fun = median, aes(label = paste0("", sprintf("%1.2f", ..y../1.075))), color = "black",  size = 4, 
               position = position_dodge(width = 0.75)) + xlab("Model") +
  ylab("Median Absolute Estimation Error") + theme_bw() + 
  scale_x_discrete(labels = rev(c("Ext #1: PoMatchCount + ScoreDiff", "Ext #1: SpctCndtsCount + ScoreDiff", 
                                  "Ext #1: ScoreDiff", "Baseline"))) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  scale_fill_brewer(palette = "Oranges",
                    guide = guide_legend(reverse=T)) + 
  ggtitle("Fungi")

csfp <- ggplot(data = csf_combined_fdrerr_df, 
               aes(x = ModelDetailed, y = MedianAbsErr, fill = `Library Size`)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() +
  stat_summary(fun.data = function(x){return(c(y = median(x)*1.075, label = length(x)))}, 
               geom = "text", fun = median, aes(label = paste0("", sprintf("%1.2f", ..y../1.075))), color = "black",  size = 4, 
               position = position_dodge(width = 0.75)) + xlab("Model") +
  ylab("Median Absolute Estimation Error") + theme_bw() + 
  scale_x_discrete(labels = rev(c("Ext #1: PoMatchCount + ScoreDiff", "Ext #1: SpctCndtsCount + ScoreDiff", 
                                  "Ext #1: ScoreDiff", "Baseline"))) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  scale_fill_brewer(palette = "Oranges",
                    guide = guide_legend(reverse=T)) + 
  ggtitle("Human CSF")


bloodp <- ggplot(data = blood_combined_fdrerr_df, 
                 aes(x = ModelDetailed, y = MedianAbsErr, fill = `Library Size`)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() +
  stat_summary(fun.data = function(x){return(c(y = median(x)*1.075, label = length(x)))}, 
               geom = "text", fun = median, aes(label = paste0("", sprintf("%1.2f", ..y../1.075))), color = "black",  size = 4, 
               position = position_dodge(width = 0.75)) + xlab("Model") +
  ylab("Median Absolute Estimation Error") + theme_bw() + 
  scale_x_discrete(labels = rev(c("Ext #1: PoMatchCount + ScoreDiff", "Ext #1: SpctCndtsCount + ScoreDiff", 
                                  "Ext #1: ScoreDiff", "Baseline"))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  scale_fill_brewer(palette = "Oranges",
                    guide = guide_legend(reverse=T)) + 
  ggtitle("Human Blood Plasma")

urinep <- ggplot(data = urine_combined_fdrerr_df, 
                 aes(x = ModelDetailed, y = MedianAbsErr, fill = `Library Size`)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() +
  stat_summary(fun.data = function(x){return(c(y = median(x)*1.075, label = length(x)))}, 
               geom = "text", fun = median, aes(label = paste0("", sprintf("%1.2f", ..y../1.075))), color = "black",  size = 4, 
               position = position_dodge(width = 0.75)) + xlab("Model") +
  ylab("Median Absolute Estimation Error") + theme_bw() + 
  scale_x_discrete(labels = rev(c("Ext #1: PoMatchCount + ScoreDiff", "Ext #1: SpctCndtsCount + ScoreDiff", 
                                  "Ext #1: ScoreDiff", "Baseline"))) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  scale_fill_brewer(palette = "Oranges",
                    guide = guide_legend(reverse=T)) + 
  ggtitle("Human Urine")


standardsp <- ggplot(data = standards_combined_fdrerr_df, 
                     aes(x = ModelDetailed, y = MedianAbsErr, fill = `Library Size`)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() +
  stat_summary(fun.data = function(x){return(c(y = median(x)*1.075, label = length(x)))}, 
               geom = "text", fun = median, aes(label = paste0("", sprintf("%1.2f", ..y../1.075))), color = "black",  size = 4, 
               position = position_dodge(width = 0.75)) + xlab("Model") +
  ylab("Median Absolute Estimation Error") + theme_bw() + 
  scale_x_discrete(labels = rev(c("Ext #1: PoMatchCount + ScoreDiff", "Ext #1: SpctCndtsCount + ScoreDiff", 
                                  "Ext #1: ScoreDiff", "Baseline"))) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  scale_fill_brewer(palette = "Oranges",
                    guide = guide_legend(reverse=T)) + 
  ggtitle("Standards")


plot <- ggpubr::ggarrange(soilp, fungip, csfp, bloodp, urinep, standardsp, ncol = 3, nrow = 2, widths = c(1.75,1,1),
                          common.legend = TRUE, legend = "right")
ggpubr::annotate_figure(plot, bottom = ggpubr::text_grob("Median Absolute Estimation Error (MAE)", size = 20, hjust = 0.25),
                        left = ggpubr::text_grob("Model", size = 20, rot = 90))
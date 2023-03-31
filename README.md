# Repository Overview
This repository contains the code for analyses presented in "Gaussian Mixture Modeling Extensions for Improved False Discovery Rate Estimation in GC-MS Metabolomics"

The code provided in this repostiory is meant for reference and is not intended for use as a pipeline. 

# Contents

Misc/ contains secondary data files used by one or more of the analysis scripts contained within R/, as well as a script for implementing the hierarchical empirical bayes model introduced by Jeong et al. (2011). The specific files contained within the Misc/ directory are:
    1. CosineCorrelationLowResGCMSDataBase.csv : This is a .csv file containing the cosin similarity scores computed between all library reference compounds used in our analyses.
    2. MetPC_modified.R : R script for fitting the hierarchical empirical bayes model introduced by Jeong et al. (2011)
    3. Total_Intensity_for_Interference.csv : This is a .csv file containing the total spectral intensities of all samples used in the study. These intensities are used in computing spectral interference
    4. ll_scoremat.rds : An R dataframe containing the same data as in CosineCorrelationLowResGCMSDataBase.csv
    5. totintensity.rds : An R dataframe containing the same data as in Total_Intensity_for_Interference.csv

Scripts/ contains all R scripts with code for replicating the analyses presented in the manuscript. The specific files contained within the Scripts/ directory are:
    1. dataset_processing_for_main_analyses.R : R script that is used to process datasets for analysis. This processing includes defining unique spectral IDs, converting missing truth annotation values to "Unknown", and computing spectrum-specific factors that are used in GMM extension models.
    2. fig2_fig3_wilcox_figsS3_S9_script.R : R script that contains code used to generate Figures 2 and 3 in the main manuscript, as well as code for computing the wilcoxon tests mentioned in the main manuscript. Also contained is code to generate Figures S3 - S9 in the supplement.  
    3. fig4_script.R : R script that contains code to generate Figure 4 in the main manuscript.
    4. figS1_script.R : R script that contains code to generate Figure S1 in the supplement.
    5. figS2_script.R : R script that contains code to generate Figure S3 in the supplement.
    6. libsize analyses.R : R script to perform library size simulation study described in main manuscript.
    7. main_analyses.R : R script to perform the main analyses in the main manuscript in which the standard GMM is compared against several model extensions.
    8. table_script.R : R script used to generate the values in Table 1 and 2
    9. template_sbatch_script_for_main_analyses.sh : Analyses were performed using high performance computing resources managed by SLURM. This is a template of the sbatch script used for running jobs using these resources. 

# Replicating Analyses

The following general steps may be implemented to replicate the analyses presented in the main text:

1. Run dataset_processing_for_main_analyses.R Follow the instructions provided within this script for defining directories and modify code with these updated paths where specified. 
2. (Optional) Assuming access to high-performance-computing (HPC) clusters, modify template_sbatch_script_for_main_analyses.sh according to the syntax of your institution/organization's workload / job submission manager. 
3. Update the directories specified in main_analyses.R / libsize_analyses.R so that they are consistent with those used in dataset_processing_for_main_analyses.R. 
4. If running analyses using HPC resources, submit the modified .sh template from (2) through your organization's workload manager. Otherwise, modify main_analyses.R / libsize_analyses.R for local use by manually specifying an integer for "myargs" on line 603 of main_analyses.R (line 681 of libsize_analyses.R). Any integer between 1 and 28 may be specified. Each choice of integer maps to a specific similarity metric.
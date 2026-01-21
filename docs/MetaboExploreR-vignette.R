## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install assistant, eval=FALSE--------------------------------------------
# source("https://raw.githubusercontent.com/Hszemray/MStargetR/master/R/install.R")
# install_MStargetR()

## ----load_library, eval=FALSE-------------------------------------------------
# library(MStargetR)

## ----project_setup, eval=FALSE------------------------------------------------
# # Create a temporary directory for the project
# project_dir <- tempdir()
# 
# # Create the raw_data subdirectory
# raw_data_dir <- file.path(project_dir, "raw_data")
# dir.create(raw_data_dir, recursive = TRUE)
# 
# # For the purpose of this vignette, we assume that the raw files are in the raw_data_dir.
# # In a real analysis, you would place your vendor-specific raw files (e.g., .wiff) here.

## ----msConvertR, eval=FALSE---------------------------------------------------
# msConvertR(input_directory = project_dir, output_directory = project_dir)

## ----mrm_example--------------------------------------------------------------
mrm_template_path <- system.file("extdata", "LGW_lipid_mrm_template_v1.tsv", package = "MStargetR")
mrm_template <- read.delim(mrm_template_path)
head(mrm_template)

## ----PeakForgeR, eval=FALSE---------------------------------------------------
# # Path to the project directory
# project_directory <- "path/to/your/project"
# 
# # List of MRM template files
# mrm_template_list <- list(system.file("extdata", "LGW_lipid_mrm_template_v1.tsv", package = "MStargetR"))
# 
# PeakForgeR(
#   user_name = "User",
#   project_directory = project_directory,
#   mrm_template_list = mrm_template_list,
#   QC_sample_label = "QC",
#   plateID_outputs = NULL
# )

## ----PeakForgeR_output_example------------------------------------------------
peakforger_report_path <- system.file("extdata", "Example_PeakForgeR_report.csv", package = "MStargetR")
peakforger_report <- read.csv(peakforger_report_path, check.names = FALSE)
head(peakforger_report)

## ----qcCheckR, eval=FALSE-----------------------------------------------------
# # In a real analysis, you would use the project directory where the PeakForgeR output is located.
# # qcCheckR can handle tsv and csv data inputs
# # See documentation ??MStargetR::qcCheckR for further information.
# 
# library(MStargetR)
# 
# #Load example mrm_template_list
#   file_path <- system.file("extdata",
#                            "LGW_lipid_mrm_template_v1.tsv",
#                            package = "MStargetR")
# 
#   sample_metadata_example <- read_tsv(file_path)
# 
# #Load example conc_guide
#   file_path <- system.file("extdata",
#                            "LGW_SIL_batch_Ultimate_2023_03_06.tsv",
#                            package = "MStargetR")
# 
#   sample_metadata_example <- read_tsv(file_path)
# 
# #Load example report file
#   file_path <- system.file("extdata",
#                            "Example_PeakForgeR_report.csv",
#                            package = "MStargetR")
# 
#   report_file <- read_csv(file_path)
# 
# #Run qcCheckR function
# qcCheckR(user_name = "user1",
#          project_directory = "path/to/project_directory",
#          mrm_template_list = list(v1 = list(
#                                     SIL_guide = path to/mrm_guide1.tsv,
#                                     conc_guide = path to/SIL_concentration_guide1.tsv),
#                                   ),
#          QC_sample_label = "qc",
#          sample_tags = c("sample","control", "qc"),
#          mv_threshold = 0.5) #default is  0.5 for 50\% missing values

## ----transition_checkR, eval= FALSE-------------------------------------------
# mrm_template_path <- system.file("extdata", "LGW_lipid_mrm_template_v1.tsv", package = "MStargetR")
# mrm_template_df <- read.delim(mrm_template_path, check.names = FALSE)
# head(mrm_template_df)
# 
# #Now lets run the function and see the output
# transition_checkR(mrm_template_df)

## ----compare_mrm_template_with_guide, eval = FALSE----------------------------
# mrm_template_path <- system.file("extdata", "LGW_lipid_mrm_template_v1.tsv", package = "MStargetR")
#   mrm_template_df <- read.delim(mrm_template_path, check.names = FALSE)
#   head(mrm_template_df)
# 
# conc_guide_path <- system.file("extdata","LGW_SIL_batch_103.tsv", package = "MStargetR")
#   conc_guide_df <- read.delim(conc_guide_path, check.names = FALSE)
#   head(conc_guide_df)
# 
# #Now lets run the function and see the output
# compare_mrm_template_with_guide(mrm_template_df, conc_guide_df)

## ----Multicohort, eval = FALSE------------------------------------------------
# # Path to the project directory
# project_directory <- "path/to/your/project"
# 
# # Create the raw_data subdirectory
# raw_data_dir <- file.path(project_dir, "raw_data")
# dir.create(raw_data_dir, recursive = TRUE)
# 
# # In a real analysis, you would place your vendor-specific raw files (e.g., .wiff and .wiff.scan) here.
# 
# #Convert vendor files to mzml
# msConvertR(input_directory = project_dir, output_directory = project_dir)
# 
# #Provide transition list paths to PeakForgeR
# #It will cycle test each transition list on a plate until the match is found
# PeakForgeR(
#   user_name = "User",
#   project_directory = project_directory,
#   mrm_template_list = list(v1 = "path to/mrm_guide1.tsv",
#                            v2 = "path to/mrm_guide2.tsv"
#                           ),
#   QC_sample_label = "LTR",
# 
# )
# 
# # Provide transition lists and their respective concentration guide paths to qcCheckR
# qcCheckR(user_name = "user",
#          project_directory = "path/to/project_directory",
#          mrm_template_list = list(v1 = list(
#                                     SIL_guide = path to/mrm_guide1.tsv,
#                                     conc_guide = path to/SIL_concentration_guide1.tsv),
#                                   v2 = list(
#                                     SIL_guide = path to/mrm_guide2.tsv,
#                                     conc_guide = path to/SIL_concentration_guide2.tsv)
#                                  ),
#          QC_sample_label = "qc",
#          sample_tags = c("sample","control", "qc"),
#          mv_threshold = 0.5) #default is  0.5 for 50\% missing values

## ----session info-------------------------------------------------------------
sessionInfo()


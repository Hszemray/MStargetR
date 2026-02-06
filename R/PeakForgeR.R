#'PeakForgeR
#' @name PeakForgeR
#' @description This function performs peak picking and integration via Skyline
#' in a Docker image. Allowing for usage across all major OS systems.
#'
#' We strongly recommend checking your mrm transition list
#' using MStargetR::TransitionCheckR prior to using it in PeakForgeR
#'
#' If the user has not used MStargetR::msConvertR to convert vendor files
#' please ensure you create a project folder containing sub folder
#' "msConvert_mzml_output" with mzml files for the project.
#' @param user_name A character string to identify user.
#' @param project_directory A path to project directory
#' @param mrm_template_list Path to MRM transition list,
#' must be in specified format. See examples and run load example
#' mrm_guide for structure. May contain more than one template for
#' multi-method projects.
#' @param QC_sample_label User specified tag to filter QC samples.
#' Character case is not sensitive
#'
#' E.g. "JANE_C5_URI_MS-LIPIDS_PLIP01_PLATE_3-PLASMA LTR_19.mzML"
#'
#' QC_sample_label = "LTR" to target files containing LTR for QC.
#' @param plateID_outputs A vector of character strings specifying plateIDs for
#' project. This parameter must only be specified by users who have not used
#' MStargetR::msConvertR..... Default is NULL
#'
#' These must match mzml files.
#' e.g. If you have two plates:
#'  - JANE_DOE_C5_URI_MS-LIPIDS_PLATE_1-PLASMA_sample_1.mzML
#'  - JANE_DOE_C5_URI_MS-LIPIDS_PLATE_2-PLASMA_sample_1.mzML
#'
#' An appropriate input would be:
#'  - plateID_outputs = c("JANE_DOE_C5_URI_MS-LIPIDS_PLATE_1",
#'                        "JANE_DOE_C5_URI_MS-LIPIDS_PLATE_2")
#' @return A curated project directory with sub folders for each plate containing Skyline exports.
#' @export
#' @examples
#' \dontrun{
#' #Load example mrm_guide
#'   file_path <- system.file("extdata", "LGW_lipid_mrm_template_v1.tsv", package = "MStargetR")
#'   example_mrm_template <- read_tsv(file_path)
#'
#' #Run PeakForgeR function
#' PeakForgeR(user_name = "Mad_max",
#'            project_directory = "USER/PATH/TO/PROJECT/DIRECTORY",
#'            mrm_template_list = list("User/path/to/user_mrm_guide_v1.tsv",
#'                                     "user/path/to/user_mrm_guide_v2.tsv"),
#'            QC_sample_label = "LTR",
#'            plateID_outputs = NULL
#'           )
#'}
#'
#' @details
#' \itemize{
#'  \item \strong{Input Validation:}
#'   \itemize{
#'    \item Validate project_directory
#'    \item Validate mrm_template_list
#'   }
#'  \item \strong{File Handling:}
#'   \itemize{
#'    \item Set plateIDs from either plate MStargetR::msConvertR or
#'          user specified plateID_outputs
#'   }
#'  \item \strong{Processing Plates:}
#'   \itemize{
#'    \item For each plateID:
#'     \itemize{
#'      \item Setup project structure
#'      \item Import mzml files
#'      \item QC optimised retention times
#'      \item QC optimised peak boundaries
#'      \item Peak picking/integration with Skyline MS through docker
#'     }
#'   }
#'  \item \strong{Final Cleanup:}
#'   \itemize{
#'    \item Archive raw files
#'    \item Message about availability of chromatograms and reports
#'   }
#' }
NULL

#' Import specific functions from packages
#' @keywords internal
#' @name PeakForgeR_import_external_functions
#' @importFrom future plan sequential multisession
#' @importFrom future.apply future_lapply
#' @importFrom parallel detectCores
NULL

PeakForgeR <- function(user_name,
                     project_directory,
                     mrm_template_list = NULL,
                     QC_sample_label = NULL,
                     plateID_outputs = NULL) {
  #Validate user
  if (!is.character(user_name) || nchar(user_name) == 0) {
    stop("Invalid user name. Please provide a valid character string.")
  }

  # Validate project_directory
  validate_project_directory(project_directory)

  # Set working directory to project directory
  setwd(project_directory)

  # Validate mrm_template_list
  validated_list <-  validate_mrm_template_list(mrm_template_list, user_name)
  ## If validate_mrm_template_list returned something, use it
  if (!is.null(validated_list)) {
    mrm_template_list <- validated_list
  }

  # Validate QC_sample_label
  if (!is.character(QC_sample_label) ||
      nchar(QC_sample_label) == 0) {
    stop(
      "Invalid QC_sample_label. Please ensure the parameter
         is of character type and length > 0"
    )
  }

  #Check install and run status of docker
  check_docker()

  # Set plateIDs
  file_paths <- list.files(project_directory)
  exclude_names <- c("raw_data", "msConvert_mzml_output", "all", "archive",
                     "error_log.txt","MStargetR_logs", "logs", "user_files")
  plateIDs <- file_paths[!basename(file_paths) %in% exclude_names]


  if (is.null(plateIDs) || length(plateIDs) == 0) {
    count_data <- list()
    for (plate in plateID_outputs) {
      mzml_files <- list.files(file.path(project_directory, plate, "data","mzml"))
      valid_count <- 0
      count_data[[plate]] <- length(mzml_files)
      valid_count <- valid_count + count_data[[plate]]
    }

    # Check for any plates with zero count
    zero_plates <- names(count_data)[unlist(count_data) == 0]

    if (length(zero_plates) > 0) {
      message(
        "The following plateIDs have zero associated mzML files:\n",
        paste(zero_plates, collapse = ", ")
      )
      stop(
        "Zero associated mzml files forbidden. Please check your input for plateID_outputs parameter."
      )
    }

    if (valid_count == length(mzml_files)) {
      message("Valid plateID_outputs have been provided")
      message(
        "plateID_outputs parameter break down:\n",
        paste(capture.output(str(count_data)), collapse = "\n")
      )
      plateIDs <- plateID_outputs
    } else {
      message(
        "Mismatch in provided plateID_outputs parameter and supplied mzML files.\n
                  Ensure plateID_outputs parameter match mzML files project/plate identifiers"
      )
      message(
        "plateID_outputs  parameter break down:\n",
        paste(capture.output(str(count_data)), collapse = "\n")
      )
      stop()
    }
  } else {
    message(
      "plateIDs gathered from existing project directory\n",
      "Plates for processing:\n ",
      paste(plateIDs, collapse = "\n")
    )
  }

  #Set failed/successful plates
  failed_plates <- c()
  successful_plates <- c()

  # Process each plate in parallel
  logs_dir <- file.path(project_directory, "MStargetR_logs")
  dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

  #Parallel settings
  future::plan(future::multisession)
  on.exit(future::plan(future::sequential), add = TRUE)
  total_cores <- parallel::detectCores()
  available_cores <- if (total_cores <= 2) 1 else total_cores - 2
  options(future.maxWorkers = available_cores)

  results <- future.apply::future_lapply(plateIDs, function(plateID) {
    start_time <- Sys.time()
    log_file <- file.path(logs_dir, paste0(plateID, "_MStargetR_log.txt"))
    library(MStargetR)

    # Helper to write a line to the log
    write_log <- function(text) {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      line <- paste0("[", timestamp, "] ", text, "\n")
      cat(enc2utf8(line), file = log_file, append = TRUE)
    }
    #Helper to write error to log
    log_error <- function(error_message, plateID) {
      write(error_message, file = log_file, append = TRUE)
    }

    tryCatch({
      write_log(paste("Starting processing for plate:", plateID))

      write_log("Step 1: Setting up project...")
      master_list <- PeakForgeR_setup_project(
        user_name,
        project_directory,
        plateID,
        mrm_template_list,
        QC_sample_label
      )
      write_log("Project setup complete.")

      write_log("Step 2: Importing mzML files...")
      master_list <- import_mzml(plateID, master_list)
      write_log("mzML import complete.")

      write_log("Step 3: Performing peak picking...")
      master_list <- peak_picking(plateID, master_list)
      write_log("Peak picking complete.")

      write_log(paste("Finished processing plate:", plateID))
      write_log("Status: SUCCESS")

      list(success = TRUE, plateID = plateID)

    }, error = function(e) {
      write_log(paste("Error during processing:", e$message))
      write_log("Status: FAILURE")

      log_error(paste("Error processing plate", plateID, ":", e$message), plateID)
      list(success = FALSE, plateID = plateID, error = e$message)
    })
  })

  # Extract successful and failed plate IDs
  successful_plates <- sapply(results, function(x) x$success && !is.null(x$plateID))
  failed_plates <- sapply(results, function(x) !x$success && !is.null(x$plateID))

  # Display summary
  message("\nProcessing complete.")

  if (any(successful_plates)) {
    message("\nPlates processed successfully:\n",
            paste(sapply(results[successful_plates], `[[`, "plateID"), collapse = "\n"))
  }

  if (any(failed_plates)) {
    message("\nPlates that failed to process:\n",
            paste(sapply(results[failed_plates], `[[`, "plateID"), collapse = "\n"))
  }

  if (!any(successful_plates)) {
    stop("All plates failed. Halting script.")
  }

  # Final cleanup and archiving
  archive_raw_files(project_directory)
  message(
    "\n Chromatograms and reports are now available per plate in ",
    paste(project_directory),
    ".\n Please run qcCheckR to calculate concentrations and QC the data"
  )
}

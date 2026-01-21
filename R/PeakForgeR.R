compute_results_anpc <- function(plateIDs, user_name, project_directory, mrm_template_list, QC_sample_label) {
  if (!(length(mrm_template_list) >= 1)) {
    stop("No MRM template provided in mrm_template_list")
  }

  #Log setup
  logs_dir <- file.path(project_directory, "MetaboExploreR_logs")
  dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

  results <- list()

  for (plateID in plateIDs) {
    message(paste("\nProcessing plateID:", plateID,"\n"))

    log_file <- file.path(logs_dir, paste0(plateID,"_Peakforger_log.txt"))

    if (!file.exists(log_file)) {
      file.create(log_file)
    }

    # Helper to write a line to the log
    write_log <- function(text) {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      line <- paste0("[", timestamp, "] ", text, "\n")
      cat(enc2utf8(line), file = log_file, append = TRUE)
    }

    # Helper to write error to log
    log_error <- function(error_message, plateID) {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      line <- paste0("[", timestamp, "] ERROR in plate ", plateID, ": ", error_message, "\n")
      cat(enc2utf8(line), file = log_file, append = TRUE)
    }

    write_log(paste0 ("Step 1: Setting up project:", plateID))
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

    ##mrm_template version selection
    if(length(master_list$templates$mrm_guides) > 1){
      ## Find all mzML files
      mzml_files <- list.files(
        path = master_list$project_details$project_dir,
        pattern = "\\.mzML$",
        recursive = TRUE,
        full.names = TRUE
      )
      ##filter for plate
      filtered_mzml <- mzml_files[grepl(master_list$project_details$plateID, mzml_files)]
      if (length(filtered_mzml) == 0) {
        stop("Error: No mzML files found in the specified directory.")
      }
      first_mzml <- filtered_mzml[1]

      ## Extract all templates except "by_plate"
      all_templates <- master_list$templates$mrm_guides
      all_templates <- all_templates[!names(all_templates) %in% "by_plate"]
      templates <- lapply(all_templates, function(x) x$mrm_guide)

      ## Run the matching function
      result <- match_templates_to_mzml(
        mrm_templates = templates,
        mzml_file = first_mzml,
        precursor_col = "Precursor Mz",
        product_col   = "Product Mz",
        name_col      = "Precursor Name",
        standardize   = TRUE,
        verbose       = TRUE
      )

      master_list$project_details$is_ver <- result$best_match_name
    }

    message(paste("Using method template: ",master_list$project_details$is_ver))
    write_log(paste("Using method template: ",master_list$project_details$is_ver))

    write_log("Step 3: Optimising retention times")
    plate_idx <- plateID
    master_list$templates$mrm_guides$by_plate[[plate_idx]] <- optimise_retention_times(master_list, plate_idx)
    write_log("Optimisation complete")

    write_log("Step 4 : Performing MSUT peak picking...")
    message ("Picking peaks and integrating...")
    tryCatch({
      mzml_dir <- file.path(project_directory, plateID, "data", "mzml")
      mzml_files <- list.files(mzml_dir, pattern="\\.mzML$", ignore.case=TRUE, full.names=TRUE)
      if (length(mzml_files) == 0L) stop("No mzML files found for plate")

      template_path <- mzml_files[[1]]
      template_bin <- readBin(template_path, "raw", n = file.info(template_path)$size)
      template_file <- msut::parse_mzml(template_bin)
      template_df <- msut::bin_to_df(template_file)
      chroms <- template_df$Ok$run$chromatograms
      targets_path <-  master_list$templates$mrm_guides[[master_list$project_details$is_ver]]$mrm_guide
      targets <- master_list$templates$mrm_guides[[master_list$project_details$is_ver]]$mrm_guide
      transitions <- build_transitions_list(chroms, targets)

      total_cores <- parallel::detectCores()
      available_cores <- if (total_cores <= 2) 1 else total_cores - 2

      plate_rows <- vector("list", length(mzml_files))
      for (i in seq_along(mzml_files)) {
        mzml_path <- mzml_files[[i]]
        bin <- readBin(mzml_path, "raw", n = file.info(mzml_path)$size)
        file <- msut::parse_mzml(bin)
        peaks <- msut::get_peaks_from_chrom(
          file,
          transitions,
          auto_baseline = TRUE,
          auto_noise = FALSE,
          cores = available_cores
        )
        if (!is.data.frame(peaks)) peaks <- as.data.frame(peaks)
        peaks$name <- tools::file_path_sans_ext(basename(mzml_path))
        peaks$plate <- plateID
        peaks <- peaks[, c("name", "plate", setdiff(names(peaks), c("name","plate")))]
        plate_rows[[i]] <- peaks
      }

      plate_table <- if (length(plate_rows)) do.call(rbind, plate_rows) else data.frame()
      out_path <- file.path(project_directory, plateID, "data", "PeakForgeR", paste0(plateID, ".tsv"))
      readr::write_tsv(plate_table, out_path)
      write_log(paste("Finished processing plate:", plateID))
      message(paste("Finished processing plate:", plateID))
      results[[plateID]] <- list(success = TRUE, plateID = plateID)
    }, error = function(e) {
      write_log(sprintf("anpc method: plate %s failed: %s", plateID, conditionMessage(e)))
      message(sprintf("anpc method: plate %s failed: %s", plateID, conditionMessage(e)))
      results[[plateID]] <- list(success = FALSE, plateID = plateID, error = e$message)
    })
    write_log("Peak picking complete.")
  }
  return(results)
}

compute_results_skyline <- function(plateIDs, user_name, project_directory, mrm_template_list, QC_sample_label) {
  #Check install and run status of docker
  check_docker()

  #Parallel settings
  future::plan(future::multisession)
  on.exit(future::plan(future::sequential), add = TRUE)
  total_cores <- parallel::detectCores()
  available_cores <- if (total_cores <= 2) 1 else total_cores - 2
  options(future.maxWorkers = available_cores)

  #Log setup
  logs_dir <- file.path(project_directory, "MetaboExploreR_logs")
  dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

  #Process
  future.apply::future_lapply(plateIDs, function(plateID) {
    start_time <- Sys.time()
    library(MetaboExploreR)

      log_file <- file.path(logs_dir, paste0(plateID,"_Peakforger_log.txt"))

      if (!file.exists(log_file)) {
        file.create(log_file)
      }

      # Helper to write a line to the log
      write_log <- function(text) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        line <- paste0("[", timestamp, "] ", text, "\n")
        cat(enc2utf8(line), file = log_file, append = TRUE)
      }

      # Helper to write error to log
      log_error <- function(error_message, plateID) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        line <- paste0("[", timestamp, "] ERROR in plate ", plateID, ": ", error_message, "\n")
        cat(enc2utf8(line), file = log_file, append = TRUE)
      }

    tryCatch({
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

      results[[plateID]] <- list(success = TRUE, plateID = plateID)
    }, error = function(e) {
      write_log(paste("Error during processing:", e$message))
      write_log("Status: FAILURE")

      log_error(paste("Error processing plate", plateID, ":", e$message), plateID)
      results[[plateID]] <- list(success = FALSE, plateID = plateID, error = e$message)
    })
  })
  return(results)
}

#'PeakForgeR
#'
#' @description This function performs peak picking and integration via Skyline
#' in a Docker image. Allowing for usage across all major OS systems.
#'
#' We strongly recommend checking your mrm transition list
#' using MetaboExploreR::TransitionCheckR prior to using it in PeakForgeR
#'
#' If the user has not used MetaboExploreR::msConvertR to convert vendor files
#' please ensure you create a project folder containing sub folder
#' "msConvert_mzml_output" with mzml files for the project.
#' @export
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
#' MetaboExploreR::msConvertR..... Default is NULL
#'
#' These must match mzml files.
#' e.g. If you have two plates:
#'  - JANE_DOE_C5_URI_MS-LIPIDS_PLATE_1-PLASMA_sample_1.mzML
#'  - JANE_DOE_C5_URI_MS-LIPIDS_PLATE_2-PLASMA_sample_1.mzML
#'
#' An appropriate input would be:
#'  - plateID_outputs = c("JANE_DOE_C5_URI_MS-LIPIDS_PLATE_1",
#'                        "JANE_DOE_C5_URI_MS-LIPIDS_PLATE_2")
#' @param method Allows user to select the peak picking method to use.
#' e.g. "skyline", "anpc", or "openms"
#' NOTE IF USING "openms" please install conda using the following method below"
#'  one-time per machine
#'  run "peakforger_setup_openms(envname="r-openms-env", python_version="3.10", backend="conda", install=TRUE)"
#'  restart R session and then run the pipeline

#' @return A curated project directory with sub folders for each plate containing Skyline exports.
#' @examples
#' \dontrun{
#' #Load example mrm_guide
#'   file_path <- system.file("extdata", "LGW_lipid_mrm_template_v1.tsv", package = "MetaboExploreR")
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
#'    \item Set plateIDs from either plate MetaboExploreR::msConvertR or
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
#' @importFrom future plan sequential multisession
#' @importFrom future.apply future_lapply
#' @importFrom parallel detectCores
PeakForgeR <- function(user_name,
                       project_directory,
                       mrm_template_list = NULL,
                       QC_sample_label = NULL,
                       plateID_outputs = NULL,
                       method = c("skyline","openms", "anpc")) {
  method <- match.arg(method)
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

  # Set plateIDs
  file_paths <- list.files(project_directory)
  exclude_names <- c("raw_data", "msConvert_mzml_output", "all", "archive",
                     "error_log.txt","MetaboExploreR_logs", "logs", "user_files")
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

  if (method == "skyline") {
    results <- compute_results_skyline(plateIDs, user_name, project_directory, mrm_template_list, QC_sample_label)
  } else if (method == "anpc") {
    results <- compute_results_anpc(plateIDs, user_name, project_directory, mrm_template_list, QC_sample_label)
  } else if (method == "openms") {
    openms_use_backend("conda", envname = "r-openms-env")
    results <- compute_results_openms(plateIDs, user_name, project_directory, mrm_template_list, QC_sample_label, template_vote_n_files = 1L, cfg = PF_CFG)
  }else {
    stop("Choose one of the valid processing methods")
  }


  # Separate successful and failed plates
  success_flags <- sapply(results, function(x) x$success)
  plate_ids <- sapply(results, function(x) x$plateID)
  successful_plates <- plate_ids[success_flags]
  failed_plates     <- plate_ids[!success_flags]

  successful_plates
  failed_plates

  message("\nProcessing complete.")

  if (length(successful_plates) > 0) {
    message("\nPlates processed successfully:\n",
            paste(successful_plates, collapse = "\n"))
  } else {
    message("\nNo plates processed successfully.")
  }

  if (length(failed_plates) > 0) {
    message("\nPlates that failed to process:\n",
            paste(failed_plates, collapse = "\n"))
  } else {
    message("\nNo plates failed.")
  }

  if (length(successful_plates) == 0) {
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

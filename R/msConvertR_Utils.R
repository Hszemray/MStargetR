# msConvertR_Utils ----
# sub functions for msConvertR function


#' Import specific functions from packages
#' @keywords internal
#' @name msConvertR_import_external_functions
#' @importFrom stringr str_subset

NULL
#'
#mzML File Conversion Functions ----
#' Validate input Directory
#'
#' This function checks if the `input_directory` parameter is a single string and if the specified directory exists.
#' @keywords internal
#' @param input_directory A character string representing the path to the project directory.
#' @return TRUE if the validation is successful, otherwise an error is thrown.
#' @examples
#' \dontrun{
#' validate_project_directory("input_directory")
#' }
validate_input_directory <- function(input_directory) {
  # Check if input_directory is a single string
  if (!is.character(input_directory) ||
      length(input_directory) != 1) {
    stop("input_directory must be a single string.")
  }

  # Check if the specified directory exists
  if (!dir.exists(input_directory)) {
    stop("The specified project directory does not exist.")
  }

  # Return TRUE if validation is successful
  message(paste("Accessing project directory ", input_directory))
}

###Primary Function----
#' msConvertR_mzml_conversion
#'
#' This function converts raw vendor files to mzML format using ProteoWizard's msconvert tool, restructures directories, and updates the script log.
#' @keywords internal
#' @param input_directory Directory path for project folder
#' @param output_directory Directory path for project folder if different from
#' input directory.
#' @param plateIDs vector of vendor files names to be converted
#' @param vendor_extension_patterns character string of vendor file extensions.
#' @return Converted mzml files.
#' @examples
#' \dontrun{
#' msConvertR_mzml_conversion(input_directory, output_directory, plateIDs)
#' }
msConvertR_mzml_conversion <- function(input_directory,
                                       output_directory,
                                       plateIDs,
                                       vendor_extension_patterns) {
  msConvertR_set_working_directory(input_directory)
  msConvertR_setup_project_directories(output_directory, plateIDs)
  commands <- msConvertR_construct_command_for_terminal(input_directory, output_directory, plateIDs)
  msConvertR_execute_command(commands, output_directory, plateIDs)
  msConvertR_restructure_directory(output_directory, plateIDs, vendor_extension_patterns)
}

###Sub Functions----

#' msConvertR_setup_project_directories
#'
#' This function sets up project directories for each plate ID.
#' @keywords internal
#' @param output_directory Output directory for mzml files.
#' @param plateIDs vector of vendor files names being converted
#' @return None. The function sets up directories.
#' @examples
#' \dontrun{
#' msConvertR_setup_project_directories(output_directory, plateIDs)
#' }
msConvertR_setup_project_directories <- function(output_directory, plateIDs) {
  for (plateID in plateIDs) {
    dir.create(file.path(output_directory,"user_files"), showWarnings = FALSE)
    base_path <- file.path(output_directory, plateID)
    dir.create(base_path, showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(base_path, "data"),
               showWarnings = FALSE,
               recursive = TRUE)
    dir.create(
      file.path(base_path, "data", "mzml"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    dir.create(
      file.path(base_path, "data", "rda"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    dir.create(
      file.path(base_path, "data", "PeakForgeR"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    dir.create(
      file.path(base_path, "data", "raw_data"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    dir.create(
      file.path(base_path, "data", "batch_correction"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    dir.create(
      file.path(base_path, "html_report"),
      showWarnings = FALSE,
      recursive = TRUE
    )
  }
}


#' msConvertR_set_working_directory
#'
#' This function sets the working directory to the project directory.
#' @keywords internal
#' @param output_directory Directory path for the project folder.
#' @return None. The function sets the working directory.
#' @examples
#' \dontrun{
#' msConvertR_set_working_directory("path/to/output_directory")
#' }
msConvertR_set_working_directory <- function(output_directory) {
  setwd(output_directory)
}

#' msConvertR_construct_commands_for_terminal
#'
#' This function constructs the command for terminal to convert files to mzML format.
#' @keywords internal
#' @param input_directory path to input directory containing vendor files
#' @param output_directory path to output directory.
#' @param plateIDs The names of vendor files to convert
#' @return The constructed command string.
#' @examples
#' \dontrun{
#' command <- msConvertR_construct_command_for_terminal(path/to/input/directory,
#'                                                      "path/to/output_directory")
#' }
msConvertR_construct_command_for_terminal <- function(input_directory, output_directory, plateIDs) {

  commands <- lapply(plateIDs, function(plateID) {
    input_path <- normalizePath(file.path(input_directory, "raw_data"), mustWork = FALSE)
    output_dir <- normalizePath(file.path(output_directory, plateID, "data", "mzml"), mustWork = FALSE)

    docker_image <- "proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:skyline_daily_25.1.1.270-67f3e15"
    container_data_path <- "/data"
    output_data_path <- "/output"

    files <- list.files(input_path)
    file_name <- files[grepl(plateID, files)]
    file_name <- file_name[!grepl(".scan", file_name)]

    sprintf(
      'docker run --rm -v "%s:%s" -v "%s:%s" %s wine msconvert %s -o %s',
      input_path,
      container_data_path,
      output_dir,
      output_data_path,
      docker_image,
      file.path(container_data_path, file_name),
      output_data_path
    )
  })

  return(commands)
}

#' msConvertR_execute_command
#'
#' This function executes the command to convert files to mzML format.
#' @keywords internal
#' @param commands The command string to execute for each plate.
#' @return None. The function executes the command.
#' @examples
#' \dontrun{
#' msConvertR_execute_command(command)
#' }
msConvertR_execute_command <- function(commands, output_directory, plateIDs) {
  message("Converting vendor files:\n", paste(plateIDs,collapse = "\n"))

  logs_dir <- file.path(output_directory, "MStargetR_logs")
  dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

  future::plan(future::multisession)
  on.exit(future::plan(future::sequential), add = TRUE)
  total_cores <- parallel::detectCores()
  available_cores <- if (total_cores <= 2) 1 else total_cores - 2
  options(future.maxWorkers = available_cores)

  # Start conversion tasks as futures
  futures <- lapply(seq_along(commands), function(i) {
    future::future({
      plateID <- plateIDs[i]
      cmd <- commands[[i]]
      log_file <- file.path(logs_dir, paste0(plateID, "_MStargetR_log.txt"))

      start_time <- Sys.time()

      output <- system(cmd, intern = TRUE)
      exit_status <- attr(output, "status")
      success <- is.null(exit_status) || exit_status == 0

      writeLines(enc2utf8(c(
        paste("Start time:", start_time),
        paste("Command:", cmd),
        "Output:",
        output,
        paste("End time:", Sys.time()),
        paste("Status:", if (success) "SUCCESS" else "FAILURE")
      )), con = log_file)

      message(sprintf("Finished conversion for %s — %s", plateID, if (success) "SUCCESS" else "FAILURE"))
      success
    })
  })

  # Start looping progress bar
  pb <- txtProgressBar(min = 0, max = 100, style = 3)
  i <- 0
  repeat {
    i <- (i + 10) %% 110
    setTxtProgressBar(pb, ifelse(i == 0, 100, i))
    Sys.sleep(3)

    # Check if all futures are resolved
    if (all(sapply(futures, future::resolved))) {
      break
    }
  }
  close(pb)
  # Collect results
  results <- lapply(futures, future::value)
  names(results) <- plateIDs

  return(results)
}

#' msConvertR_restructure_directory
#'
#' This function restructures the directory by moving raw_data and mzML files to correct locations.
#' @keywords internal
#' @param output_directory Output directory where the mzML files will be stored.
#' @param plateIDs filenames for plates being converted with no extension.
#' @param vendor_extension_patterns vector of file extensions for vendor files
#' @return The updated master list object with mzML file paths.
#' @examples
#' \dontrun{
#' master_list <- msConvertR_restructure_directory(output_directory,
#'                                                 plateIDs,
#'                                                 vendor_extension_patterns)
#' }
msConvertR_restructure_directory <- function(output_directory,
                                             plateIDs,
                                             vendor_extension_patterns) {
  for (plateID in plateIDs) {
    # Define key paths
    project_dir <- output_directory
    raw_data_dir <- file.path(output_directory, "raw_data")
    mzml_output_dir <- file.path(output_directory, "msConvert_mzml_output")
    plate_data_dir <- file.path(output_directory, plateID, "data")
    raw_data_dest <- file.path(plate_data_dir, "raw_data")
    mzml_dest <- file.path(plate_data_dir, "mzml")

    # Create destination directories if they don't exist
    dir.create(raw_data_dest,
               recursive = TRUE,
               showWarnings = FALSE)
    dir.create(mzml_dest, recursive = TRUE, showWarnings = FALSE)

    # Find raw data files matching plateID
    raw_files <- list.files(path = raw_data_dir,
                            pattern = vendor_extension_patterns,
                            full.names = TRUE)
    matched_raw_files <- str_subset(raw_files, plateID)

    # Copy individual raw files (non-directories)
    file_info <- file.info(matched_raw_files)
    raw_file_paths <- matched_raw_files[!file_info$isdir]
    file.copy(from = raw_file_paths,
              to = raw_data_dest,
              recursive = FALSE)

    # Copy entire .d directories
    d_dirs <- matched_raw_files[file_info$isdir]
    for (dir in d_dirs) {
      dest_dir <- file.path(raw_data_dest)
      dir.create(dest_dir,
                 recursive = TRUE,
                 showWarnings = FALSE)
      file.copy(from = dir,
                to = dest_dir,
                recursive = TRUE)
    }

    # Filter and copy mzML files
    mzml_files <- list.files(path = mzml_output_dir,
                             pattern = "\\.mzML$",
                             full.names = TRUE)
    mzml_files <- mzml_files[!grepl("-COND\\d+_|-BLANK\\d+_|-BLANK_\\d+|-ISTDs_\\d+", mzml_files, ignore.case = TRUE)]
    matched_mzml <- str_subset(mzml_files, plateID)
    file.copy(from = matched_mzml,
              to = mzml_dest,
              recursive = FALSE)

  }
}

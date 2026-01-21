# config.R

#Update Script Log Functions----

##Primary Function----
#' Update Script Log
#'
#' This function updates the script log in the `master_list` object by capturing the current time, calculating the runtime for the current section, and creating a message for the log.
#' @keywords internal
#' @param master_list A list containing project details and script log information.
#' @param section_name A string representing the name of the current section.
#' @param previous_section_name A string representing the name of the previous section.
#' @param next_section_name A string representing the name of the next section.
#' @return The updated `master_list` object with the new log information.
#' @examples
#' \dontrun{
#' master_list <- list(
#'                 project_details = list(
#'                   script_log = list(
#'                     timestamps = list(
#'                       start_time = Sys.time()
#'                     ),
#'                     runtimes = list(),
#'                     messages = list()
#'                    )
#'                  )
#'                )
#'
#' update_script_log(master_list, "section_1", "start_time", "section_2")
#' }
update_script_log <- function(master_list,
                              section_name,
                              previous_section_name,
                              next_section_name) {
  validate_previous_section(master_list, previous_section_name)
  master_list <- capture_current_time(master_list, section_name)
  master_list <- calculate_runtime(master_list, section_name, previous_section_name)
  master_list <- calculate_total_runtime(master_list, section_name)
  master_list <- create_message(master_list, section_name, next_section_name)
  master_list <- print_message(master_list, section_name)
  return(master_list)
}

##Seconday Functions----
#' Validate Previous Section
#'
#' This function validates the previous section name.
#' @keywords internal
#' @param master_list A list containing project details and script log information.
#' @param previous_section_name A string representing the name of the previous section.
#'
#' @return updated master list

validate_previous_section <- function(master_list, previous_section_name) {
  if (!previous_section_name %in% names(master_list$project_details$script_log$timestamps)) {
    stop("Invalid previous_section_name")
  }
}

#' Capture Current Time
#'
#' This function captures the current time for the given section.
#' @keywords internal
#' @param master_list A list containing project details and script log information.
#' @param section_name A string representing the name of the current section.
#'
#' @return updated master list
capture_current_time <- function(master_list, section_name) {
  master_list$project_details$script_log$timestamps[[section_name]] <- Sys.time()
  return(master_list)
}

#' Calculate Runtime
#'
#' This function calculates the runtime for the given section.
#' @keywords internal
#' @param master_list A list containing project details and script log information.
#' @param section_name A string representing the name of the current section.
#' @param previous_section_name A string representing the name of the previous section.
#'
#' @return None

calculate_runtime <- function(master_list,
                              section_name,
                              previous_section_name) {
  master_list$project_details$script_log$runtimes[[section_name]] <- difftime(
    master_list$project_details$script_log$timestamps[[section_name]],
    master_list$project_details$script_log$timestamps[[previous_section_name]],
    units = "mins"
  )
  return(master_list)
}

#' Calculate Total Runtime
#'
#' This function calculates the total runtime from the start.
#' @keywords internal
#' @param master_list A list containing project details and script log information.
#' @param section_name A string representing the name of the current section.
#'
#' @return None
calculate_total_runtime <- function(master_list, section_name) {
  master_list$project_details$script_log$runtimes$total_runtime <- difftime(
    master_list$project_details$script_log$timestamps[[section_name]],
    master_list$project_details$script_log$timestamps$start_time,
    units = "mins"
  )
  return(master_list)
}

#' Create Message
#'
#' This function creates a message for the log.
#' @keywords internal
#' @param master_list A list containing project details and script log information.
#' @param section_name A string representing the name of the current section.
#' @param next_section_name A string representing the name of the next section.
#'
#' @return None

create_message <- function(master_list,
                           section_name,
                           next_section_name) {
  master_list$project_details$script_log$messages[[section_name]] <- paste0(
    "\n",
    toupper(gsub("_", " ", section_name)),
    " complete!",
    "\n\n Section runtime: ",
    signif(
      as.numeric(master_list$project_details$script_log$runtimes[[section_name]]),
      digits = 3
    ),
    " minutes",
    "\n\n Total runtime: ",
    signif(
      as.numeric(
        master_list$project_details$script_log$runtimes$total_runtime
      ),
      digits = 3
    ),
    " minutes",
    "\n",
    "\nInitialising: ",
    toupper(gsub("_", " ", next_section_name)),
    "...."
  )
  return(master_list)
}

#' Print Message
#'
#' This function prints the message for the log.
#' @keywords internal
#' @param master_list A list containing project details and script log information.
#' @param section_name A string representing the name of the current section.
#'
#' @return None

print_message <- function(master_list, section_name) {
  msg <- master_list$project_details$script_log$messages[[section_name]]
  if (is.null(msg)) {
    stop("No message found for the given section.")
  }
  message(msg)
  return(master_list)
}



#Validate Parameter Functions----
#' Validate Project Directory
#'
#' This function checks if the `project_directory` parameter is a single string and if the specified directory exists.
#' @keywords internal
#' @param project_directory A character string representing the path to the project directory.
#' @return TRUE if the validation is successful, otherwise an error is thrown.
#' @examples
#' \dontrun{
#' validate_project_directory("path/to/project_directory")
#' }
validate_project_directory <- function(project_directory) {
  # Check if project_directory is a single string
  if (!is.character(project_directory) ||
      length(project_directory) != 1) {
    stop("project_directory must be a single string.")
  }

  # Check if the specified directory exists
  if (!dir.exists(project_directory)) {
    stop("The specified project directory does not exist.")
  }

  # Return TRUE if validation is successful
  message(paste("Accessing project directory ", project_directory))
}

#' validate_master_list_project_directory
#' This function validates the existence of the project directory specified in the master list.
#' @keywords internal
#' @param master_list A list containing project details.
#' @return Stops execution if the project directory does not exist.
#' @examples
#' \dontrun{
#' validate_master_list_project_directory(master_list)
#' }
validate_master_list_project_directory <- function(master_list) {
  if (!dir.exists(master_list$project_details$project_dir)) {
    stop(
      paste(
        "Project directory does not exist:",
        master_list$project_details$project_dir
      )
    )
  }
}




#' Validate MRM Template List
#'
#' This function checks if the `mrm_template_list` parameter is valid and contains required columns.
#' @keywords internal
#' @param mrm_template_list A list of character strings or a named list of data frames representing MRM templates.
#' @param user_name A character string identifying the user.
#' @return NULL or an ANPC mrm_template_list if mrm_template_list is NULL and user_name is ANPC
#' @examples
#' \dontrun{
#' validate_mrm_template_list(list("path/to/template1.csv", "path/to/template2.csv"), "user")
#' }
validate_mrm_template_list <- function(mrm_template_list, user_name) {
  required_columns <- c(
    "Molecule List Name",
    "Precursor Name",
    "Precursor Mz",
    "Precursor Charge",
    "Product Mz",
    "Product Charge",
    "Explicit Retention Time",
    "Explicit Retention Time Window",
    "Note",
    "control_chart"
  )

  if (user_name == "ANPC") {
    if (is.null(mrm_template_list)) {
      mrm_template_list <- list(
        v1 = system.file("extdata", "LGW_lipid_mrm_template_v1.tsv", package = "MStargetR"),
        v2 = system.file("extdata", "LGW_lipid_mrm_template_v2.tsv", package = "MStargetR"),
        v4 = system.file("extdata", "LGW_lipid_mrm_template_v4.tsv", package = "MStargetR")
      )
      message("mrm_template validation complete")
      return(mrm_template_list)
    } else {
      message("mrm_template validation complete")
      return(NULL)
    }
  }

  if (is.null(mrm_template_list)) {
    stop("Please provide a valid mrm_template_list.")
  }

  if (!is.list(mrm_template_list)) {
    stop("mrm_template_list must be a list.")
  }

  for (version in names(mrm_template_list)) {
    version_list <- mrm_template_list[[version]]

    if (!is.data.frame(version_list)) {
      stop(
        paste(
          "Each version in mrm_template_list must be a data frame. Problem with:",
          version
        )
      )
    }

    missing_columns <- setdiff(required_columns, colnames(version_list))
    if (length(missing_columns) > 0) {
      stop(paste(
        "Missing required columns in version",
        version,
        ":\n",
        paste(missing_columns, collapse = "\n")
      ))
    }
  }
  message("mrm_template validation complete")
  return(NULL)
}

#' Log Error to File
#'
#' This function logs error messages to a file named `error_log.txt`.
#' @keywords internal
#' @param error_message A character string representing the error message to be logged.
#' @return None. The function writes the error message to the log file.
#' @examples
#' \dontrun{
#' log_error("An error occurred while processing the data.")
#' }
log_error <- function(error_message, plateID) {
  log_file <- file.path(project_directory, "MStargetR_logs", paste0(plateID, "_MStargetR_log.txt"))
  write(error_message, file = log_file, append = TRUE)
}


#' Validate qcCheckR mrm template list
#'
#' This function validates the mrm_template_list list by checking the column headers and ensuring there are no NA or NULL values in the SIL_guide and conc_guide files.
#' @keywords internal
#' @param master_list A list containing all project details and data
#' @return TRUE if validation passes. Stops execution if validation fails.
#' @examples
#' \dontrun{
#' validate_qcCheckR_mrm_template_list(master_list)
#' }
validate_qcCheckR_mrm_template_list <- function(master_list) {
  mrm_template_list <- master_list$templates$mrm_guides

  if (!is.list(mrm_template_list)) {
    stop("mrm_template_list must be a list.")
  }

  for (version in names(mrm_template_list)) {
    version_list <- mrm_template_list[[version]]

    # Check that version is a list
    if (!is.list(version_list)) {
      stop(paste(
        "Each version in mrm_template_list must be a list. Problem with:",
        version
      ))
    }

    for (guide in c("SIL_guide", "conc_guide")) {
      if (!guide %in% names(version_list)) {
        stop(paste("Missing", guide, "in version", version))
      }

      # Check required columns
      if (guide == "SIL_guide") {
        required_columns <- c(
          "Molecule List Name",
          "Precursor Name",
          "Precursor Mz",
          "Precursor Charge",
          "Product Mz",
          "Product Charge",
          "Explicit Retention Time",
          "Explicit Retention Time Window",
          "Note",
          "control_chart"
        )

        data <- mrm_template_list[[version]][[guide]]

        if (!all(required_columns %in% colnames(data))) {
          #store missing columns
          missing_columns <- setdiff(required_columns, colnames(data))
          stop(paste(
            guide,
            "for version",
            version,
            "\n Missing required columns: ",
            paste(missing_columns, collapse = "\n")
          ))
        }


        # Define the columns you want to check
        check_cols <- c(
          "Molecule List Name",
          "Precursor Name",
          "Precursor Mz",
          "Precursor Charge",
          "Product Mz",
          "Product Charge",
          "Explicit Retention Time",
          "Explicit Retention Time Window",
          "control_chart"
        )

        # Subset the data to only include the filtered check columns
        data_to_check <- data[, check_cols, drop = FALSE]

        # Check for NA or NULL values
        if (any(is.na(data_to_check))) {
          stop(paste("NA values found in", guide, "for version", version))
        }

        transition_result <- transition_checkR(data)
        if (is.data.frame(transition_result)) {
          stop(paste("Non-unique transitions found in version", version, "\n",
                     paste(capture.output(print(transition_result)), collapse = "\n")))
        }

      } else if (guide == "conc_guide") {
        required_columns <- c("concentration_factor", "SIL_name")

        data <- mrm_template_list[[version]][[guide]]

        if (!all(required_columns %in% colnames(data))) {
          missing_columns <- setdiff(required_columns, colnames(data))
          stop(paste(
            guide,"for version",version,
            "\n Missing required columns: ",
            paste(missing_columns, collapse = "\n")
          ))
        }
        sil_guide <- version_list[["SIL_guide"]]
        compare_result <- compare_mrm_template_with_guide(sil_guide, data)
        if (is.character(compare_result)) {
          stop(paste("Unmatched Note values in version", version, "\n",
                     paste(compare_result, collapse = "\n")))
        }
      }
    }
  }

  message("Validation passed: mrm_template_list structure and contents are valid.")
  return(TRUE)
}

# Check Docker----
#'
#' Function to check Docker installation, daemon, and containers
#' @keywords internal
#' @return Returns current status of docker
#' @examples
#' \dontrun{
#' check_docker()
#' }
# Function to check Docker installation, daemon status, and container execution
check_docker <- function() {
  # Check if Docker is installed
  docker_installed <- tryCatch({
    system("docker --version",
           intern = TRUE,
           ignore.stderr = TRUE)
    TRUE
  }, error = function(e) {
    FALSE
  })

  if (!docker_installed) {
    message("Docker is NOT installed. Please install Docker to proceed")
    message("Directing you to Docker website for install")
    Sys.sleep(5)
    browseURL("https://www.docker.com/products/docker-desktop/")
    stop("Execution halted due to missing Docker installation")
    return(invisible(FALSE))
  }

  # Docker is installed now check it's running
  docker_container_status <- system("docker run hello-world")

  if (docker_container_status == 0) {
       message("Pulling proteowizard docker...  ")
    proteowizard_status <- system("docker pull proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:3.0.25114-e35aac0")
    if (proteowizard_status == 0) {
      message("Successfully pulled proteowizard docker!")
    } else{
      stop("Awwww snap an error occured during pull!")
    }

  } else{
    stop(
      "\n!!!Execution halted!!!
            \nDocker is installed, but NOT running.
            \nPlease open Docker Application.
            \nIf docker is open please restart the application."
    )
  }
}

# Validate Raw files----
#' Validate Raw files
#'
#' This function checks project directories contains vendor files.
#' @keywords internal
#' @param input_directory directory path for vendor file locations
#' @return validated paths and returns message on outcome of check.
#' @examples
#' \dontrun{
#' all_file_paths <- validate_file_types(input_directory)
#' }
validate_file_types <- function(input_directory) {
  file_path <- file.path(input_directory, "raw_data")
  files <- list.files(path = file_path, full.names = TRUE)

  validated_files <- c()
  invalid_files <- c()

  for (file in files) {
    if (grepl("\\.wiff$", file)) {
      scan_file <- paste0(file, ".scan")
      if (file.exists(scan_file)) {
        validated_files <- c(validated_files, file)
      } else {
        message("Missing .wiff.scan for: ", basename(file))
        invalid_files <- c(invalid_files, file)
      }
    # } else if (grepl("\\.wiff2$", file)) {
    #   validated_files <- c(validated_files, file)
    } else if (grepl("\\.raw$", file, ignore.case = TRUE)) {
      validated_files <- c(validated_files, file)
    } else if (grepl("\\.d$", file) && dir.exists(file)) {
      validated_files <- c(validated_files, file)
    } else if (grepl("\\.baf$", file)) {
      validated_files <- c(validated_files, file)
    } else if (grepl("\\.fid$", file)) {
      validated_files <- c(validated_files, file)
    } else if (grepl("\\.yep$", file)) {
      validated_files <- c(validated_files, file)
    } else if (grepl("\\.tsf$", file)) {
      validated_files <- c(validated_files, file)
    } else if (grepl("\\.tdf$", file)) {
      validated_files <- c(validated_files, file)
    } else if (grepl("\\.mbi$", file)) {
      validated_files <- c(validated_files, file)
    } else if (grepl("\\.qgd$", file) ||
               grepl("\\.qgb$", file) || grepl("\\.qgm$", file)) {
      validated_files <- c(validated_files, file)
    } else if (grepl("\\.lcd$", file) ||
               grepl("\\.lcdproj$", file)) {
      validated_files <- c(validated_files, file)
    } else if (grepl("\\.uep$", file) ||
               grepl("\\.sdf$", file) || grepl("\\.dat$", file)) {
      validated_files <- c(validated_files, file)
    } else if (grepl("\\.wcf$", file) ||
               grepl("\\.wproj$", file) || grepl("\\.wdata$", file)) {
      validated_files <- c(validated_files, file)
    } else {
      message("Unsupported file type found: ", basename(file))
      invalid_files <- c(invalid_files, file)
    }
  }


  if (length(validated_files) > 0) {
    message("Returning validated files for processing:\n",
            paste(basename(validated_files), collapse = "\n"))
    return(validated_files)
  }

  if (length(invalid_files) > 0) {
    message(
      "Removed following unsupported files:\n",
      paste(invalid_wiff_files, collapse = "\n")
    )
  }

}

# special character replacement for mrm_templates----
#' replace_precursor_symbols
#'
#' This function replaces forward or backwards slashes in 'Precursor Name' while preserving the original naming convention
#' This is due to skyline cmd being unable to handle the special character
#' @keywords internal
#' @param mrm_template dataframe of transitions (mrm_template) for SkylineR or qcCheckR
#' @return Updated mrm_template with special characters replaced in 'Precursor Name' and 'Note', while the original names are preserved for the columns in original_col
#' @examples
#' \dontrun{
#' replace_precursor_symbols(mrm_template, columns = c("Precursor Name", "Note"))
#' }
replace_precursor_symbols <- function(mrm_template, columns = c("Precursor Name", "Note")) {
  for (col in columns) {
    original_col <- paste0("original_", gsub(" ", "_", col))
    mrm_template[[original_col]] <- mrm_template[[col]]

    # Replace / and \ with underscores
    mrm_template[[col]] <- gsub("[/\\\\]", "_", mrm_template[[col]])
  }
  return(mrm_template)
}

# Functions for tests ----

check_dir_exists <- function(path)
  dir.exists(path)
create_dir <- function(path)
  dir.create(path, recursive = TRUE)

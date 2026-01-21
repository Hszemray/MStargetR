#qcCheckR internal functions ----

#' Import specific functions from packages
#' @name qcCheckR_import_external_functions
#' @keywords internal
#' @importFrom utils browseURL capture.output
#' @importFrom readr read_csv write_csv read_tsv
#' @importFrom dplyr vars mutate_all mutate bind_rows bind_cols filter select rename relocate arrange contains intersect pull left_join right_join any_of all_of across distinct rowwise c_across ungroup group_by case_when everything row_number summarise
#' @importFrom purrr map set_names
#' @importFrom plotly ggplotly
#' @importFrom ggplot2 ggplot aes theme element_text labs geom_vline geom_hline geom_point theme_bw scale_shape_manual scale_color_manual scale_size_manual guides guide_legend facet_wrap scale_fill_manual ylab geom_text
#' @importFrom tibble tibble add_column as_tibble is_tibble column_to_rownames rownames_to_column
#' @importFrom tidyr replace_na pivot_wider
#' @importFrom stringr str_extract str_detect
#' @importFrom magrittr %>%
#' @importFrom stats median setNames na.omit sd
#' @importFrom data.table :=
#' @importFrom tidyselect where
#' @importFrom openxlsx write.xlsx
#' @importFrom tools file_ext
#' @importFrom viridis viridis
#' @importFrom statTarget shiftCor
#' @importFrom ropls opls
NULL

#.----
#Setup Project Functions----

###Primary Function----
#' qcCheckR_setup_project
#'
#' This function sets up the project by initialising the master list, setting up project directories, and updating the script log.
#' @keywords internal
#' @param user_name Character string representing the user name for the project.
#' @param project_directory Directory path for the project folder containing the wiff folder and .wiff and .wiff.scan files for each plate.
#' @param mrm_template_list List of lists for mrm_guides.
#' @param QC_sample_label Key for filtering QC samples from sample list.
#' @param sample_tags Vector of character strings to pull sample types for names
#' @param mv_threshold threshold for missing value filter. default is 50%.
#' @return The updated `master_list` object with the project setup details.
#' @examples
#' \dontrun{
#' qcCheckR_setup_project(user_name, project_directory, QC_sample_label, sample_tags, mv_threshold)
#' }
qcCheckR_setup_project <- function(user_name,
                                   project_directory,
                                   mrm_template_list,
                                   QC_sample_label,
                                   sample_tags,
                                   mv_threshold) {
  validate_project_directory(project_directory)
  master_list <- initialise_master_list()
  master_list <- store_environment_details(master_list)
  master_list <- qcCheckR_set_project_details(
    master_list,
    user_name,
    project_directory,
    QC_sample_label,
    sample_tags,
    mv_threshold
  )
  master_list <- qcCheckR_read_mrm_guides(master_list, mrm_template_list)
  qcCheckR_setup_project_directories(master_list)
  master_list <- qcCheckR_import_PeakForgeR_reports(master_list)
  master_list <- find_method_version(master_list)
  master_list <- update_script_log(master_list,
                                   "project_setup",
                                   "start_time",
                                   "data_preparation")
  return(master_list)
}

###Sub Functions----
#' qcCheckR_set_project_details
#'
#' This function sets project details in the master list.
#' @keywords internal
#' @param master_list The master list object.
#' @param user_name Character string representing the user name for the project.
#' @param project_directory Directory path for the project folder.
#' @param QC_sample_label Key for filtering QC samples from sample list.
#' @param sample_tags Character vector of sample tags to filter sample types from file_names.
#' @param mv_threshold Numeric value for the missing value sample threshold.
#' @return The updated master list object with project details.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_set_project_details(master_list,
#'                                             user_name,
#'                                             project_directory,
#'                                             QC_sample_label,
#'                                             sample_tags,
#'                                             mv_threshold)
#' }
qcCheckR_set_project_details <- function(master_list,
                                         user_name,
                                         project_directory,
                                         QC_sample_label,
                                         sample_tags,
                                         mv_threshold) {
  master_list$project_details$project_dir <- project_directory
  master_list$project_details$user_name <- user_name
  master_list$project_details$project_name <- stringr::str_extract(master_list$project_details$project_dir, "[^/]*$")
  master_list$project_details$plateIDs <- c()
  master_list$project_details$qc_type <- QC_sample_label
  master_list$project_details$script_log$timestamps$start_time <- Sys.time()
  master_list$project_details$mv_sample_threshold <- mv_threshold
  #Set sample tags for ANPC
  if (is.null(master_list$project_details$sample_tags) &&
      master_list$project_details$user_name == "ANPC") {
    master_list$project_details$sample_tags <- c("pqc",
                                                 "qc",
                                                 "vltr",
                                                 "sltr",
                                                 "ltr",
                                                 "blank",
                                                 "istds",
                                                 "cond",
                                                 "sample")
  } else {
    master_list$project_details$sample_tags <- sample_tags
  }
  return(master_list)
}


#' qcCheckR_read_mrm_guides
#'
#' This function reads MRM guides from user-supplied paths in the mrm_template_list.
#' @keywords internal
#' @param master_list The master list object.
#' @param mrm_template_list List of MRM guide file paths for tsv.
#' @return The updated master list object with MRM guides.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_read_mrm_guides(master_list, mrm_template_list)
#' }
qcCheckR_read_mrm_guides <- function(master_list, mrm_template_list) {
  if (master_list$project_details$user_name == "ANPC" &&
      is.null(mrm_template_list)) {

    read_and_clean_sil <- function(path) {
      guide <- readr::read_tsv(path, show_col_types = FALSE)
      replace_precursor_symbols(guide, columns = c("Precursor Name", "Note"))
    }

    read_and_clean_conc <- function(path) {
      guide <- readr::read_tsv(path, show_col_types = FALSE)
      replace_precursor_symbols(guide, columns = c("SIL_name"))
    }

    master_list$templates$mrm_guides <- list(
      v1 = list(
        SIL_guide  = read_and_clean_sil(
          system.file("extdata", "LGW_lipid_mrm_template_v1.tsv", package = "MStargetR")
        ),
        conc_guide = read_and_clean_conc(
          system.file("extdata", "LGW_SIL_batch_103.tsv", package = "MStargetR")
        )
      ),
      v2 = list(
        SIL_guide  = read_and_clean_sil(
          system.file("extdata", "LGW_lipid_mrm_template_v2.tsv", package = "MStargetR")
        ),
        conc_guide = read_and_clean_conc(
          system.file(
            "extdata",
            "LGW_SIL_batch_Ultimate_2023_03_06.tsv",
            package = "MStargetR"
          )
        )
      ),
      v3 = list(
        SIL_guide  = read_and_clean_sil(
          system.file("extdata", "LGW_lipid_mrm_template_v3.tsv", package = "MStargetR")
        ),
        conc_guide = read_and_clean_conc(
          system.file(
            "extdata",
            "LGW_SIL_batch_Ultimate_2023_03_06.tsv",
            package = "MStargetR"
          )
        )
      ),
      v4 = list(
        SIL_guide  = read_and_clean_sil(
          system.file("extdata", "LGW_lipid_mrm_template_v4.tsv", package = "MStargetR")
        ),
        conc_guide = read_and_clean_conc(
          system.file("extdata", "v4_ISTD_conc_updated.tsv", package = "MStargetR")
        )
      )
    )

  } else {
    for (version in names(mrm_template_list)) {
      guide_paths <- mrm_template_list[[version]]

      read_guide <- function(path, columns) {
        ext <- tools::file_ext(path)
        guide <- if (ext == "tsv") {
          readr::read_tsv(path, show_col_types = FALSE)
        } else if (ext == "csv") {
          readr::read_csv(path, show_col_types = FALSE)
        } else {
          stop(paste("Unsupported file type:", ext))
        }
        replace_precursor_symbols(guide, columns = columns)
      }

      master_list$templates$mrm_guides[[version]]$SIL_guide  <- read_guide(guide_paths$SIL_guide,
                                                                           columns = c("Precursor Name", "Note"))
      master_list$templates$mrm_guides[[version]]$conc_guide <- read_guide(guide_paths$conc_guide, columns = c("SIL_name"))
    }
  }

  validate_qcCheckR_mrm_template_list(master_list)

  return(master_list)
}

#' qcCheckR_setup_project_directories
#'
#' This function sets up the project directories for the master list.
#' @keywords internal
#' @param master_list The master list object.
#' @return The updated master list object with project directories set up.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_setup_project_directories(master_list)
#' }
qcCheckR_setup_project_directories <- function(master_list) {
  #set qcCheckR write path
  qcCheckR_path <- file.path(master_list$project_details$project_dir, "all")
  # Create project directories if they do not exist
  if (!dir.exists(qcCheckR_path)) {
    dir.create(qcCheckR_path, recursive = TRUE, showWarnings = FALSE)
  }

  subdirs <- c("data", "html_report", "xlsx_report")
  for (subdir in subdirs) {
    dir_path <- file.path(qcCheckR_path, subdir)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, showWarnings = FALSE)
    }
  }
  if (!dir.exists(qcCheckR_path)) {
    stop("Failed to create qcCheckR directory at: ", qcCheckR_path)
  } else {
    message("qcCheckR directory set up at: ", qcCheckR_path)
  }
}

#' qcCheckR_import_PeakForgeR_reports
#'
#' This function imports PeakForgeR reports from the project directory and stores them in the master list.
#' @keywords internal
#' @param master_list The master list object.
#' @return The updated master list object with PeakForgeR reports imported.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_import_PeakForgeR_reports(master_list)
#' }
qcCheckR_import_PeakForgeR_reports <- function(master_list) {
  PeakForgeR_report_files <- list.files(
    master_list$project_details$project_dir,
    pattern = "_PeakForgeR_.*\\.(csv|tsv)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = FALSE,
    include.dirs = TRUE,
    no.. = FALSE
  )

  if (length(PeakForgeR_report_files) == 0) {
    stop(
      "No report files found in the specified project directory. Please ensure reports are in .csv or .tsv format and named with '_PeakForgeR_' pattern."
    )
  }

  for (file in PeakForgeR_report_files) {

    if (.Platform$OS.type == "windows" && nchar(file) > 260) {
      path <- dirname(file)
      short_link <- "C:/PeakForgeR"
      if (dir.exists(short_link)) {
        unlink(short_link, recursive = TRUE, force = TRUE)
      }
      cmd <- paste("cmd", "/c", "mklink", "/J", shQuote(short_link), shQuote(path))
      system(cmd, intern = FALSE)

      file_to_read <- file.path(short_link, basename(file))
    } else {
      file_to_read <- file
    }

    file_ext <- tools::file_ext(file_to_read)
    if (file_ext == "csv") {
      PeakForgeR_report <- readr::read_csv(file_to_read, show_col_types = FALSE)
    } else if (file_ext == "tsv") {
      PeakForgeR_report <- readr::read_tsv(file_to_read, show_col_types = FALSE)
    } else {
      warning(paste("Unsupported file type:", file))
      next
    }

    if (exists("short_link") && dir.exists(short_link)) {
      unlink(short_link, recursive = TRUE, force = TRUE)
    }

    file_name <- basename(file_to_read) %>%
      sub("\\.(csv|tsv)$", "", .) %>%
      sub("_PeakForgeR_", "_", .)

    colnames(PeakForgeR_report) <- gsub("\\.", "", colnames(PeakForgeR_report))
    colnames(PeakForgeR_report) <- gsub(" ", "", colnames(PeakForgeR_report))

    PeakForgeR_report <- PeakForgeR_report[!is.na(PeakForgeR_report$FileName), ]
    PeakForgeR_report <- PeakForgeR_report[PeakForgeR_report$Area != 0 &
                                             PeakForgeR_report$Height != 0, ]
    PeakForgeR_report <- PeakForgeR_report[!apply(is.na(PeakForgeR_report), 1, all), ]
    PeakForgeR_report <- PeakForgeR_report[
      !grepl("(?i)\\bCOND\\b|\\bBLANK\\b|\\bISTDs\\b", PeakForgeR_report$FileName, perl = TRUE),
    ]
    master_list$data$PeakForgeRReport[[file_name]] <- PeakForgeR_report
    master_list$project_details$plateIDs <- c(master_list$project_details$plateIDs, file_name)
  }

  return(master_list)
}

#' Find method versions
#'
#' This function finds the method versions for each plate in the master list.
#' It matches internal standards to the SIL guide and stores the matching SIL_guide.
#' @keywords internal
#' @param master_list A list containing porject details and data.
#' @return The updated `master_list` with method versions found for each plate.
find_method_version <- function(master_list) {
  for (plate_id in names(master_list$data$PeakForgeRReport)) {
    #Gather unique SIL molecule names from the PeakForgeR report
    report_SILS <-  master_list$data$PeakForgeRReport[[plate_id]] %>%
      dplyr::select(dplyr::contains("MoleculeName")) %>%
      dplyr::filter(grepl("SIL", MoleculeName, ignore.case = TRUE)) %>%
      unique()

    #Check for SILs in the SIL guides
    for (version in names(master_list$templates$mrm_guides)) {
      sil_guide <- master_list$templates$mrm_guides[[version]]$SIL_guide
      sil_targets <- sil_guide %>%
        dplyr::select(dplyr::contains("Precursor Name")) %>%
        dplyr::filter(grepl("SIL", `Precursor Name`, ignore.case = TRUE)) %>%
        unique()

      #Check if all SILs in the report are present in the SIL guide
      if (all(report_SILS$MoleculeName %in% sil_targets$`Precursor Name`)) {
        master_list$project_details$plate_method_versions[[plate_id]] <- version
        #message found version
        message("Method version found for plate ", plate_id, ": ", version)
        break
      }
    }
  }
  # If no version is found, stop script
  if (length(master_list$project_details$plate_method_versions) == 0) {
    stop(
      "No method version found for the plates in the master list. Please check the SIL guides and report data."
    )
  } else {
    #report unfound version
    if (length(master_list$project_details$plate_method_versions) < length(master_list$project_details$plateIDs)) {
      missing_plates <- setdiff(
        master_list$project_details$plateIDs,
        names(master_list$project_details$plate_method_versions)
      )
      message(
        "No method version found for the following plates: ",
        paste(missing_plates, collapse = ", "),
        ". \n Please check the SIL guides and report data."
      )
    }
  }
  return(master_list)
}


#.----

#Phase 1: Data preparation ----

##Transpose PeakForgeR Report Data Functions ----
###Primary Function ----

#' Transpose PeakForgeR Report Data
#'
#' This function transposes PeakForgeR report data for each plate in the master list.
#' It reshapes the data, cleans sample names, converts values to numeric, and stores the result.
#' @param master_list A list containing project details and data.
#' @keywords internal
#' @return The updated `master_list` object with transposed peak area data.
#'
#' @examples
#' \dontrun{
#' qcCheckR_transpose_data(master_list)
#' }
qcCheckR_transpose_data <- function(master_list) {
  validate_master_list(master_list)
  master_list$data$peakArea$transposed <- list()

  for (plate_id in master_list$project_details$plateIDs) {
    message("\nTransposing plate: ", plate_id)
    matching_name <- find_matching_report(master_list$data$PeakForgeRReport, plate_id)

    if (length(matching_name) == 1) {
      tryCatch({
        transposed <- transpose_plate_data(master_list$data$PeakForgeRReport[[matching_name]])
        master_list$data$peakArea$transposed[[plate_id]] <- transposed
        message("\nSuccessfully transposed plate: ", plate_id)
      }, error = function(e) {
        message("\nError transposing plate: ", plate_id)
        message("\nError message: ", e$message)
      })
    } else if (length(matching_name) > 1) {
      message("\nMultiple matching PeakForgeR reports found for plate: ",
              plate_id)
    } else {
      message("\nNo matching PeakForgeR report found for plate: ", plate_id)
    }
  }

  return(master_list)
}

###Sub Functions ----
#' Validate Master List Structure
#'
#' This function checks if the master list has the required structure and data.
#' @keywords internal
#' @param master_list The master list object to validate.
#' @return NULL if the structure is valid, otherwise throws an error.
validate_master_list <- function(master_list) {
  if (!is.list(master_list) ||
      !is.list(master_list$data$PeakForgeRReport)) {
    stop("Invalid master_list format.")
  }
}

#' Find Matching PeakForgeR Report Name
#'
#' This function finds the matching PeakForgeR report name based on the plate ID.
#' @keywords internal
#' @param report_list A list of PeakForgeR reports.
#' @param plate_id The ID of the plate to match.
#' @return The name of the matching report, or NULL if not found.
find_matching_report <- function(report_list, plate_id) {
  names(report_list)[grepl(plate_id, names(report_list))]
}

#' Transpose and Clean Plate Data
#'
#' This function transposes the plate data from the PeakForgeR report, reshaping it into a wide format.
#' It cleans the sample names by removing the file extension and converts area values to numeric.
#' @keywords internal
#' @param data The PeakForgeR report data frame for a specific plate.
#' @return A transposed tibble with sample names as rows and molecule names as columns.
transpose_plate_data <- function(data) {
  transposed <- data %>%
    tidyr::pivot_wider(
      id_cols = FileName,
      names_from = MoleculeName,
      values_from = Area,
      names_glue = "{MoleculeName}"
    ) %>%
    dplyr::rename(sample_name = FileName)

  transposed$sample_name <- sub("\\.mzML$", "", transposed$sample_name)
  transposed[, -1] <- sapply(transposed[, -1], as.numeric)
  tibble::as_tibble(transposed)
  transposed <- transposed %>%
    replace(is.na(.), 0) %>%
    dplyr::select(tidyselect::where(~ !all(. == 0)))



}

#.----
##Sort and QC Check Data Functions ----
### Primary Function ----
#' Sort and QC Check Data
#'
#' This function sorts the transposed peak area data by run order and performs QC checks.
#' It assigns sample types, validates QC coverage, and sets the appropriate QC type for the project.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with sorted data and QC check results.
#'
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_sort_data(master_list)
#' }
qcCheckR_sort_data <- function(master_list) {
  master_list$project_details$run_orders <- list()
  master_list$data$peakArea$sorted <- list()

  for (plate_id in names(master_list$data$peakArea$transposed)) {
    run_order <- extract_run_order(master_list$data$PeakForgeRReport[[plate_id]], plate_id)
    run_order <- assign_sample_type(master_list$project_details$sample_tags, run_order)
    validate_qc_types(run_order, master_list$project_details$sample_tags)
    sorted_data <- sort_and_filter_data(
      run_order,
      master_list$data$peakArea$transposed[[plate_id]],
      master_list$project_details$sample_tags
    )
    sorted_data$sample_ID <- extract_sample_id(sorted_data$sample_name)
    sorted_data <- sorted_data %>% dplyr::relocate(sample_ID, .after = sample_name)
    master_list$project_details$run_orders[[plate_id]] <- run_order
    master_list$data$peakArea$sorted[[plate_id]] <- sorted_data
    rm(run_order, sorted_data)
  }

  master_list <- assess_qc_coverage(master_list)
  master_list <- set_project_qc_type(master_list)
  master_list <- finalise_sorted_data(master_list)

  return(master_list)
}

###Sub Functions ----
#' Extract Run Order from PeakForgeR Report
#'
#' This function extracts the run order data from a PeakForgeR report, filtering by plate ID.
#' @keywords internal
#' @param report The PeakForgeR report data frame.
#' @param plate_id The ID of the plate to filter by.
#' @return A data frame containing the sample names, timestamps, and other relevant information.
extract_run_order <- function(report, plate_id) {
  extracted_data <- report %>%
    dplyr::select(dplyr::contains(c("FileName", "AcquiredTime"))) %>%
    dplyr::distinct(FileName, .keep_all = TRUE) %>%
    dplyr::mutate(FileName = sub(".mzML", "", FileName)) %>%
    dplyr::rename(sample_name = FileName, sample_timestamp = AcquiredTime) %>%
    dplyr::mutate(sample_timestamp = as.POSIXct(
      sample_timestamp,
      tryFormats = c(
        "%d/%m/%Y %H:%M",               # e.g. 9/02/2022 16:00
        "%m/%d/%Y %H:%M:%S",            # e.g. 09/27/2024 10:41:28
        # ISO 8601 formats
        "%Y-%m-%dT%H:%M:%SZ",           # e.g. 2021-03-13T18:12:31Z (UTC)
        "%Y-%m-%dT%H:%M:%S",            # e.g. 2021-03-13T18:12:31
        "%Y-%m-%dT%H:%M:%S%z",          # e.g. 2021-03-13T18:12:31+0800
        # Standard formats
        "%Y-%m-%d %H:%M:%S",            # e.g. 2021-03-13 18:12:31
        "%Y/%m/%d %H:%M:%S",            # e.g. 2021/03/13 18:12:31
        "%d/%m/%Y %H:%M:%S",            # e.g. 13/03/2021 18:12:31
        "%d-%m-%Y %H:%M:%S",            # e.g. 13-03-2021 18:12:31
        # AM/PM formats
        "%d/%m/%Y %I:%M %p",            # e.g. 13/03/2021 6:12 PM
        "%Y-%m-%d %I:%M %p",            # e.g. 2021-03-13 6:12 PM
        # Long formats
        "%B %d, %Y %H:%M",              # e.g. March 13, 2021 18:12
        "%b %d, %Y %I:%M %p",           # e.g. Mar 13, 2021 6:12 PM
        # Date-only formats
        "%Y-%m-%d",                     # e.g. 2021-03-13
        "%d/%m/%Y",                     # e.g. 13/03/2021
        "%B %d, %Y"                     # e.g. March 13, 2021
      ))) %>%

    # Now filter and arrange safely
    dplyr::filter(!is.na(sample_timestamp)) %>%
    dplyr::arrange(sample_timestamp) %>%

    dplyr::mutate(
      sample_plate_id = plate_id,
      sample_plate_order = dplyr::row_number(),
      sample_matrix = dplyr::case_when(
        stringr::str_detect(sample_name, "(?i)_SER_") ~ "SER",
        stringr::str_detect(sample_name, "(?i)_PLA_") ~ "PLA",
        stringr::str_detect(sample_name, "(?i)_URI_") ~ "URI",
        TRUE ~ NA_character_
      )
    )

  return(extracted_data)
}


#' Assign Sample Type
#'
#' This function assigns sample types based on the sample names and user-defined QC types.
#' @keywords internal
#' @param sample_tags A character vector passed from master_list$project_details$sample_tags, containing QC types.
#' @param run_order A data frame containing the run order information.
assign_sample_type <- function(sample_tags, run_order) {
  if (is.null(sample_tags) || length(sample_tags) < 1) {
    stop("sample_tags must be a character vector with at least one QC type.")
  }


  qc_case_expr <- paste(purrr::map_chr(
    sample_tags,
    ~ paste0("stringr::str_detect(sample_name, '(?i)", .x, "') ~ '", .x, "'")
  ), collapse = ",\n")
  qc_case_expr <- paste(qc_case_expr, "TRUE ~ 'sample'", sep = ",\n")

  data <- run_order %>%
    dplyr::mutate(sample_type = eval(parse(
      text = paste0("dplyr::case_when(", qc_case_expr, ")")
    )))

  return(data)
}

#' Validate QC Types
#'
#' This function validates the sample types in the run order against the provided sample tags.
#' @keywords internal
#' @param run_order A data frame containing the run order information.
#' @param sample_tags A character vector of QC types to validate against.
#' @return NULL if validation passes, otherwise throws an error.
validate_qc_types <- function(run_order, sample_tags) {
  invalid <- run_order %>%
    dplyr::filter(sample_type != "sample" & !sample_type %in% sample_tags)
  if (nrow(invalid) > 0) {
    print(invalid)
    stop("Invalid QC sample_type detected.")
  }
  if (all(run_order$sample_type == "sample")) {
    stop("No QC types were identified.")
  }
}

#' Sort and Filter Data
#'
#' This function sorts the run order data and filters it based on sample types.
#' @keywords internal
#' @param run_order A data frame containing the run order information.
#' @param transposed_data A data frame containing the transposed peak area data.
#' @param sample_tags A character vector of sample tags to filter by.
#' @return A data frame containing the sorted and filtered data.
sort_and_filter_data <- function(run_order, transposed_data, sample_tags) {
  data <- run_order %>%
    dplyr::left_join(transposed_data, by = "sample_name") %>%
    dplyr::arrange(sample_timestamp) %>%
    dplyr::mutate(sample_run_index = dplyr::row_number())

  return(data)
}

#' extract_sample_id
#'
#' This function pulls sampleID from file name.
#' @keywords internal
#' @param filenames sample_name column of .data
#' @return sample_ID column with unique sample identifiers
extract_sample_id <- function(filenames) {
  if (length(filenames) == 1)
    return(filenames)
  # Tokenise using _, -, .
  token_lists <- strsplit(filenames, "[-_.]")

  # Flatten all tokens to count frequency
  all_tokens <- unlist(token_lists)
  token_freq <- table(all_tokens)

  # Identify tokens that appear in all filenames
  common_tokens <- names(token_freq[token_freq >= length(filenames)])

  # Remove common tokens from each filename
  sample_ids <- lapply(token_lists, function(tokens) {
    unique_tokens <- tokens[!tokens %in% common_tokens]
    paste(unique_tokens, collapse = "_")
  })

  return(unlist(sample_ids))
}

#' Assess QC Coverage
#'
#' This function assesses the QC coverage for each plate in the master list.
#' It checks the ratio of QC samples to total samples and determines if the QC passed or failed.
#' @keywords internal
#' @param master_list A list containing project details and sorted peak area data.
#' @return The updated master list with QC coverage results.
assess_qc_coverage <- function(master_list) {
  master_list$project_details$qc_passed <- list()
  master_list$project_details$global_qc_pass <- list()

  for (plate_id in names(master_list$data$peakArea$sorted)) {
    plate_data <- master_list$data$peakArea$sorted[[plate_id]]
    qc_types <- plate_data %>%
      dplyr::filter(sample_type != "sample") %>%
      dplyr::distinct(sample_type) %>%
      dplyr::pull()

    for (qc in qc_types) {
      total <- nrow(plate_data)
      count <- sum(tolower(plate_data$sample_type) == tolower(qc))
      ratio <- count / total
      status <- if (ratio < 8 / 120 || count < 2)
        "fail"
      else
        "pass"
      master_list$project_details$qc_passed[[plate_id]][[qc]] <- status
      master_list$project_details$global_qc_pass[[qc]] <- "pass"
    }
  }

  for (plate_id in names(master_list$data$peakArea$sorted)) {
    for (qc in names(master_list$project_details$qc_passed[[plate_id]])) {
      if (master_list$project_details$qc_passed[[plate_id]][[qc]] == "fail") {
        master_list$project_details$global_qc_pass[[qc]] <- "fail"
      }
    }
  }

  return(master_list)
}

#' Set Project QC Type
#'
#' This function sets the QC type for the project based on user-specified QC types and global QC pass status.
#' @keywords internal
#' @param master_list A list containing project details and QC pass status.
#' @return The updated master list with the project QC type set.
set_project_qc_type <- function(master_list) {
  user_qc <- tolower(master_list$project_details$qc_type)
  global_pass <- master_list$project_details$global_qc_pass

  if (!is.null(global_pass[[user_qc]]) &&
      global_pass[[user_qc]] == "pass") {
    message("qcCheckeR has set QC to user-specified type: ", user_qc)
    master_list$project_details$qc_type <- user_qc
  } else {
    warning("User-specified QC type ", user_qc, " did not pass QC checks.")
    passed_qcs <- names(global_pass)[global_pass == "pass"]
    if (length(passed_qcs) > 0) {
      qc_counts <- sapply(passed_qcs, function(qc) {
        sum(sapply(master_list$project_details$qc_passed, function(plate) {
          if (!is.null(plate[[qc]]) && plate[[qc]] == "pass")
            1
          else
            0
        }))
      })
      best_qc <- names(which.max(qc_counts))
      master_list$project_details$qc_type <- best_qc
      message("qcCheckeR has set QC to alternative type: ", best_qc)
    } else {
      stop(
        "No QC types passed. Stopping script.
            \n Please check sample_tags and filenames are correct."
      )
    }
  }

  return(master_list)
}

#' Finalise Sorted Data
#'
#' This function finalises the sorted data by adding sample type factors, reversing the order of sample types, and setting the sample data source.
#' @keywords internal
#' @param master_list A list containing project details and sorted peak area data.
#' @return The updated master list with finalised sorted data.
finalise_sorted_data <- function(master_list) {
  for (plate_id in names(master_list$data$peakArea$sorted)) {
    sorted <- master_list$data$peakArea$sorted[[plate_id]]
    qc_type <- master_list$project_details$qc_type

    sorted <- sorted %>%
      dplyr::mutate(
        sample_type_factor = factor(sample_type, levels = unique(
          c("sample", qc_type, "ltr", "pqc", "vltr", "sltr")
        ), ordered = TRUE),
        sample_type_factor_rev = factor(
          sample_type_factor,
          levels = rev(levels(sample_type_factor)),
          ordered = TRUE
        ),
        sample_type = ifelse(tolower(sample_type) == tolower(qc_type), "qc", "sample"),
        sample_data_source = ".peakArea"
      ) %>%
      dplyr::relocate(sample_type_factor,
               sample_type_factor_rev,
               sample_data_source,
               .after = sample_type)


    master_list$data$peakArea$sorted[[plate_id]] <- sorted
  }

  all_samples <- dplyr::bind_rows(master_list$data$peakArea$sorted) %>%
    dplyr::arrange(sample_timestamp) %>%
    dplyr::mutate(sample_run_index = dplyr::row_number()) %>%
    dplyr::relocate(sample_run_index, .before = sample_name)

  master_list$data$peakArea$sorted <- split(all_samples, all_samples$sample_plate_id)


  master_list$data$peakArea$sorted <- lapply(
    master_list$data$peakArea$sorted,
    function(df) {
      df %>% dplyr::select(
        which(
          !sapply(df, function(col) all(is.na(col))) | grepl("sample", names(df))
        )
      )
    }
  )



  return(master_list)
}

#.----
##Impute Missing Data ----
###Primary Function ----
#' Impute Missing Data
#'
#' This function imputes missing and zero values in the `master_list` data using the minimum intensity of each feature in the batch divided by 2.
#' It handles infinite and NaN values, applies imputation, and merges metadata back into the result.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with imputed peak area data.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_impute_data(master_list)
#' }
qcCheckR_impute_data <- function(master_list) {
  master_list$data$peakArea$imputed <- list()

  for (plate_id in names(master_list$data$peakArea$sorted)) {
    sorted_data <- master_list$data$peakArea$sorted[[plate_id]]
    impute_matrix <- prepare_imputation_matrix(sorted_data)
    imputed <- apply_lgw_imputation(impute_matrix)
    imputed <- merge_metadata(sorted_data, imputed)
    master_list$data$peakArea$imputed[[plate_id]] <- imputed
  }

  return(master_list)
}

###Sub Functions ----
#' Prepare Imputation Matrix
#'
#' This function prepares the imputation matrix from the sorted peak area data.
#' It removes sample name columns, converts the data to a matrix, and replaces problematic values (zeros, infinite, NaN).
#' @keywords internal
#' @param data A tibble containing the sorted peak area data.
#' @return A tibble ready for imputation, with problematic values replaced.
prepare_imputation_matrix <- function(data) {
  rownames(data) <- NULL
  prepped_data <- data %>%
    tibble::column_to_rownames("sample_name") %>%
    dplyr::select(-dplyr::contains("sample")) %>%
    replace_problematic_values()
  return(prepped_data)
}

#' Replace Problematic Values
#'
#' This function replaces zeros, infinite, and NaN values in a matrix with NA.
#' @keywords internal
#' @param mat A numeric matrix.
#' @return A matrix with zeros replaced by NA, and infinite and NaN values also replaced by NA.
replace_problematic_values <- function(mat) {
  mat <- as.matrix(mat)
  mat[mat == 0] <- NA
  mat[is.infinite(mat)] <- NA
  mat[is.nan(mat)] <- NA
  return(mat)
}

#' LGW Impute
#'
#' This function performs imputation on a tibble by replacing zero values with half the minimum non-zero value in each column.
#' It calculates the minimum non-zero value for each column, divides it by 2, and replaces zero values with this calculated value.
#' @keywords internal
#' @param x A tibble containing numeric data.
#' @return A tibble with zero values replaced by half the minimum non-zero value in each column.
lgw_impute <- function(x) {
  # Validate parameter
  if (!tibble::is_tibble(x)) {
    tibble::as_tibble(x)
  }

  if (!all(sapply(x, is.numeric))) {
    x <- x %>% dplyr::mutate(dplyr::across(
      tidyselect::where(is.character) &
        !matches("sample_name"),
      ~ suppressWarnings(as.numeric(.))
    ))

  }

  # Calculate min/2 values for each column
  min_half_values <- purrr::map(.x = x, .f = ~ {
    non_zero_values <- .x[.x > 0]
    if (length(non_zero_values) == 0 ||
        all(is.na(non_zero_values))) {
      return(NA)
    }

    min(non_zero_values, na.rm = TRUE) / 2
  })

  # Replace NAs with min/2 values
  x %>%
    dplyr::mutate_all(~ replace(., . == 0, NA)) %>%
    tidyr::replace_na(replace = min_half_values)
}

#' Apply LGW Imputation
#'
#' A wrapper to apply lgw_impute to the dataframe
#' @keywords internal
#' @param mat dataframe or matrix for imputation to be applied
apply_lgw_imputation <- function(mat) {
  mat %>%
    as.data.frame() %>%
    lgw_impute() %>%
    tibble::rownames_to_column("sample_name") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ ifelse(is.infinite(.), 1, .)))
}

#'Merge_Metadata
#'
#' This function merges metadata from the original data with the imputed data.
#' It selects sample name columns from the original data and performs a left join with the imputed data.
#' It also adds a column indicating the data source of the imputed data.
#' @keywords internal
#' @param original_data A tibble containing the original peak area data with sample names.
#' @param imputed_data A tibble containing the imputed peak area data with sample names.
#' @return A tibble with merged metadata and imputed data, including a sample data source column.
merge_metadata <- function(original_data, imputed_data) {
  metadata <- dplyr::select(original_data, dplyr::contains("sample"))
  merged <- dplyr::left_join(metadata, imputed_data, by = "sample_name")
  merged$sample_data_source <- ".peakAreaImputed"
  return(merged)
}

#.----
##Calculate Response and Concentration ----
###Primary Function ----
#' Calculate Response and Concentration
#'
#' This function calculates the response ratio and concentration for each sample in the `master_list` data.
#' It uses SIL internal standards and template guides to compute values for both sorted and imputed data.
#' @keywords internal
#' @param master_list A list containing project details, peak area data, and SIL templates.
#' @return The updated `master_list` with calculated response and concentration data.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_calculate_response_concentration(master_list)
#' }
qcCheckR_calculate_response_concentration <- function(master_list) {
  master_list$data$response <- list()
  master_list$data$concentration <- list()

  batches <- names(master_list$data$peakArea$imputed)

  for (plate_id in batches) {
    template_version <- master_list$project_details$plate_method_versions[[plate_id]]
    master_list$project_details$is_ver <- template_version
    master_list$templates[["Plate SIL version"]][[plate_id]] <- template_version

    for (data_type in c("sorted", "imputed")) {
      master_list <- calculate_plate_response_concentration(master_list, plate_id, data_type, template_version)
    }
  }

  master_list <- harmonise_lipid_columns(master_list)
  return(master_list)
}

###Sub Functions ----

#' Calculate Response and Concentration for a Plate
#'
#' This function calculates the response and concentration for a specific plate in the master list.
#' It retrieves the peak area data, identifies SIL columns, and processes each SIL target to compute response ratios and concentrations.
#' @keywords internal
#' @param master_list A list containing project details, peak area data, and SIL templates.
#' @param plate_id The ID of the plate to process.
#' @param data_type The type of data to process (either "sorted" or "imputed").
#' @param template_version The version of the template to use for processing.
#' @return The updated master list with calculated response and concentration data for the specified plate.
calculate_plate_response_concentration <- function(master_list,
                                                   plate_id,
                                                   data_type,
                                                   template_version) {
  data <- master_list$data$peakArea[[data_type]][[plate_id]]
  sil_notes <- master_list$templates$mrm_guides[[template_version]]$SIL_guide$Note
  sil_cols <- dplyr::intersect(colnames(dplyr::select(data, dplyr::contains("SIL"))), sil_notes)

  # Reorder columns: non-SIL first, then SIL
  data <- dplyr::bind_cols(dplyr::select(data, -dplyr::contains("SIL")), dplyr::select(data, dplyr::all_of(sil_cols)))
  master_list$data$peakArea[[data_type]][[plate_id]] <- data

  # Initialise response and concentration lists
  master_list$data$response[[data_type]][[plate_id]] <- dplyr::select(data, dplyr::contains("sample"))
  master_list$data$concentration[[data_type]][[plate_id]] <- dplyr::select(data, dplyr::contains("sample"))

  for (sil in sil_cols) {
    master_list <- process_sil_target(master_list, plate_id, data_type, template_version, sil)
  }

  # Tag data sources
  master_list$data$response[[data_type]][[plate_id]]$sample_data_source <- paste0(".response.", data_type)
  master_list$data$concentration[[data_type]][[plate_id]]$sample_data_source <- paste0("concentration.", data_type)

  return(master_list)
}

#' Process SIL Target
#'
#' This function processes a specific SIL target by calculating the response ratio and concentration for each sample.
#' It retrieves the precursor names from the SIL guide, calculates the response ratio by dividing the peak area by the SIL value, and computes the concentration using the concentration factor from the template.
#' @keywords internal
#' @param master_list A list containing project details, peak area data, and SIL templates.
#' @param plate_id The ID of the plate to process.
#' @param data_type The type of data to process (either "sorted" or "imputed").
#' @param template_version The version of the template to use for processing.
#' @param sil The name of the SIL target to process.
#' @return The updated master list with calculated response and concentration data for the specified SIL target.
process_sil_target <- function(master_list,
                               plate_id,
                               data_type,
                               template_version,
                               sil) {
  sil_guide <- master_list$templates$mrm_guides[[template_version]]$SIL_guide
  precursors <- sil_guide %>% dplyr::filter(Note == sil) %>% dplyr::pull(`Precursor Name`)


  if (length(precursors) == 0)
    return(master_list)

  data <- master_list$data$peakArea[[data_type]][[plate_id]]
  sil_values <- data[[sil]]
  target_data <- dplyr::select(data, sample_name, dplyr::any_of(precursors)) %>%
                  tibble::column_to_rownames("sample_name")

  if (ncol(target_data) == 0)
    return(master_list)

  # Calculate response
  response <- as.matrix(target_data / sil_values)
  response[is.na(response) | is.infinite(response)] <- 0
  response_df <- as.data.frame(response) %>% tibble::rownames_to_column("sample_name")

  master_list$data$response[[data_type]][[plate_id]] <- dplyr::left_join(master_list$data$response[[data_type]][[plate_id]], response_df, by = "sample_name")

  # Calculate concentration
  conc_factor <- master_list$templates$mrm_guides[[template_version]]$conc_guide %>%
    dplyr::filter(SIL_name == sil) %>%
    dplyr::pull(concentration_factor)


  if (length(conc_factor) == 1) {
    concentration <- as.matrix(response * conc_factor)
    concentration[is.na(concentration) |
                    is.infinite(concentration)] <- 0
    concentration_df <- as.data.frame(concentration) %>% tibble::rownames_to_column("sample_name")

    master_list$data$concentration[[data_type]][[plate_id]] <- dplyr::left_join(master_list$data$concentration[[data_type]][[plate_id]], concentration_df, by = "sample_name")
  }

  return(master_list)
}

#' Harmonise Lipid Columns
#'
#' This function harmonises the lipid columns across response and concentration data types in the master list.
#' It ensures that all lipid subtypes have the same columns by selecting only the common lipids across all plates.
#' @keywords internal
#' @param master_list A list containing project details and response/concentration data.
#' @return The updated master list with harmonised lipid columns.
harmonise_lipid_columns <- function(master_list) {
  for (data_type in c("response", "concentration")) {
    for (subtype in names(master_list$data[[data_type]])) {
      common_lipids <- Reduce(dplyr::intersect, lapply(master_list$data[[data_type]][[subtype]], colnames))
      for (plate_id in names(master_list$data[[data_type]][[subtype]])) {
        master_list$data[[data_type]][[subtype]][[plate_id]] <- dplyr::select(master_list$data[[data_type]][[subtype]][[plate_id]], dplyr::all_of(common_lipids))
      }
    }
  }
  return(master_list)
}

#.----
##Batch Correction and Signal Drift Adjustment ----
### Primary Function ----
#' Batch Correction and Signal Drift Adjustment
#'
#' This function performs batch correction and signal drift adjustment on the concentration data in `master_list` using the `statTarget` package.
#' It prepares phenotype and profile files, runs `statTarget::shiftCor`, and integrates corrected data back into the master list.
#' @keywords internal
#' @param master_list A list containing project details, concentration data, and metadata.
#' @return The updated `master_list` with corrected concentration and peak area data.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_statTarget_batch_correction(master_list)
#' }
qcCheckR_statTarget_batch_correction <- function(master_list) {
  FUNC_list <- initialise_statTarget_environment(master_list)
  FUNC_list <- prepare_statTarget_files(FUNC_list)
  FUNC_list <- run_statTarget_shiftCor(FUNC_list, master_list)
  master_list <- integrate_corrected_data(master_list, FUNC_list)
  master_list <- update_script_log(master_list,
                                   "data_preparation",
                                   "project_setup",
                                   "data_filtering")
  return(master_list)
}

### Sub Functions ----
#' Initialise statTarget Environment
#'
#' This function initialises the environment for `statTarget` batch correction.
#' It creates necessary directories, sets up the master data, and flags failed QC injections.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return A list containing the project directory, master data, and metabolite list for `statTarget`.
initialise_statTarget_environment <- function(master_list) {
  dirs <- c("data",
            "data/batch_correction",
            "data/rda",
            "xlsx_report",
            "html_report")
  lapply(dirs, function(d) {
    dir_path <- file.path(master_list$project_details$project_dir, "all", d)
    if (!check_dir_exists(dir_path))
      create_dir(dir_path)
  })


  master_list$project_details$statTarget_qc_type <- master_list$project_details$qc_type
  FUNC_list <- list()
  FUNC_list$project_dir <- file.path(master_list$project_details$project_dir,
                                     "all",
                                     "data",
                                     "batch/correction")
  FUNC_list$master_data <- dplyr::bind_rows(master_list$data$concentration$imputed)

  FUNC_list$master_data$sample_type <- "sample"
  FUNC_list$master_data$sample_type[tolower(FUNC_list$master_data$sample_type_factor) == tolower(master_list$project_details$statTarget_qc_type)] <- "qc"

  FUNC_list$metabolite_list <- FUNC_list$master_data %>%
    dplyr::select(-dplyr::contains("sample")) %>%
    names()

  FUNC_list <- flag_failed_qc_injections(FUNC_list)
  return(FUNC_list)
}

#' Flag Failed QC Injections
#'
#' This function flags failed QC injections by checking the signal intensity of QC samples.
#' It identifies QC samples with low signal intensity and marks them as "sample" in the master data.
#' @keywords internal
#' @param FUNC_list A list containing the master data and metabolite list.
#' @return The updated `FUNC_list` with flagged failed QC injections.
flag_failed_qc_injections <- function(FUNC_list) {
  qc_fail <- lapply(unique(FUNC_list$master_data$sample_plate_id), function(batch) {
    qc_data <- FUNC_list$master_data %>%
      dplyr::filter(sample_type == "qc", sample_plate_id == batch) %>%
      dplyr::select(-dplyr::contains("sample"), -dplyr::contains("SIL"))

    low_signal <- rowSums(qc_data) < stats::median(rowSums(qc_data)) * 0.1

    # Get the sample names of the low signal QC samples
    failed_samples <- FUNC_list$master_data %>%
      dplyr::filter(sample_type == "qc", sample_plate_id == batch) %>%
      dplyr::pull(sample_name)

    failed_samples[low_signal]
  })

  qc_fail <- unlist(qc_fail)

  FUNC_list$master_data$sample_type[FUNC_list$master_data$sample_name %in% qc_fail] <- "sample"

  return(FUNC_list)
}


#' Prepare statTarget Files
#'
#' This function prepares the phenotype and profile files required for `statTarget` batch correction.
#' It creates a phenotype file with sample metadata and a profile file with metabolite data.
#' @keywords internal
#' @param FUNC_list A list containing the project directory, master data, and metabolite list.
#' @return The updated `FUNC_list` with created phenotype and profile files.
prepare_statTarget_files <- function(FUNC_list) {
  FUNC_list <- create_pheno_file(FUNC_list)
  FUNC_list <- create_profile_file(FUNC_list)
  return(FUNC_list)
}

#' Create Pheno File
#'
#' This function creates a phenotype file from the master data.
#' It selects relevant columns, renames them, and formats the sample IDs and classes.
#' @keywords internal
#' @param FUNC_list A list containing the master data and project directory.
#' @return The updated `FUNC_list` with the created phenotype file.
create_pheno_file <- function(FUNC_list) {
  pheno <- FUNC_list$master_data %>%
    dplyr::select(sample_name, sample_plate_id, sample_type, sample_run_index) %>%
    dplyr::rename(batch = sample_plate_id,
           class = sample_type,
           order = sample_run_index) %>%
    dplyr::arrange(order)

  # Initialise QC-ordered template
  FUNC_list$PhenoFile$template_qc_order <- NULL
  qc_idx <- NULL

  # Get unique batch IDs
  batch_ids <- unique(pheno$batch)

  # Ensure QC is first and last in each batch
  for (batch_id in batch_ids) {
    batch_data <- pheno %>% dplyr::filter(batch == batch_id)
    qc_positions <- which(batch_data$class == "qc")

    # Move first QC to position 1 if needed
    if (length(qc_positions) > 0 && qc_positions[1] > 1) {
      first_qc <- batch_data[qc_positions[1], ]
      batch_data <- batch_data[-qc_positions[1], ]
      batch_data <- dplyr::bind_rows(first_qc, batch_data)
    }

    # Move last QC to final position if needed
    qc_positions <- which(batch_data$class == "qc")
    if (length(qc_positions) > 0 &&
        qc_positions[length(qc_positions)] < nrow(batch_data)) {
      last_qc <- batch_data[qc_positions[length(qc_positions)], ]
      batch_data <- batch_data[-qc_positions[length(qc_positions)], ]
      batch_data <- dplyr::bind_rows(batch_data, last_qc)
    }

    # Store updated batch
    qc_idx <- c(qc_idx, qc_positions)
    FUNC_list$PhenoFile$template_qc_order <- dplyr::bind_rows(FUNC_list$PhenoFile$template_qc_order, batch_data)
  }

  # Assign sample names
  FUNC_list$PhenoFile$template_qc_order <- FUNC_list$PhenoFile$template_qc_order %>%
    dplyr::group_by(class) %>%
    dplyr::mutate(sample = dplyr::case_when(
      class == "qc" ~ paste0("QC", dplyr::row_number()),
      class == "sample" ~ paste0("sample", dplyr::row_number())
    )) %>%
    dplyr::ungroup() %>%
    dplyr::relocate(sample, .after = sample_name)

  # Final formatting
  FUNC_list$PhenoFile$template_sample_id <- FUNC_list$PhenoFile$template_qc_order
  FUNC_list$PhenoFile$template_sample_id$class[FUNC_list$PhenoFile$template_sample_id$class == "qc"] <- NA
  FUNC_list$PhenoFile$template_sample_id$order <- seq_len(nrow(FUNC_list$PhenoFile$template_sample_id))
  FUNC_list$PhenoFile$template_sample_id$batch <- as.numeric(factor(FUNC_list$PhenoFile$template_sample_id$batch))

  #Final output file for statTarget
  Output <-  FUNC_list$PhenoFile$template_sample_id %>%
    dplyr::select(sample, batch, class, order)

  # Write to CSV
  output_dir <- file.path(FUNC_list$project_dir,
                          paste0(Sys.Date(), "_signal_correction_results"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  readr::write_csv(Output, file.path(output_dir, "PhenoFile.csv"))

  return(FUNC_list)
}

#' Create Profile File
#'
#' This function creates a profile file from the master data.
#' It selects the sample names and metabolite data, renames columns, and formats the data into a matrix.
#' @keywords internal
#' @param FUNC_list A list containing the master data, metabolite list, and project directory.
#' @param FUNC_list A list containing the master data, metabolite list, and project directory.
#' @return The updated `FUNC_list` with the created profile file.
create_profile_file <- function(FUNC_list) {

  if (length(FUNC_list$metabolite_list) == 0) {
    stop("metabolite_list cannot be empty")
  }

  profile <- FUNC_list$master_data %>%
    dplyr::select(sample_name, dplyr::all_of(FUNC_list$metabolite_list))

  ordered <- FUNC_list$PhenoFile$template_sample_id %>%
    dplyr::select(sample, sample_name) %>%
    dplyr::left_join(profile, by = "sample_name") %>%
    dplyr::select(-sample_name)

  profile_matrix <- tibble::as_tibble(cbind(nms = names(ordered), t(ordered))) %>%
    stats::setNames(.[1, ]) %>%
    dplyr::rename(name = sample) %>%
    dplyr::filter(name != "sample") %>%
    dplyr::mutate(dplyr::across(-name, as.numeric))

  metabolite_map <- profile_matrix %>%
    dplyr::select(name) %>%
    dplyr::filter(!stringr::str_detect(name, "SIL")) %>%
    dplyr::mutate(metabolite_code = paste0("M", dplyr::row_number()))

  profile_matrix <- dplyr::left_join(metabolite_map, profile_matrix, by = "name") %>%
    dplyr::select(-name) %>%
    dplyr::rename(name = metabolite_code)

  FUNC_list$ProfileFile <- list(ProfileFile = profile_matrix, metabolite_list = metabolite_map)

  output_dir <- file.path(FUNC_list$project_dir,
                          paste0(Sys.Date(), "_signal_correction_results"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  readr::write_csv(profile_matrix, file.path(output_dir, "ProfileFile.csv"))
  return(FUNC_list)
}

#' Run statTarget shiftCor
#'
#' This function runs the `statTarget::shiftCor` function to perform signal drift correction on the prepared phenotype and profile files.
#' It specifies the correction parameters and reads the corrected data from the output file.
#' @keywords internal
#' @param FUNC_list A list containing the project directory, phenotype file, and profile file.
#' @param master_list A list containing all project details and data.
#' @return The updated `FUNC_list` with corrected data and adjusted QC means.
run_statTarget_shiftCor <- function(FUNC_list, master_list) {
  samPeno <- file.path(
    FUNC_list$project_dir,
    paste0(Sys.Date(), "_signal_correction_results/PhenoFile.csv")
  )
  samFile <- file.path(
    FUNC_list$project_dir,
    paste0(Sys.Date(), "_signal_correction_results/ProfileFile.csv")
  )

  cur_wd <- getwd()
  setwd(file.path(
    FUNC_list$project_dir,
    paste0(Sys.Date(), "_signal_correction_results")
  ))

  statTarget::shiftCor(
    samPeno = samPeno,
    samFile = samFile,
    Frule = 0,
    ntree = 500,
    MLmethod = 'QCRFSC',
    imputeM = "minHalf",
    plot = FALSE,
    coCV = 10000
  )

  setwd(cur_wd)
  rm(cur_wd)

  corrected <- readr::read_csv(file.path(
    FUNC_list$project_dir,
    paste0(
      Sys.Date(),
      "_signal_correction_results/statTarget/shiftCor/After_shiftCor/shift_all_cor.csv"
    )
  ), show_col_types = FALSE)

  corrected <- clean_statTarget_output(corrected)
  FUNC_list$corrected_data <- list(data = corrected)
  FUNC_list <- transpose_and_merge_corrected(FUNC_list)
  FUNC_list <- adjust_qc_means(FUNC_list, master_list)
  return(FUNC_list)
}

#' Clean statTarget Output
#'
#' This function cleans the output data from `statTarget` by filtering out unwanted rows and renaming columns.
#' It handles different column structures based on the presence of specific sample columns.
#' @keywords internal
#' @param data A tibble containing the output data from `statTarget`.
#' @return A cleaned tibble with renamed columns and numeric data types.
clean_statTarget_output <- function(data) {
  if ("sample1" %in% colnames(data)) {
    data <- data %>%
      dplyr::filter(sample != "class") %>%
      dplyr::rename(name = sample) %>%
      dplyr::mutate(dplyr::across(-name, as.numeric))
  } else if ("M1" %in% colnames(data)) {
    data <- data %>%
      t() %>%
      data.frame() %>%
      tibble::rownames_to_column() %>%
      stats::setNames(.[1, ]) %>%
      dplyr::filter(sample != "class" & sample != "sample") %>%
      dplyr::rename(name = sample) %>%
      dplyr::mutate(dplyr::across(-name, as.numeric))
  }
  return(data)
}


#' Transpose and Merge Corrected Data
#' This function transposes the corrected data and merges it with the metabolite list from the profile file.
#' It renames columns, filters out unwanted rows, and adds sample metadata.
#' @keywords internal
#' @param FUNC_list A list containing the corrected data and profile file.
#' @return The updated `FUNC_list` with transposed and merged corrected data.
transpose_and_merge_corrected <- function(FUNC_list) {
  transposed <- FUNC_list$ProfileFile$metabolite_list %>%
    dplyr::rename(lipid = metabolite_code) %>%
    dplyr::right_join(FUNC_list$corrected_data$data %>% dplyr::rename(lipid = name), by = "lipid") %>%
    dplyr::select(-lipid) %>%
    as.matrix() %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column() %>%
    stats::setNames(.[1, ]) %>%
    dplyr::filter(name != "name") %>%
    dplyr::mutate(dplyr::across(-name, as.numeric)) %>%
    dplyr::rename(sample = name) %>%
    dplyr::left_join(FUNC_list$PhenoFile$template_sample_id, by = "sample") %>%
    dplyr::left_join(FUNC_list$master_data %>% dplyr::select(dplyr::contains("sample")), by = "sample_name") %>%
    dplyr::select(-sample, -batch, -class, -order)

  #reorder columns to follow FUNC_list$master_data
  transposed <- transposed %>%
    dplyr::select(
      sample_run_index,
      sample_name,
      sample_timestamp,
      sample_plate_id,
      sample_plate_order,
      sample_matrix,
      sample_type,
      sample_type_factor,
      sample_type_factor_rev,
      sample_data_source,
      dplyr::everything()
    )

  FUNC_list$corrected_data$data_transposed <- transposed
  return(FUNC_list)
}

#' Adjust QC Means
#'
#' This function adjusts the means of QC samples in the corrected data by calculating the correction ratio based on original and corrected means.
#' It applies the correction ratio to the corrected data and updates the sample type for QC samples.
#' @keywords internal
#' @param FUNC_list A list containing the master data and corrected data.
#' @param master_list A list containing the project details and data.
#' @return The updated `FUNC_list` with adjusted QC means in the corrected data.
adjust_qc_means <- function(FUNC_list, master_list) {


  required_cols <- setdiff(names(FUNC_list$master_data), "sample_type")
  if (length(required_cols) == 0) {
    stop("No metabolite columns found in master_data.")
  }

  original_means <- FUNC_list$master_data %>%
    dplyr::filter(sample_type == "qc") %>%
    dplyr::select(-dplyr::contains("sample")) %>%
    colMeans()

  if (any(original_means == 0)) {
    warning("Original mean contains zero(s), which may lead to invalid correction ratios.")
  }

  corrected_means <- FUNC_list$corrected_data$data_transposed %>%
    dplyr::filter(sample_type == "qc") %>%
    dplyr::select(-dplyr::contains("sample")) %>%
    colMeans()

  qc_means <- tibble(
    metabolite = names(original_means),
    original_mean = original_means,
    corrected_mean = corrected_means,
    correction_ratio = corrected_mean / original_mean
  )

  adjusted <- FUNC_list$corrected_data$data_transposed
  for (met in qc_means$metabolite) {
    adjusted[[met]] <- adjusted[[met]] / qc_means$correction_ratio[qc_means$metabolite == met]
  }

  qc_type <- unique(dplyr::bind_rows(master_list$data$concentration$imputed)$sample_type_factor[dplyr::bind_rows(master_list$data$concentration$imputed)$sample_type == "qc"]) %>% as.character()

  adjusted$sample_type <- ifelse(adjusted$sample_type_factor == qc_type, "qc", "sample")
  FUNC_list$corrected_data$data_qc_mean_adjusted <- adjusted
  return(FUNC_list)
}


#' Integrate Corrected Data into Master List
#'
#' This function integrates the corrected data from `statTarget` into the master list.
#' It splits the corrected data by sample plate ID, updates the sample data source, and processes the peak area and concentration data.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param FUNC_list A list containing the corrected data from `statTarget`.
#' @return The updated `master_list` with integrated corrected data.
integrate_corrected_data <- function(master_list, FUNC_list) {
  if (!"sample_plate_id" %in% colnames(FUNC_list$corrected_data$data_qc_mean_adjusted)) {
    stop("Missing required column: sample_plate_id")
  }

  corrected <- FUNC_list$corrected_data$data_qc_mean_adjusted %>%
    dplyr::relocate(sample_ID, .after = sample_name)
  master_list$data$concentration$corrected <- split(corrected, corrected$sample_plate_id)

  for (batch in names(master_list$data$concentration$corrected)) {
    master_list$data$concentration$corrected[[batch]]$sample_data_source <- ".peakAreaCorrected"
  }

  master_list$data$peakArea$statTargetProcessed <- list()
  master_list$data$concentration$statTargetProcessed <- list()
  for (batch in unique(corrected$sample_plate_id)) {
    master_list$data$peakArea$statTargetProcessed[[batch]] <- corrected %>%
      dplyr::filter(sample_plate_id == batch) %>%
      dplyr::relocate(sample_ID, .after = sample_name)
    master_list$data$concentration$statTargetProcessed[[batch]] <- corrected %>%
      dplyr::filter(sample_plate_id == batch) %>%
      dplyr::relocate(sample_ID, .after = sample_name)
    master_list$data$peakArea$statTargetProcessed[[batch]]$sample_data_source <- "concentration.statTarget"
  }

  return(master_list)
}


#. ----
#Phase 2: Data Filtering ----
# mising value data filtering is performed on peakArea data from the imported user report.
# First the filtering is performed per sample to identify failed samples (e.g. injection, extraction, preparation errors, or if sample is missing from well) it will return a high % of missing values
# Once failed samples have been identified - the filtering then identifies lipids that have >50% missing values
# note: missing also refers to <limit of detection [<LOD]. This refers to instances of peak areas that are <5000 counts, as PeakForgeR will sometimes integrate noise giving a small value.

## Set QC Type for Filtering ----
###Primary Function ----
#' Set QC Type for Filtering
#'
#' Determines and sets the QC type for filtering based on the global QC pass status in the `master_list`.
#' If no viable QC type is found, the function stops execution and prints a detailed error message.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with the QC type set and filters initialized.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_set_qc(master_list)
#' }
qcCheckR_set_qc <- function(master_list) {
  qc_type <- determine_qc_type(master_list)

  if (qc_type == "unknown") {
    stop_with_qc_error(
      project_name = master_list$project_details$project_name,
      global_qc_pass = master_list$project_details$global_qc_pass,
      plate_qc_passed = master_list$project_details$qc_passed
    )
  } else {
    notify_qc_type(qc_type)
  }

  master_list$project_details$qc_type <- qc_type
  master_list$filters <- list()

  return(master_list)
}

###Sub Functions ----

#' Determine QC Type
#'
#' This function determines the QC type based on the global QC pass status.
#' It checks if there are multiple QC types that have passed.
#' Assesses if a secondary QC sample is available such as a pooled Quality controls.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return A string indicating the QC type ("pqc", "ltr", or "unknown").
determine_qc_type <- function(master_list) {
  global_qc_pass <- master_list$project_details$global_qc_pass
  passed_qc <- names(global_qc_pass[global_qc_pass == "pass"])
  user_supplied_qc <- master_list$project_details$qc_type

  if (is.null(passed_qc) || length(passed_qc) == 0) {
    message("No viable QC type found for filtering. Please check your QC_sample_label is correct.")
    return("unknown")

  } else if (length(passed_qc) == 1) {
    message(passed_qc[1], " is the only viable QC.")
    return(passed_qc[1])

  } else if (length(passed_qc) > 1) {
    if (user_supplied_qc %in% passed_qc) {
      message("Multiple valid QC types found. Reverting to default user choice.")
      return(master_list$project_details$qc_type)
    } else {
      message("Multiple QC types found, but user-supplied QC is not among them.")
      message("Valid QC types: ", paste(passed_qc, collapse = ", "))
      message("Please enter a valid QC_sample_label from the list above and rerun")
      return("unknown")
    }
  }

}




#' Notify QC Type
#'
#' This function prints a message indicating the QC type that has been set for filtering.
#' @keywords internal
#' @param qc_type A string indicating the QC type (e.g. "pqc" or "ltr").
#' @return NULL
notify_qc_type <- function(qc_type) {
  message <- paste("qcCheckeR has set filtering QC to:", qc_type)
  message(message, "\n")
}

#' Stop with QC Error
#'
#' This function stops the script execution and prints a detailed error message if no viable QC type is found.
#' It includes information about the global QC pass status and plate QC assessment.
#' @keywords internal
#' @param project_name A string containing the project name.
#' @param global_qc_pass A list containing the global QC pass status.
#' @param plate_qc_passed A list containing the plate QC assessment.
#' @return Stops the script execution with an error message.
stop_with_qc_error <- function(project_name,
                               global_qc_pass,
                               plate_qc_passed) {
  global_qc_pass_str <- utils::capture.output(str(global_qc_pass))
  plate_qc_passed_str <- utils::capture.output(str(plate_qc_passed))

  error_message <- paste0(
    "STOPPING SCRIPT\n",
    "There are no viable qc samples for RSD filtering: ",
    project_name,
    ".\n",
    "Please refer to the details below:\n",
    "Global QC Assessment:\n",
    paste(global_qc_pass_str, collapse = "\n"),
    "\n",
    "Plate QC Assessment:\n",
    paste(plate_qc_passed_str, collapse = "\n")
  )

  stop(error_message)
}


#. ----
## Sample Filter ----
###Primary Function ----
#' Sample Filter
#'
#' Flags samples based on missing values and summed signal intensity in the `master_list` data.
#' Applies thresholds to identify low-quality samples and aggregates results across plates.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with sample filter flags and failed sample lists.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_sample_filter(master_list)
#' }
qcCheckR_sample_filter <- function(master_list) {
  master_list$filters$samples.missingValues <- list()
  master_list$filters$failed_samples <- list()

  for (idx_batch in names(master_list$data$peakArea$sorted)) {
    sample_data <- master_list$data$peakArea$sorted[[idx_batch]]
    sample_meta <- sample_data %>%
      dplyr::select(sample_run_index,
             sample_name,
             sample_plate_id,
             sample_type_factor)

    #Extract lipid data
    lipid_data <- sample_data %>%
      dplyr::select(!dplyr::contains("sample")) %>%
      dplyr::select(!dplyr::contains("SIL")) %>%
      as.matrix()
    # Extract SIL data
    sil_data <- sample_data %>%
      dplyr::select(!dplyr::contains("sample")) %>%
      dplyr::select(dplyr::contains("SIL")) %>%
      as.matrix()

    #Flags for missing values and signal intensity
    flags <- list(
      #Summed signal intensities
      summed.lipid.signal = rowSums(lipid_data, na.rm = TRUE),
      summed.SIL.Int.Std.signal = rowSums(sil_data, na.rm = TRUE),
      #Missing /zero values <5000
      missing.lipid = rowSums(lipid_data < 5000, na.rm = TRUE),
      missing.SIL = rowSums(sil_data < 5000, na.rm = TRUE),
      #NA values
      na.lipid = rowSums(is.na(lipid_data), na.rm = TRUE),
      na.SIL = rowSums(is.na(sil_data), na.rm = TRUE),
      #NaN values
      nan.lipid = rowSums(is.nan(lipid_data), na.rm = TRUE),
      nan.SIL = rowSums(is.nan(sil_data), na.rm = TRUE),
      #Infinite values
      inf.lipid = rowSums(is.infinite(lipid_data), na.rm = TRUE),
      inf.SIL = rowSums(is.infinite(sil_data), na.rm = TRUE)
    )

    #Join metadata with flags
    flags <- cbind(sample_meta, flags)

    flags$totalMissingValues.lipid <- rowSums(do.call(cbind, flags[c("missing.lipid", "na.lipid", "nan.lipid", "inf.lipid")]))
    flags$totalMissingValues.SIL <- rowSums(do.call(cbind, flags[c("missing.SIL", "na.SIL", "nan.SIL", "inf.SIL")]))

    flags$sample.lipid.intensity.flag <- as.integer(flags$summed.lipid.signal < stats::median(flags$summed.lipid.signal) * 0.20)
    flags$sample.SIL.Int.Std.intensity.flag <- as.integer(
      flags$summed.SIL.Int.Std.signal < stats::median(flags$summed.SIL.Int.Std.signal) * 0.33
    )

    lipid_threshold <- ncol(lipid_data) * (master_list$project_details$mv_sample_threshold / 100)
    sil_threshold <- ncol(sil_data) * 0.33

    flags$sample.missing.value.flag <- as.integer(
      flags$totalMissingValues.lipid > lipid_threshold |
        flags$totalMissingValues.SIL > sil_threshold
    )

    flags$sample.flag <- as.integer(rowSums(as.matrix(flags[c(
      "sample.lipid.intensity.flag",
      "sample.SIL.Int.Std.intensity.flag",
      "sample.missing.value.flag"
    )]), na.rm = TRUE) > 0)

    flags$failed_samples <- ifelse(flags$sample.flag == 1, flags$sample_name, NA)

    master_list$filters$samples.missingValues[[idx_batch]] <- as.data.frame(flags)
    master_list$filters$failed_samples[[idx_batch]] <- flags$failed_samples
  }

  master_list$filters$samples.missingValues <- do.call(rbind, master_list$filters$samples.missingValues)
  master_list$filters$failed_samples <- master_list$filters$samples.missingValues %>%
                                        filter(!is.na(failed_samples)) %>%
                                        pull(sample_name)

  return(master_list)
}

#.----
##SIL Internal Standard Filter ----
###Primary Function ----
#' SIL Internal Standard Filter
#'
#' Filters SIL internal standards based on missing values in the `master_list` data.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with SIL internal standard filter flags.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_sil_IntStd_filter(master_list)
#' }
qcCheckR_sil_IntStd_filter <- function(master_list) {
  master_list$filters$sil.intStd.missingValues <- list()
  master_list$filters$sil.intStd.missingValues$summary <- initialise_sil_summary()

  for (idx_batch in names(master_list$data$peakArea$sorted)) {
    sil_flags <- calculate_sil_flags_per_plate(master_list, idx_batch)
    master_list$filters$sil.intStd.missingValues[[idx_batch]] <- sil_flags
    master_list$filters$sil.intStd.missingValues$summary <- rbind(master_list$filters$sil.intStd.missingValues$summary,
                                                                  sil_flags)
  }

  master_list <- calculate_sil_flags_per_version(master_list)

  master_list$filters$failed_sil.intStds <- unique(unlist(master_list$filters$failed_sil.intStds, use.names = FALSE))

  return(master_list)
}


### Sub Functions ----
#' Initialise SIL Summary Data Frame
#'
#' This function initializes an empty data frame to store SIL summary statistics.
#' It includes columns for lipid names, template versions, plate IDs, and various flags related to SIL internal standards.
#' @keywords internal
#' @return An empty data frame with specified columns for SIL summary statistics.
initialise_sil_summary <- function() {
  data.frame(
    lipid = character(),
    template_version = character(),
    plateID = character(),
    peakArea_5000_LOD = numeric(),
    naValues = numeric(),
    nanValues = numeric(),
    infValues = numeric(),
    totalMissingValues = numeric(),
    flag_SIL_intStd_Plate = numeric(),
    stringsAsFactors = FALSE
  )
}

#' Calculate SIL Flags per Plate
#'
#' This function calculates flags for SIL internal standards on a per-plate basis.
#' It computes the number of peak areas below a threshold, counts missing values, and flags plates with excessive missing values.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param idx_batch The index of the batch (plate) to process.
#' @return A tibble containing SIL flags for each lipid, including counts of peak areas below a threshold, missing values, and flags for excessive missing values.
calculate_sil_flags_per_plate <- function(master_list, idx_batch) {
  sil_names <- master_list$data$peakArea$sorted[[idx_batch]] %>%
    dplyr::select(dplyr::contains("SIL")) %>%
    names()

  sil_matrix <- master_list$data$peakArea$sorted[[idx_batch]] %>%
    dplyr::filter(
      sample_name %in% dplyr::filter(master_list$filters$samples.missingValues, sample.flag == 0)[["sample_name"]]
    ) %>%
    dplyr::select(dplyr::contains("SIL")) %>%
    as.matrix()

  flags <- tibble(
    lipid = sil_names,
    peakArea_5000_LOD = colSums(sil_matrix < 5000, na.rm = TRUE),
    naValues = colSums(is.na(sil_matrix), na.rm = TRUE),
    nanValues = colSums(is.nan(sil_matrix), na.rm = TRUE),
    infValues = colSums(is.infinite(sil_matrix), na.rm = TRUE)
  )

  flags$totalMissingValues <- rowSums(flags %>% dplyr::select(-lipid), na.rm = TRUE)

  #Flag if SIL is missing in greater than 5% of samples
  valid_sample_count <- nrow(sil_matrix)
  flags$flag_SIL_intStd_Plate <- as.integer(flags$totalMissingValues > (valid_sample_count * 0.05))

  flags$template_version <- master_list$templates$`Plate SIL version`[[idx_batch]]
  flags$plateID <- idx_batch

  return(flags)
}

#' Calculate SIL Flags per Version
#'
#' This function calculates flags for SIL internal standards across different
#' template versions. It aggregates SIL flags from all plates for each version,
#' counts missing values, and flags versions with excessive missing values.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with SIL flags calculated for each version.
calculate_sil_flags_per_version <- function(master_list) {
  master_list$filters$sil.intStd.missingValues$PROJECT.flag.SIL.intStd <- list()
  master_list$filters$failed_sil.intStds <- list()

  if (!"template_version" %in% colnames(master_list$filters$sil.intStd.missingValues$summary)) {
    stop("Missing required column: template_version")
  }

  for (version in unique(master_list$filters$sil.intStd.missingValues$summary$template_version)) {
    version_data <- master_list$filters$sil.intStd.missingValues$summary %>%
      dplyr::filter(template_version == version)

    lipid_list <- unique(version_data$lipid)
    plate_list <- unique(version_data$plateID)

    sil_matrix <- master_list$data$peakArea$sorted %>%
      dplyr::bind_rows() %>%
      dplyr::filter(sample_plate_id %in% plate_list) %>%
      dplyr::select(tidyselect::where(~ !all(is.na(.)))) %>%
      dplyr::filter(
        sample_name %in% dplyr::filter(
          master_list$filters$samples.missingValues,
          sample.flag == 0
        )[["sample_name"]]
      ) %>%
      dplyr::select(dplyr::contains("SIL")) %>%
      as.matrix()

    version_flags <- tibble::tibble(
      lipid = colnames(sil_matrix),
      peakArea_5000_LOD = colSums(sil_matrix < 5000, na.rm = TRUE),
      naValues = colSums(is.na(sil_matrix), na.rm = TRUE),
      nanValues = colSums(is.nan(sil_matrix), na.rm = TRUE),
      infValues = colSums(is.infinite(sil_matrix), na.rm = TRUE)
    )

    version_flags$totalMissingValues <- rowSums(version_flags %>% dplyr::select(-lipid), na.rm = TRUE)
    master_list$filters$sil.intStd.missingValues$allPlates[[version]] <- version_flags

    valid_sample_count <- nrow(sil_matrix)
    master_list$filters$sil.intStd.missingValues$PROJECT.flag.SIL.intStd[[version]] <- as.integer(version_flags$totalMissingValues > (valid_sample_count * 0.05))

    master_list$filters$failed_sil.intStds[[version]] <- version_flags$lipid[master_list$filters$sil.intStd.missingValues$PROJECT.flag.SIL.intStd[[version]] == 1]
  }

  return(master_list)
}


# . ----

## Lipid Filter ----
###Primary Function ----
#' Lipid Filter
#'
#' Filters lipids based on missing values across plates and template versions.
#' Flags lipids with more than 50% missing values and compiles a list of failed lipids.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with lipid filter flags and failed lipid list.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_lipid_filter(master_list)
#' }
qcCheckR_lipid_filter <- function(master_list) {
  master_list <- initialise_lipid_filter(master_list)

  for (idx_batch in names(master_list$data$peakArea$sorted)) {
    lipid_data <- get_lipid_data(master_list, idx_batch)
    lipid_flags <- calculate_lipid_flags(master_list, idx_batch, lipid_data)

    master_list$filters$lipid.missingValues[[idx_batch]] <- lipid_flags
    master_list$filters$lipid.missingValues$summary <- rbind(master_list$filters$lipid.missingValues$summary,
                                                             lipid_flags)
  }

  master_list <- process_lipid_versions(master_list)

  master_list$filters$failed_lipids <- unique(unlist(master_list$filters$failed_lipids, use.names = FALSE))

  return(master_list)
}

### Sub Functions ----

#' Initialise lipid filter structure
#'
#' This function initialses the structure for lipid filtering in the `master_list`.
#' It creates a list to store lipid missing values, a summary data frame, and lists for project flags and failed lipids.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with initialised lipid filter structure.
initialise_lipid_filter <- function(master_list) {
  master_list$filters$lipid.missingValues <- list()
  master_list$filters$lipid.missingValues$summary <- data.frame(
    lipid = character(),
    silFilter.flag.Lipid = numeric(),
    peakArea_5000_LOD = numeric(),
    naValues = numeric(),
    nanValues = numeric(),
    infValues = numeric(),
    totalMissingValues = numeric(),
    flag.Lipid.Plate = numeric(),
    template_version = character(),
    plateID = character(),
    stringsAsFactors = FALSE
  )
  master_list$filters$lipid.missingValues$PROJECT.flag.lipid <- list()
  master_list$filters$failed_lipids <- list()

  return(master_list)
}

#' Extract lipid matrix for a batch
#'
#' This function extracts the lipid data matrix for a specific batch from the `master_list`.
#' It filters out samples with missing values and selects relevant columns.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param idx_batch The index of the batch (plate) to process.

get_lipid_data <- function(master_list, idx_batch) {
  valid_samples <- master_list$filters$samples.missingValues %>%
    dplyr::filter(sample.flag == 0) %>%
    dplyr::pull(sample_name)

  lipid_data <- master_list$data$peakArea$sorted[[idx_batch]] %>%
    dplyr::filter(sample_name %in% valid_samples) %>%
    dplyr::select(!dplyr::contains("sample") & !dplyr::contains("SIL")) %>%
    as.matrix()

  return(lipid_data)
}


#' Calculate lipid flags for a batch
#'
#' This function calculates flags for lipids based on their signal intensity and missing values in a specific batch.
#' It checks for SIL internal standards, counts peak areas below a threshold, and flags plates with excessive missing values.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param idx_batch The index of the batch (plate) to process.
#' @param lipid_matrix A matrix containing lipid data for the batch.
#' @return A tibble containing lipid flags, including counts of peak areas below a threshold, missing values, and flags for excessive missing values.
calculate_lipid_flags <- function(master_list, idx_batch, lipid_matrix) {

  lipid_names <- colnames(lipid_matrix)
  SIL_version <- master_list$templates$`Plate SIL version`[[idx_batch]]
  failed_sil <- master_list$filters$failed_sil.intStds

  if (ncol(lipid_matrix) == 0) {
    return(tibble::tibble(
      lipid = character(),
      silFilter.flag.Lipid = numeric(),
      peakArea_5000_LOD = numeric(),
      naValues = numeric(),
      nanValues = numeric(),
      infValues = numeric(),
      totalMissingValues = numeric(),
      flag.Lipid.Plate = numeric(),
      template_version = SIL_version,
      plateID = idx_batch
    )[0, ])
  }

  if (is.null(SIL_version)) {
    stop("Missing template version for batch: ", idx_batch)
  }

  if (is.null(master_list$templates$mrm_guides[[SIL_version]]$SIL_guide)) {
    sil_flag <- rep(0, length(lipid_names))
  } else {
    sil_flag <- as.integer(
      lipid_names %in% dplyr::filter(
        master_list$templates$mrm_guides[[SIL_version]]$SIL_guide,
        Note %in% failed_sil
      )[["Precursor Name"]]
    )
  }

  lipid_flags <- tibble::tibble(
    lipid = lipid_names,
    silFilter.flag.Lipid = sil_flag,
    peakArea_5000_LOD = colSums(lipid_matrix < 5000, na.rm = TRUE),
    naValues = colSums(is.na(lipid_matrix), na.rm = TRUE),
    nanValues = colSums(is.nan(lipid_matrix), na.rm = TRUE),
    infValues = colSums(is.infinite(lipid_matrix), na.rm = TRUE),
    template_version = SIL_version,
    plateID = idx_batch
  )
  lipid_flags$totalMissingValues <- rowSums(lipid_flags %>% dplyr::select(-c(lipid, template_version, plateID)), na.rm = TRUE)
  lipid_flags$flag.Lipid.Plate <- as.integer(lipid_flags$totalMissingValues > (
    nrow(lipid_matrix) * (master_list$project_details$mv_sample_threshold / 100)
  ) |
    lipid_flags$silFilter.flag.Lipid == 1)

  return(lipid_flags)
}

#' Process lipid flags across template versions
#'
#' This function processes lipid flags across different template versions.
#' It aggregates lipid flags from all plates for each version, counts missing values, and flags versions with excessive missing values.
#' @keywords internal
#' @param master_list A list containing project details and data.
process_lipid_versions <- function(master_list) {
  for (version in unique(master_list$filters$lipid.missingValues$summary$template_version)) {
    version_data <- master_list$filters$lipid.missingValues$summary %>%
      dplyr::filter(template_version == version)

    lipid_list <- unique(version_data$lipid)
    plate_list <- unique(version_data$plateID)

    lipid_matrix <- master_list$data$peakArea$sorted %>%
      dplyr::bind_rows() %>%
      dplyr::filter(sample_plate_id %in% plate_list) %>%
      dplyr::select(tidyselect::where(~ !all(is.na(.)))) %>%
      dplyr::filter(
        sample_name %in% dplyr::filter(
          master_list$filters$samples.missingValues,
          sample.flag == 0
        )[["sample_name"]]
      ) %>%
      dplyr::select(!dplyr::contains("sample") & !dplyr::contains("SIL")) %>%
      as.matrix()

    version_flags <- tibble::tibble(
      lipid = colnames(lipid_matrix),
      peakArea_5000_LOD = colSums(lipid_matrix < 5000, na.rm = TRUE),
      naValues = colSums(is.na(lipid_matrix), na.rm = TRUE),
      nanValues = colSums(is.nan(lipid_matrix), na.rm = TRUE),
      infValues = colSums(is.infinite(lipid_matrix), na.rm = TRUE)
    )

    version_flags$totalMissingValues <- rowSums(version_flags %>% dplyr::select(-lipid), na.rm = TRUE)
    master_list$filters$lipid.missingValues$allPlates[[version]] <- version_flags

    valid_sample_count <- nrow(lipid_matrix)
    master_list$filters$lipid.missingValues$PROJECT.flag.lipid[[version]] <- as.integer(version_flags$totalMissingValues > (
      valid_sample_count * (master_list$project_details$mv_sample_threshold / 100)
    ))

    master_list$filters$failed_lipids[[version]] <- version_flags$lipid[master_list$filters$lipid.missingValues$PROJECT.flag.lipid[[version]] == 1]
  }
  return(master_list)
}



#. ----

##RSD Filter ----
###Primary Function ----
#' RSD Filter
#'
#' Filters features per plate with a %RSD > 30% based on the precision of measurement in the `master_list` data.
#' Applies filtering to peakArea, concentration, and statTarget concentration data sources.
#' @keywords internal
#' @param master_list Master list from previous functions.
#' @return The updated `master_list` with RSD filter flags.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_RSD_filter(master_list)
#' }
qcCheckR_RSD_filter <- function(master_list) {
  master_list$filters$rsd <- tibble::tibble()

  # Apply RSD filtering for each data source
  master_list$filters$rsd <- dplyr::bind_rows(
    calculate_rsd(master_list, "peakArea", master_list$data$peakArea$imputed),
    calculate_rsd(master_list, "peakArea", list(
      allBatches = dplyr::bind_rows(master_list$data$peakArea$imputed)
    )),
    calculate_rsd(
      master_list,
      "concentration",
      master_list$data$concentration$imputed
    ),
    calculate_rsd(master_list, "concentration", list(
      allBatches = dplyr::bind_rows(master_list$data$concentration$imputed)
    )),
    calculate_rsd(
      master_list,
      "concentration[statTarget]",
      master_list$data$concentration$statTargetProcessed
    ),
    calculate_rsd(
      master_list,
      "concentration[statTarget]",
      list(
        allBatches = dplyr::bind_rows(master_list$data$concentration$statTargetProcessed)
      )
    )
  )

  # Clean and format RSD table
  master_list$filters$rsd <- master_list$filters$rsd %>%
    dplyr::rename(dataSource = V1, dataBatch = V2) %>%
    dplyr::mutate(dplyr::across(!dplyr::contains("data"), as.numeric)) %>%
    dplyr::mutate(dplyr::across(!dplyr::contains("data"), round, 2))

  # Update script log
  master_list <- update_script_log(master_list,
                                   "data_filtering",
                                   "data_preparation",
                                   "summary_report")

  return(master_list)
}

###Sub Functions ----
#'Calculate RSD for a given data source and batch list
#'
#'This function calculates the relative standard deviation (RSD) for each feature in the provided data batches.
#'It filters out failed samples and selects only QC samples, then computes the RSD values.
#' @keywords internal
#' @param master_list Master list containing project details and data.
#' @param source_name Name of the data source (e.g., "peakArea", "concentration").
#' @param data_batches A list of data batches to process.
calculate_rsd <- function(master_list, source_name, data_batches) {
  rsd_results <- list()

  for (batch_name in names(data_batches)) {
    data <- data_batches[[batch_name]] %>%
      filter(!sample_name %in% master_list$filters$failed_samples) %>%
      filter(sample_type %in% c("qc")) %>%
      select(!contains("sample")) %>%
      select(!contains("SIL"))

    if (nrow(data) == 0)
      next

    rsd_values <- (apply(data, 2, stats::sd) / apply(data, 2, mean)) * 100
    rsd_results[[batch_name]] <- rbind(c(source_name, batch_name, rsd_values)) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(dplyr::across(names(rsd_values), as.numeric))
  }

  bind_rows(rsd_results)
}

#. ----

#Phase 3: Summary Report----
###Primary Function ----
#' Summary Report
#'
#' Generates a summary report for the `master_list` data,
#' including metrics for cohorts, matrix types, sample counts, lipid targets,
#' SIL versions, missing value filter flags, and RSD percentages.
#' @keywords internal
#' @param master_list Master list generated by previous functions.
#' @return The updated `master_list` with the summary report.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_summary_report(master_list)
#' }
qcCheckR_summary_report <- function(master_list) {
  #Create sample metrics
  sample_tags <- dplyr::bind_rows(master_list$data$peakArea$sorted) %>%
    dplyr::select(sample_type_factor) %>%
    unique() %>%
    dplyr::pull(sample_type_factor) %>%
    as.character() %>%
    .[!grepl("sample", .)]

  sample_metrics <- c()
  for (tag in sample_tags) {
    sample_metrics <- c(sample_metrics, paste0(tag, "Samples"))
  }

  metrics <- c(
    "MatrixType",
    "totalSamples",
    "studySamples",
    paste0(sample_metrics),
    "lipidTargets",
    "matchedLipidTargets",
    "SIL.version",
    "SIL.IntStds",
    "missingValueFilterFlags[samples]",
    "missingValueFilterFlags[SIL.IS]",
    "missingValueFilterFlags[lipidTargets]",
    "rsd<30%[peakArea]",
    "rsd<20%[peakArea]",
    "rsd<10%[peakArea]",
    "rsd<30%[concentration]",
    "rsd<20%[concentration]",
    "rsd<10%[concentration]",
    "rsd<30%[concentration.statTarget]",
    "rsd<20%[concentration.statTarget]",
    "rsd<10%[concentration.statTarget]"
  )

  master_list$summary_tables <- list()
  master_list$summary_tables$projectOverview <- tibble::tibble(metric = metrics)

  # Per-plate summary
  for (idx_batch in names(master_list$data$peakArea$sorted)) {
    plate_summary <- generate_plate_summary(master_list, idx_batch, metrics, sample_tags)
    master_list$summary_tables$projectOverview <- dplyr::left_join(master_list$summary_tables$projectOverview,
                                                            plate_summary,
                                                            by = "metric")
  }

  # Inter-plate summary
  inter_plate_summary <- generate_inter_plate_summary(master_list, metrics, sample_tags)
  master_list$summary_tables$projectOverview <- dplyr::left_join(master_list$summary_tables$projectOverview,
                                                          inter_plate_summary,
                                                          by = "metric")

  # Update script log
  master_list <- update_script_log(master_list,
                                   "summary_report",
                                   "data_filtering",
                                   "plot_generation")

  return(master_list)
}

###Sub Functions ----
#' Generate Plate Summary
#'
#' This function generates a summary for a specific plate in the `master_list`.
#' It includes metrics such as matrix type, sample counts, lipid targets,
#' SIL versions, missing value filter flags, and RSD percentages.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param idx_batch The index of the batch (plate) to process.
#' @param metrics tibble of the metrics to be included in the summary.
#' @param sample_tags A vector of sample type tags to be included in the summary.
#' @return A tibble containing the summary metrics for the specified plate.
generate_plate_summary <- function(master_list,
                                   idx_batch,
                                   metrics,
                                   sample_tags) {
  data <- master_list$data$peakArea$sorted[[idx_batch]]
  lipid_data <- data %>% dplyr::select(-dplyr::contains("sample"), -dplyr::contains("SIL"))
  sil_data <- data %>% dplyr::select(dplyr::contains("SIL"))

  #Create RSD metrics function
  rsd_metrics <- function(source) {
    rsd <- master_list$filters$rsd %>%
      dplyr::filter(dataBatch == idx_batch, dataSource == source) %>%
      dplyr::select(-dplyr::contains("data")) %>%
      dplyr::select(-dplyr::any_of(master_list$filters$failed_lipids))
    c(sum(rsd < 30, na.rm = TRUE),
      sum(rsd < 20, na.rm = TRUE),
      sum(rsd < 10, na.rm = TRUE))
  }

  plate_summary <- tibble::tibble(
    metric = metrics,
    !!idx_batch := c(
      unique(data$sample_matrix),
      nrow(data),
      nrow(dplyr::filter(data, sample_type_factor == "sample")),
      sapply(sample_tags, function(tag)
        nrow(dplyr::filter(
          data, sample_type_factor == tag
        ))),
      ncol(lipid_data),
      ncol(lipid_data),
      master_list$templates$`Plate SIL version`[[idx_batch]],
      ncol(sil_data),
      nrow(
        dplyr::filter(
          master_list$filters$samples.missingValues,
          sample_plate_id == idx_batch,
          sample.flag == 1
        )
      ),
      sum(master_list$filters$sil.intStd.missingValues[[idx_batch]][["flag_SIL_intStd_Plate"]] == 1),
      sum(master_list$filters$lipid.missingValues[[idx_batch]][["flag.Lipid.Plate"]] == 1),
      rsd_metrics("peakArea"),
      rsd_metrics("concentration"),
      rsd_metrics("concentration[statTarget]")
    ) %>% unlist()
  )

  return(plate_summary)
}

#' Generate Inter-Plate Summary
#'
#' This function generates a summary of all plates in the `master_list`.
#' It aggregates data across all plates, including metrics such as matrix type,
#' sample counts, lipid targets, SIL versions, missing value filter flags,
#' and RSD percentages.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param metrics tibble of the metrics to be included in the summary.
#' @param sample_tags A vector of sample type tags to be included in the summary.
#' @return A tibble containing the inter-plate summary metrics.
generate_inter_plate_summary <- function(master_list, metrics, sample_tags) {
  data <- dplyr::bind_rows(master_list$data$peakArea$sorted)
  lipid_data <- data %>% dplyr::select(-dplyr::contains("sample"), -dplyr::contains("SIL"))
  sil_data <- data %>% dplyr::select(dplyr::contains("SIL"))

  rsd_metrics <- function(source) {
    rsd <- master_list$filters$rsd %>%
      dplyr::filter(dataBatch == "allBatches", dataSource == source) %>%
      dplyr::select(-dplyr::contains("data")) %>%
      dplyr::select(-dplyr::any_of(master_list$filters$failed_lipids))
    c(sum(rsd < 30, na.rm = TRUE),
      sum(rsd < 20, na.rm = TRUE),
      sum(rsd < 10, na.rm = TRUE))
  }

  interplate_summary <- tibble::tibble(
    metric = metrics,
    all_plates = c(
      paste(unique(data$sample_matrix), collapse = ","),
      nrow(data),
      nrow(dplyr::filter(data, sample_type_factor == "sample")),
      sapply(sample_tags, function(tag)
        nrow(dplyr::filter(
          data, sample_type_factor == tag
        ))),
      ncol(lipid_data),
      ncol(dplyr::bind_rows(master_list$data$concentration$statTargetProcessed) %>%
             dplyr::select(-dplyr::contains("sample"))),
      paste(
        unique(master_list$templates$`Plate SIL version`),
        collapse = ","
      ),
      ncol(sil_data),
      nrow(
        dplyr::filter(master_list$filters$samples.missingValues, sample.flag == 1)
      ),
      length(master_list$filters$failed_sil.intStds),
      length(master_list$filters$failed_lipids),
      rsd_metrics("peakArea"),
      rsd_metrics("concentration"),
      rsd_metrics("concentration[statTarget]")
    ) %>% unlist()
  )
  return(interplate_summary)
}

#. ----

#Phase 4: Plot Generation----
## Plot Options ----
###Primary Function ----
#' Set Plot Options
#'
#' Sets the plot color, fill, shape, and size options for the `master_list` data based on sample types and QC type.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with plot options set.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_plot_options(master_list)
#' }
qcCheckR_plot_options <- function(master_list) {

  sample_types <- master_list$project_details$sample_tags
  colors <- viridis::viridis(n = length(sample_types))
  master_list$project_details$plot_fill <- stats::setNames(colors, sample_types)

  master_list$project_details$plot_colour <- purrr::set_names(rep("black", length(sample_types)), sample_types)

  master_list$project_details$plot_shape <- purrr::set_names(rep(21, length(sample_types)), sample_types)

  master_list$project_details$plot_size <- purrr::set_names(rep(2, length(sample_types)), sample_types)

  qc_type <- tolower(master_list$project_details$qc_type)

  if (qc_type %in% names(master_list$project_details$plot_colour)) {
    master_list$project_details$plot_colour[[qc_type]] <- "red"
    master_list$project_details$plot_shape[[qc_type]] <- 23
    master_list$project_details$plot_size[[qc_type]] <- 3
  }

  return(master_list)
}

#.----

## PCA Analysis and Plotting ----
### Primary Function ----
#' PCA Analysis and Plotting
#'
#' Performs PCA analysis on the `master_list` data and generates PCA plots using `ggplot2`.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with PCA models, scores, and ggplot objects.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_PCA(master_list)
#' }
qcCheckR_PCA <- function(master_list) {
  master_list$pca <- list(models = list(),
                          scores = list(),
                          plot = list())

  for (source in c("peakArea", "concentration")) {
    master_list <- run_pca_model(master_list, source, preprocessed = FALSE)
    master_list <- run_pca_model(master_list, source, preprocessed = TRUE)
  }

  for (fill_var in c("sample_type_factor", "sample_plate_id")) {
    master_list$pca$plot[[fill_var]] <- generate_pca_ggplot(master_list, fill_var)
  }

  return(master_list)
}

###Sub Functions ----
#' Run PCA Model
#'
#' This function runs a PCA model on the specified data source from the `master_list`.
#' It preprocesses the data by filtering out failed samples and high RSD lipids if `preprocessed` is set to TRUE.
#' It then performs PCA using the `ropls` package and stores the model and scores in the `master_list`.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param source The data source to run PCA on (e.g., "peakArea", "concentration", "concentration.statTarget").
#' @param preprocessed Logical indicating whether to preprocess the data (default is FALSE).
#' @return The updated `master_list` with PCA models and scores.
run_pca_model <- function(master_list, source, preprocessed = FALSE) {
  # Determine the data key based on the source
  if (source == "peakArea") {
    data_keys <- "imputed"
  } else if (source == "concentration") {
    data_keys <- c("imputed", "statTargetProcessed")
  }

  for (data_key in data_keys) {
    #validate source and data_key
    if (!source %in% names(master_list$data) ||
        !data_key %in% names(master_list$data[[source]])) {
      stop(paste(
        "Data source",
        source,
        "or data key",
        data_key,
        "not found in master_list."
      ))
    }

    #Here we remove the cols that are not present in all three dataframes
    cols_to_keep <- colnames(master_list$data$concentration$sorted[[1]])
    data <- bind_rows(master_list$data[[source]][[data_key]]) %>%
      select(all_of(cols_to_keep))




    if (preprocessed) {
      data <- data %>%
        dplyr::filter(!sample_name %in% master_list$filters$failed_samples)

      failed_lipids <- master_list$filters$failed_lipids
      high_rsd_lipids <- master_list$filters$rsd %>%
        dplyr::filter(dataSource == source, dataBatch == "allBatches") %>%
        dplyr::select(-dplyr::contains("data")) %>%
        dplyr::summarise(dplyr::across(dplyr::everything(), ~ ifelse(. > 30, TRUE, FALSE))) %>%
        dplyr::select(tidyselect::where(~ any(., na.rm = TRUE))) %>%
        names()

      data <- data %>%
        dplyr::select(-dplyr::any_of(c(failed_lipids, high_rsd_lipids)))
    }

    pca_matrix <- data %>%
      tibble::column_to_rownames("sample_name") %>%
      dplyr::select(-dplyr::contains("sample"), -dplyr::contains("SIL")) %>%
      as.matrix()

    model_name <- if (preprocessed)
      paste0(source, ".", data_key, ".preProcessed")
    else
      paste0(source, ".", data_key)

    master_list$pca$models[[model_name]] <- ropls::opls(
      x = pca_matrix,
      y = NULL,
      crossvalI = 1,
      predI = 3,
      algoC = "nipals",
      log10L = FALSE,
      scale = "pareto",
      plotSubC = NA,
      fig.pdfC = "none"
    )

    scores <- master_list$pca$models[[model_name]]@scoreMN %>%
      tibble::as_tibble() %>%
      dplyr::rename(PC1 = p1, PC2 = p2, PC3 = p3) %>%
      tibble::add_column(sample_name = rownames(pca_matrix), .before = 1) %>%
      dplyr::left_join(data %>% dplyr::select(dplyr::contains("sample")), by = "sample_name")

    scores$sample_data_source <- paste0(
      if (preprocessed)
        paste0(source, ".", data_key, ".preProcessed")
      else
        paste0(source, ".", data_key),
      ": s=",
      nrow(pca_matrix),
      ", l=",
      ncol(pca_matrix)
    )

    master_list$pca$scores[[model_name]] <- scores
  }

  return(master_list)
}

#' Generate PCA ggplot
#'
#' This function generates a PCA plot using `ggplot2` for the PCA scores stored in the `master_list`.
#' It allows for coloring the points by a specified variable (e.g., sample type or plate ID).
#' @keywords internal
#' @param master_list A list containing project details and PCA scores.
#' @param fill_var The variable to color the points by (e.g., "sample_type_factor", "sample_plate_id").
#' @return A `ggplot` object representing the PCA scores.
generate_pca_ggplot <- function(master_list, fill_var) {
  all_scores <- dplyr::bind_rows(master_list$pca$scores)

  #Factor and arrange to ensure consistent ordering in the plot
  all_scores$source_prefix <- sub(":.*", "", all_scores$sample_data_source) %>%
    factor(
      levels = c(
        "peakArea.imputed",
        "peakArea.imputed.preProcessed",
        "concentration.imputed",
        "concentration.imputed.preProcessed",
        "concentration.statTargetProcessed",
        "concentration.statTargetProcessed.preProcessed"
      )
    )
  all_scores <- all_scores %>% dplyr::arrange(source_prefix)
  all_scores$source_suffix <- sub(".*: ", "", all_scores$sample_data_source)
  all_scores$facet_label <- paste0(all_scores$source_prefix, ":", all_scores$source_suffix)
  all_scores$facet_label <- factor(all_scores$facet_label, levels = unique(all_scores$facet_label))

  plot <- plotly::ggplotly(
    ggplot2::ggplot(
      data = all_scores,
      ggplot2::aes(
        x = PC1,
        y = PC2,
        group = sample_name,
        fill = .data[[fill_var]],
        color = sample_type_factor,
        shape = sample_type_factor,
        size = sample_type_factor
      )
    ) +
      ggplot2::geom_vline(xintercept = 0, colour = "darkgrey") +
      ggplot2::geom_hline(yintercept = 0, color = "darkgrey") +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::scale_shape_manual(values = master_list$project_details$plot_shape) +
      ggplot2::scale_color_manual(values = master_list$project_details$plot_colour) +
      ggplot2::scale_size_manual(values = master_list$project_details$plot_size) +
      ggplot2::guides(
        shape = "none",
        size = "none",
        color = "none",
        fill = ggplot2::guide_legend(title = fill_var)
      ) +
      ggplot2::facet_wrap(
        facets = vars(facet_label),
        scales = "free",
        ncol = 2,
        nrow = 3
      ) +
      ggplot2::labs(
        title = paste0(
          "PCA scores; coloured by ",
          fill_var,
          "; ",
          master_list$project_details$project_name
        )
      )
  )
  return(plot)
}

#.----

## Run Order Plots ----
###Primary Function ----
#' Run Order Plots
#'
#' Generates run order plots for the `master_list` data, showing PCA scores versus run order using `ggplot2`.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with run order plots stored in `master_list$pca$scoresRunOrder`.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_run_order_plots(master_list)
#' }
qcCheckR_run_order_plots <- function(master_list) {
  #Get plate boundaries and labels
  boundary_data <- dplyr::bind_rows(master_list$pca$scores) %>%
    dplyr::select(sample_run_index, sample_plate_id) %>%
    dplyr::distinct()

  boundaries <- c(0.5)
  labels <- c()

  for (plate in unique(boundary_data$sample_plate_id)) {
    plate_data <- boundary_data %>% dplyr::filter(sample_plate_id == plate)
    boundaries <- c(
      boundaries,
      min(plate_data$sample_run_index) - 0.5,
      max(plate_data$sample_run_index) + 0.5
    )
    labels <- c(labels, stats::median(plate_data$sample_run_index))
  }

  boundaries <- unique(boundaries)

  run_order_labels <- tibble::tibble(
    sample_data_source = "plate_boundaries",
    sample_plate_id = rep(unique(boundary_data$sample_plate_id), each = 1),
    sample_run_index = labels,
    PC1 = 0,
    PC2 = 0,
    PC3 = 0,
    value = 0
  )

  master_list$pca$scoresRunOrder <- list()

  for (pc in c("PC1", "PC2", "PC3")) {
    master_list$pca$scoresRunOrder[[pc]] <- plot_run_order(
      scores = dplyr::bind_rows(master_list$pca$scores),
      pc = pc,
      boundaries = boundaries,
      annotations = run_order_labels,
      plot_settings = master_list$project_details
    )
  }

  return(master_list)
}

### Sub Functions ----
#' Plot Run Order
#'
#' This function generates a run order plot for PCA scores using `ggplot2`.
#' It plots the PCA scores against the sample run index, with vertical lines indicating plate boundaries and annotations for plate IDs.
#' @keywords internal
#' @param scores A tibble containing PCA scores and sample information.
#' @param pc The principal component to plot (e.g., "PC1", "PC2", "PC3").
#' @param boundaries A vector of boundaries for the plates.
#' @param annotations A tibble containing annotations for the plates.
#' @param plot_settings A list containing plot settings such as colors, shapes, and sizes.
#' @return A `ggplot` object representing the run order plot.
plot_run_order <- function(scores,
                           pc,
                           boundaries,
                           annotations,
                           plot_settings) {
  plotly::ggplotly(
    ggplot2::ggplot(
      scores,
      ggplot2::aes(
        x = sample_run_index,
        y = .data[[pc]],
        group = sample_name,
        fill = sample_type_factor,
        color = sample_type_factor,
        shape = sample_type_factor,
        size = sample_type_factor
      )
    ) +
      ggplot2::geom_vline(xintercept = boundaries, linetype = "dashed") +
      ggplot2::geom_point() +
      ggplot2::geom_text(
        data = annotations,
        ggplot2::aes(x = sample_run_index, y = .data[[pc]], label = sample_plate_id),
        inherit.aes = FALSE,
        color = "black"
      ) +
      ggplot2::theme_bw() +
      ggplot2::scale_shape_manual(values = plot_settings$plot_shape) +
      ggplot2::scale_fill_manual(values = plot_settings$plot_fill) +
      ggplot2::scale_color_manual(values = plot_settings$plot_colour) +
      ggplot2::scale_size_manual(values = plot_settings$plot_size) +
      ggplot2::ylab(pc) +
      ggplot2::guides(
        shape = "none",
        size = "none",
        color = "none",
        fill = ggplot2::guide_legend(title = "sample_type_factor")
      ) +
      ggplot2::facet_wrap(
        facets = vars(sample_data_source),
        ncol = 1,
        scales = "free_y"
      ) +
      ggplot2::labs(
        title = paste0(
          pc,
          "; run order (x) vs PCA scores (y); ",
          plot_settings$project_name
        )
      ) +
      ggplot2::theme(text = ggplot2::element_text(size = 12))
  )
}

#.----

## Target Control Charts ----
### Primary Function ----
#' Target Control Charts
#'
#' Generates control charts for the `master_list` data, showing the values of target metabolites and SIL internal standards across different sample types and data sources.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with control charts stored in `master_list$control_charts`.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_target_control_charts(master_list)
#' }
qcCheckR_target_control_charts <- function(master_list) {
  master_list$control_charts <- list()

  common_control_charts <- get_common_control_metabolites(master_list)
  plate_boundaries <- get_plate_boundaries(master_list)
  annotate_label <- get_plate_annotations(master_list)

  for (metabolite in common_control_charts) {
    master_list$control_charts[[metabolite]] <- plot_control_chart(
      master_list = master_list,
      metabolite = metabolite,
      plate_boundaries = plate_boundaries,
      annotate_label = annotate_label
    )
  }

  master_list <- update_script_log(master_list,
                                   "plot_generation",
                                   "summary_report",
                                   "data_exports")

  return(master_list)
}

###Sub Functions ----
#' Get Common Control Metabolites
#' This function retrieves the common control metabolites across different SIL versions in the `master_list`.
#' It intersects the precursor names of control charts from each SIL version's guide.
#' @keywords internal
#' @param master_list A list containing project details and templates.
#' @return A vector of common control metabolites.
get_common_control_metabolites <- function(master_list) {
  SIL_versions <- unique(unlist(master_list$templates[["Plate SIL version"]]))
  precursor_lists <- list()
  for (ver in SIL_versions) {
    sil_guide <- master_list$templates$mrm_guides[[ver]]$SIL_guide
    precursor_lists[[ver]] <- sil_guide$`Precursor Name`[sil_guide$control_chart == TRUE]
  }
  Reduce(dplyr::intersect, precursor_lists)
}

#' Get Plate Boundaries
#'
#' This function retrieves the boundaries for each plate in the `master_list` data.
#' It calculates the minimum and maximum run indices for each plate and creates a list of boundaries and labels for plotting.
#' @keywords internal
#' @param master_list A list containing project details and PCA scores.
#' @return A vector of unique boundaries for the plates.
get_plate_boundaries <- function(master_list) {
  boundary_data <- dplyr::bind_rows(master_list$pca$scores) %>%
    dplyr::select(sample_run_index, sample_plate_id) %>%
    dplyr::distinct()

  boundaries <- c(0.5)
  for (plate in unique(boundary_data$sample_plate_id)) {
    plate_data <- boundary_data %>% dplyr::filter(sample_plate_id == plate)
    boundaries <- c(
      boundaries,
      min(plate_data$sample_run_index) - 0.5,
      max(plate_data$sample_run_index) + 0.5
    )
  }
  unique(boundaries)
}

#' Get Plate Annotations
#'
#' This function retrieves the median run index for each plate in the `master_list` data.
#' It creates a tibble with sample data source, plate ID, run index, and placeholder values for PCA components and value.
#' @keywords internal
#' @param master_list A list containing project details and PCA scores.
#' @return A tibble containing plate annotations with median run indices.
get_plate_annotations <- function(master_list) {
  boundary_data <- dplyr::bind_rows(master_list$pca$scores) %>%
    dplyr::select(sample_run_index, sample_plate_id) %>%
    dplyr::distinct()

  coords <- sapply(unique(boundary_data$sample_plate_id), function(plate) {
    stats::median(boundary_data$sample_run_index[boundary_data$sample_plate_id == plate])
  })

  tibble::tibble(
    sample_data_source = "x.plateID",
    sample_plate_id = names(coords),
    sample_run_index = coords,
    PC1 = 0,
    PC2 = 0,
    PC3 = 0,
    value = 0
  )
}

#' Plot Control Chart
#'
#' This function generates a control chart for a specific metabolite in the `master_list`.
#' It combines data from peak area, SIL peak area, concentration, and stat target concentration,
#' and plots the values against the sample run index.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param metabolite The metabolite to plot in the control chart.
#' @param plate_boundaries A vector of plate boundaries for vertical lines in the plot.
#' @param annotate_label A tibble containing annotations for the plates.
#' @return A `ggplot` object representing the control chart for the specified metabolite.
plot_control_chart <- function(master_list,
                               metabolite,
                               plate_boundaries,
                               annotate_label) {
  sil_versions <- unique(unlist(master_list$templates[["Plate SIL version"]]))
  sil_notes <- unique(unlist(lapply(sil_versions, function(ver) {
    guide <- master_list$templates$mrm_guides[[ver]]$SIL_guide
    guide$Note[guide$`Precursor Name` == metabolite &
                 guide$control_chart == TRUE]
  })))

  data_combined <- dplyr::bind_rows(
    dplyr::bind_rows(master_list$data$peakArea$imputed) %>%
      dplyr::select(dplyr::contains("sample"), dplyr::any_of(metabolite)) %>%
      dplyr::mutate(sample_data_source = ".peakArea"),

    dplyr::bind_rows(master_list$data$peakArea$imputed) %>%
      dplyr::select(dplyr::contains("sample"), dplyr::any_of(sil_notes)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(SIL = paste(stats::na.omit(
        dplyr::c_across(dplyr::contains("SIL"))
      ), collapse = " ")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(SIL = as.numeric(SIL)) %>%
      dplyr::rename(!!metabolite := SIL) %>%
      dplyr::select(-dplyr::contains("SIL")) %>%
      dplyr::mutate(sample_data_source = ".SIL.peakArea"),

    dplyr::bind_rows(master_list$data$concentration$imputed) %>%
      dplyr::select(dplyr::contains("sample"), dplyr::any_of(metabolite)) %>%
      dplyr::mutate(sample_data_source = "concentration.preprocessed"),

    dplyr::bind_rows(master_list$data$concentration$statTargetProcessed) %>%
      dplyr::select(dplyr::contains("sample"), dplyr::any_of(metabolite)) %>%
      dplyr::mutate(sample_data_source = "statTargetConcentration.preprocessed")
  ) %>%
    dplyr::rename(value = !!metabolite)


  plotly::ggplotly(
    ggplot2::ggplot(
      data_combined,
      ggplot2::aes(
        x = sample_run_index,
        y = value,
        group = sample_name,
        fill = sample_type_factor,
        color = sample_type_factor,
        shape = sample_type_factor,
        size = sample_type_factor,
      )
    ) +
      ggplot2::geom_vline(xintercept = plate_boundaries, linetype = "dashed") +
      ggplot2::geom_point() +
      ggplot2::geom_text(
        data = annotate_label,
        ggplot2::aes(x = sample_run_index, y = value, label = sample_plate_id),
        inherit.aes = FALSE,
        color = "black"
      ) +
      ggplot2::theme_bw() +
      ggplot2::scale_shape_manual(values = master_list$project_details$plot_shape) +
      ggplot2::scale_fill_manual(values = master_list$project_details$plot_fill) +
      ggplot2::scale_color_manual(values = master_list$project_details$plot_colour) +
      ggplot2::scale_size_manual(values = master_list$project_details$plot_size) +
      ggplot2::ylab(metabolite) +
      ggplot2::guides(
        shape = "none",
        size = "none",
        color = "none",
        fill = ggplot2::guide_legend(title = "sample_type")
      ) +
      ggplot2::facet_wrap(
        facets = vars(sample_data_source),
        ncol = 1,
        scales = "free_y"
      ) +
      ggplot2::labs(
        title = paste0(
          metabolite,
          "; control chart; ",
          master_list$project_details$project_name
        )
      )
  )
}

#. ----

#Phase 5: Exports----
###Primary Function ----
#' Export All Project Outputs
#'
#' Exports the `master_list` to XLSX, HTML, and RDA formats, including summary tables, QC metrics, and processed data.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with exported files.
#' @examples
#' \dontrun{
#' master_list <- qcCheckR_export_all(master_list)
#' }
qcCheckR_export_all <- function(master_list) {
  master_list <- export_xlsx_file(master_list)
  master_list <- export_html_report(master_list)
  master_list <- export_master_list_rda(master_list)
  return(master_list)
}

### Sub Functions ----
#' Export XLSX File
#'
#' This function exports the `master_list` data to an XLSX file, including user
#' guide, QC metrics, and processed data.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with the XLSX file exported.
export_xlsx_file <- function(master_list) {
  master_list$summary_tables$odsAreaOverview <- create_user_guide(master_list)

  output_path <- file.path(
    master_list$project_details$project_dir,
    "all",
    "xlsx_report",
    paste0(
      Sys.Date(),
      "_",
      master_list$project_details$user_name,
      "_",
      master_list$project_details$project_name,
      "_",
      "_lipidData_qcCheckeR.xlsx"
    )
  )

  openxlsx::write.xlsx(
    x = list(
      "userGuide" = master_list$summary_tables$odsAreaOverview,
      "QC.platePerformance" = master_list$summary_tables$projectOverview,
      "QC.sampleMV" = master_list$filters$samples.missingValues,
      "QC.lipidsMV" = dplyr::bind_rows(master_list$filters$lipid.missingValues$allPlates),
      "QC.lipidQcRsd" = format_rsd_table(master_list),
      "DATA.peakArea" = dplyr::bind_rows(master_list$data$peakArea$sorted) %>%
        dplyr::select(-dplyr::contains("SIL"))%>%
        dplyr::mutate(dplyr::across(
          .cols = where(is.numeric) & !matches("sample", ignore.case = TRUE),
          .fns = ~ ifelse(is.na(.), NA,
                          ifelse(. < 1, signif(., 3), round(., 2)))
        )),
      "DATA.silPeakArea" = dplyr::bind_rows(master_list$data$peakArea$sorted) %>%
        dplyr::select(dplyr::contains("sample") | dplyr::contains("SIL"))%>%
        dplyr::mutate(dplyr::across(
          .cols = where(is.numeric) & !matches("sample", ignore.case = TRUE),
          .fns = ~ ifelse(is.na(.), NA,
                          ifelse(. < 1, signif(., 3), round(., 2)))
        )),
      "DATA.all.concentration" = dplyr::bind_rows(master_list$data$concentration$sorted)%>%
        dplyr::mutate(dplyr::across(
          .cols = where(is.numeric) & !matches("sample", ignore.case = TRUE),
          .fns = ~ ifelse(is.na(.), NA,
                          ifelse(. < 1, signif(., 3), round(., 2)))
        )),
      "DATA.preProcessed.concentration" = filter_concentration(master_list, "concentration")%>%
        dplyr::mutate(dplyr::across(
          .cols = where(is.numeric) & !matches("sample", ignore.case = TRUE),
          .fns = ~ ifelse(is.na(.), NA,
                          ifelse(. < 1, signif(., 3), round(., 2)))
        )),
      "DATA.all.concentration.S.T." = dplyr::bind_rows(master_list$data$concentration$statTargetProcessed)%>%
        dplyr::mutate(dplyr::across(
          .cols = where(is.numeric) & !matches("sample", ignore.case = TRUE),
          .fns = ~ ifelse(is.na(.), NA,
                          ifelse(. < 1, signif(., 3), round(., 2)))
        )),
      "DATA.preProcessed.conc.S.T." = filter_concentration(master_list, "concentration[statTarget]")%>%
      dplyr::mutate(dplyr::across(
        .cols = where(is.numeric) & !matches("sample", ignore.case = TRUE),
        .fns = ~ ifelse(is.na(.), NA,
                        ifelse(. < 1, signif(., 3), round(., 2)))
      ))),
    file = output_path,
    overwrite = TRUE
  )

  return(master_list)
}

#' Export HTML Report
#' This function exports the `master_list` data to an HTML report using a predefined R Markdown template.
#' It renders the report and opens it in the default web browser.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with the HTML report exported.
export_html_report <- function(master_list) {
  output_file <- file.path(
    master_list$project_details$project_dir,
    "all",
    "html_report",
    paste0(
      Sys.Date(),
      "_",
      master_list$project_details$user_name,
      "_",
      master_list$project_details$project_name,
      "_",
      "_lipidExploreR_qcCheckeR_report.html"
    )
  )


  control_chart_code <- character()

  for (i in names(master_list$control_charts)) {
    text <- paste0(
      "#### ",
      i,
      "\n",
      "```{r ",i,", echo=FALSE, message=FALSE, warning=FALSE, fig.width=14, fig.height=7}\n",
      "master_list$control_charts[['",
      i,
      "']]\n",
      "```\n\n"
    )
    control_chart_code <- c(control_chart_code, text)
  }


  template_path <- system.file("templates",
                               "qcCheckR_report_template_HS_V1.Rmd",
                               package = "MStargetR")
  template_content <- readLines(template_path)


  filled_template <- gsub(
    "control_charts_custom_code_placeholder",
    paste(control_chart_code, collapse = "\n"),
    paste(template_content, collapse = "\n")
  )

  temp_file <- tempfile(fileext = ".Rmd")
  writeLines(filled_template, temp_file)

  rmarkdown::render(
    input = temp_file,
    output_format = "html_document",
    output_dir = dirname(output_file),
    output_file = basename(output_file)
  )

  utils::browseURL(output_file)
  unlink(temp_file)

  return(master_list)
}

#' Export Master List as RDA File
#'
#' This function exports the `master_list` to an RDA file, which can be used for future analysis or sharing.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return The updated `master_list` with the RDA file exported.
export_master_list_rda <- function(master_list) {
  output_file <- file.path(
    master_list$project_details$project_dir,
    "all/data/rda",
    paste0(
      Sys.Date(),
      "_",
      master_list$project_details$project_name,
      "_",
      "_qcCheckR.rda"
    )
  )
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  callr::r_bg(function(master_list, output_file) {
    save(master_list, file = output_file, compress = FALSE)
  },args = list(master_list = master_list, output_file = output_file))

  master_list <- update_script_log(
    master_list,
    "data_exports",
    "plot_generation",
    "Script Complete \n\n\n Thank you for choosing MStargetR"
  )

  return(master_list)
}

#' Create User Guide
#'
#' This function creates a user guide for the `master_list` project, summarizing key project details and metrics.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return A tibble containing the user guide with key metrics and descriptions.
create_user_guide <- function(master_list) {
  #Create sample metrics
  sample_tags <- dplyr::bind_rows(master_list$data$peakArea$sorted) %>%
    dplyr::select(sample_type_factor) %>%
    unique() %>%
    dplyr::pull(sample_type_factor) %>%
    as.character() %>%
    .[!grepl("sample", .)]

  user_guide <- tibble::tibble(
    key = c(
      "projectName",
      "user",
      "qcType_preProcessing|filtering",
      "SIL_Int.Std_version",
      "total.StudyPlates",
      "total.Samples",
      "studySamples",
      paste0(sample_tags),
      "total.LipidFeatures",
      "total.MatchedLipidFeatures",
      "total.SIL.Int.Stds",
      "",
      "TAB_DESCRIPTION:",
      "QC.platePerformance",
      "QC.samplesMV",
      "QC.lipidsMV",
      "QC.lipidQcRsd",
      "DATA.lipidPeakArea",
      "DATA.silPeakArea",
      "DATA.all.concentration",
      "DATA.preProcessed.concentration",
      "DATA.all.concentration.S.T.",
      "DATA.preProcessed.conc.S.T."
    ),
    value = c(
      master_list$project_details$project_name,
      master_list$project_details$user_name,
      master_list$project_details$qc_type,
      paste(
        unique(master_list$templates$`Plate SIL version`),
        collapse = ","
      ),
      length(unique(
        dplyr::bind_rows(master_list$data$peakArea$sorted)$sample_plate_id
      )),
      nrow(dplyr::bind_rows(master_list$data$peakArea$sorted)),
      nrow(dplyr::filter(
        dplyr::bind_rows(master_list$data$peakArea$sorted),
        sample_type_factor == "sample"
      )),
      sapply(sample_tags, function(tag)
        nrow(
          dplyr::filter(
            dplyr::bind_rows(master_list$data$peakArea$sorted),
            sample_type_factor == tag
          )
        )),
      ncol(
        dplyr::bind_rows(master_list$data$peakArea$sorted) %>% dplyr::select(-dplyr::contains("sample"), -dplyr::contains("SIL"))
      ),
      ncol(
        dplyr::bind_rows(master_list$data$concentration$statTargetProcessed) %>% dplyr::select(-dplyr::contains("sample"))
      ),
      ncol(
        dplyr::bind_rows(master_list$data$peakArea$sorted) %>% dplyr::select(dplyr::contains("SIL"))
      ),
      "",
      "",
      "overview of project quality performance (per plate)",
      "detailed overview of sample quality (missing values)",
      "detailed overview of lipid quality (missing values)",
      "detailed overview of lipid quality (% RSD in QC samples)",
      "lipid target peak area integrals (PeakForgeR)",
      "stable isotope labelled internal standard peak area integrals (PeakForgeR)",
      "peakArea >> SIL ratio >> concentration factor adjusted",
      "peakArea >> imputed >> SIL ratio >> concentration factor adjusted >> filtered",
      "peakArea >> imputed >> SIL ratio >> concentration factor adjusted >> statTarget correction",
      "same as above, but filtered"
    )
  )
  return(user_guide)
}

#' Format RSD Table
#'
#' This function formats the RSD table from the `master_list` filters.
#' It rounds the RSD values, adds a data column, and transposes the table for better readability.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @return A tibble containing the formatted RSD table with rounded values and transposed structure.
format_rsd_table <- function(master_list) {
  master_list$filters$rsd %>%
    dplyr::mutate(dplyr::across(!dplyr::contains("data"), round, 2)) %>%
    tibble::add_column(data = paste0(.$dataSource, ".", .$dataBatch),
               .before = 1) %>%
    dplyr::select(-dataSource, -dataBatch) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    stats::setNames(.[1, ]) %>%
    dplyr::filter(data != "data") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(dplyr::across(!dplyr::contains("data"), as.numeric))
}

#' Filter Concentration Data
#'
#' This function filters the concentration data from the `master_list` based on the specified source.
#' It removes failed samples and lipids, and applies RSD filters based on the specified source.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param source The data source to filter (e.g., "concentration").
#' @return A tibble containing the filtered concentration data.
filter_concentration <- function(master_list, source) {
  data <- dplyr::bind_rows(master_list$data$concentration[[if (source == "concentration")
    "imputed"
    else
      "statTargetProcessed"]])
  rsd_filter <- master_list$filters$rsd %>%
    dplyr::filter(dataSource == source, dataBatch == "allBatches") %>%
    dplyr::select(-dplyr::contains("data")) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ ifelse(. > 30, TRUE, FALSE))) %>%
    dplyr::select(tidyselect::where(~ any(., na.rm = TRUE))) %>%
    names()

  data %>%
    dplyr::filter(!sample_name %in% master_list$filters$failed_samples) %>%
    dplyr::select(-dplyr::any_of(master_list$filters$failed_lipids)) %>%
    dplyr::select(-dplyr::any_of(rsd_filter))
}

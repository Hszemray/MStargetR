#' Utils_Global.R

#' Import specific functions from packages
#' @keywords internal
#' @name PeakForgeR_import_external_functions
#' @importFrom utils installed.packages sessionInfo data str flush.console setTxtProgressBar txtProgressBar
#' @importFrom readr read_csv write_csv read_tsv write_tsv
#' @importFrom dplyr bind_rows bind_cols filter select setequal rename arrange contains intersect
#' @importFrom tibble tibble add_column
#' @importFrom stringr str_extract str_sub
#' @importFrom mzR openMSfile chromatogramHeader chromatograms
#' @importFrom janitor clean_names
#' @importFrom magrittr %>%
#' @importFrom stats median setNames
NULL

#PeakForgeR internal functions----

#.----
#Setup Project Functions----

###Primary Function----
#' PeakForgeR_setup_project
#'
#' This function sets up the project by initialising the master list,
#' setting up project directories, and updating the script log.
#' @keywords internal
#' @param user_name string identifying user
#' @param project_directory Directory path for the project folder containing
#' the wiff folder and wiff files.
#' @param plateID Plate ID for the current plate.
#' @param mrm_template_list List of MRM guides.
#' @param QC_sample_label Key for filtering QC samples from sample list.
#' @return The updated `master_list` object with the project setup details.
#' @examples
#' \dontrun{
#' PeakForgeR_setup_project("path/to/project_directory", "plateID",
#'mrm_template_list, "QC_sample_label")
#' }
PeakForgeR_setup_project <- function(user_name,
                                  project_directory,
                                  plateID,
                                  mrm_template_list,
                                  QC_sample_label) {
  master_list <- initialise_master_list()
  master_list <- store_environment_details(master_list)
  master_list <- set_project_details(master_list, user_name, project_directory, plateID, QC_sample_label)
  master_list <- read_mrm_guides(master_list, mrm_template_list)
  setup_project_directories(master_list)
  master_list <- update_script_log(master_list, "project_setup", "start_time", "mzR_mzml_import")
  return(master_list)
}

###Sub Functions----
#' initialise_master_list
#'
#' This function initialises the master list with default values.
#' @keywords internal
#' @return A list representing the initialised master list.
#' @examples
#' \dontrun{
#' master_list <- initialise_master_list()
#' }
initialise_master_list <- function() {
  master_list <- list()
  master_list$environment <- list()
  master_list$environment$user_functions <- list()
  master_list$templates <- list()
  master_list$templates$mrm_guides <- list()
  master_list$project_details <- list()
  master_list$data <- list()
  master_list$summary_tables <- list()
  master_list$process_lists <- list()
  return(master_list)
}

#' store_environment_details
#'
#' This function stores environment details in the master list.
#' @keywords internal
#' @param master_list The master list object.
#' @return The updated master list object with environment details.
#' @examples
#' \dontrun{
#' master_list <- store_environment_details(master_list)
#' }
store_environment_details <- function(master_list) {

  master_list$environment$r_version <- utils::sessionInfo()$R.version$version.string
  master_list$environment$base_packages <- utils::sessionInfo()$basePkgs

  installed <- utils::installed.packages()
  attached <- names(utils::sessionInfo()$otherPkgs)
  valid_pkgs <- attached[attached %in% rownames(installed)]
  versions <- installed[valid_pkgs, "Version"]
  master_list$environment$user_packages <- paste0(valid_pkgs, ": ", versions)

  return(master_list)
}

#' set_project_details
#'
#' This function sets project details in the master list.
#' @keywords internal
#' @param master_list The master list object.
#' @param user_name string specifying user
#' @param project_directory Directory path for the project folder.
#' @param plateID Plate ID for the current plate.
#' @param QC_sample_label Key for filtering QC samples from sample list.
#' @return The updated master list object with project details.
#' @examples
#' \dontrun{
#' master_list <- set_project_details(master_list,
#'                                    "John Smith"
#'                                    "path/to/project_directory",
#'                                    "plateID",
#'                                    "QC_sample_label")
#' }
set_project_details <- function(master_list,
                                user_name,
                                project_directory,
                                plateID,
                                QC_sample_label) {
  master_list$project_details$project_dir <- project_directory
  master_list$project_details$lipidExploreR_version <- "Automated"
  master_list$project_details$user_name <- user_name
  master_list$project_details$project_name <- stringr::str_extract(master_list$project_details$project_dir, "[^/]*$")
  master_list$project_details$plateID <- plateID
  master_list$project_details$qc_type <- QC_sample_label
  master_list$project_details$script_log$timestamps$start_time <- Sys.time()
  return(master_list)
}

#' read_mrm_guides
#'
#' This function reads MRM guides from user-supplied paths in the mrm_template_list.
#' Capable of reading .tsv or .csv files
#' @keywords internal
#' @param master_list The master list object.
#' @param mrm_template_list List of MRM guide file paths.
#' @return The updated master list object with MRM guides.
#' @examples
#' \dontrun{
#' master_list <- read_mrm_guides(master_list, mrm_template_list)
#' }

read_mrm_guides <- function(master_list, mrm_template_list) {
  names(mrm_template_list) <- sapply(mrm_template_list, basename)

  for (version in names(mrm_template_list)) {
    file_path <- mrm_template_list[[version]]

    # Determine file extension and read accordingly
    if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
      guide_data <- readr::read_csv(file_path,
                                    show_col_types = FALSE,
                                    name_repair = "minimal")
    } else if (grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
      guide_data <- readr::read_tsv(file_path,
                                    show_col_types = FALSE ,
                                    name_repair = "minimal")
    } else {
      stop(
        paste(
          "Unsupported file format for:",
          file_path,
          ". \n Please ensure mrm_templates are .csv or .tsv"
        )
      )
    }

    # validate mrm_template_list
    ##Check column  names and length
    mandatory_cols <- c(
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
    guide_cols <- colnames(guide_data)
    matching_cols  <- dplyr::intersect(mandatory_cols, guide_cols)

    # report if all cols are matching
    if (length(matching_cols) != length(mandatory_cols)) {
      missing_cols <- setdiff(mandatory_cols, guide_cols)
      stop(paste(
        version,
        ": Missing mandatory columns: ",
        paste(missing_cols, collapse = ", ")
      ))
    }

    ## Check there are no NA in all columns
    ## Except "Note" if "Precursor Name" contains "SIL" then expect an NA
    sil_rows <- grepl("SIL", guide_data[["Precursor Name"]], ignore.case = TRUE)
    cols_to_check <- setdiff(mandatory_cols, "Note")
    for (col in cols_to_check) {
      if (any(is.na(guide_data[[col]]))) {
        stop(paste(
          version,
          ": Column",
          col,
          "contains NA values, which are not allowed."
        ))
      }
    }
    note_na <- is.na(guide_data[["Note"]])
    invalid_note_na <- note_na & !sil_rows
    if (any(invalid_note_na)) {
      stop(version,
          ": NA values in 'Note' column
        are only allowed for rows where 'Precursor Name' contains 'SIL'."
      )
    }

    full_message <- paste(version,
                          ": All required columns are validated and contain no unexpected NA values.")
    message(full_message)

    guide_data <- replace_precursor_symbols(guide_data, columns = c("Precursor Name", "Note"))

    transition_check_result <- transition_checkR(guide_data)
    if (!is.null(transition_check_result)) {
      stop(paste(version, ": Non-unique transitions detected. Please amend mrm_template prior rerunning"))
      print(transition_check_result)
    }

    #Store valdiated mrm_template in master_list
    master_list$templates$mrm_guides[[version]]$mrm_guide <- guide_data
  }

  return(master_list)
}


#' setup_project_directories
#'
#' This function sets up project directories for each plate ID.
#' @keywords internal
#' @param master_list The master list object.
#' @return None. The function sets up directories.
#' @examples
#' \dontrun{
#' setup_project_directories(master_list)
#' }
setup_project_directories <- function(master_list) {
  for (plate_ID in master_list$project_details$plateID) {
    base_path <- file.path(master_list$project_details$project_dir, plate_ID)
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

#.----

#mzR Mrm FindR Functions-----

###Primary Function----
#' mzR_mrm_findR
#'
#' This function processes mzML files to find peak apex and boundaries using QC mzR data, updates the mrm guide, and provides peak boundary information.
#' @keywords internal
#' @param FUNC_mzR List from master_list containing mzR object for each sample, mzR_header, mzR_chromatogram parsed internally.
#' @param FUNC_mrm_guide Tibble of mrm details parsed internally. See run mrm_template_guide for example.
#' @param FUNC_OPTION_qc_type QC type used in the experiment passed internally.
#' @return A list containing updated mrm guide and peak boundary information.
#' @examples
#' \dontrun{
#' mzR_mrm_findR(FUNC_mzR, FUNC_mrm_guide, FUNC_OPTION_qc_type)
#' }
mzR_mrm_findR <- function(FUNC_mzR,
                          FUNC_mrm_guide,
                          FUNC_OPTION_qc_type) {
  validate_mzR_parameters(FUNC_mzR, FUNC_mrm_guide, FUNC_OPTION_qc_type)
  mzML_filelist <- get_mzML_filelist(FUNC_mzR)
  mzML_filelist_qc <- filter_mzML_filelist_qc(mzML_filelist, FUNC_OPTION_qc_type)
  FUNC_tibble <- process_files(FUNC_mzR, FUNC_mrm_guide, mzML_filelist_qc)
  FUNC_output <- create_output(FUNC_tibble, FUNC_mrm_guide, mzML_filelist)
  return(FUNC_output)
}

###Sub Functions----

#' validate_mzR_parameters
#'
#' This function validates parameters mzR_mrmfindR. If parameters fail validation script stops
#' @keywords internal
#' @param FUNC_mzR List from master_list containing mzR object for each sample, mzR_header, mzR_chromatogram parsed internally.
#' @param FUNC_mrm_guide Tibble of mrm details parsed internally. See run mrm_template_guide for example.
#' @param FUNC_OPTION_qc_type QC type used in the experiment (LTR, PQC, none) parsed internally.
#' @return NULL
#' @examples
#' \dontrun{
#' validate_mzR_parameters(FUNC_mzR, FUNC_mrm_guide, FUNC_OPTION_qc_type)
#' }
validate_mzR_parameters <- function(FUNC_mzR,
                                    FUNC_mrm_guide,
                                    FUNC_OPTION_qc_type) {
  if (!is.list(FUNC_mzR) ||
      !is.data.frame(FUNC_mrm_guide) ||
      !is.character(FUNC_OPTION_qc_type)) {
    stop("Invalid input parameters.")
  }
}

#' get_mzML_filelist
#'
#' This function generates the mzML_filelist from FUNC_mzR
#' @keywords internal
#' @param FUNC_mzR List from master_list containing mzR object for each sample, mzR_header, mzR_chromatogram parsed internally.
#' @return A list containing mzML file information from plate.
#' @examples
#' \dontrun{
#' get_mzML_filelist(FUNC_mzR)
#' }
get_mzML_filelist <- function(FUNC_mzR) {
  mzML_filelist <- NULL
  for (idx_plate in names(FUNC_mzR)) {
    mzML_filelist <- c(mzML_filelist, names(FUNC_mzR[[idx_plate]]))
  }
  return(mzML_filelist)
}

#' filter_mzML_filelist_qc
#'
#' This function filters the mzML_filelist for user selected quality control samples
#' @keywords internal
#' @param mzML_filelist A list containing mzML file information from plate.
#' @param FUNC_OPTION_qc_type QC type used in the experiment parsed internally.
#' @return A list containing quality control sample data.
#' @examples
#' \dontrun{
#' filter_mzML_filelist_qc(mzML_filelist,FUNC_OPTION_qc_type)
#' }
filter_mzML_filelist_qc <- function(mzML_filelist, FUNC_OPTION_qc_type) {
  mzML_filelist_qc <- mzML_filelist[grep(FUNC_OPTION_qc_type, mzML_filelist)]
  return(mzML_filelist_qc)
}

#' process_files
#'
#' This function process each mzml to produce a tibble of mrm data for all samples on plate.
#' @keywords internal
#' @param FUNC_mzR List from master_list containing mzR object for each sample, mzR_header, mzR_chromatogram parsed internally.
#' @param FUNC_mrm_guide Tibble of mrm details parsed internally. See run mrm_template_guide for example.
#' @param mzML_filelist_qc A list containing mzML file information from plate for only quality control samples.
#' @return A tibble containing mrm data for all samples on the plate.
#' @examples
#' \dontrun{
#' process_files(FUNC_mzR, FUNC_mrm_guide, mzML_filelist_qc)
#' }
process_files <- function(FUNC_mzR,
                          FUNC_mrm_guide,
                          mzML_filelist_qc) {
  FUNC_tibble <- list()

  for (idx_mzML in mzML_filelist_qc) {
    FUNC_tibble[[idx_mzML]] <- tibble::tibble()

    for (idx_plate in names(FUNC_mzR)) {
      if (length(grep(idx_mzML, names(FUNC_mzR[[idx_plate]]))) == 1) {
        result <- tryCatch({
          process_mrm_transitions(FUNC_mzR, FUNC_mrm_guide, idx_plate, idx_mzML)
        }, error = function(e) {
          message(sprintf("Removing %s from optimisation due to lack of datapoints during aquisition: \n Error:%s ",
                          idx_mzML, e$message))
          return(NULL)
        })

        if (!is.null(result)) {
          FUNC_tibble[[idx_mzML]] <- result
        }
      }
    }
  }

  FUNC_tibble <- dplyr::bind_rows(FUNC_tibble) %>%
    dplyr::filter(lipid_class != "no match") %>%
    dplyr::filter(lipid_class != "multiple match")

  return(FUNC_tibble)
}



#' process_mrm_transitions
#'
#' This function processes mrm_transition for each sample on a plate.
#' @keywords internal
#' @param FUNC_mzR List from master_list containing mzR object for each sample, mzR_header, mzR_chromatogram parsed internally.
#' @param FUNC_mrm_guide Tibble of mrm details parsed internally. See run mrm_template_guide for example.
#' @param idx_plate A string naming the current plate ID being processed. Passed from process_files function.
#' @param idx_mzML A string naming the current mzML file being processed. Passed from process_files function.
#' @return A tibble of mrm data for the specific idx_plate/idx_mzML.
#' @examples
#' \dontrun{
#' process_mrm_transitions(FUNC_mzR, FUNC_mrm_guide, idx_plate, idx_mzML)
#' }
process_mrm_transitions <- function(FUNC_mzR,
                                    FUNC_mrm_guide,
                                    idx_plate,
                                    idx_mzML) {
  FUNC_tibble <- tibble::tibble()
  for (idx_mrm in 3:nrow(FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_header)) {
    if (nrow(FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]]) > 0) {
      precursor_mz <- FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_header$precursorIsolationWindowTargetMZ[idx_mrm]
      product_mz <- FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_header$productIsolationWindowTargetMZ[idx_mrm]
      baseline_value <- calculate_baseline(FUNC_mzR, idx_plate, idx_mzML, idx_mrm)
      peak_apex_idx <- find_peak_apex_idx(FUNC_mzR, idx_plate, idx_mzML, idx_mrm)
      peak_start_idx <- find_peak_start_idx(FUNC_mzR,
                                            idx_plate,
                                            idx_mzML,
                                            idx_mrm,
                                            peak_apex_idx,
                                            baseline_value)
      peak_end_idx <- find_peak_end_idx(FUNC_mzR,
                                        idx_plate,
                                        idx_mzML,
                                        idx_mrm,
                                        peak_apex_idx,
                                        baseline_value)
      mzml_rt_apex <- FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]]$rtime[peak_apex_idx] %>% round(2)
      mzml_rt_start <- FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]]$rtime[peak_start_idx] %>% round(2)
      mzml_rt_end <- FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]]$rtime[peak_end_idx] %>% round(2)
      lipid_info <- find_lipid_info(
        FUNC_mrm_guide,
        precursor_mz,
        product_mz,
        mzml_rt_apex,
        FUNC_mzR,
        idx_plate,
        idx_mzML,
        idx_mrm
      )
      FUNC_tibble <- FUNC_tibble %>% dplyr::bind_rows(
        dplyr::bind_cols(
          "mzml" = idx_mzML,
          "lipid_class" = lipid_info$class,
          "lipid" = lipid_info$name,
          "precursor_mz" = precursor_mz,
          "product_mz" = product_mz,
          "peak_apex" = mzml_rt_apex,
          "peak_start" = mzml_rt_start,
          "peak_end" = mzml_rt_end
        )
      )
    }
  }
  return(FUNC_tibble)
}

#' calculate_baseline
#'
#' This function calculates baseline for a single transition.
#' @keywords internal
#' @param FUNC_mzR List from master_list containing mzR object for each sample, mzR_header, mzR_chromatogram parsed internally.
#' @param idx_plate A string naming the current plate ID being processed. Passed from process_files function.
#' @param idx_mzML A string naming the current mzML file being processed. Passed from process_files function.
#' @param idx_mrm A string identifying the current mrm. Passed from process_files function.
#' @return A numeric value for baseline.
#' @examples
#' \dontrun{
#' calculate_baseline(FUNC_mzR, idx_plate, idx_mzML, idx_mrm)
#' }
calculate_baseline <- function(FUNC_mzR, idx_plate, idx_mzML, idx_mrm) {
  baseline_value <- FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]][, 2] %>% stats::median()
  if (baseline_value < 1) {
    baseline_value <- FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]][, 2] %>% base::mean()
  }
  return(baseline_value)
}

#' find_peak_apex_idx
#'
#' This function finds peak apex for a single transition.
#' @keywords internal
#' @param FUNC_mzR List from master_list containing mzR object for each sample, mzR_header, mzR_chromatogram parsed internally.
#' @param idx_plate A string naming the current plate ID being processed. Passed from process_files function.
#' @param idx_mzML A string naming the current mzML file being processed. Passed from process_files function.
#' @param idx_mrm A string identifying the current mrm. Passed from process_files function.
#' @return A numeric value for peak apex index.
#' @examples
#' \dontrun{
#' find_peak_apex_idx(FUNC_mzR, idx_plate, idx_mzML, idx_mrm)
#' }
find_peak_apex_idx <- function(FUNC_mzR, idx_plate, idx_mzML, idx_mrm) {
  peak_apex_idx <- which.max(FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]][, 2])
  if (peak_apex_idx < 5 ||
      peak_apex_idx > (length(FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]][, 2]) -
                       5)) {
    peak_apex_idx <- which.max((FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]][, 2])[5:(length(FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]][, 2]) -
                                                                                                         5)]) + 4
  }
  return(peak_apex_idx)
}

#' find_peak_start_idx
#'
#' This function finds peak start for a single transition.
#' @keywords internal
#' @param FUNC_mzR List from master_list containing mzR object for each sample, mzR_header, mzR_chromatogram parsed internally.
#' @param idx_plate A string naming the current plate ID being processed. Passed from process_files function.
#' @param idx_mzML A string naming the current mzML file being processed. Passed from process_files function.
#' @param idx_mrm A string identifying the current mrm. Passed from process_files function.
#' @param peak_apex_idx Index of peak apex passed from process_files function.
#' @param baseline_value Numeric value indicating transitions baseline. Passed from process_files function.
#' @return A numeric value for peak start index.
#' @examples
#' \dontrun{
#' find_peak_start_idx(FUNC_mzR, idx_plate, idx_mzML, idx_mrm, peak_apex_idx, baseline_value)
#' }
find_peak_start_idx <- function(FUNC_mzR,
                                idx_plate,
                                idx_mzML,
                                idx_mrm,
                                peak_apex_idx,
                                baseline_value) {
  baseline_idx <- which((FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]][, 2]) < baseline_value)
  peak_start_idx <- baseline_idx[which(baseline_idx < peak_apex_idx)][(length(which(baseline_idx < peak_apex_idx))) -
                                                                        3]
  if (length(peak_start_idx) == 0 || peak_start_idx < 1) {
    peak_start_idx <- 1
  }
  return(peak_start_idx)
}

#' find_peak_end_idx
#'
#' This function finds peak end for a single transition.
#' @keywords internal
#' @param FUNC_mzR List from master_list containing mzR object for each sample, mzR_header, mzR_chromatogram parsed internally.
#' @param idx_plate A string naming the current plate ID being processed. Passed from process_files function.
#' @param idx_mzML A string naming the current mzML file being processed. Passed from process_files function.
#' @param idx_mrm A string identifying the current mrm. Passed from process_files function.
#' @param peak_apex_idx Index of peak apex passed from process_files function.
#' @param baseline_value Numeric value indicating transitions baseline. Passed from process_files function.
#' @return A numeric value for peak end index.
#' @examples
#' \dontrun{
#' find_peak_end_idx(FUNC_mzR, idx_plate, idx_mzML, idx_mrm, peak_apex_idx, baseline_value)
#' }
find_peak_end_idx <- function(FUNC_mzR,
                              idx_plate,
                              idx_mzML,
                              idx_mrm,
                              peak_apex_idx,
                              baseline_value) {
  baseline_idx <- which((FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]][, 2]) < baseline_value)
  peak_end_idx <- baseline_idx[which(baseline_idx > peak_apex_idx)][3]
  if (length(peak_end_idx) == 0 ||
      peak_end_idx > length(FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]][, 2]) ||
      is.na(peak_end_idx)) {
    peak_end_idx <- length(FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]][, 2])
  }
  return(peak_end_idx)
}


#' find_lipid_info
#'
#' This function finds and matches lipids to mrm data.
#' @keywords internal
#' @param FUNC_mrm_guide Tibble of mrm details parsed internally. See run mrm_template_guide for example.
#' @param precursor_mz A numeric value for precursor mass to charge ratio.
#' @param product_mz A numeric value for product mass to charge ratio.
#' @param mzml_rt_apex A numeric value for retention time apex.
#' @param FUNC_mzR List from master_list containing mzR object for each sample, mzR_header, mzR_chromatogram parsed internally.
#' @param idx_plate A string naming the current plate ID being processed. Passed from process_files function.
#' @param idx_mzML A string naming the current mzML file being processed. Passed from process_files function.
#' @param idx_mrm A string identifying the current mrm. Passed from process_files function.
#' @return A list contain lipid class and lipid species information.
#' @examples
#' \dontrun{
#' find_lipid_info(FUNC_mrm_guide, precursor_mz, product_mz,
#'                 mzml_rt_apex, FUNC_mzR, idx_plate, idx_mzML, idx_mrm)
#' }
find_lipid_info <- function(FUNC_mrm_guide,
                            precursor_mz,
                            product_mz,
                            mzml_rt_apex,
                            FUNC_mzR,
                            idx_plate,
                            idx_mzML,
                            idx_mrm) {
  lipid_idx <- which(
    FUNC_mrm_guide$precursor_mz == precursor_mz &
      FUNC_mrm_guide$product_mz == product_mz
  )
  if (length(lipid_idx) != 1) {
    lipid_idx <- which(
      FUNC_mrm_guide$precursor_mz > (precursor_mz - 0.25) &
        FUNC_mrm_guide$precursor_mz < (precursor_mz + 0.25) &
        FUNC_mrm_guide$product_mz > (product_mz - 0.25) &
        FUNC_mrm_guide$product_mz < (product_mz + 0.25) &
        FUNC_mrm_guide$explicit_retention_time > (min(FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]]$rtime) -
                                                    0.1) &
        FUNC_mrm_guide$explicit_retention_time < (max(FUNC_mzR[[idx_plate]][[idx_mzML]]$mzR_chromatogram[[idx_mrm]]$rtime) +
                                                    0.1)
    )
  }
  if (length(lipid_idx) > 1) {
    lipid_class <- "multiple match"
    lipid <- "multiple match"
  } else if (length(lipid_idx) == 0) {
    lipid_class <- "no match"
    lipid <- "no match"
  } else if (length(lipid_idx == 1)) {
    lipid_class <- FUNC_mrm_guide$molecule_list_name[lipid_idx]
    lipid <- FUNC_mrm_guide$precursor_name[lipid_idx]
  }
  return(list(class = lipid_class, name = lipid))
}


#' create_output
#'
#' This function updated mrm_guide and peak_boundary_update for later use in PeakForgeR
#' @keywords internal
#' @param FUNC_tibble A tibble of mrm details per sample per lipid species.
#' @param FUNC_mrm_guide Tibble of mrm details parsed internally. See run mrm_template_guide for example.
#' @param mzML_filelist A list containing mzML file information from plate.
#' @return A list containing two tibbles for mrm_guide_updated and peak_boundary_update.
#' @examples
#' \dontrun{
#' create_output(FUNC_tibble, FUNC_mrm_guide, mzML_filelist)
#' }
create_output <- function(FUNC_tibble,
                          FUNC_mrm_guide,
                          mzML_filelist) {
  if (nrow(FUNC_tibble) == 0) {
    return(
      list(
        mrm_guide_updated = tibble::tibble(),
        peak_boundary_update = tibble::tibble()
      )
    )
  }

  FUNC_output <- list()
  FUNC_output$mrm_guide_updated <- tibble::tibble()
  FUNC_output$peak_boundary_update <- tibble::tibble()
  for (idx_lipid in unique(FUNC_tibble$lipid)) {
    FUNC_output$mrm_guide_updated <- dplyr::bind_rows(
      FUNC_output$mrm_guide_updated,
      dplyr::bind_cols(
        "Molecule List Name" = (FUNC_tibble %>% dplyr::filter(lipid == idx_lipid))[["lipid_class"]] %>% unique(),
        "Precursor Name" = idx_lipid,
        "Precursor Mz" = (
          FUNC_mrm_guide %>% dplyr::filter(precursor_name == idx_lipid)
        )[["precursor_mz"]],
        "Precursor Charge" = (
          FUNC_mrm_guide %>% dplyr::filter(precursor_name == idx_lipid)
        )[["precursor_charge"]],
        "Product Mz" = (
          FUNC_mrm_guide %>% dplyr::filter(precursor_name == idx_lipid)
        )[["product_mz"]],
        "Product Charge" = (
          FUNC_mrm_guide %>% dplyr::filter(precursor_name == idx_lipid)
        )[["product_charge"]],
        "Explicit Retention Time" = (FUNC_tibble %>% dplyr::filter(lipid == idx_lipid))[["peak_apex"]] %>% stats::median(),
        "Explicit Retention Time Window" = (
          FUNC_mrm_guide %>% dplyr::filter(precursor_name == idx_lipid)
        )[["explicit_retention_time_window"]],
        "Note" = (
          FUNC_mrm_guide %>% dplyr::filter(precursor_name == idx_lipid)
        )[["note"]]
      )
    )

    FUNC_output$peak_boundary_update <- dplyr::bind_rows(
      FUNC_output$peak_boundary_update,
      dplyr::bind_cols(
        "FileName" = mzML_filelist,
        "FullPeptideName" = rep(idx_lipid, length(mzML_filelist)),
        "MinStartTime" = rep(((FUNC_tibble %>% dplyr::filter(lipid == idx_lipid))[["peak_start"]] %>% summary()
        )[["1st Qu."]], length(mzML_filelist)),
        "MaxEndTime" = rep(((FUNC_tibble %>% dplyr::filter(lipid == idx_lipid))[["peak_end"]] %>% summary()
        )[["3rd Qu."]], length(mzML_filelist))
      )
    )
  }
  FUNC_output$mrm_guide_updated <- FUNC_output$mrm_guide_updated %>% dplyr::arrange(`Precursor Name`)
  FUNC_output$peak_boundary_update <- FUNC_output$peak_boundary_update %>% dplyr::arrange(`FullPeptideName`)
  return(FUNC_output)
}


#.----
#Import mzML Files Functions ----

###Primary Function----
#' import_mzml
#'
#' This function imports mzML files for each plate using the mzR package, extracts relevant information, and updates the script log.
#' @keywords internal
#' @param plateID Plate ID for the current plate.
#' @param master_list Master list generated internally.
#' @return The updated `master_list` object with the mzML import details.
#' @examples
#' \dontrun{
#' import_mzml("plateID", master_list)
#' }
import_mzml <- function(plateID, master_list) {
  validate_master_list_project_directory(master_list)
  mzml_filelist <- initialise_mzml_filelist(master_list)
  master_list <- process_plates(master_list, mzml_filelist)
  master_list <- update_script_log(
    master_list,
    "mzR_mzml_import",
    "project_setup",
    "peak_picking_and_integration"
  )
  return(master_list)
}

###Sub Functions----

#' initialise_mzml_filelist
#'
#' This function initializes a list of mzML files for each plate in the master list, excluding files with specific patterns.
#' @keywords internal
#' @param master_list A list containing project details, including plate IDs and project directory.
#' @return A list of mzML files for each plate, excluding files with "COND", "Blank", or "ISTDs" in their names.
#' @examples
#' \dontrun{
#' initialise_mzml_filelist(master_list)
#' }
initialise_mzml_filelist <- function(master_list) {
  mzml_filelist <- list()

  for (idx_plate in master_list$project_details$plateID) {
    mzml_filelist[[idx_plate]] <- list.files(
      file.path(
        master_list$project_details$project_dir,
        idx_plate,
        "data",
        "mzml"
      ),
      pattern = ".mzML",
      full.names = FALSE
    )
    mzml_filelist[[idx_plate]] <- mzml_filelist[[idx_plate]][
      !grepl("(?i)\\bCOND\\b|\\bBLANK\\b|\\bISTDs\\b", mzml_filelist[[idx_plate]], perl = TRUE)
    ]  }
  return(mzml_filelist)
}


#' process_plates
#'
#' This function processes mzML files for each plate in the master list, extracting relevant information and updating the master list.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param mzml_filelist A list of mzML files for each plate.
#' @return The updated `master_list` object with processed mzML data.
#' @examples
#' \dontrun{
#' process_plates(master_list, mzml_filelist)
#' }
process_plates <- function(master_list, mzml_filelist) {
  for (idx_plate in master_list$project_details$plateID) {
    master_list$data[[idx_plate]]$mzR <- list()
    long_path <- file.path(master_list$project_details$project_dir,
                           idx_plate,
                           "data",
                           "mzml")
    if (length(long_path) > 260 & .Platform$OS.type == "windows") {
      cmd <- paste0('cmd /c mklink /J "C:\\mzml_short" "', long_path, '"')
      system(cmd,intern = FALSE, ignore.stdout = TRUE)
      mzml_path <- "C:/mzml_short"
    } else {
      mzml_path <- file.path(master_list$project_details$project_dir,
                             idx_plate,
                             "data",
                             "mzml")
    }

    for (idx_mzML in mzml_filelist[[idx_plate]]) {
      master_list$data[[idx_plate]]$mzR[[idx_mzML]] <- list()
      master_list$data[[idx_plate]]$mzR[[idx_mzML]]$mzR_object <- mzR::openMSfile(filename = file.path(mzml_path, idx_mzML))
      master_list$data[[idx_plate]]$mzR[[idx_mzML]]$mzR_header <- mzR::chromatogramHeader(master_list$data[[idx_plate]]$mzR[[idx_mzML]]$mzR_object)
      master_list$data[[idx_plate]]$mzR[[idx_mzML]]$mzR_chromatogram <- mzR::chromatograms(master_list$data[[idx_plate]]$mzR[[idx_mzML]]$mzR_object)
      master_list$data[[idx_plate]]$mzR[[idx_mzML]]$mzR_timestamp <- master_list$data[[idx_plate]]$mzR[[idx_mzML]]$mzR_object@backend$getRunStartTimeStamp()
    }
    master_list$data$global_timestamp[[idx_plate]] <- extract_timestamp(master_list$data[[idx_plate]]$mzR[[idx_mzML]]$mzR_timestamp)
    master_list$project_details$mzml_sample_list[[idx_plate]] <- update_sample_list(master_list, idx_plate)

    if (length(long_path) > 260 & .Platform$OS.type == "windows") {
      unlink("C:/mzml_short", recursive = TRUE)
    }
  }
  return(master_list)
}


#' extract_timestamp
#'
#' This function extracts the year from a timestamp string and converts it to a numeric value.
#' @keywords internal
#' @param timestamp A string representing the timestamp.
#' @return A numeric value representing the year extracted from the timestamp.
#' @examples
#' \dontrun{
#' extract_timestamp("2025-06-17T13:31:58Z")
#' }
extract_timestamp <- function(timestamp) {
  return(timestamp %>% stringr::str_sub(., 1, 4) %>% as.numeric())
}


#' update_sample_list
#'
#' This function updates the sample list for a given plate in the master list.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param idx_plate The index of the plate to update the sample list for.
#' @return The updated sample list for the specified plate.
#' @examples
#' \dontrun{
#' update_sample_list(master_list, idx_plate)
#' }
update_sample_list <- function(master_list, idx_plate) {
  if (is.null(master_list$project_details$mzml_sample_list[[idx_plate]])) {
    master_list$project_details$mzml_sample_list[[idx_plate]] <- character()
  }
  return(c(
    master_list$project_details$mzml_sample_list[[idx_plate]],
    names(master_list$data[[idx_plate]]$mzR)
  ))
}

#.----
#Peak Picking Functions ----

###Primary Function----
#' peak_picking
#'
#' This function processes mzML files for each plate, optimizes retention times, updates peak boundaries, and checks for SIL internal standards.
#' @keywords internal
#' @param plateID Plate ID for the current plate.
#' @param master_list Master list generated internally.
#' @return The updated `master_list` object with peak picking details.
#' @examples
#' \dontrun{
#' peak_picking(plateID, master_list)
#' }
peak_picking <- function(plateID, master_list) {
  validate_master_list_project_directory(master_list)
  plate_idx <- plateID
  message(paste("Processing plate:", plateID))

  sil_found <- FALSE

  versions <- names(master_list$templates$mrm_guides)
  versions <- setdiff(versions, "by_plate")

  #mrm_template selection
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

    likely_version <- result$best_match_name
    if (!is.na(likely_version) && likely_version %in% versions) {
      versions <- c(likely_version, setdiff(versions, likely_version))
    }
  }

  for (i in seq_along(versions)) {
    version <- versions[i]

    if (sil_found)
      break

    full_msg <- paste("Starting peak picking and integration using version:",
                      version)
    border <- paste(rep("=", nchar(full_msg) + 4), collapse = "")
    message("\n", border)
    message("= ", full_msg, " =")
    message(border, "\n")

    version_success <- tryCatch({
      master_list$project_details$is_ver <- version
      master_list$summary_tables$project_summary <- create_summary_table(master_list, plate_idx)
      master_list$templates$mrm_guides$by_plate[[plate_idx]] <- optimise_retention_times(master_list, plate_idx)

      export_files(master_list, plate_idx)
      PeakForgeR_command <- execute_PeakForgeR_command(master_list, plate_idx)

      system_success <- tryCatch({
        output_file <- file.path(
          master_list$project_details$project_dir,
          plate_idx,
          "data",
          "PeakForgeR",
          "CMD_output.txt"
        )

        run_system_command(PeakForgeR_command, output_file)
        TRUE

      }, error = function(e) {
        cat("System command failed in version",
            version,
            ":",
            e$message,
            "\n")
        flush.console()
        FALSE
      })

      if (!system_success) {
        message("Skipping version due to system command failure.\n")
        return(FALSE)
      }

      master_list$data$PeakForgeR_report[[plate_idx]] <- reimport_PeakForgeR_file(master_list, plate_idx)
      sil_result <- check_sil_standards(master_list, plate_idx, version)

      if (isTRUE(sil_result)) {
        sil_found <- TRUE
        save_plate_data(master_list, plate_idx)
        message("SIL matched to version:", version, "\n")
        message("PeakForgeR data saved for plate:", plate_idx, "\n")
      } else {
        message("No SIL standards detected with version: ",
                version,
                "- trying next version\n")
      }

      TRUE
    }, error = function(e) {
      message("\nError during processing of version ",
              version,
              ":\n",
              e$message,
              "\n")
      message("--- End of error for version ", version, "---\n")
      flush.console()
      FALSE
    })

    if (!version_success)
      next
  }

  if (!sil_found) {
    stop(
      paste(
        "No SIL internal standards detected in plate",
        plate_idx,
        "after trying all method versions. Please ensure your mrm_guide is the transitions used in the project!"
      )
    )
  }

  setwd(master_list$project_details$project_dir)
  master_list <- update_script_log(
    master_list,
    "peak_picking_and_integration",
    "mzR_mzml_import",
    "next_plate_for_processing"
  )
  return(master_list)
}


###Sub Functions----

#' create_summary_table
#'
#' This function creates a summary table for a given plate in the master list, including various project details.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param plate_idx The index of the plate to create the summary table for.
#' @return A tibble containing the summary table with project details and their values.
#' @examples
#' \dontrun{
#' create_summary_table(master_list, plate_idx)
#' }
create_summary_table <- function(master_list, plate_idx) {
  Temp_list <- master_list$project_details[c(
    "project_dir",
    "lipidExploreR_version",
    "user_name",
    "project_name",
    "qc_type",
    "plateID",
    "is_ver"
  )]
  Temp_list$plateID <- paste(plate_idx)
  project_summary <- tibble::tibble(unlist(Temp_list)) %>%
    tibble::add_column(
      "Project detail" = c(
        "local directory",
        "lipidExploreR version",
        "user initials",
        "project name",
        "project QC",
        "plateID",
        "int. std. version"
      ),
      .before = 1
    )
  project_summary <- stats::setNames(project_summary, c("Project detail", "value"))
  return(project_summary)
}


#' optimise_retention_times
#'
#' This function optimises retention times for each plate in the master list using the mzR_mrm_findR function and updates the MRM guide.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param plate_idx A vector of plate indices to optimise retention times for.
#' @return A list containing the optimised retention times and updated MRM guide for each plate.
#' @examples
#' \dontrun{
#' optimise_retention_times(master_list, plate_idx)
#' }
optimise_retention_times <- function(master_list, plate_idx) {
  by_plate <- list()
  idx <- plate_idx
  result <- mzR_mrm_findR(
    FUNC_mzR = master_list$data[[idx]],
    FUNC_mrm_guide = master_list$templates$mrm_guides[[master_list$project_details$is_ver]]$mrm_guide %>% janitor::clean_names(),
    FUNC_OPTION_qc_type = master_list$project_details$qc_type
  ) %>% append(master_list$templates$mrm_guides[[master_list$project_details$is_ver]])

  #Replace 0 retention times with default
  zero_rt_indices <- which(result$mrm_guide_updated$`Explicit Retention Time` == 0)
  result[["mrm_guide_updated"]][["Explicit Retention Time"]][zero_rt_indices] <-  result[["mrm_guide"]][["Explicit Retention Time"]][zero_rt_indices]

  by_plate[[plate_idx]] <- result
  by_plate[[plate_idx]]$mrm_guide_updated <- stats::setNames(by_plate[[plate_idx]]$mrm_guide_updated, names(by_plate[[plate_idx]]$mrm_guide))


  message("Successfully optimised retention times!")
  return(by_plate)
}


#' export_files
#'
#' This function exports various files related to the project for a given plate, including updated MRM guides, peak boundaries, and default templates.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param plate_idx The index of the plate to export files for.
#' @return Exports CSV, SKY, and TSV files to the specified project directory.
#' @examples
#' \dontrun{
#' export_files(master_list, plate_idx)
#' }
export_files <- function(master_list, plate_idx) {
  long_path <- file.path(master_list$project_details$project_dir,
                         plate_idx,
                         "data",
                         "PeakForgeR")
  if (length(long_path) > 260 & .Platform$OS.type == "windows") {
    cmd <- paste0('cmd /c mklink /J "C:\\PeakForgeR_short" "', long_path, '"')
    system(cmd, intern = FALSE, ignore.stdout = TRUE)
    PeakForgeR_path <- "C:/PeakForgeR_short"
  } else {
    PeakForgeR_path <- long_path
  }

  readr::write_csv(
    x = master_list$templates$mrm_guides$by_plate[[plate_idx]][[plate_idx]]$mrm_guide_update,
    file = file.path(
      PeakForgeR_path,
      paste0(Sys.Date(), "_RT_update_", plate_idx, ".csv")
    )
  )

  readr::write_csv(
    x = master_list$templates$mrm_guides$by_plate[[plate_idx]][[plate_idx]]$peak_boundary_update,
    file = file.path(
      PeakForgeR_path,
      paste0(Sys.Date(), "_peak_boundary_update_", plate_idx, ".csv")
    )
  )

  file.copy(
    from = system.file("templates", "default_skyline_file.sky", package = "MetaboExploreR"),
    to = paste0(PeakForgeR_path)
  )

  file.rename(
    from = file.path(PeakForgeR_path, "default_skyline_file.sky"),
    to = file.path(PeakForgeR_path, paste0(Sys.Date(), "_", plate_idx, ".sky"))
  )
  file.copy(
    from = system.file("templates", "default_csv.csv", package = "MetaboExploreR"),
    to = file.path(PeakForgeR_path)
  )

  file.rename(
    from = file.path(PeakForgeR_path, "default_csv.csv"),
    to = file.path(
      PeakForgeR_path,
      paste0(Sys.Date(), "_PeakForgeR_", plate_idx, ".csv")
    )
  )
  file.copy(
    from = system.file("templates", "default_tsv.tsv", package = "MetaboExploreR"),
    to = file.path(PeakForgeR_path)
  )

  file.rename(
    from = file.path(PeakForgeR_path , "default_tsv.tsv"),
    to = file.path(
      PeakForgeR_path,
      paste0(Sys.Date(), "_", plate_idx, "_chromatograms.tsv")
    )
  )
  file.copy(
    from = system.file(
      "templates",
      "YYYY-MM-DD_PeakForgeR_project_name.skyr",
      package = "MetaboExploreR"
    ),
    to = file.path(PeakForgeR_path)
  )

  if (length(long_path) > 260 & .Platform$OS.type == "windows") {
    unlink("C:/PeakForgeR_short", recursive = TRUE)
  }
}


#' execute_PeakForgeR_command
#'
#' This function executes a Skyline system command to process mzML files and generate various reports for a given plate.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param plate_idx The index of the plate to execute the Skyline command for.
#' @return Executes the Skyline command and generates reports and chromatogram files saving to project directory.
#' @examples
#' \dontrun{
#' execute_PeakForgeR_command(master_list, plate_idx)
#' }
execute_PeakForgeR_command <- function(master_list, plate_idx) {
  # Get absolute path to data directory
  data_dir <- normalizePath(file.path(master_list$project_details$project_dir, plate_idx, "data"))

  # Build filenames
  date_str <- Sys.Date()
  base_path <- file.path("PeakForgeR")

  in_file <- file.path(base_path, paste0(date_str, "_", plate_idx, ".sky"))
  import_transition_list <- file.path(base_path,
                                      paste0(date_str, "_RT_update_", plate_idx, ".csv"))
  mzml_path <- "mzml"
  import_peak_boundaries <- file.path(base_path,
                                      paste0(date_str, "_peak_boundary_update_", plate_idx, ".csv"))
  out_file <- file.path(base_path, paste0(date_str, "_", plate_idx, ".sky"))
  report_name <- "YYYY-MM-DD_PeakForgeR_project_name"
  report_file <- file.path(base_path,
                           paste0(date_str, "_PeakForgeR_", plate_idx, ".csv"))
  report_template <- file.path(base_path, "YYYY-MM-DD_PeakForgeR_project_name.skyr")
  chromatogram_file <- file.path(base_path,
                                 paste0(date_str, "_", plate_idx, "_chromatograms.tsv"))

  #Add the following to the command if missing targets to adjust mz tolerance when matching template to mrm_guide
  #--instrument-method-mz-tolerance=0.0006 \

  docker_command <- sprintf(
    'docker run --rm -v "%s:/data" proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:latest wine SkylineCmd \
      --dir=/data \
      --in=%s \
      --instrument-method-mz-tolerance=0.055 \
      --import-transition-list=%s \
      --import-all=%s \
      --import-peak-boundaries=%s \
      --save-settings --overwrite \
      --out=%s \
      --report-conflict-resolution=overwrite \
      --report-name=%s \
      --report-file=%s \
      --report-add=%s \
      --report-format=csv \
      --report-invariant \
      --chromatogram-file=%s \
      --chromatogram-precursors \
      --chromatogram-products \
      --chromatogram-base-peaks \
      --chromatogram-tics',
    data_dir,
    #-v "%s:/data"
    shQuote(in_file),
    #--in=%s
    shQuote(import_transition_list),
    #--import-transition-list=%s
    shQuote(mzml_path),
    #--import-all=%s
    shQuote(import_peak_boundaries),
    #--import-peak-boundaries=%s
    shQuote(out_file),
    #--out=%s
    shQuote(report_name),
    #--report-name=%s
    shQuote(report_file),
    # --report-file=%s
    shQuote(report_template),
    # --report-add=%s
    shQuote(chromatogram_file) #--chromatogram-file=%s
  )

  return(docker_command)
}


#' run_system_command
#'
#' This function wraps the system command to run docker skyline allowing for
#' unit testing. It captures all output (stdout and stderr) and writes it to a .txt file.
#' @keywords internal
#' @param PeakForgeR_command developed skyline docker command
#' @param output_file optional path to save the command output
#' @examples
#' \dontrun{
#' run_system_command("docker run skyline", "output.txt")
#' }
run_system_command <- function(PeakForgeR_command, output_file) {
  # Redirect both stdout and stderr
  full_command <- paste0(PeakForgeR_command, " 2>&1")

  result <- tryCatch({
    suppressWarnings(system(full_command, intern = TRUE))
  }, error = function(e) {
    message("Error occurred: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("Warning occurred: ", w$message)
    return(NULL)
  })

  if (is.null(result)) {
    message("The system command failed, but the error was suppressed.")
  } else {
    if (!is.null(output_file)) {
      writeLines(result, con = output_file)
    }
    return(result)
  }
}

#' reimport_PeakForgeR_file
#'
#' This function reimports a PeakForgeR file for a given plate, converts specific columns to numeric, and cleans the column names.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param plate_idx The index of the plate to reimport the Skyline file for.
#' @return A data frame containing the reimported PeakForgeR data with cleaned column names.
#' @examples
#' \dontrun{
#' reimport_PeakForgeR_file(master_list, plate_idx)
#' }
reimport_PeakForgeR_file <- function(master_list, plate_idx) {
  long_path <- file.path(master_list$project_details$project_dir,
                         plate_idx,
                         "data",
                         "PeakForgeR")
  if (length(long_path) > 260 & .Platform$OS.type == "windows") {
    cmd <- paste0('cmd /c mklink /J "C:\\PeakForgeR_short" "', long_path, '"')
    system(cmd, intern = FALSE, ignore.stdout = TRUE)
    PeakForgeR_path <- "C:/PeakForgeR_short"
  } else {
    PeakForgeR_path <- file.path(master_list$project_details$project_dir,
                              plate_idx,
                              "data",
                              "PeakForgeR")
  }

  PeakForgeR_data <- suppressWarnings(readr::read_csv(
    file = list.files(
      PeakForgeR_path,
      pattern = paste0("_PeakForgeR_", plate_idx),
      full.names = TRUE
    ),
    show_col_types = FALSE
  ))
  cols_to_convert <- c(
    "PrecursorMz",
    "ProductMz",
    "RetentionTime",
    "StartTime",
    "EndTime",
    "Area",
    "Height"
  )
  suppressWarnings(PeakForgeR_data[cols_to_convert] <- lapply(PeakForgeR_data[cols_to_convert], as.numeric))
  PeakForgeR_data <- janitor::clean_names(PeakForgeR_data)

  if (length(long_path) > 260 & .Platform$OS.type == "windows") {
    unlink("C:/PeakForgeR_short", recursive = TRUE)
  }

  return(PeakForgeR_data)
}

#' check_sil_standards
#'
#' This function checks if the SIL standards on a given plate match the SIL standards for a specified mrm_guide version.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param plate_idx The index of the plate to check the SIL standards for.
#' @param current_version The version of the MRM guide to compare against.
#' @return A logical value indicating whether the SIL standards on the plate match the SIL standards for the specified version.
#' @examples
#' \dontrun{
#' check_sil_standards(master_list, plate_idx, current_version)
#' }
check_sil_standards <- function(master_list, plate_idx, current_version) {
  # Extract SILs from the plate
  sil_on_plate <- master_list$data$PeakForgeR_report[[plate_idx]] %>%
    dplyr::filter(grepl("SIL", molecule_name, ignore.case = TRUE)) %>%
    dplyr::pull(molecule_name) %>%
    unique()

  template_versions <- setdiff(names(master_list$templates$mrm_guides), "by_plate")

  sil_by_version <- lapply(template_versions, function(version) {
    master_list$templates$mrm_guides[[version]]$mrm_guide %>%
      dplyr::filter(grepl("SIL", `Precursor Name`, ignore.case = TRUE)) %>%
      dplyr::pull(`Precursor Name`) %>%
      unique()
  })
  names(sil_by_version) <- template_versions

  current_sils <- sil_by_version[[current_version]]

  other_versions <- setdiff(template_versions, current_version)
  other_sils <- unlist(sil_by_version[other_versions], use.names = FALSE)

  # Compute SILs unique to current version
  unique_sils_current <- setdiff(current_sils, other_sils)

  # Match plate SILs to unique SILs and keep if 90% match
  matching_sils <- intersect(sil_on_plate, unique_sils_current)
  match_ratio <- length(matching_sils) / length(unique_sils_current)

  sil_found <- match_ratio >= 0.90

  return(sil_found)
}

#' save_plate_data
#'
#' This function saves the master list data for a given plate to an RDA file.
#' @keywords internal
#' @param master_list A list containing project details and data.
#' @param plate_idx The index of the plate to save the data for.
#' @return Saves the master list data to an RDA file in the specified project directory.
#' @examples
#' \dontrun{
#' save_plate_data(master_list, plate_idx)
#' }
save_plate_data <- function(master_list, plate_idx) {
  long_path <- file.path(master_list$project_details$project_dir,
                         plate_idx,
                         "data",
                         "rda")
  if (length(long_path) > 260 & .Platform$OS.type == "windows") {
    cmd <- paste0('cmd /c mklink /J "C:\\rda" "', long_path, '"')
    system(cmd, intern = FALSE, ignore.stdout = TRUE)
    save_path <- "C:/rda"
  } else {
    # Use full path directly on macOS/Linux
    save_path <- file.path(master_list$project_details$project_dir,
                           plate_idx,
                           "data",
                           "rda")
  }


  callr::r_bg(function(master_list, save_path, plate_idx) {
    save(
      master_list,
      file = file.path(
        save_path,
        paste0(Sys.Date(), "_", master_list$project_details$project_name, "_",
               plate_idx, "_PeakForgeR.rda")
      ),
      compress = FALSE
    )}, args = list(master_list = master_list, save_path = save_path, plate_idx = plate_idx)
  )


  if (length(long_path) > 260 & .Platform$OS.type == "windows") {
    unlink("C:/rda", recursive = TRUE)
  }
}
#.----
#Archive Raw Files ----

###Primary Function----
#' Archive Raw Files
#'
#' This function moves raw files (wiff and mzML) to an archive directory after processing is complete.
#'
#' @param project_directory Path to the directory for the project parsed from PeakForgeR.
#' @return None. The function performs the archive operation and a message upon successful completion.
#' @examples
#' \dontrun{
#' archive_raw_files("path/to/project_directory")
#' }
archive_raw_files <- function(project_directory) {
  validate_project_directory(project_directory)
  archive_files(project_directory, "MetaboExploreR_logs")
  archive_files(project_directory, "raw_data")
  #archive_files(project_directory, "msConvert_mzml_output")
  message("\n PeakForgeR is now finished running all plates :)")
}

###Sub Functions----

#' Archive Files
#'
#' This function moves a specified folder to the archive directory within the project directory.
#' @keywords internal
#' @param project_directory The directory of the project.
#' @param folder_name The name of the folder to be archived.
#' @return Moves the specified folder to the archive directory.
#' @examples
#' \dontrun{
#' archive_files(project_directory, folder_name)
#' }
archive_files <- function(project_directory, folder_name) {
  move_folder(
    file.path(project_directory, folder_name),
    file.path(project_directory, "archive")
  )
}

#Move Folder Functions----

###Primary Function----
#' move_folder
#' This function moves a folder from the source directory to the destination directory, waits until files are not in use, and then deletes the source directory.
#' @keywords internal
#' @param source_dir Directory path for the folder to copy.
#' @param dest_dir Directory path for the folder to be copied to.
#' @param max_wait Maximum wait time (in seconds) for the system prior to deleting the moved folder.
#' @param max_retries Maximum number of attempts to try move/delete prior to error. Default 30
#' @return None. The function performs the move operation and a message upon successful completion.
#' @examples
#' \dontrun{
#' move_folder(source_dir = "path/to/source",
#'             dest_dir = "path/to/destination",
#'             max_wait = 60,
#'             max_retries = 30)
#' }
move_folder <- function(source_dir,
                        dest_dir,
                        max_wait = 60,
                        max_retries = 30) {
  validate_directories(source_dir, dest_dir)
  files_to_copy <- copy_files(source_dir, dest_dir)
  wait_until_files_free(files_to_copy, max_wait, max_retries)
  delete_source_directory(source_dir)
  return(TRUE)
}

###Sub Functions----

#' validate_directories
#' This function validates the existence of the source directory and creates the destination directory if it does not exist.
#' @keywords internal
#' @param source_dir Directory path for the folder to copy.
#' @param dest_dir Directory path for the folder to be copied to.
#' @return None. The function performs directory validation.
#' @examples
#' \dontrun{
#' validate_directories(source_dir = "path/to/source", dest_dir = "path/to/destination")
#' }
validate_directories <- function(source_dir, dest_dir) {
  if (!dir.exists(source_dir)) {
    message(paste("Source directory does not exist:", source_dir))
  }
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
  }
}


  #' copy_files
  #' This function copies files from the source directory to the destination directory.
  #' @keywords internal
  #' @param source_dir Directory path for the folder to copy.
  #' @param dest_dir Directory path for the folder to be copied to.
  #' @return A vector of file paths that were copied.
  #' @examples
  #' \dontrun{
  #' copy_files(source_dir = "path/to/source", dest_dir = "path/to/destination")
  #' }
  copy_files <- function(source_dir, dest_dir) {
    files_to_copy <- list.files(source_dir, full.names = TRUE)

    if (!dir.exists(dest_dir)) {
      message(paste("Destination is not a directory:", dest_dir))
    }

    success <- file.copy(files_to_copy, dest_dir, recursive = TRUE)
    if (!all(success)) {
      message("Some files failed to copy.")
    }
    return(files_to_copy)
  }

  #' wait_until_files_free
  #' This function waits until files are not in use by attempting to rename them.
  #' @keywords internal
  #' @param files_to_copy A vector of file paths to check.
  #' @param max_wait Maximum wait time (in seconds) for the system prior to deleting the moved folder. Default is 60
  #' @param max_retries Maximum number of attempts to try move/delete prior to error. Default is 30.
  #' @return None. The function waits until files are free.
  #' @examples
  #' \dontrun{
  #' wait_until_files_free(files_to_copy, max_wait = 60, max_retries = 30)
  #' }
  wait_until_files_free <- function(files_to_copy,
                                    max_wait = 60,
                                    max_retries = 30) {
    for (test_file in files_to_copy) {
      retries <- 0
      start_time <- Sys.time()
      while (TRUE) {
        if (!file.exists(test_file))
          break

        test_rename <- try(file.rename(test_file, paste0(test_file, ".tmp")), silent = TRUE)
        if (!inherits(test_rename, "try-error") && test_rename) {
          file.rename(paste0(test_file, ".tmp"), test_file)
          break
        }

        if (as.numeric(Sys.time() - start_time, units = "secs") > max_wait ||
            retries >= max_retries) {
          stop(paste("File still in use after waiting:", test_file))
        }

        Sys.sleep(2)
        retries <- retries + 1
      }
    }
  }

  #' delete_source_directory
  #' This function deletes the source directory after files have been moved.
  #' @keywords internal
  #' @param source_dir Directory path for the folder to delete.
  #' @return None. The function deletes the source directory.
  #' @examples
  #' \dontrun{
  #' delete_source_directory(source_dir = "path/to/source")
  #' }
  delete_source_directory <- function(source_dir) {
    unlink(source_dir, recursive = TRUE, force = TRUE)
    if (dir.exists(source_dir)) {
      message(paste("Failed to delete directory:", source_dir))
    } else {
      message(paste("Successfully moved and deleted:", source_dir))
    }
  }

  extract_transition_values <- function(ids, tag) {
    TRANSITION_REGEX <- "([+-]?(?:\\d+(?:\\.\\d+)?|\\.\\d+)(?:[eE][+-]?\\d+)?)"
    rx <- paste0("\\b", tag, "=", TRANSITION_REGEX)
    matches <- regmatches(ids, regexec(rx, ids, perl = TRUE))
    as.numeric(vapply(matches, function(m) if (length(m) >= 2) m[2] else NA_character_, ""))
  }

  build_transitions_list <- function(chroms, targets) {
    stopifnot(is.data.frame(chroms), "id" %in% names(chroms))
    stopifnot(is.data.frame(targets))
    count_decimals <- function(x) {
      s <- sub("e.*", "", format(x, digits = 22, scientific = FALSE, trim = TRUE))
      frac <- sub("^-?\\d+\\.?", "", s)
      nchar(sub("0+$", "", frac))
    }
    rounded_equal <- function(a, vec) {
      da <- count_decimals(a)
      dv <- count_decimals(vec)
      dec <- pmin.int(da, dv)
      round(vec, dec) == round(a, dec)
    }
    ids_chr <- as.character(chroms$id)
    q1_vals <- extract_transition_values(ids_chr, "Q1")
    q3_vals <- extract_transition_values(ids_chr, "Q3")
    n_targets <- nrow(targets)
    idx1 <- rep.int(NA_integer_, n_targets)
    idx0 <- rep.int(NA_integer_, n_targets)
    match_id <- rep(NA_character_, n_targets)
    explicit_rt <- suppressWarnings(as.numeric(targets$`Explicit Retention Time`))
    explicit_window <- suppressWarnings(as.numeric(targets$`Explicit Retention Time Window`))
    for (i in seq_len(n_targets)) {
      target_q1 <- suppressWarnings(as.numeric(targets$`Precursor Mz`[i]))
      target_q3 <- suppressWarnings(as.numeric(targets$`Product Mz`[i]))
      if (!is.finite(target_q1) || !is.finite(target_q3)) next
      have_both <- !is.na(q1_vals) & !is.na(q3_vals)
      q1_match <- have_both & rounded_equal(target_q1, q1_vals)
      q3_match <- have_both & rounded_equal(target_q3, q3_vals)
      hits <- which(q1_match & q3_match)
      if (length(hits)) {
        # idx1[i] <- hits[1L]
        idx0[i] <- hits[1L] - 1L
        match_id[i] <- ids_chr[hits[1L]]
      }
    }
    res <- data.frame(
      idx = idx0,
      # idx1 = idx1,
      id = match_id,
      rt = explicit_rt,
      window = explicit_window,
      stringsAsFactors = FALSE
    )
    res <- res[order(res$idx, na.last = TRUE), , drop = FALSE]
    rownames(res) <- NULL
    res
  }

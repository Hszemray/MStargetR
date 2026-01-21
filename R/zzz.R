#Ignore column names used in dplyr and tibble pipes

#' @importFrom utils globalVariables
utils::globalVariables(
  c(
    ".",
    ".data",
    "AcquiredTime",
    "Area",
    "failed_samples",
    "FileName",
    "FullPeptideName",
    "MoleculeName",
    "molecule_name",
    "Note",
    "PC1",
    "PC2",
    "Precursor",
    "Precursor Name",
    "Name",
    "SIL",
    "SIL_name",
    "V1",
    "V2",
    "batch",
    "concentration_factor",
    "corrected_mean",
    "dataBatch",
    "dataSource",
    "facet_label",
    "invalid_wiff_files",
    "lipid",
    "lipid_class",
    "matches",
    "metabolite_code",
    "name",
    "original_mean",
    "p1",
    "p2",
    "p3",
    "plateID",
    "precursor_name",
    "project_directory",
    "sample.flag",
    "sample_data_source",
    "sample_ID",
    "sample_matrix",
    "sample_name",
    "sample_plate_id",
    "sample_plate_order",
    "sample_run_index",
    "sample_timestamp",
    "sample_type",
    "sample_type_factor",
    "sample_type_factor_rev",
    "source_prefix",
    "template_version",
    "value"
  )
)

.onAttach <- function(libname, pkgname) {
  # Don't run in non-interactive sessions (e.g., R CMD check, CI, tests)
  if (!interactive()) return()

  # Don't run during install or test environments
  if (Sys.getenv("R_METABOEXPORER_SKIP_ONATTACH", "false") == "true") return()

  user <- NULL

  # Try prompting with rstudioapi
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    user <- tryCatch(
      rstudioapi::showPrompt("Login", "Enter your username:", ""),
      error = function(e) NULL
    )
  }

  # Fallback to tcltk
  if (is.null(user) || identical(user, "")) {
    if (requireNamespace("tcltk", quietly = TRUE)) {
      try({
        tt <- tcltk::tktoplevel()
        tcltk::tkwm.title(tt, "Enter Username")

        user_var <- tcltk::tclVar("")
        entry <- tcltk::tkentry(tt, textvariable = user_var)
        tcltk::tkpack(tcltk::tklabel(tt, text = "Username:"), entry)

        done <- tcltk::tclVar(0)
        onOK <- function()
          tcltk::tclvalue(done) <- 1
        ok_button <- tcltk::tkbutton(tt, text = "OK", command = onOK)
        tcltk::tkpack(ok_button)

        tcltk::tkwait.variable(done)
        user <- tcltk::tclvalue(user_var)
        tcltk::tkdestroy(tt)
      }, silent = TRUE)
    }
  }

  # If user is ANPC, show the R Markdown
  if (identical(user, "ANPC")) {
    rmd_file <- system.file("rmd", "workflow_HS.Rmd", package = pkgname)

    if (nzchar(rmd_file)) {
      if (requireNamespace("rstudioapi", quietly = TRUE) &&
          rstudioapi::isAvailable()) {
        try(rstudioapi::navigateToFile(rmd_file), silent = TRUE)
      } else if (requireNamespace("rmarkdown", quietly = TRUE)) {
        out_html <- tempfile(fileext = ".html")
        try({
          rmarkdown::render(rmd_file, output_file = out_html, quiet = TRUE)
          browseURL(out_html)
        }, silent = TRUE)
      }
    }
  }else{
    rmd_file <- system.file("rmd", "workflow_generic.Rmd", package = pkgname)

    if (nzchar(rmd_file)) {
      if (requireNamespace("rstudioapi", quietly = TRUE) &&
          rstudioapi::isAvailable()) {
        try(rstudioapi::navigateToFile(rmd_file), silent = TRUE)
      } else if (requireNamespace("rmarkdown", quietly = TRUE)) {
        out_html <- tempfile(fileext = ".html")
        try({
          rmarkdown::render(rmd_file, output_file = out_html, quiet = TRUE)
          browseURL(out_html)
        }, silent = TRUE)
      }
    }

  }
}


#openms---
.openms_state <- new.env(parent = emptyenv())
.openms_state$om <- NULL
.openms_state$configured <- FALSE

.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible())

  # Package-best practice: declare Python deps so reticulate can provision them.
  reticulate::py_require("numpy")
  reticulate::py_require("pyopenms")

  .openms_state$configured <- TRUE
  invisible(TRUE)
}


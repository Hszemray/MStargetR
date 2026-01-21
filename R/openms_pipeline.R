
# PeakForgeR / MetaboExploreR - Production OpenMS (pyOpenMS) Backend (Single Script)

# Supports:
#   - SRM/MRM chromatogram mzML (0 spectra, many chromatograms) -> uses existing traces
#   - PRM/DIA/full scan mzML (spectra present) -> extracts chromatograms using
#     OpenSWATH-style ChromatogramExtractor (pyOpenMS 3.5+ API)
#
# Output TSV schema:
# FileName, MoleculeListName, MoleculeName, PrecursorMz, ProductMz, RetentionTime,
# StartTime, EndTime, Area, Height, AcquiredTime
#


# ------------------------------ Defaults -------------------------------------
PF_CFG <- list(
  envname = "r-openms-env",
  python_version = "3.10",
  backend = "conda",            # conda recommended
  auto_install_env = FALSE,     # if TRUE, attempts to install conda env + pyopenms
  mz_tol_q1 = 0.02,             # Da tolerance for Q1 matching (SRM)
  mz_tol_q3 = 0.02,             # Da tolerance for Q3 matching (SRM)
  mz_extraction_window = 0.05,  # Da window for spectra->chrom extraction (PRM/DIA)
  ppm = FALSE,
  rt_extraction_window = 0,     # 0 means "no RT restriction" in our wrapper
  filter = "tophat",
  max_log_mb = 10,              # rotate logs at ~N MB
  debug = FALSE,
  parallel_files = TRUE,     # enable parallel over mzML files
  n_workers = max(1L, parallel::detectCores() - 1L),
  notify = TRUE,                         # show progress bar + per-file messages
  notify_handler = "auto"                # "auto" | "rstudio" | "txt" | "progress"

)

# --------------------------- Small Utilities --------------------------------
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

assert_has_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install it with install.packages('%s').", pkg, pkg),
         call. = FALSE)
  }
  invisible(TRUE)
}

assert_has_function <- function(fname) {
  if (!exists(fname, mode = "function", inherits = TRUE)) {
    stop(sprintf(
      "Required function '%s' was not found in the R session.\n",
      fname
    ), call. = FALSE)
  }
  invisible(TRUE)
}

normalize_path_safe <- function(path) normalizePath(path, winslash = "/", mustWork = FALSE)

make_logger <- function(log_file, max_mb = 10) {
  dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)

  rotate_if_needed <- function() {
    if (!file.exists(log_file)) return(invisible(TRUE))
    sz <- file.info(log_file)$size
    if (is.na(sz)) return(invisible(TRUE))
    if (sz > max_mb * 1024^2) {
      ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
      file.rename(log_file, paste0(log_file, ".", ts, ".bak"))
    }
    invisible(TRUE)
  }

  rotate_if_needed()
  con <- file(log_file, open = "a", encoding = "UTF-8")
  force(con)

  structure(list(
    write = function(level = "INFO", msg) {
      rotate_if_needed()
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      writeLines(sprintf("[%s] %s: %s", timestamp, level, enc2utf8(msg)), con = con)
      flush(con)
    },
    close = function() try(close(con), silent = TRUE)
  ), class = "pf_logger")
}

# ---------------------- Reticulate + pyOpenMS env ----------------------------
.openms_state <- local({
  e <- new.env(parent = emptyenv())
  e$om <- NULL
  e$configured <- FALSE
  e
})

openms_use_backend <- function(backend = c("auto", "conda", "virtualenv"),
                               envname = "r-openms-env") {
  assert_has_package("reticulate")
  backend <- match.arg(backend)
  options(PeakForgeR.openms_backend = backend)
  options(PeakForgeR.openms_envname = envname)
  .openms_state$om <- NULL
  .openms_state$configured <- FALSE
  invisible(TRUE)
}

install_openms_env <- function(envname = "r-openms-env", python_version = "3.10") {
  assert_has_package("reticulate")
  reticulate::install_miniconda()

  envs <- reticulate::conda_list()
  if (!envname %in% envs$name) {
    reticulate::conda_create(envname, python_version = python_version)
  }

  # Stability: numpy via conda, pyopenms via pip
  reticulate::conda_install(envname, "pip", pip = FALSE)
  reticulate::conda_install(envname, "numpy", pip = FALSE)
  reticulate::conda_install(envname, "pyopenms", pip = TRUE)

  invisible(TRUE)
}

openms_maybe_use_env <- function() {
  assert_has_package("reticulate")
  backend <- getOption("PeakForgeR.openms_backend", default = "auto")
  envname <- getOption("PeakForgeR.openms_envname", default = "r-openms-env")
  if (identical(backend, "conda")) {
    reticulate::use_condaenv(envname, required = TRUE)
  } else if (identical(backend, "virtualenv")) {
    reticulate::use_virtualenv(envname, required = TRUE)
  }
  invisible(TRUE)
}

openms_get <- function() {
  assert_has_package("reticulate")
  if (!is.null(.openms_state$om)) return(.openms_state$om)

  openms_maybe_use_env()

  if (!isTRUE(.openms_state$configured)) {
    if ("py_require" %in% getNamespaceExports("reticulate")) {
      reticulate::py_require("numpy")
      reticulate::py_require("pyopenms")
    }
    .openms_state$configured <- TRUE
  }

  om <- try(reticulate::import("pyopenms", delay_load = FALSE), silent = TRUE)
  if (inherits(om, "try-error")) {
    cfg <- try(reticulate::py_config(), silent = TRUE)
    stop(
      paste0(
        "Failed to import 'pyopenms' via reticulate.\n\n",
        "Likely causes:\n",
        "  • Python already initialized to a different interpreter in this R session\n",
        "  • pyopenms not installed in selected env\n\n",
        "Fix (recommended): Restart R, then run install_openms_env(), then retry.\n\n",
        "Diagnostics:\n",
        "  py_config(): ", if (!inherits(cfg, "try-error")) cfg$python else "<py_config failed>"
      ),
      call. = FALSE
    )
  }

  .openms_state$om <- om
  om
}

peakforger_setup_openms <- function(envname = PF_CFG$envname,
                                    python_version = PF_CFG$python_version,
                                    backend = c("conda", "auto", "virtualenv"),
                                    install = FALSE) {
  backend <- match.arg(backend)
  assert_has_package("reticulate")
  openms_use_backend(backend, envname = envname)
  if (isTRUE(install) && backend == "conda") {
    install_openms_env(envname = envname, python_version = python_version)
  }
  openms_get()
  invisible(TRUE)
}

# ------------------------------ mzML helpers --------------------------------
infer_mzml_time_unit <- function(mzml_path) {
  if (requireNamespace("xml2", quietly = TRUE)) {
    doc <- xml2::read_xml(mzml_path, options = c("NOBLANKS", "RECOVER"))
    nodes <- xml2::xml_find_all(doc, ".//cvParam[@accession='MS:1000595']")
    if (length(nodes)) {
      unit <- tolower(xml2::xml_attr(nodes[[1]], "unitName") %||% "")
      if (unit %in% c("minute", "minutes", "min")) return("min")
      if (unit %in% c("second", "seconds", "sec", "s")) return("sec")
      if (nzchar(unit)) return(unit)
    }
  }
  NA_character_
}

get_mzml_acquired_time <- function(mzml_path,
                                   out_format = "%m/%d/%Y %H:%M:%S",
                                   tz_out = "UTC",
                                   fallback_to_mtime = TRUE) {

  iso <- NA_character_

  if (requireNamespace("xml2", quietly = TRUE)) {
    doc <- xml2::read_xml(mzml_path, options = c("NOBLANKS", "RECOVER"))
    run <- xml2::xml_find_first(doc, ".//run")
    if (!inherits(run, "xml_missing")) {
      iso <- xml2::xml_attr(run, "startTimeStamp")
      if (is.na(iso) || !nzchar(iso)) iso <- NA_character_
    }
  } else {
    txt <- readLines(mzml_path, n = 5000, warn = FALSE)
    hit <- grep("startTimeStamp=", txt, value = TRUE)
    if (length(hit)) {
      m <- regmatches(hit[1], regexec('startTimeStamp="([^"]+)"', hit[1]))[[1]]
      if (length(m) >= 2) iso <- m[2]
    }
  }

  if (is.na(iso)) {
    if (fallback_to_mtime) {
      mt <- file.info(mzml_path)$mtime
      if (!is.na(mt)) return(format(mt, out_format, tz = tz_out))
    }
    return(NA_character_)
  }

  iso_norm <- iso
  iso_norm <- sub("Z$", "+0000", iso_norm)
  iso_norm <- sub("([+-]\\d\\d):(\\d\\d)$", "\\1\\2", iso_norm)

  t <- suppressWarnings(as.POSIXct(strptime(iso_norm, "%Y-%m-%dT%H:%M:%OS%z", tz = "UTC")))
  if (is.na(t)) {
    iso_no_tz <- sub("([+-]\\d{4})$", "", iso_norm)
    t <- suppressWarnings(as.POSIXct(strptime(iso_no_tz, "%Y-%m-%dT%H:%M:%OS", tz = "UTC")))
  }

  if (is.na(t)) {
    if (fallback_to_mtime) {
      mt <- file.info(mzml_path)$mtime
      if (!is.na(mt)) return(format(mt, out_format, tz = tz_out))
    }
    return(NA_character_)
  }

  format(t, out_format, tz = tz_out)
}

# ------------------ Transition table standardisation -------------------------
standardize_one_transition_table <- function(df) {
  df <- as.data.frame(df)

  required <- c(
    "Molecule List Name", "Precursor Name", "Precursor Mz", "Product Mz",
    "Explicit Retention Time", "Explicit Retention Time Window"
  )
  miss <- setdiff(required, names(df))
  if (length(miss)) {
    stop("Transition table missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  rt_min  <- suppressWarnings(as.numeric(df[["Explicit Retention Time"]]))
  win_min <- suppressWarnings(as.numeric(df[["Explicit Retention Time Window"]]))

  out <- data.frame(
    molecule_list_name = as.character(df[["Molecule List Name"]]),
    molecule_name      = as.character(df[["Precursor Name"]]),
    precursor_mz       = suppressWarnings(as.numeric(df[["Precursor Mz"]])),
    product_mz         = suppressWarnings(as.numeric(df[["Product Mz"]])),
    rt_min             = rt_min,
    rt_start_min       = ifelse(is.finite(rt_min) & is.finite(win_min), rt_min - win_min, NA_real_),
    rt_end_min         = ifelse(is.finite(rt_min) & is.finite(win_min), rt_min + win_min, NA_real_),
    stringsAsFactors   = FALSE
  )

  out$transition_name <- paste(
    out$molecule_list_name,
    out$molecule_name,
    sprintf("Q1%.4f", out$precursor_mz),
    sprintf("Q3%.4f", out$product_mz),
    sep = "|"
  )

  out <- out[is.finite(out$precursor_mz) & is.finite(out$product_mz), , drop = FALSE]
  if (!nrow(out)) stop("No valid transitions found after standardisation.", call. = FALSE)

  out
}

pick_guide_table <- function(inner, prefer = c("mrm_guide_updated", "mrm_guide")) {
  prefer <- match.arg(prefer)
  if (!is.list(inner)) stop("Expected an inner list with mrm_guide_updated / mrm_guide.", call. = FALSE)
  if (!is.null(inner[[prefer]])) return(inner[[prefer]])
  fallback <- setdiff(c("mrm_guide_updated", "mrm_guide"), prefer)
  if (!is.null(inner[[fallback]])) return(inner[[fallback]])
  stop("No guide table found (expected mrm_guide_updated or mrm_guide).", call. = FALSE)
}

standardize_transition_table_lipidyzer <- function(x, prefer = c("mrm_guide_updated", "mrm_guide")) {
  prefer <- match.arg(prefer)

  if (is.data.frame(x)) return(standardize_one_transition_table(x))

  if (is.list(x) && any(c("mrm_guide_updated", "mrm_guide") %in% names(x))) {
    df <- pick_guide_table(x, prefer = prefer)
    return(standardize_one_transition_table(df))
  }

  if (is.list(x)) {
    dfs <- lapply(x, function(inner) {
      df <- pick_guide_table(inner, prefer = prefer)
      standardize_one_transition_table(df)
    })
    dfs <- Filter(function(z) is.data.frame(z) && nrow(z) > 0, dfs)
    if (!length(dfs)) stop("No valid transitions found in list input.", call. = FALSE)
    out <- do.call(rbind, dfs)
    rownames(out) <- NULL
    return(out)
  }

  stop("Unsupported input type for x.", call. = FALSE)
}

sort_transitions_openms_like <- function(trans_df) {
  trans_df <- as.data.frame(trans_df)
  o <- order(trans_df$precursor_mz, trans_df$product_mz,
             trans_df$molecule_list_name, trans_df$molecule_name,
             trans_df$transition_name, na.last = TRUE)
  trans_df[o, , drop = FALSE]
}

# -------------------- Template selection (INCLUDED) --------------------------
select_method_template_for_plate <- function(master_list, plateID, project_directory,
                                             prefer_n_files = 1L,
                                             qc_label = NULL,
                                             verbose = TRUE) {
  # Requires match_templates_to_mzml from PeakForgeR ecosystem
  assert_has_function("match_templates_to_mzml")
  assert_has_package("tools")

  mzml_dir <- normalize_path_safe(file.path(project_directory, plateID, "data", "mzml"))
  mzml_files <- list.files(mzml_dir, pattern = "\\.mzML(\\.gz)?$", ignore.case = TRUE, full.names = TRUE)
  if (!length(mzml_files)) stop("No mzML files found in: ", mzml_dir, call. = FALSE)

  if (!is.null(qc_label) && nzchar(qc_label)) {
    non_qc <- mzml_files[!grepl(qc_label, basename(mzml_files), ignore.case = TRUE)]
    if (length(non_qc)) mzml_files <- non_qc
  }

  all_templates <- master_list$templates$mrm_guides
  all_templates <- all_templates[!names(all_templates) %in% "by_plate"]

  if (length(all_templates) < 1) {
    stop("No method templates found in master_list$templates$mrm_guides (excluding by_plate).",
         call. = FALSE)
  }

  templates <- lapply(all_templates, function(x) x$mrm_guide)

  prefer_n_files <- max(1L, min(as.integer(prefer_n_files), length(mzml_files)))
  candidates <- mzml_files[seq_len(prefer_n_files)]

  matches <- lapply(candidates, function(f) {
    match_templates_to_mzml(
      mrm_templates = templates,
      mzml_file = f,
      precursor_col = "Precursor Mz",
      product_col   = "Product Mz",
      name_col      = "Precursor Name",
      standardize   = TRUE,
      verbose       = verbose
    )
  })

  best_names <- vapply(matches, function(x) x$best_match_name %||% NA_character_, character(1))
  best_names <- best_names[!is.na(best_names) & nzchar(best_names)]
  if (!length(best_names)) stop("Template matching failed: no best_match_name.", call. = FALSE)

  tab <- sort(table(best_names), decreasing = TRUE)
  best <- names(tab)[1]
  master_list$project_details$is_ver <- best

  if (length(unique(best_names)) > 1) {
    warning(sprintf(
      "Template selection varied across %d files: %s. Using majority vote: %s",
      prefer_n_files, paste(unique(best_names), collapse = ", "), best
    ), call. = FALSE)
  }

  master_list
}

# -------------------- TargetedExperiment builder -----------------------------
set_compound_id_safe <- function(comp, cid) {
  assert_has_package("reticulate")
  for (m in c("setId", "setID", "setIdentifier")) {
    if (reticulate::py_has_attr(comp, m)) {
      ok <- tryCatch({ comp[[m]](cid); TRUE }, error = function(e) FALSE)
      if (ok) return(invisible(TRUE))
    }
  }
  if (reticulate::py_has_attr(comp, "id")) {
    comp$id <- cid
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

set_transition_ids_safe <- function(tr, id) {
  assert_has_package("reticulate")
  if (reticulate::py_has_attr(tr, "setNativeID")) try(tr$setNativeID(as.character(id)), silent = TRUE)
  if (reticulate::py_has_attr(tr, "setName"))     try(tr$setName(as.character(id)), silent = TRUE)
  if (reticulate::py_has_attr(tr, "setId"))       try(tr$setId(as.character(id)), silent = TRUE)
  else if (reticulate::py_has_attr(tr, "id"))     try({ tr$id <- as.character(id) }, silent = TRUE)
  invisible(TRUE)
}

set_transition_mzs_safe <- function(tr, q1, q3) {
  assert_has_package("reticulate")
  q1 <- as.numeric(q1); q3 <- as.numeric(q3)

  if (reticulate::py_has_attr(tr, "setPrecursorMZ")) tr$setPrecursorMZ(q1)
  else if (reticulate::py_has_attr(tr, "setPrecursorMz")) tr$setPrecursorMz(q1)
  else stop("No supported precursor m/z setter found.", call. = FALSE)

  if (reticulate::py_has_attr(tr, "setProductMZ")) tr$setProductMZ(q3)
  else if (reticulate::py_has_attr(tr, "setProductMz")) tr$setProductMz(q3)
  else stop("No supported product m/z setter found.", call. = FALSE)

  invisible(TRUE)
}

build_targeted_experiment_pf <- function(trans_df_sorted) {
  assert_has_package("reticulate")
  om <- openms_get()
  trans_df_sorted <- as.data.frame(trans_df_sorted)

  te <- om$TargetedExperiment()

  comp_key <- paste(trans_df_sorted$molecule_list_name, trans_df_sorted$molecule_name, sep = "||")
  uniq_keys <- unique(comp_key)
  comp_id_map <- setNames(vector("character", length(uniq_keys)), uniq_keys)

  for (k in uniq_keys) {
    parts <- strsplit(k, "\\|\\|", fixed = FALSE)[[1]]
    mol_list <- parts[1]; mol_name <- parts[2]
    comp <- om$Compound()
    cid  <- paste0(mol_list, "|", mol_name)
    set_compound_id_safe(comp, cid)
    te$addCompound(comp)
    comp_id_map[[k]] <- cid
  }

  for (i in seq_len(nrow(trans_df_sorted))) {
    row <- trans_df_sorted[i, , drop = FALSE]
    tid <- as.character(row$transition_name[[1]])
    tr <- om$ReactionMonitoringTransition()
    set_transition_ids_safe(tr, tid)
    set_transition_mzs_safe(tr, row$precursor_mz[[1]], row$product_mz[[1]])

    ck <- paste(row$molecule_list_name[[1]], row$molecule_name[[1]], sep = "||")
    if (reticulate::py_has_attr(tr, "setCompoundRef")) {
      tr$setCompoundRef(comp_id_map[[ck]])
    } else {
      stop("Missing setCompoundRef() in this pyopenms build.", call. = FALSE)
    }
    te$addTransition(tr)
  }

  te
}

# ---------------------- Fast transition index --------------------------------
build_transition_index <- function(trans_df_sorted, tol_q1, tol_q3) {
  trans_df_sorted <- as.data.frame(trans_df_sorted)
  b1 <- as.integer(round(trans_df_sorted$precursor_mz / tol_q1))
  b3 <- as.integer(round(trans_df_sorted$product_mz   / tol_q3))
  key <- paste(b1, b3, sep = ":")
  idx <- split(seq_len(nrow(trans_df_sorted)), key)
  list(idx = idx, tol_q1 = tol_q1, tol_q3 = tol_q3)
}

find_transition_match <- function(index, trans_df_sorted, q1, q3) {
  tol_q1 <- index$tol_q1; tol_q3 <- index$tol_q3
  b1 <- as.integer(round(q1 / tol_q1))
  b3 <- as.integer(round(q3 / tol_q3))
  keys <- as.vector(outer((b1-1):(b1+1), (b3-1):(b3+1),
                          FUN = function(x, y) paste(x, y, sep=":")))
  cand <- unlist(index$idx[keys], use.names = FALSE)
  if (!length(cand)) return(NA_integer_)
  sub <- trans_df_sorted[cand, , drop = FALSE]
  d1 <- abs(sub$precursor_mz - q1); d3 <- abs(sub$product_mz - q3)
  ok <- which(d1 <= tol_q1 & d3 <= tol_q3)
  if (!length(ok)) return(NA_integer_)
  cand[ ok[which.min(d1[ok] + d3[ok])] ]
}

# -------------------- SRM parsing + filtering --------------------------------
parse_q1q3_from_nativeid <- function(nid) {
  nid <- as.character(nid %||% "")
  q1m <- regmatches(nid, regexec("\\bQ1=([0-9]+(?:\\.[0-9]+)?)", nid))[[1]]
  q3m <- regmatches(nid, regexec("\\bQ3=([0-9]+(?:\\.[0-9]+)?)", nid))[[1]]
  q1 <- if (length(q1m) >= 2) as.numeric(q1m[2]) else NA_real_
  q3 <- if (length(q3m) >= 2) as.numeric(q3m[2]) else NA_real_
  list(q1=q1, q3=q3)
}

is_transition_like_chrom <- function(ch) {
  nid <- tryCatch(as.character(ch$getNativeID()), error = function(e) "")
  if (toupper(nid) == "TIC") return(FALSE)
  q1 <- suppressWarnings(as.numeric(ch$getPrecursor()$getMZ()))
  q3 <- suppressWarnings(as.numeric(ch$getProduct()$getMZ()))
  prod_mz <- suppressWarnings(as.numeric(ch$getMZ()))
  (is.finite(q1) && q1 > 0 && is.finite(q3) && q3 > 0) ||
    (is.finite(prod_mz) && prod_mz > 0) ||
    grepl("\\bQ1=", nid)
}

# -------------------- Hybrid chrom acquisition --------------------------------
extract_all_chromatograms <- function(mzml_path,
                                      targeted_exp,
                                      mz_extraction_window,
                                      ppm,
                                      rt_extraction_window,
                                      filter,
                                      log = NULL) {
  assert_has_package("reticulate")
  om <- openms_get()

  exp <- om$MSExperiment()
  om$MzMLFile()$load(mzml_path, exp)

  ns <- exp$getNrSpectra()
  nc <- exp$getNrChromatograms()
  if (!is.null(log) && isTRUE(getOption("PeakForgeR.debug", FALSE))) {
    log$write("DEBUG", sprintf("Loaded mzML: spectra=%d chromatograms=%d", ns, nc))
  }

  # SRM/MRM chromatogram mzML: use existing traces
  if (ns == 0 && nc > 0) return(exp)

  # Spectra-based mzML: OpenSWATH-style extraction
  extractor <- om$ChromatogramExtractor()
  sa <- om$SpectrumAccessOpenMS(exp)

  out_ptrs <- reticulate::py_eval("[]")
  coords   <- reticulate::py_eval("[]")

  # Avoid RT-limited extraction unless rt_extraction_window > 0
  rt_win_for_prepare <- if (is.numeric(rt_extraction_window) && rt_extraction_window > 0) {
    as.numeric(rt_extraction_window)
  } else {
    -1.0
  }

  om$ChromatogramExtractor$prepare_coordinates(
    out_ptrs, coords, targeted_exp,
    rt_win_for_prepare,
    FALSE, 0L
  )

  extractor$extractChromatograms(
    sa, out_ptrs, coords,
    as.numeric(mz_extraction_window),
    as.logical(ppm),
    0.0,
    filter
  )

  helper <- om$OpenSwathDataAccessHelper()
  ms_chroms <- vector("list", reticulate::py_len(out_ptrs))
  for (i in seq_along(ms_chroms)) {
    ch <- om$MSChromatogram()
    helper$convertToOpenMSChromatogram(out_ptrs[[i]], ch)
    ms_chroms[[i]] <- ch
  }

  out_exp <- om$MSExperiment()
  out_exp$setChromatograms(ms_chroms)
  out_exp
}

# ----------------------- PeakPicker float arrays ------------------------------



extract_float_arrays_safe <- function(picked_chrom, ch, n_expected, rt_scale) {
  assert_has_package("reticulate")
  na_out <- list(area = rep(NA_real_, n_expected),
                 left = rep(NA_real_, n_expected),
                 right = rep(NA_real_, n_expected))

  # Helper: robust numeric conversion for FloatDataArray / IntegerDataArray
  to_num <- function(arr) {
    # numpy path
    if (reticulate::py_has_attr(arr, "to_numpy")) {
      v_np <- try(arr$to_numpy(), silent = TRUE)
      if (!inherits(v_np, "try-error")) {
        v <- try(reticulate::py_to_r(v_np), silent = TRUE)
        if (!inherits(v, "try-error") && is.numeric(v)) return(as.numeric(v))
      }
    }
    # direct conversion
    v <- try(reticulate::py_to_r(arr), silent = TRUE)
    if (!inherits(v, "try-error") && is.numeric(v)) return(as.numeric(v))

    # sequence iteration fallback
    if (reticulate::py_has_attr(arr, "__len__")) {
      n <- reticulate::py_len(arr)
      out <- numeric(n)
      for (i in seq_len(n)) {
        xi <- try(arr[[i]], silent = TRUE)
        out[i] <- suppressWarnings(as.numeric(reticulate::py_to_r(xi)))
      }
      return(out)
    }
    NULL
  }

  get_names <- function(arrs) {
    vapply(arrs, function(a) {
      nm <- try(a$getName(), silent = TRUE)
      if (inherits(nm, "try-error")) NA_character_ else tolower(as.character(nm))
    }, character(1))
  }

  # Attempt to pull arrays from picked chromatogram
  fda <- picked_chrom$getFloatDataArrays()
  ida <- picked_chrom$getIntegerDataArrays()

  f_names <- if (length(fda)) get_names(fda) else character(0)
  i_names <- if (length(ida)) get_names(ida) else character(0)

  # Selectors
  idx_area_f  <- if (length(fda)) which(grepl("area|abundance|integrated|integral", f_names)) else integer(0)
  idx_left_f  <- if (length(fda)) which(grepl("left|rt_left|start", f_names)) else integer(0)
  idx_right_f <- if (length(fda)) which(grepl("right|rt_right|end", f_names)) else integer(0)

  idx_left_i  <- if (length(ida)) which(grepl("left|lbound|start", i_names)) else integer(0)
  idx_right_i <- if (length(ida)) which(grepl("right|rbound|end", i_names)) else integer(0)

  # Extract peaks (apex list) from picked
  pk <- picked_chrom$get_peaks()
  rt_apex <- as.numeric(pk[[1]])         # in raw units (likely seconds)
  ht_apex <- as.numeric(pk[[2]])
  n_peaks <- length(rt_apex)
  if (n_peaks != n_expected) {
    # align later by min(n_peaks, n_expected)
    n_expected <- n_peaks
  }
  rt_apex_min <- rt_apex * rt_scale

  # Raw chromatogram points (for fallback area/bounds)
  pk_raw <- ch$get_peaks()
  rt_raw <- as.numeric(pk_raw[[1]])
  int_raw <- as.numeric(pk_raw[[2]])
  rt_raw_min <- rt_raw * rt_scale

  # --- Initialize outputs
  area  <- rep(NA_real_, n_expected)
  left  <- rep(NA_real_, n_expected)
  right <- rep(NA_real_, n_expected)

  # --- First: use float arrays if present
  if (length(idx_area_f)) {
    a <- to_num(fda[[idx_area_f[1]]])
    if (!is.null(a)) area <- align_vec(a, n_expected)
  }
  if (length(idx_left_f)) {
    l <- to_num(fda[[idx_left_f[1]]])
    if (!is.null(l)) left <- align_vec(l, n_expected) * rt_scale
  }
  if (length(idx_right_f)) {
    r <- to_num(fda[[idx_right_f[1]]])
    if (!is.null(r)) right <- align_vec(r, n_expected) * rt_scale
  }

  # --- Second: if float boundaries missing, use integer arrays (indices)
  if (all(is.na(left)) && length(idx_left_i)) {
    l <- to_num(ida[[idx_left_i[1]]])
    if (!is.null(l)) {
      l <- pmin(pmax(as.integer(l), 1L), n_expected)
      # map index → apex times (or to raw times?)
      # OpenMS often stores boundaries as indices on the "picked peaks" list:
      # we convert indices to apex RT in minutes; downstream code maps indices to actual RT.
      left <- rt_apex_min[l]
    }
  }
  if (all(is.na(right)) && length(idx_right_i)) {
    r <- to_num(ida[[idx_right_i[1]]])
    if (!is.null(r)) {
      r <- pmin(pmax(as.integer(r), 1L), n_expected)
      right <- rt_apex_min[r]
    }
  }

  # --- Third: compute boundaries & area from raw chromatogram if arrays are missing
  need_bounds <- any(is.na(left)) || any(is.na(right))
  need_area   <- any(is.na(area))

  if (need_bounds || need_area) {
    # Precompute nearest raw index per apex
    apex_raw_idx <- vapply(seq_len(n_expected), function(i) {
      # locate nearest raw time to apex time
      which.min(abs(rt_raw_min - rt_apex_min[i]))
    }, integer(1))

    # Boundary finder: try height threshold then local valley fallback
    find_left_idx <- function(int, apex_idx, peak_h, frac = 0.1) {
      # threshold method
      thr <- peak_h * frac
      cand <- which(int[seq_len(apex_idx)] <= thr)
      if (length(cand)) return(max(cand))
      # valley: walk left while intensity decreases
      i <- apex_idx
      while (i > 1 && int[i-1] <= int[i]) i <- i - 1
      i
    }

    find_right_idx <- function(int, apex_idx, peak_h, frac = 0.1, nmax = length(int)) {
      thr <- peak_h * frac
      cand <- which(int[apex_idx:nmax] <= thr)
      if (length(cand)) return(apex_idx - 1 + min(cand))
      # valley: walk right while intensity decreases
      i <- apex_idx
      while (i < nmax && int[i+1] <= int[i]) i <- i + 1
      i
    }

    # trapezoidal integration
    trapz <- function(x, y) {
      if (length(x) < 2) return(0)
      sum((x[-1] - x[-length(x)]) * (y[-1] + y[-length(y)]) / 2)
    }

    for (i in seq_len(n_expected)) {
      idx <- apex_raw_idx[i]
      ph  <- int_raw[idx]

      # Fill boundaries if missing
      if (is.na(left[i]) || is.na(right[i])) {
        li <- find_left_idx(int_raw, idx, ph, frac = 0.1)
        ri <- find_right_idx(int_raw, idx, ph, frac = 0.1, nmax = length(int_raw))
        # guard ordering
        if (ri < li) { tmp <- li; li <- ri; ri <- tmp }
        # Convert to minutes
        left[i]  <- rt_raw_min[li]
        right[i] <- rt_raw_min[ri]
      }

      # Fill area if missing
      if (is.na(area[i])) {
        li <- max(1L, which.min(abs(rt_raw_min - left[i])))
        ri <- min(length(rt_raw_min), which.min(abs(rt_raw_min - right[i])))
        if (ri > li) {
          area[i] <- trapz(rt_raw_min[li:ri], int_raw[li:ri])
        } else {
          area[i] <- ph # fallback: use apex height if degenerate
        }
      }
    }
  }

  list(area = area, left = left, right = right)
}

# helper to align vector lengths
align_vec <- function(x, n) {
  if (is.null(x) || !length(x)) return(rep(NA_real_, n))
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == n) return(x)
  if (length(x) > n) return(x[seq_len(n)])
  c(x, rep(NA_real_, n - length(x)))
}



# ----------------------- Peak selection --------------------------------------

get_callable_or_attr <- function(obj, name) {
  # returns either a callable (function) or attribute object
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("reticulate required")
  if (!reticulate::py_has_attr(obj, name)) return(NULL)
  obj[[name]]
}

get_mz_from_precursor_like <- function(x) {
  if (is.null(x)) return(NA_real_)
  tryCatch({
    if (requireNamespace("reticulate", quietly = TRUE) && reticulate::py_has_attr(x, "getMZ")) {
      suppressWarnings(as.numeric(x$getMZ()))
    } else {
      suppressWarnings(as.numeric(reticulate::py_to_r(x)))
    }
  }, error = function(e) NA_real_)
}


unwrap_pyobj <- function(x) {
  # If we accidentally have list(<python object>), unwrap it
  if (is.list(x) && length(x) == 1) return(x[[1]])
  x
}

get_q1q3_safe <- function(ch) {
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("reticulate required")

  ch <- unwrap_pyobj(ch)

  if (!reticulate::is_py_object(ch)) {
    stop(sprintf(
      "get_q1q3_safe() expected a pyOpenMS chromatogram object.\nGot class: %s\nTip: use ch <- chroms[[i]] not ch <- chroms",
      paste(class(ch), collapse = ", ")
    ), call. = FALSE)
  }

  nid <- tryCatch(as.character(ch$getNativeID()), error = function(e) "")

  # BEST for your files: parse Q1/Q3 from the nativeID (you have Q1=... Q3=...)
  pq <- parse_q1q3_from_nativeid(nid)
  q1 <- pq$q1
  q3 <- pq$q3

  # Fallback: use precursor/product objects if parsing fails
  if (!is.finite(q1) || q1 <= 0) {
    q1 <- tryCatch(suppressWarnings(as.numeric(ch$getPrecursor()$getMZ())), error = function(e) NA_real_)
  }
  if (!is.finite(q3) || q3 <= 0) {
    q3 <- tryCatch(suppressWarnings(as.numeric(ch$getProduct()$getMZ())), error = function(e) NA_real_)
  }

  # Fallback: SRM sometimes stores product m/z directly in chromatogram MZ
  if (!is.finite(q3) || q3 <= 0) {
    q3 <- tryCatch(suppressWarnings(as.numeric(ch$getMZ())), error = function(e) NA_real_)
  }

  list(q1 = q1, q3 = q3, nid = nid)
}


choose_best_lenient <- function(df) {
  df <- as.data.frame(df)

  score_vec <- function(d) {
    if ("Area" %in% names(d) && any(is.finite(d$Area))) as.numeric(d$Area) else as.numeric(d$Height)
  }

  have_window <- all(c("rt_start_min", "rt_end_min") %in% names(df)) &&
    any(is.finite(df$rt_start_min)) &&
    any(is.finite(df$rt_end_min))

  pick_max_score <- function(d) {
    sc <- score_vec(d)
    ht <- if ("Height" %in% names(d)) as.numeric(d$Height) else rep(NA_real_, nrow(d))
    o <- order(-sc, -ht, na.last = TRUE)
    d[o[1], , drop = FALSE]
  }

  if (have_window) {
    inside <- is.finite(df$rt_start_min) & is.finite(df$rt_end_min) &
      df$RetentionTime >= df$rt_start_min & df$RetentionTime <= df$rt_end_min
    din <- df[inside, , drop = FALSE]
    if (nrow(din) > 0) return(pick_max_score(din))

    if (any(is.finite(df$rt_min))) {
      exp_rt <- df$rt_min[which(is.finite(df$rt_min))[1]]
      df$rt_delta <- abs(df$RetentionTime - exp_rt)
      sc <- score_vec(df)
      ht <- as.numeric(df$Height)
      o <- order(df$rt_delta, -sc, -ht, na.last = TRUE)
      return(df[o[1], , drop = FALSE])
    }
    return(pick_max_score(df))
  }

  if (any(is.finite(df$rt_min))) {
    exp_rt <- df$rt_min[which(is.finite(df$rt_min))[1]]
    df$rt_delta <- abs(df$RetentionTime - exp_rt)
    sc <- score_vec(df)
    ht <- as.numeric(df$Height)
    o <- order(df$rt_delta, -sc, -ht, na.last = TRUE)
    return(df[o[1], , drop = FALSE])
  }

  pick_max_score(df)
}



compute_results_openms <- function(
    plateIDs,
    user_name,
    project_directory,
    mrm_template_list,
    QC_sample_label,
    template_vote_n_files = 1L,
    cfg = PF_CFG
) {
  assert_has_package("readr")
  assert_has_package("tools")
  assert_has_package("reticulate")

  # Optional parallel packages (only needed if cfg$parallel_files = TRUE)
  parallel_enabled <- isTRUE(cfg$parallel_files %||% FALSE)
  if (parallel_enabled) {
    if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Parallel processing requested but packages 'future' and 'future.apply' are not installed.", call. = FALSE)
    }
  }
  n_workers <- as.integer(cfg$n_workers %||% max(1L, parallel::detectCores() - 1L))

  # Optional notifications (progress bar + per-file messages)
  notify_enabled <- isTRUE(cfg$notify %||% TRUE)
  if (notify_enabled && !requireNamespace("progressr", quietly = TRUE)) {
    warning("Notifications requested but package 'progressr' is not installed. Proceeding without progress UI.")
    notify_enabled <- FALSE
  }

  # Choose a handler for progressr
  if (notify_enabled) {
    handler <- tolower(cfg$notify_handler %||% "auto")
    if (handler == "rstudio") {
      progressr::handlers(global = TRUE, progressr::handler_rstudio())
    } else if (handler == "txt") {
      progressr::handlers(global = TRUE, progressr::handler_txtprogressbar())
    } else if (handler == "progress") {
      progressr::handlers(global = TRUE, progressr::handler_progress())
    } else {
      # auto: pick the best available
      progressr::handlers(global = TRUE)
    }
  }

  project_directory <- normalize_path_safe(project_directory)

  # Setup OpenMS backend once in the parent (workers will also init their own)
  peakforger_setup_openms(
    envname = cfg$envname,
    python_version = cfg$python_version,
    backend = cfg$backend,
    install = isTRUE(cfg$auto_install_env)
  )

  # allow verbose debug
  options(PeakForgeR.debug = isTRUE(cfg$debug))

  logs_dir <- file.path(project_directory, "MetaboExploreR_logs")
  dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

  results <- list()

  # ------------------------- Worker function (per file) -------------------------
  process_one_mzml <- function(
    mzml_path,
    project_directory,
    trans_df_sorted,
    tindex,
    cfg,
    mzml_all = NULL,
    p = NULL  # progressor (optional)
  ) {
    # Each worker must init pyOpenMS independently
    peakforger_setup_openms(
      envname = cfg$envname,
      python_version = cfg$python_version,
      backend = cfg$backend,
      install = FALSE  # do NOT install envs in workers
    )
    om <- openms_get()

    # Rebuild TargetedExperiment inside worker (safe for reticulate)
    targeted_exp <- build_targeted_experiment_pf(trans_df_sorted)

    file_stat <- list(
      file = basename(mzml_path), spectra = NA_integer_,
      chromatograms = NA_integer_, transition_like = 0L,
      matched = 0L, peaks = 0L
    )

    # Resolve index for nice "[i/m]" messages
    i <- NA_integer_
    m <- NA_integer_
    if (!is.null(mzml_all)) {
      i <- match(mzml_path, mzml_all)
      m <- length(mzml_all)
    }

    res <- tryCatch({
      # Notify START (amount = 0 means "update message; do not advance progress")
      if (!is.null(p)) p(sprintf("⇢ [%s/%s] Start: %s", i %||% "?", m %||% "?", basename(mzml_path)), amount = 0)

      exp <- om$MSExperiment()
      om$MzMLFile()$load(mzml_path, exp)
      file_stat$spectra <- exp$getNrSpectra()
      file_stat$chromatograms <- exp$getNrChromatograms()

      # Extract (or reuse SRM)
      chrom_exp <- extract_all_chromatograms(
        mzml_path = mzml_path,
        targeted_exp = targeted_exp,
        mz_extraction_window = cfg$mz_extraction_window,
        ppm = cfg$ppm,
        rt_extraction_window = cfg$rt_extraction_window,
        filter = cfg$filter,
        log = NULL # avoid shared logger; messages are returned
      )

      chroms_all <- chrom_exp$getChromatograms()
      if (!length(chroms_all)) {
        # Advance progress (file done, even though empty)
        if (!is.null(p)) p(sprintf("⚠ [%s/%s] No chromatograms: %s", i %||% "?", m %||% "?", basename(mzml_path)))
        return(list(best_out = NULL, file_stat = file_stat, messages = paste0("WARN: No chromatograms: ", basename(mzml_path))))
      }

      chroms <- Filter(is_transition_like_chrom, chroms_all)
      file_stat$transition_like <- length(chroms)
      if (!length(chroms)) {
        if (!is.null(p)) p(sprintf("⚠ [%s/%s] No transition-like: %s", i %||% "?", m %||% "?", basename(mzml_path)))
        return(list(best_out = NULL, file_stat = file_stat, messages = paste0("WARN: No transition-like chromatograms: ", basename(mzml_path))))
      }

      # RT scaling (minutes)
      mzml_rt_unit <- infer_mzml_time_unit(mzml_path)
      rt_scale <- 1
      if (identical(mzml_rt_unit, "sec")) rt_scale <- 1/60
      if (is.na(mzml_rt_unit)) {
        probe <- try(as.numeric(chroms[[1]]$getMaxRT()), silent = TRUE)
        if (!inherits(probe, "try-error") && is.finite(probe) && probe > 2000) rt_scale <- 1/60
      }

      picker <- om$PeakPickerChromatogram()
      rows_for_file <- list()
      matched_ct <- 0L
      peak_ct <- 0L

      sample_name <- tools::file_path_sans_ext(basename(mzml_path))
      acquired_time <- get_mzml_acquired_time(mzml_path)

      for (ii in seq_along(chroms)) {
        ch <- chroms[[ii]]
        qq <- get_q1q3_safe(ch)
        q1 <- qq$q1; q3 <- qq$q3
        if (!(is.finite(q1) && q1 > 0 && is.finite(q3) && q3 > 0)) next

        mi <- find_transition_match(tindex, trans_df_sorted, q1, q3)
        if (!is.finite(mi)) next
        meta <- trans_df_sorted[mi, , drop = FALSE]
        matched_ct <- matched_ct + 1L

        picked <- om$MSChromatogram()
        picker$pickChromatogram(ch, picked)
        pk <- picked$get_peaks()
        rt <- as.numeric(pk[[1]])
        intens <- as.numeric(pk[[2]])
        if (!length(rt)) next

        rt_min <- rt * rt_scale
        n_peaks <- length(rt_min)
        peak_ct <- peak_ct + n_peaks

        # Updated extractor that understands IntegratedIntensity and left/right widths
        fa <- extract_float_arrays_safe(
          picked_chrom = picked,
          n_expected   = n_peaks,
          rt_scale     = rt_scale,
          ch           = ch
        )
        area  <- fa$area
        left  <- fa$left
        right <- fa$right

        df <- data.frame(
          FileName         = sample_name,
          MoleculeListName = meta$molecule_list_name[[1]],
          MoleculeName     = meta$molecule_name[[1]],
          PrecursorMz      = meta$precursor_mz[[1]],
          ProductMz        = meta$product_mz[[1]],
          RetentionTime    = rt_min,
          StartTime        = left,
          EndTime          = right,
          Area             = area,
          Height           = intens,
          rt_min           = meta$rt_min[[1]],
          rt_start_min     = meta$rt_start_min[[1]],
          rt_end_min       = meta$rt_end_min[[1]],
          stringsAsFactors = FALSE
        )
        rows_for_file[[length(rows_for_file) + 1]] <- df
      }

      file_stat$matched <- matched_ct
      file_stat$peaks <- peak_ct

      all_peaks <- if (length(rows_for_file)) do.call(rbind, rows_for_file) else data.frame()
      if (!nrow(all_peaks)) {
        if (!is.null(p)) p(sprintf("⚠ [%s/%s] No peaks: %s", i %||% "?", m %||% "?", basename(mzml_path)))
        return(list(best_out = NULL, file_stat = file_stat, messages = paste0("WARN: No peaks after picking/matching: ", basename(mzml_path))))
      }

      split_by <- split(all_peaks, paste(
        all_peaks$MoleculeListName, all_peaks$MoleculeName,
        all_peaks$PrecursorMz, all_peaks$ProductMz, sep="|"
      ))
      best <- do.call(rbind, lapply(split_by, choose_best_lenient))
      rownames(best) <- NULL

      best_out <- data.frame(
        FileName         = best$FileName,
        MoleculeListName = best$MoleculeListName,
        MoleculeName     = best$MoleculeName,
        PrecursorMz      = best$PrecursorMz,
        ProductMz        = best$ProductMz,
        RetentionTime    = best$RetentionTime,
        StartTime        = best$StartTime,
        EndTime          = best$EndTime,
        Area             = best$Area,
        Height           = best$Height,
        AcquiredTime     = acquired_time,
        stringsAsFactors = FALSE
      )

      # Notify DONE (advance progress by +1 for this file)
      if (!is.null(p)) p(sprintf("✔ [%s/%s] Done: %s (rows=%d)", i %||% "?", m %||% "?", basename(mzml_path), nrow(best_out)), amount = 1)

      list(
        best_out = best_out,
        file_stat = file_stat,
        messages = sprintf("OK: %s matched=%d rows=%d", basename(mzml_path), matched_ct, nrow(best_out))
      )
    }, error = function(e) {
      if (!is.null(p)) p(sprintf("✖ [%s/%s] ERROR: %s :: %s", i %||% "?", m %||% "?", basename(mzml_path), conditionMessage(e)), amount = 1)
      list(best_out = NULL, file_stat = file_stat, messages = paste("ERROR:", conditionMessage(e)))
    })

    res
  }


  for (plateID in plateIDs) {
    log_file <- file.path(logs_dir, paste0(plateID, "_PeakForgeR_log.txt"))
    log <- make_logger(log_file, max_mb = cfg$max_log_mb)

    plate_res <- list(success = FALSE, plateID = plateID, out_path = NA_character_,
                      n_rows = 0, file_stats = list(), errors = character())

    # Use a local plan and restore at the end to avoid leaking global plan changes
    old_plan <- NULL
    if (parallel_enabled) {
      old_plan <- future::plan()
      future::plan(future::multisession, workers = max(1L, n_workers))
    }

    tryCatch({
      log$write("INFO", paste0("=== Plate start: ", plateID, " ==="))

      # Step 1/2: project + mzML listing
      master_list <- PeakForgeR_setup_project(user_name, project_directory, plateID, mrm_template_list, QC_sample_label)
      master_list <- import_mzml(plateID, master_list)

      mzml_dir <- normalize_path_safe(file.path(project_directory, plateID, "data", "mzml"))
      mzml_files <- list.files(mzml_dir, pattern = "\\.mzML(\\.gz)?$", ignore.case = TRUE, full.names = TRUE)
      if (!length(mzml_files)) stop("No mzML files found for plate: ", plateID, call. = FALSE)

      # Step 2b: template selection
      template_names <- setdiff(names(master_list$templates$mrm_guides), "by_plate")
      if (length(template_names) > 1) {
        master_list <- select_method_template_for_plate(
          master_list = master_list,
          plateID = plateID,
          project_directory = project_directory,
          prefer_n_files = template_vote_n_files,
          qc_label = QC_sample_label,
          verbose = cfg$debug
        )
      } else if (length(template_names) == 1) {
        master_list$project_details$is_ver <- template_names[1]
      } else {
        stop("No method templates found in master_list$templates$mrm_guides.", call. = FALSE)
      }
      log$write("INFO", paste0("Template: ", master_list$project_details$is_ver))

      # Step 3: RT optimisation per plate
      plate_idx <- plateID
      if (is.null(master_list$templates$mrm_guides$by_plate) || !is.list(master_list$templates$mrm_guides$by_plate)) {
        master_list$templates$mrm_guides$by_plate <- list()
      }
      master_list$templates$mrm_guides$by_plate[[plate_idx]] <- optimise_retention_times(master_list, plate_idx)

      targets_raw <- master_list$templates$mrm_guides$by_plate[[plate_idx]]
      if (is.null(targets_raw)) {
        targets_raw <- master_list$templates$mrm_guides[[master_list$project_details$is_ver]]$mrm_guide
      }

      # Standardize transitions and build TE spec
      trans_df <- standardize_transition_table_lipidyzer(targets_raw)
      trans_df <- trans_df[!duplicated(trans_df$transition_name), , drop = FALSE]
      trans_df_sorted <- sort_transitions_openms_like(trans_df)

      # Build fast transition index (Q1/Q3 matching)
      tindex <- build_transition_index(trans_df_sorted, cfg$mz_tol_q1, cfg$mz_tol_q3)

      out_dir <- file.path(project_directory, plateID, "data", "PeakForgeR")
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      out_path <- file.path(out_dir, paste0(plateID, "_PeakForgeR_OpenMS.tsv"))

      # ------------------- Map over mzML files (with progress UI) -------------------
      do_map <- function() {
        if (parallel_enabled) {
          future.apply::future_lapply(
            mzml_files,
            FUN = process_one_mzml,
            project_directory = project_directory,
            trans_df_sorted   = trans_df_sorted,
            tindex            = tindex,
            cfg               = cfg,
            mzml_all          = mzml_files,
            p                 = if (notify_enabled) p else NULL,
            future.seed       = TRUE
          )
        } else {
          lapply(
            mzml_files,
            FUN = process_one_mzml,
            project_directory = project_directory,
            trans_df_sorted   = trans_df_sorted,
            tindex            = tindex,
            cfg               = cfg,
            mzml_all          = mzml_files,
            p                 = if (notify_enabled) p else NULL
          )
        }
      }

      if (notify_enabled) {
        # Progressor with one step per file; we update "start" with amount=0 and "finish" with amount=1
        p <- progressr::progressor(steps = length(mzml_files))
        mapped <- progressr::with_progress({
          do_map()
        })
      } else {
        mapped <- do_map()
      }

      # Collect results + write messages to plate log (sequentially)
      all_rows <- list()
      for (k in seq_along(mapped)) {
        m <- mapped[[k]]
        if (!is.null(m$messages) && nzchar(m$messages)) log$write("INFO", m$messages)
        if (!is.null(m$best_out)) all_rows[[length(all_rows) + 1]] <- m$best_out
        plate_res$file_stats[[basename(mzml_files[k])]] <- m$file_stat
      }

      plate_table <- if (length(all_rows)) do.call(rbind, all_rows) else data.frame()

      required_cols <- c("FileName","MoleculeListName","MoleculeName",
                         "PrecursorMz","ProductMz",
                         "RetentionTime","StartTime","EndTime",
                         "Area","Height","AcquiredTime")

      if (!nrow(plate_table)) {
        plate_table <- setNames(
          as.data.frame(replicate(length(required_cols), character(0), simplify = FALSE),
                        stringsAsFactors = FALSE),
          required_cols
        )
      } else {
        num_cols <- c("PrecursorMz","ProductMz","RetentionTime","StartTime","EndTime","Area","Height")
        for (cc in intersect(num_cols, names(plate_table))) plate_table[[cc]] <- as.numeric(plate_table[[cc]])
        plate_table <- plate_table[, required_cols, drop = FALSE]
      }

      readr::write_tsv(plate_table, out_path)
      log$write("INFO", paste0("Wrote: ", out_path, " rows=", nrow(plate_table)))

      plate_res$success <- TRUE
      plate_res$out_path <- out_path
      plate_res$n_rows <- nrow(plate_table)

    }, error = function(e) {
      plate_res$errors <- c(plate_res$errors, conditionMessage(e))
      log$write("ERROR", conditionMessage(e))
    }, finally = {
      log$write("INFO", paste0("=== Plate end: ", plateID, " success=", plate_res$success, " ==="))
      log$close()
      # Restore previous future plan if we changed it
      if (parallel_enabled) {
        if (!is.null(old_plan)) {
          future::plan(old_plan)
        } else {
          future::plan(future::sequential)
        }
      }
      results[[plateID]] <- plate_res
    })
  }

  results
}


#' Discover Visium sample metadata on disk
#'
#' Automatically scans a root directory for 10x Genomics Visium
#' CellRanger outputs and constructs a per-sample metadata table
#' describing file locations and sample annotations.
#'
#' Each sample directory is expected to contain:
#' \itemize{
#'   \item filtered_feature_bc_matrix.h5
#'   \item tissue_positions*.csv
#'   \item scalefactors_json.json
#'   \item tissue_hires_image.png
#'   \item tissue_lowres_image.png
#' }
#'
#' @param root_dir Directory containing per-sample Visium output folders.
#' @param sample_dirs Optional character vector of sample subdirectories
#'   (relative to \code{root_dir} or full paths). If NULL, all first-level
#'   subdirectories of \code{root_dir} are used.
#' @param sample_ids Optional character vector of full sample IDs.
#'   Defaults to directory names.
#' @param sample_labels Optional vector of colloquial sample labels.
#'   Defaults to \code{sample_ids}.
#' @param sample_groups Optional vector of grouping labels.
#'   Defaults to \code{sample_labels}.
#' @param barcode_suffix Optional vector of barcode suffixes (for example "_1").
#'   Defaults to sequential suffixes.
#'
#' @return A tibble with one row per sample and columns:
#' \describe{
#'   \item{sample_id}{Full sample identifier}
#'   \item{sample_label}{Colloquial sample label}
#'   \item{sample_group}{Grouping label}
#'   \item{h5_file}{Path to filtered_feature_bc_matrix.h5}
#'   \item{positions_csv}{Path to tissue_positions CSV}
#'   \item{scalefactors_json}{Path to scalefactors JSON}
#'   \item{hires_image}{Path to high-resolution tissue image}
#'   \item{lowres_image}{Path to low-resolution tissue image}
#'   \item{barcode_suffix}{Barcode suffix used for this sample}
#' }
#'
#' @export
discover_visium_sample_table <- function(
    root_dir,
    sample_dirs    = NULL,
    sample_ids     = NULL,
    sample_labels  = NULL,
    sample_groups  = NULL,
    barcode_suffix = NULL
) {

  if (!dir.exists(root_dir)) {
    stop("root_dir does not exist: ", root_dir)
  }

  ## ------------------------------------------------------------
  ## Determine sample directories
  ## ------------------------------------------------------------
  if (is.null(sample_dirs)) {
    sample_paths <- list.dirs(
      path       = root_dir,
      full.names = TRUE,
      recursive  = FALSE
    )
  } else {
    sample_paths <- ifelse(
      dir.exists(sample_dirs),
      sample_dirs,
      file.path(root_dir, sample_dirs)
    )
  }

  if (length(sample_paths) == 0) {
    stop("No sample directories found under: ", root_dir)
  }

  n <- length(sample_paths)

  ## ------------------------------------------------------------
  ## Sample identifiers
  ## ------------------------------------------------------------
  if (is.null(sample_ids)) {
    sample_ids <- basename(sample_paths)
  } else if (length(sample_ids) != n) {
    stop("sample_ids must have length equal to number of samples")
  }

  if (is.null(sample_labels)) {
    sample_labels <- sample_ids
  } else if (length(sample_labels) != n) {
    stop("sample_labels must have length equal to number of samples")
  }

  if (is.null(sample_groups)) {
    sample_groups <- sample_labels
  } else if (length(sample_groups) != n) {
    stop("sample_groups must have length equal to number of samples")
  }

  if (is.null(barcode_suffix)) {
    barcode_suffix <- paste0("_", seq_len(n))
  } else if (length(barcode_suffix) != n) {
    stop("barcode_suffix must have length equal to number of samples")
  }

  ## ------------------------------------------------------------
  ## Discover per-sample files
  ## ------------------------------------------------------------
  h5_files          <- character(n)
  positions_files   <- character(n)
  scalefactors_json <- character(n)
  hires_imgs        <- character(n)
  lowres_imgs       <- character(n)

  for (i in seq_len(n)) {
    sp <- sample_paths[i]

    h5_files[i] <- list.files(
      sp,
      pattern     = "filtered_feature_bc_matrix\\.h5$",
      full.names = TRUE,
      recursive  = TRUE
    )[1]

    positions_files[i] <- list.files(
      sp,
      pattern     = "tissue_positions.*\\.csv$",
      full.names = TRUE,
      recursive  = TRUE
    )[1]

    scalefactors_json[i] <- list.files(
      sp,
      pattern     = "scalefactors_json\\.json$",
      full.names = TRUE,
      recursive  = TRUE
    )[1]

    hires_imgs[i] <- list.files(
      sp,
      pattern     = "tissue_hires_image\\.png$",
      full.names = TRUE,
      recursive  = TRUE
    )[1]

    lowres_imgs[i] <- list.files(
      sp,
      pattern     = "tissue_lowres_image\\.png$",
      full.names = TRUE,
      recursive  = TRUE
    )[1]
  }

  ## ------------------------------------------------------------
  ## Assemble result table
  ## ------------------------------------------------------------
  tibble::tibble(
    sample_id         = sample_ids,
    sample_label      = sample_labels,
    sample_group      = sample_groups,
    h5_file           = h5_files,
    positions_csv     = positions_files,
    scalefactors_json = scalefactors_json,
    hires_image       = hires_imgs,
    lowres_image      = lowres_imgs,
    barcode_suffix    = barcode_suffix
  )
}

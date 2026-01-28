#' Build Visium image pixel data frames aligned to spatial coordinates
#'
#' Reads high- and low-resolution H&E images for each Visium sample
#' and converts them into pixel-level data frames aligned to spatial
#' spot coordinates in a monocle3 cell_data_set.
#'
#' @param cds A monocle3 cell_data_set with spatial coordinates in colData.
#' @param sample_table A data.frame produced by discover_visium_sample_table().
#' @param img_buffer Pixel buffer around spots used to flag image regions
#'   containing tissue.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{highres}{Pixel-level data.frame for high-resolution images}
#'     \item{lowres}{Pixel-level data.frame for low-resolution images}
#'   }
#'
#' @export
build_visium_images <- function(
    cds,
    sample_table,
    img_buffer = 5
) {

  required_cols <- c(
    "x_highres", "y_highres",
    "x_lowres", "y_lowres",
    "sample"
  )

  if (!all(required_cols %in% colnames(SummarizedExperiment::colData(cds)))) {
    stop(
      "cds must contain spatial coordinates in colData: ",
      paste(required_cols, collapse = ", "),
      ". Did you run attach_spatial_to_cds()?"
    )
  }

  df_spots <- as.data.frame(SummarizedExperiment::colData(cds))
  df_spots$sample <- as.character(df_spots$sample)

  hi_list <- list()
  lo_list <- list()

  for (i in seq_len(nrow(sample_table))) {

    sample_id <- as.character(sample_table$sample_id[i])
    df_samp   <- df_spots[df_spots$sample == sample_id, , drop = FALSE]

    if (nrow(df_samp) == 0) next

    ## Bounding boxes for spot coordinates
    xh_rng <- range(df_samp$x_highres, na.rm = TRUE)
    yh_rng <- range(df_samp$y_highres, na.rm = TRUE)
    xl_rng <- range(df_samp$x_lowres,  na.rm = TRUE)
    yl_rng <- range(df_samp$y_lowres,  na.rm = TRUE)

    ## ---------------- High resolution image ----------------
    hi_file <- sample_table$hires_image[i]
    if (!is.na(hi_file) && file.exists(hi_file)) {

      img_hi <- imager::load.image(hi_file)

      df_hi <- as.data.frame(img_hi, wide = "c")
      df_hi$rgb.val <- rgb(df_hi$c.1, df_hi$c.2, df_hi$c.3)

      df_hi$sample       <- sample_id
      df_hi$sample_label <- sample_table$sample_label[i]
      df_hi$sample_group <- sample_table$sample_group[i]
      df_hi$resolution   <- "highres"

      df_hi$contains_spots <- ifelse(
        df_hi$x >= (xh_rng[1] - img_buffer) &
          df_hi$x <= (xh_rng[2] + img_buffer) &
          df_hi$y >= (yh_rng[1] - img_buffer) &
          df_hi$y <= (yh_rng[2] + img_buffer),
        "yes", "no"
      )

      hi_list[[length(hi_list) + 1]] <- df_hi
    }

    ## ---------------- Low resolution image ----------------
    lo_file <- sample_table$lowres_image[i]
    if (!is.na(lo_file) && file.exists(lo_file)) {

      img_lo <- imager::load.image(lo_file)

      df_lo <- as.data.frame(img_lo, wide = "c")
      df_lo$rgb.val <- rgb(df_lo$c.1, df_lo$c.2, df_lo$c.3)

      df_lo$sample       <- sample_id
      df_lo$sample_label <- sample_table$sample_label[i]
      df_lo$sample_group <- sample_table$sample_group[i]
      df_lo$resolution   <- "lowres"

      df_lo$contains_spots <- ifelse(
        df_lo$x >= (xl_rng[1] - img_buffer) &
          df_lo$x <= (xl_rng[2] + img_buffer) &
          df_lo$y >= (yl_rng[1] - img_buffer) &
          df_lo$y <= (yl_rng[2] + img_buffer),
        "yes", "no"
      )

      lo_list[[length(lo_list) + 1]] <- df_lo
    }
  }

  list(
    highres = if (length(hi_list) > 0) dplyr::bind_rows(hi_list) else NULL,
    lowres  = if (length(lo_list) > 0) dplyr::bind_rows(lo_list) else NULL
  )
}


compute_spatial_xy_ratio <- function(
    cds,
    sample_id,
    use_fullres = TRUE
) {

  df <- as.data.frame(SummarizedExperiment::colData(cds))
  df <- df[df$sample == sample_id, , drop = FALSE]

  if (nrow(df) == 0) {
    stop("No spots found for sample_id: ", sample_id)
  }

  if (use_fullres) {
    required <- c("pxl_col_in_fullres", "pxl_row_in_fullres")
  } else {
    required <- c("x_highres", "y_highres")
  }

  if (!all(required %in% colnames(df))) {
    stop(
      "Required spatial columns not found in colData(cds): ",
      paste(required, collapse = ", ")
    )
  }

  if (use_fullres) {
    x_rng <- range(df$pxl_col_in_fullres, na.rm = TRUE)
    y_rng <- range(df$pxl_row_in_fullres, na.rm = TRUE)
  } else {
    x_rng <- range(df$x_highres, na.rm = TRUE)
    y_rng <- range(df$y_highres, na.rm = TRUE)
  }

  (x_rng[2] - x_rng[1]) / (y_rng[2] - y_rng[1])
}

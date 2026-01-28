#' Read Visium spatial coordinates and scale factors
#'
#' Reads the tissue_positions CSV and scalefactors_json.json
#' files produced by CellRanger Visium and computes high- and
#' low-resolution pixel coordinates.
#'
#' @param positions_csv Path to tissue_positions*.csv
#' @param scalefactors_json Path to scalefactors_json.json
#' @param barcode_suffix Optional barcode suffix (e.g. "_1")
#' @param sample_id Optional sample identifier to attach
#'
#' @return A data.frame with spatial coordinates per barcode
#' @export
read_visium_positions <- function(
    positions_csv,
    scalefactors_json,
    barcode_suffix = NULL,
    sample_id      = NULL
) {

  if (!file.exists(positions_csv)) {
    stop("positions_csv does not exist: ", positions_csv)
  }
  if (!file.exists(scalefactors_json)) {
    stop("scalefactors_json does not exist: ", scalefactors_json)
  }

  dt <- data.table::fread(positions_csv)

  if (!"barcode" %in% colnames(dt)) {
    stop("positions_csv must contain a 'barcode' column")
  }

  if (!is.null(barcode_suffix)) {
    dt$barcode <- paste0(dt$barcode, barcode_suffix)
  }

  scale_factors <- jsonlite::fromJSON(scalefactors_json)

  required_sf <- c("tissue_hires_scalef", "tissue_lowres_scalef")
  if (!all(required_sf %in% names(scale_factors))) {
    stop(
      "scalefactors_json is missing required fields: ",
      paste(required_sf, collapse = ", ")
    )
  }

  dt$y_highres <- dt$pxl_row_in_fullres * scale_factors$tissue_hires_scalef
  dt$x_highres <- dt$pxl_col_in_fullres * scale_factors$tissue_hires_scalef
  dt$y_lowres  <- dt$pxl_row_in_fullres * scale_factors$tissue_lowres_scalef
  dt$x_lowres  <- dt$pxl_col_in_fullres * scale_factors$tissue_lowres_scalef

  if (!is.null(sample_id)) {
    dt$sample <- sample_id
  }

  as.data.frame(dt)
}

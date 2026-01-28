#' Plot a spatial signature overlay for a single Visium sample
#'
#' Overlays a numeric signature score stored in colData(cds) onto a
#' Visium tissue image at high or low resolution.
#'
#' @param cds A monocle3 cell_data_set containing spatial coordinates
#'   and the signature column.
#' @param images List returned by build_visium_images().
#' @param signature_col Name of numeric column in colData(cds).
#' @param sample_id Sample identifier matching colData(cds)$sample.
#' @param resolution Image resolution: "highres" or "lowres".
#' @param point_size Size of spatial points.
#' @param point_alpha Transparency of spatial points.
#' @param img_buffer_only_spots Logical; restrict image to region
#'   containing spots.
#' @param palette Color palette type.
#' @param custom_cols Custom colors if palette = "custom".
#' @param viridis_opt Viridis palette option.
#' @param diverging_midpoint Midpoint for diverging color scale.
#'
#' @return A ggplot object.
#'
#' @export
plot_signature_spatial_sample <- function(
    cds,
    images,
    signature_col,
    sample_id,
    resolution  = c("highres", "lowres"),
    point_size  = 0.5,
    point_alpha = 0.8,
    img_buffer_only_spots = TRUE,
    palette      = c("diverging", "sequential", "viridis", "custom"),
    custom_cols  = NULL,
    viridis_opt  = "magma",
    diverging_midpoint = 0
) {

  resolution <- match.arg(resolution)
  palette    <- match.arg(palette)

  ## -------------------- sanity checks --------------------
  cd <- SummarizedExperiment::colData(cds)

  if (!signature_col %in% colnames(cd)) {
    stop("signature_col not found in colData(cds): ", signature_col)
  }

  if (!is.numeric(cd[[signature_col]])) {
    stop("signature_col must be numeric: ", signature_col)
  }

  df_spots <- as.data.frame(cd)
  df_spots <- df_spots[df_spots$sample == sample_id, , drop = FALSE]

  if (nrow(df_spots) == 0) {
    stop("No spots found for sample_id: ", sample_id)
  }

  ## -------------------- coordinates --------------------
  if (resolution == "highres") {
    x_col <- "x_highres"
    y_col <- "y_highres"
  } else {
    x_col <- "x_lowres"
    y_col <- "y_lowres"
  }

  if (!all(c(x_col, y_col) %in% colnames(df_spots))) {
    stop("Spatial coordinates missing from colData(cds).")
  }

  ## -------------------- image pixels --------------------
  image_df <- images[[resolution]]
  if (is.null(image_df)) {
    stop("No images available for resolution: ", resolution)
  }

  df_img <- image_df[image_df$sample == sample_id, , drop = FALSE]
  if (nrow(df_img) == 0) {
    stop("No image pixels found for sample_id: ", sample_id)
  }

  if (img_buffer_only_spots && "contains_spots" %in% colnames(df_img)) {
    df_img <- df_img[df_img$contains_spots == "yes", , drop = FALSE]
  }

  ## -------------------- aspect ratio --------------------
  xy_ratio <- compute_spatial_xy_ratio(
    cds       = cds,
    sample_id = sample_id,
    use_fullres = TRUE
  )

  ## -------------------- color scale --------------------
  color_scale <- switch(
    palette,
    "sequential" = ggplot2::scale_colour_gradient(
      low = "white", high = "red", na.value = "grey80"
    ),
    "diverging" = ggplot2::scale_colour_gradient2(
      low = "blue",
      mid = "grey90",
      high = "red",
      midpoint = diverging_midpoint,
      na.value = "grey80"
    ),
    "viridis" = ggplot2::scale_colour_viridis_c(
      option = viridis_opt,
      na.value = "grey80"
    ),
    "custom" = {
      if (is.null(custom_cols) || length(custom_cols) < 2) {
        stop("custom_cols must contain at least two colors.")
      }
      ggplot2::scale_colour_gradientn(
        colors = custom_cols,
        na.value = "grey80"
      )
    }
  )

  sample_label <- unique(df_spots$sample_label)

  ## -------------------- plot --------------------
  ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = df_img,
      ggplot2::aes(x = x, y = y, fill = rgb.val)
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_reverse() +
    ggplot2::geom_point(
      data = df_spots,
      ggplot2::aes(
        x = .data[[x_col]],
        y = .data[[y_col]],
        colour = .data[[signature_col]]
      ),
      size  = point_size,
      alpha = point_alpha
    ) +
    color_scale +
    ggplot2::theme_void() +
    ggplot2::theme(
      aspect.ratio = 1 / xy_ratio,
      legend.text  = ggplot2::element_text(size = 10, color = "black"),
      legend.title = ggplot2::element_text(size = 10, color = "black")
    ) +
    ggplot2::labs(
      title  = paste0("Sample: ", sample_label, " (", sample_id, ")"),
      colour = signature_col
    )
}

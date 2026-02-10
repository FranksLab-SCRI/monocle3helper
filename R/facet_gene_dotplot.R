#' Faceted gene expression dot plot grouped by an arbitrary metadata variable
#'
#' This function generates a faceted dot plot summarizing gene expression
#' across cell groups and an arbitrary cell-level metadata variable. For each
#' gene, the mean (log-transformed) expression and the percentage of cells
#' expressing the gene are computed within each combination of cell group and
#' metadata category. Dot color represents mean expression, and dot size
#' represents the fraction of cells expressing the gene.
#'
#' The function is designed for single-cell expression data stored in a
#' \code{SummarizedExperiment}-compatible object (for example, a Monocle3
#' \code{cell_data_set}) and supports flexible grouping by any column in
#' \code{colData(cds)}.
#'
#' @param cds A \code{SummarizedExperiment}-compatible object containing
#'   single-cell expression data.
#' @param markers A character vector of gene IDs or gene short names to plot.
#' @param variable A character string specifying the column in
#'   \code{colData(cds)} to use for grouping cells along the x-axis
#'   (for example, \code{"Timepoint"} or \code{"Condition"}).
#' @param group_cells_by A character string specifying the column in
#'   \code{colData(cds)} used to group cells along the y-axis.
#' @param lower_threshold Numeric value indicating the minimum expression
#'   required for a gene to be considered expressed when calculating the
#'   percentage of cells.
#' @param max.size Maximum point size used to scale the percentage of cells
#'   expressing a gene.
#' @param pseudocount Numeric pseudocount added prior to log transformation
#'   when computing mean expression.
#' @param scale_max Upper bound for clamping mean expression values used for
#'   color scaling.
#' @param scale_min Lower bound for clamping mean expression values used for
#'   color scaling.
#'
#' @return A \code{ggplot} object showing faceted dot plots for each gene.
#'
#' @examples
#' \dontrun{
#' facet_gene_dotplot_group_by_variable(
#'   cds,
#'   markers = c("TNNT2", "MYH6"),
#'   variable = "Timepoint",
#'   group_cells_by = "clusters"
#' )
#' }
#'
#' @export
facet_gene_dotplot_group_by_variable <- function(
    cds,
    markers,
    variable = "Timepoint",
    group_cells_by = "clusters",
    lower_threshold = 0,
    max.size = 10,
    pseudocount = 1,
    scale_max = 3,
    scale_min = -3
) {

  ## ----------------------------
  ## Gene IDs
  ## ----------------------------
  gene_ids <- as.data.frame(SummarizedExperiment::rowData(cds)) |>
    tibble::rownames_to_column("gene_id") |>
    dplyr::filter(
      gene_id %in% markers |
        gene_short_name %in% markers
    ) |>
    dplyr::pull(gene_id)

  ## ----------------------------
  ## Expression matrix (long)
  ## ----------------------------
  exprs_mat <- as.matrix(SummarizedExperiment::assay(cds)[gene_ids, , drop = FALSE])

  exprs_df <- as.data.frame(t(exprs_mat)) |>
    tibble::rownames_to_column("Cell") |>
    tidyr::pivot_longer(
      cols = -Cell,
      names_to = "Gene",
      values_to = "Expression"
    )

  ## ----------------------------
  ## Metadata
  ## ----------------------------
  meta <- as.data.frame(SummarizedExperiment::colData(cds))

  exprs_df[[variable]] <- meta[exprs_df$Cell, variable]
  exprs_df$Group       <- meta[exprs_df$Cell, group_cells_by]

  exprs_df <- dplyr::filter(
    exprs_df,
    !is.na(.data[[variable]]),
    !is.na(Group)
  )

  ## ----------------------------
  ## Summarize
  ## ----------------------------
  ExpVal <- exprs_df |>
    dplyr::group_by(
      Group,
      Gene,
      .data[[variable]]
    ) |>
    dplyr::summarize(
      mean = mean(log(Expression + pseudocount)),
      percentage = sum(Expression > lower_threshold) / dplyr::n(),
      .groups = "drop"
    )

  ## Clamp
  ExpVal$mean <- pmax(pmin(ExpVal$mean, scale_max), scale_min)

  ## Map gene IDs → symbols
  ExpVal$Gene <- SummarizedExperiment::rowData(cds)[ExpVal$Gene, "gene_short_name"]

  ## ----------------------------
  ## Plot
  ## ----------------------------
  ggplot2::ggplot(
    ExpVal,
    ggplot2::aes(
      y = Group,
      x = .data[[variable]]
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        colour = mean,
        size = percentage
      )
    ) +
    viridis::scale_color_viridis(
      name = paste0("log(mean + ", pseudocount, ")")
    ) +
    ggplot2::scale_size(
      name = "Percentage",
      range = c(0, max.size)
    ) +
    ggplot2::facet_wrap(~Gene) +
    ggplot2::theme_bw() +
    ggplot2::xlab(variable) +
    ggplot2::ylab(group_cells_by)
}

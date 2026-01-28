#' Compute an aggregate gene signature score
#'
#' Computes a scaled aggregate score for a gene set using
#' size-factor–normalized expression.
#'
#' @param cds monocle3 cell_data_set
#' @param genes Character vector of gene symbols
#' @param col_name Name of output column
#'
#' @return CDS with aggregate score added to colData
#' @export
aggregate_gene_score <- function(
    cds,
    genes,
    col_name = "gene_signature_score"
) {

  if ("gene_short_name" %in% colnames(SummarizedExperiment::rowData(cds))) {
    gene_names <- SummarizedExperiment::rowData(cds)$gene_short_name
  } else {
    gene_names <- rownames(cds)
  }

  gene_names <- as.character(gene_names)
  idx <- gene_names %in% genes

  if (!any(idx)) {
    SummarizedExperiment::colData(cds)[[col_name]] <- NA_real_
    return(cds)
  }

  counts_mat <- monocle3::counts(cds)[idx, , drop = FALSE]
  size_factors <- SummarizedExperiment::colData(cds)$Size_Factor

  expr <- t(t(counts_mat) / size_factors)
  expr <- log1p(expr)
  expr <- t(stats::scale(t(expr)))
  expr <- expr[rowSums(is.na(expr)) == 0, , drop = FALSE]

  SummarizedExperiment::colData(cds)[[col_name]] <- colMeans(expr)

  cds
}

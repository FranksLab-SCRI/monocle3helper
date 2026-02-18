#' Convert Gene Identifiers Using g:Profiler Orthology Mapping
#'
#' Internal helper that maps gene identifiers between species using
#' the g:Profiler \code{gorth()} service. This function provides robust
#' ortholog mapping without relying on Ensembl BioMart mirrors.
#'
#' @param genes Character vector of gene symbols.
#' @param from_organism Source organism (g:Profiler format, e.g. "hsapiens").
#' @param to_organism Target organism (g:Profiler format, e.g. "mmusculus").
#'
#' @return Character vector of mapped gene symbols in the target species.
#'
#' @keywords internal
.convert_gene_list_gorth <- function(
    genes,
    from_organism,
    to_organism
) {

  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop(
      "Package 'gprofiler2' is required for ortholog mapping. ",
      "Install it with install.packages('gprofiler2')."
    )
  }

  # Remove NA and duplicates
  genes <- unique(na.omit(genes))

  if (length(genes) == 0) {
    return(character())
  }

  # Query g:Profiler
  res <- try(
    gprofiler2::gorth(
      query = genes,
      source_organism = from_organism,
      target_organism = to_organism
    ),
    silent = TRUE
  )

  if (inherits(res, "try-error") || is.null(res)) {
    warning("g:Profiler ortholog mapping failed.")
    return(character())
  }

  # Extract target genes
  unique(res$ortholog_name)
}

#' Convert a Human Gene List to Mouse Orthologs
#'
#' Uses the g:Profiler \code{gorth()} service to map human gene symbols
#' (HGNC) to mouse gene symbols (MGI).
#'
#' This approach is more robust and reproducible than Ensembl BioMart
#' for cross-species transcriptomic analyses.
#'
#' @param genes Character vector of human gene symbols (HGNC).
#'
#' @return Character vector of mouse gene symbols (MGI).
#'
#' @export
convert_human_to_mouse <- function(genes) {

  .convert_gene_list_gorth(
    genes = genes,
    from_organism = "hsapiens",
    to_organism = "mmusculus"
  )
}


#' Convert a Mouse Gene List to Human Orthologs
#'
#' Uses the g:Profiler \code{gorth()} service to map mouse gene symbols
#' (MGI) to human gene symbols (HGNC).
#'
#' @param genes Character vector of mouse gene symbols (MGI).
#'
#' @return Character vector of human gene symbols (HGNC).
#'
#' @export
convert_mouse_to_human <- function(genes) {

  .convert_gene_list_gorth(
    genes = genes,
    from_organism = "mmusculus",
    to_organism = "hsapiens"
  )
}


.convert_gene_list_biomart <- function(
    genes,
    from_dataset,
    from_attr,
    to_dataset,
    to_attr
) {

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(
      "Package 'biomaRt' is required for gene conversion. ",
      "Install it with BiocManager::install('biomaRt')."
    )
  }

  mart_from <- biomaRt::useMart(
    "ensembl",
    dataset = from_dataset
  )

  mart_to <- biomaRt::useMart(
    "ensembl",
    dataset = to_dataset
  )

  res <- biomaRt::getLDS(
    attributes = from_attr,
    filters = from_attr,
    values = genes,
    mart = mart_from,
    attributesL = to_attr,
    martL = mart_to,
    uniqueRows = TRUE
  )

  unique(res[[to_attr]])
}


#' Convert a human gene list to mouse orthologs
#'
#' Uses Ensembl BioMart to map human HGNC symbols to mouse MGI symbols.
#'
#' @param genes Character vector of human gene symbols (HGNC)
#'
#' @return Character vector of mouse gene symbols (MGI)
#' @export
convert_human_to_mouse <- function(genes) {

  .convert_gene_list_biomart(
    genes = genes,
    from_dataset = "hsapiens_gene_ensembl",
    from_attr = "hgnc_symbol",
    to_dataset = "mmusculus_gene_ensembl",
    to_attr = "mgi_symbol"
  )
}


#' Convert a mouse gene list to human orthologs
#'
#' Uses Ensembl BioMart to map mouse MGI symbols to human HGNC symbols.
#'
#' @param genes Character vector of mouse gene symbols (MGI)
#'
#' @return Character vector of human gene symbols (HGNC)
#' @export
convert_mouse_to_human <- function(genes) {

  .convert_gene_list_biomart(
    genes = genes,
    from_dataset = "mmusculus_gene_ensembl",
    from_attr = "mgi_symbol",
    to_dataset = "hsapiens_gene_ensembl",
    to_attr = "hgnc_symbol"
  )
}

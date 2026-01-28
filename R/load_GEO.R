#' Discover GEO matrix samples in a directory
#'
#' @param folder_path Path containing *_matrix.mtx.gz files
#'
#' @return Character vector of sample prefixes
#' @export
find_geo_samples <- function(folder_path) {
  files <- list.files(
    folder_path,
    pattern = "_matrix.mtx.gz$",
    full.names = FALSE
  )

  sub("_matrix.mtx.gz$", "", files)
}


#' Load a single GEO matrix into a monocle3 CDS
#'
#' @param sample_name Sample prefix (no suffix)
#' @param folder_path Path containing GEO files
#'
#' @return monocle3 cell_data_set
#' @export
load_geo_sample <- function(sample_name, folder_path) {

  cds <- monocle3::load_mm_data(
    mat_path = file.path(folder_path, paste0(sample_name, "_matrix.mtx.gz")),
    feature_anno_path = file.path(folder_path, paste0(sample_name, "_features.tsv.gz")),
    cell_anno_path = file.path(folder_path, paste0(sample_name, "_barcodes.tsv.gz"))
  )

  colnames(SummarizedExperiment::fData(cds)) <- "gene_short_name"
  SummarizedExperiment::fData(cds)$id <- rownames(SummarizedExperiment::fData(cds))

  annotate_sample_metadata(cds, sample_name)
}

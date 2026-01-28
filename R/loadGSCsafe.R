#' Safely load gene set collections from common formats
#'
#' Reads gene set collections (GSCs) from GMT, SIF, or data.frame inputs.
#' Returns a standardized object of class \code{GSC}.
#'
#' @param file Character file path (for GMT or SIF) or a data.frame
#' @param type Input type. One of "auto", "gmt", "sif", or "data.frame"
#' @param addInfo Optional two-column data.frame with additional gene set metadata
#' @param sep Field separator for GMT files (default tab)
#' @param encoding File encoding for GMT files
#'
#' @return An object of class \code{GSC} with elements:
#' \describe{
#'   \item{gsc}{Named list of gene sets}
#'   \item{addInfo}{Optional data.frame of additional information}
#' }
#'
#' @export
loadGSCSafe <- function(
    file,
    type = "auto",
    addInfo = NULL,
    sep = "\t",
    encoding = "latin1"
) {

  addUserInfo <- if (is.null(addInfo)) "skip" else "yes"

  type <- try(
    match.arg(type, c("auto", "gmt", "sif", "data.frame")),
    silent = TRUE
  )
  if (inherits(type, "try-error")) {
    stop("argument 'type' set to unknown value")
  }

  ## ----------------------------
  ## Auto-detect type
  ## ----------------------------
  if (type == "auto") {
    if (is.character(file)) {
      ext <- tolower(tail(strsplit(file, "\\.")[[1]], 1))
      if (!ext %in% c("gmt", "sif")) {
        stop("Unsupported file extension: ", ext)
      }
      type <- ext
    } else {
      type <- "data.frame"
    }
  }

  ## ----------------------------
  ## GMT
  ## ----------------------------
  if (type == "gmt") {

    con <- base::file(file, encoding = encoding)
    on.exit(close(con), add = TRUE)

    gsc <- list()
    addInfo_local <- NULL

    while (length(line <- scan(
      con,
      nlines = 1,
      what = "character",
      quiet = TRUE,
      sep = sep
    )) > 0) {

      if (addUserInfo == "skip") {
        addInfo_local <- rbind(addInfo_local, line[1:2])
      }

      genes <- unique(line[-c(1, 2)])
      genes <- genes[genes != "" & !is.na(genes)]
      gsc[[line[1]]] <- genes
    }

    gsc <- gsc[!duplicated(names(gsc))]

    if (addUserInfo == "skip" && !is.null(addInfo_local)) {
      addInfo <- unique(addInfo_local)
    }
  }

  ## ----------------------------
  ## SIF
  ## ----------------------------
  else if (type == "sif") {

    gsc_df <- try(
      utils::read.delim(
        file,
        header = FALSE,
        quote = "",
        stringsAsFactors = FALSE
      ),
      silent = TRUE
    )

    if (inherits(gsc_df, "try-error") || ncol(gsc_df) != 3) {
      stop("SIF file must contain exactly three columns")
    }

    if (addUserInfo == "skip") {
      addInfo <- gsc_df[, 1:2]
    }

    gsc_df <- unique(gsc_df[, c(3, 1)])
    geneSets <- unique(gsc_df[, 2])

    gsc <- lapply(geneSets, function(gs) {
      gsc_df[gsc_df[, 2] == gs, 1]
    })
    names(gsc) <- geneSets
  }

  ## ----------------------------
  ## data.frame
  ## ----------------------------
  else if (type == "data.frame") {

    gsc_df <- as.data.frame(file, stringsAsFactors = FALSE)
    if (ncol(gsc_df) != 2) {
      stop("data.frame input must contain exactly two columns")
    }

    gsc_df <- unique(gsc_df)
    geneSets <- unique(gsc_df[, 2])

    gsc <- lapply(geneSets, function(gs) {
      gsc_df[gsc_df[, 2] == gs, 1]
    })
    names(gsc) <- geneSets
  }

  ## ----------------------------
  ## Additional info handling
  ## ----------------------------
  if (!is.null(addInfo)) {
    addInfo <- as.data.frame(addInfo, stringsAsFactors = FALSE)
    if (ncol(addInfo) != 2) {
      stop("addInfo must contain exactly two columns")
    }
    addInfo <- unique(addInfo[addInfo[, 1] %in% names(gsc), ])
  }

  res <- list(gsc = gsc, addInfo = addInfo)
  class(res) <- "GSC"
  res
}

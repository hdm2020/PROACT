#' Example dataset for the PROACT workflow
#'
#' A small example dataset for demonstrating the core PROACT workflow.
#' It contains a peak-by-motif occupancy matrix, selected accessible peaks,
#' expressed TFs, target TFs, and a precomputed TF prior matrix.
#'
#' @format A list with five elements:
#' \describe{
#'   \item{peak_motif}{A data frame containing a `peakID` column and binary motif occupancy columns.}
#'   \item{da_peaks}{A character vector of selected accessible peak IDs.}
#'   \item{expressed_tfs}{A character vector of TFs retained after expression filtering.}
#'   \item{target_tfs}{A character vector of TFs used as response TFs in PROACT.}
#'   \item{prior_mat}{A symmetric TF prior matrix used by `run_PROACT()`.}
#' }
#'
#' @usage data(PROACT_example)
#'
#' @examples
#' data(PROACT_example)
#' names(PROACT_example)
#'
"PROACT_example"

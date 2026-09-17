# 1. Load motif annotation
load_motif_annotation <- function(motif_anno_file = NULL){
    if(is.null(motif_anno_file)){
        motif_anno_file <- system.file("extdata", "metadata.tsv", package = "PROACT")
        if(!nzchar(motif_anno_file)) stop("Bundled motif annotation 'inst/extdata/metadata.tsv' was not found.")
    }
    if(!file.exists(motif_anno_file)) stop("Motif annotation file does not exist: ", motif_anno_file)
    motif_anno <- read.table(motif_anno_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    required_cols <- c("motif_id", "cluster", "tf_name")
    missing_cols <- setdiff(required_cols, colnames(motif_anno))
    if(length(missing_cols) > 0) stop("motif annotation is missing columns: ", paste(missing_cols, collapse = ", "))
    motif_anno <- motif_anno[, required_cols]
    name_map <- setNames(motif_anno$tf_name, motif_anno$motif_id)
    list(motif_anno = motif_anno, name_map = name_map)
}

# 2. Remove redundant motifs
remove_redundant_motifs <- function(peak_motif_matrix, remove_motifs = NULL){
    if(is.null(remove_motifs)) remove_motifs <- DEFAULT_REMOVE_MOTIFS
    keep_cols <- setdiff(colnames(peak_motif_matrix), remove_motifs)
    peak_motif_matrix[, keep_cols, drop = FALSE]
}

# 3. Read peak-by-motif matrix
load_peak_motif_matrix <- function(pm_matrix_file){
    read.table(pm_matrix_file, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

# 4. Read differential peaks
load_da_peaks <- function(da_peaks_file, ident){
    da <- read.csv(da_peaks_file, check.names = FALSE, stringsAsFactors = FALSE)
    da <- na.omit(da[, c("Peak", ident)])
    unique(da$Peak)
}

# 5. Subset peak-by-motif matrix by differential peaks
subset_peak_motif_by_da <- function(peak_motif_matrix, da_peaks){
    peak_motif_matrix[peak_motif_matrix$peakID %in% da_peaks, , drop = FALSE]
}

# 6. Map motif IDs to TF names
convert_motif_to_tf <- function(peak_motif_matrix, name_map){
    motif_ids <- colnames(peak_motif_matrix)
    mapped <- motif_ids %in% names(name_map)
    if(any(!mapped)){
        warning(sum(!mapped), " motif columns are absent from motif annotation and were removed.")
        peak_motif_matrix <- peak_motif_matrix[, mapped, drop = FALSE]
        motif_ids <- colnames(peak_motif_matrix)
    }
    colnames(peak_motif_matrix) <- unname(name_map[motif_ids])
    peak_motif_matrix
}

# 7. Filter TFs by minimum number of binding peaks
filter_low_binding_tfs <- function(peak_motif_matrix, min_sites = 5){
    counts <- colSums(peak_motif_matrix != 0)
    keep <- names(counts[counts >= min_sites])
    peak_motif_matrix[, keep, drop = FALSE]
}

# 8. Filter TFs by minimum binding proportion
filter_low_binding_percent <- function(peak_motif_matrix, min_percent = 0.005){
    threshold <- nrow(peak_motif_matrix) * min_percent
    counts <- colSums(peak_motif_matrix != 0)
    keep <- names(counts[counts >= threshold])
    peak_motif_matrix[, keep, drop = FALSE]
}

# 9. Read expressed TFs
load_expressed_tfs <- function(expr_file, ident){
    expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
    expr$TF <- rownames(expr)
    expr <- expr[!is.na(expr[[ident]]), ]
    expr$TF
}

# 10. Public API
#' Prepare a binary peak-by-TF matrix
#'
#' Convert a peak-by-motif matrix into the binary peak-by-TF matrix used by
#' PROACT. The function restricts the input to selected peaks, removes redundant
#' motifs, maps motif IDs to TF names, collapses multiple motifs for the same TF,
#' binarizes the matrix, and optionally filters TFs by binding frequency and
#' expression.
#'
#' @param peak_motif_mat A data frame or matrix containing a `peakID` column and
#'   motif columns. Non-zero motif entries are interpreted as motif presence.
#' @param da_peaks Character vector of peak IDs to retain.
#' @param motif_anno_file Optional path to a tab-delimited motif annotation file with
#'   columns `motif_id`, `cluster`, and `tf_name`. If `NULL`, the bundled
#'   `extdata/metadata.tsv` file is used.
#' @param filter_expr Logical. If `TRUE`, retain only TFs listed in `expr_tfs`.
#' @param expr_tfs Character vector of expressed TFs. Required when
#'   `filter_expr = TRUE`.
#' @param min_sites Minimum number of peaks containing a TF motif when
#'   `min_percent` is `NULL`.
#' @param min_percent Optional minimum fraction of peaks containing a TF motif.
#'   When supplied, this criterion is used instead of `min_sites`.
#' @param remove_motifs Character vector of motif IDs to exclude. If `NULL`,
#'   `DEFAULT_REMOVE_MOTIFS` is used. Use `character(0)` to disable motif removal.
#'
#' @return A list with elements `X` (binary peak-by-TF numeric matrix), `Peaks`,
#'   `TFs`, `nPeaks`, `nTFs`, and `expressed_tfs`.
#'
#' @details Motif columns not found in `motif_anno_file` are removed with a
#'   warning. Binarization is applied after duplicated TF columns are collapsed,
#'   ensuring that the final matrix contains only 0/1 values.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' peak_tf <- prepare_peak_tf_matrix(
#'     peak_motif_mat = peak_motif_mat,
#'     da_peaks = da_peaks
#' )
#' }
prepare_peak_tf_matrix <- function(
    peak_motif_mat,
    da_peaks,
    motif_anno_file = NULL,
    filter_expr = FALSE,
    expr_tfs = NULL,
    min_sites = 5,
    min_percent = NULL,
    remove_motifs = NULL
){
    if(is.null(remove_motifs)) remove_motifs <- DEFAULT_REMOVE_MOTIFS
    if(!"peakID" %in% colnames(peak_motif_mat)) stop("peak_motif_mat must contain a 'peakID' column.")
    expressed_tfs <- NULL

    message("========================================")
    message("Preparing Peak x TF matrix")
    message("========================================")

    # 1. Motif annotation
    motif_info <- load_motif_annotation(motif_anno_file = motif_anno_file)
    name_map <- motif_info$name_map

    # 2. Subset peaks
    df <- subset_peak_motif_by_da(peak_motif_matrix = peak_motif_mat, da_peaks = da_peaks)
    if(nrow(df) == 0) stop("No peaks remain after peak filtering.")
    message("After peak filtering: ", nrow(df), " peaks")

    # 3. Remove redundant motifs
    n_removed <- sum(colnames(df) %in% remove_motifs)
    df <- remove_redundant_motifs(peak_motif_matrix = df, remove_motifs = remove_motifs)
    message("Removed ", n_removed, " redundant motif columns")

    # 4. Save peak IDs
    peak_ids <- df$peakID
    df$peakID <- NULL

    # 5. Motif -> TF
    df <- convert_motif_to_tf(peak_motif_matrix = df, name_map = name_map)
    if(ncol(df) == 0) stop("No mapped TF motif columns remain.")

    # 6. Collapse duplicated TFs
    if(any(duplicated(colnames(df)))){
        message("Collapsing duplicated TF motifs ...")
        uniq_tf <- unique(colnames(df))
        df2 <- lapply(uniq_tf, function(tf){
            idx <- which(colnames(df) == tf)
            if(length(idx) == 1) df[, idx] else rowSums(df[, idx, drop = FALSE])
        })
        df <- as.data.frame(df2, check.names = FALSE)
        colnames(df) <- uniq_tf
    }

    # Binarize after duplicated motifs are collapsed
    df[df != 0] <- 1

    # 7. Remove low-binding TFs
    if(!is.null(min_percent)){
        if(min_percent < 0 || min_percent > 1) stop("min_percent must be between 0 and 1.")
        df <- filter_low_binding_percent(peak_motif_matrix = df, min_percent = min_percent)
        message("Using motif frequency threshold = ", min_percent)
    } else {
        if(min_sites < 0) stop("min_sites must be non-negative.")
        df <- filter_low_binding_tfs(peak_motif_matrix = df, min_sites = min_sites)
        message("Using minimum motif count = ", min_sites)
    }
    if(ncol(df) == 0) stop("No TFs remain after motif frequency filtering.")
    message(ncol(df), " TFs retained")

    # 8. Expression filtering
    if(filter_expr){
        if(is.null(expr_tfs)) stop("expr_tfs must be provided when filter_expr = TRUE.")
        expressed_tfs <- intersect(expr_tfs, colnames(df))
        df <- df[, expressed_tfs, drop = FALSE]
        if(ncol(df) == 0) stop("No TFs remain after expression filtering.")
        message(ncol(df), " TFs retained after expression filtering")
    }

    # 9. Convert to numeric matrix
    X <- as.matrix(df)
    storage.mode(X) <- "numeric"
    rownames(X) <- peak_ids
    message("Final matrix: ", nrow(X), " peaks x ", ncol(X), " TFs")

    list(
        X = X,
        Peaks = peak_ids,
        TFs = colnames(X),
        nPeaks = nrow(X),
        nTFs = ncol(X),
        expressed_tfs = expressed_tfs
    )
}


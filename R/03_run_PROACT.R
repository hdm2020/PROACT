# 1. Build priority blocks
build_priority_blocks <- function(
    target_tf,
    prior_mat,
    motif_anno_file = NULL,
    remove_motifs = NULL
){
    if(is.null(remove_motifs)) remove_motifs <- DEFAULT_REMOVE_MOTIFS

    prior_vec <- prior_mat[, target_tf]
    prior_vec <- prior_vec[names(prior_vec) != target_tf]

    # Remove predictors in the same motif cluster as the response TF
    motif_info <- load_motif_annotation(motif_anno_file = motif_anno_file)
    motif_cluster <- motif_info$motif_anno
    motif_cluster <- motif_cluster[!(motif_cluster$motif_id %in% remove_motifs), ]

    if(target_tf %in% motif_cluster$tf_name){
        clst <- unique(motif_cluster$cluster[motif_cluster$tf_name == target_tf])
        same_cluster_tf <- unique(motif_cluster$tf_name[motif_cluster$cluster %in% clst])
        prior_vec <- prior_vec[setdiff(names(prior_vec), same_cluster_tf)]
    } else {
        message("Target TF not in motif annotation, skipping cluster removal.")
    }

    message(target_tf, ": ", length(prior_vec), " TFs as predictors")

    block1 <- which(prior_vec == 1)
    block2 <- which(prior_vec == 2)
    block3 <- which(prior_vec == 0)

    # prioritylasso requires at least two predictors per block
    if(length(block1) >= 2 && length(block2) >= 2){
        blocks <- list(block1, block2, block3)
        message(paste(target_tf, length(block1), length(block2), length(block3), sep = ":"))
    } else {
        block12 <- c(block1, block2)
        if(length(block12) >= 2){
            blocks <- list(block12, block3)
            message(paste(target_tf, length(block12), length(block3), sep = ":"))
        } else {
            block123 <- c(block12, block3)
            blocks <- list(block123)
            message(paste(target_tf, length(block123), sep = ":"))
        }
    }

    # Remove empty blocks; merge a one-predictor lower-priority block upward
    blocks <- blocks[lengths(blocks) > 0]
    if(length(blocks) > 1 && any(lengths(blocks) < 2)){
        for(i in rev(seq_along(blocks)[-1])){
            if(length(blocks[[i]]) < 2){
                blocks[[i - 1]] <- c(blocks[[i - 1]], blocks[[i]])
                blocks[[i]] <- NULL
            }
        }
    }

    if(length(blocks) == 0 || any(lengths(blocks) < 2)) return(list(blocks = NULL, predictors = names(prior_vec)))
    list(blocks = blocks, predictors = names(prior_vec))
}

# 2. Fit one target-TF model
run_single_tf_lasso <- function(
    target_tf,
    X,
    prior_mat,
    motif_anno_file = NULL,
    family = "binomial",
    block1.penalization = TRUE,
    remove_motifs = NULL
){
    if(is.null(remove_motifs)) remove_motifs <- DEFAULT_REMOVE_MOTIFS
    if(!target_tf %in% colnames(X)) return(NULL)

    info <- build_priority_blocks(target_tf, prior_mat, motif_anno_file = motif_anno_file, remove_motifs = remove_motifs)
    predictors <- info$predictors
    if(is.null(info$blocks) || length(predictors) < 2) return(NULL)

    X_sub <- X[, predictors, drop = FALSE]
    Y <- X[, target_tf]
    if(length(unique(Y)) < 2) return(NULL)

    set.seed(1234)
    fit <- prioritylasso(
        X = X_sub, Y = Y, family = family, type.measure = "auc", blocks = info$blocks,
        lambda.type = "lambda.min", block1.penalization = block1.penalization, standardize = FALSE
    )

    coef <- fit$coefficients
    result <- data.frame(TF1 = target_tf, TF2 = names(coef), beta = as.numeric(coef), stringsAsFactors = FALSE)
    result[result$TF2 != "(Intercept)", , drop = FALSE]
}

# 3. Public API
#' Run PROACT cooperative TF inference
#'
#' Fit prior-regularized priority LASSO models to a binary peak-by-TF matrix and
#' infer candidate cooperative TF pairs. For each target TF, motif occupancy is
#' modeled as the response and the remaining TFs are used as predictors, with
#' predictor blocks defined by the prior matrix.
#'
#' @param peak_tf_matrix Output from [prepare_peak_tf_matrix()].
#' @param prior_mat Symmetric TF prior matrix, typically produced by
#'   [build_prior()].
#' @param target_tfs Character vector of TFs to use as response TFs.
#' @param all_jaspar_tf Logical. If `TRUE`, all TFs shared by the peak-by-TF
#'   matrix and prior matrix may be predictors. If `FALSE`, the prior matrix is
#'   first restricted to `target_tfs`.
#' @param return.pos Logical. If `TRUE`, retain only associations with positive
#'   LASSO coefficients after unordered-pair deduplication.
#' @param motif_anno_file Optional path to a motif annotation file with columns
#'   `motif_id`, `cluster`, and `tf_name`. If `NULL`, the bundled
#'   `extdata/metadata.tsv` file is used.
#' @param remove_motifs Character vector of motif IDs to exclude when defining
#'   motif clusters. If `NULL`, `DEFAULT_REMOVE_MOTIFS` is used. Use
#'   `character(0)` to disable motif removal.
#' @param ncores Number of parallel workers. A serial backend is used when
#'   `ncores = 1`; socket parallelism is used on Windows and multicore
#'   parallelism on Unix-like systems.
#' @param block1.penalization Logical passed to [prioritylasso::prioritylasso()]
#'   to determine whether the first priority block is penalized.
#'
#' @return A data frame of candidate TF pairs with columns including `TF1`,
#'   `TF2`, `beta`, `Cluster1`, `Cluster2`, `prior`, and `pairs`.
#'
#' @details TFs belonging to the same motif cluster as the response TF are
#'   removed from its predictor set. When the same unordered TF pair is inferred
#'   in both directions, only the result with the largest absolute coefficient
#'   is retained. A positive coefficient represents a positive conditional motif
#'   association and should not by itself be interpreted as direct physical
#'   interaction or causal regulation.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pairs <- run_PROACT(
#'     peak_tf_matrix = peak_tf,
#'     prior_mat = prior_mat,
#'     target_tfs = c("JUN", "FOS"),
#'     ncores = 4
#' )
#' }
run_PROACT <- function(
    peak_tf_matrix,
    prior_mat,
    target_tfs,
    all_jaspar_tf = TRUE,
    return.pos = TRUE,
    motif_anno_file = NULL,
    remove_motifs = NULL,
    ncores = 4,
    block1.penalization = TRUE
){
    if(is.null(remove_motifs)) remove_motifs <- DEFAULT_REMOVE_MOTIFS
    if(!is.list(peak_tf_matrix) || is.null(peak_tf_matrix$X)) stop("peak_tf_matrix must be the output of prepare_peak_tf_matrix().")
    if(is.null(rownames(prior_mat)) || is.null(colnames(prior_mat))) stop("prior_mat must have TF row and column names.")
    if(length(ncores) != 1 || is.na(ncores) || ncores < 1) stop("ncores must be a positive integer.")
    ncores <- as.integer(ncores)

    X <- peak_tf_matrix$X

    if(!all_jaspar_tf){
        target_prior_tf <- intersect(target_tfs, colnames(prior_mat))
        prior_mat <- prior_mat[target_prior_tf, target_prior_tf, drop = FALSE]
    }

    common_tf <- intersect(colnames(X), colnames(prior_mat))
    if(length(common_tf) < 2) stop("Fewer than 2 TFs are available for priority lasso.")

    X <- X[, common_tf, drop = FALSE]
    prior_mat <- prior_mat[common_tf, common_tf, drop = FALSE]
    message(length(common_tf), " TFs used for priority lasso")

    target_tfs <- intersect(target_tfs, colnames(X))
    if(length(target_tfs) == 0) stop("No target TFs remain after filtering.")
    message(length(target_tfs), " target TFs retained for modeling")

    if(ncores == 1){
        param <- SerialParam()
    } else {
        param <- SnowParam(workers = ncores, progressbar = TRUE)
    }

    results <- bplapply(target_tfs, function(tf){
        run_single_tf_lasso(
            target_tf = tf, X = X, prior_mat = prior_mat, motif_anno_file = motif_anno_file,
            block1.penalization = block1.penalization, remove_motifs = remove_motifs
        )
    }, BPPARAM = param)

    results <- do.call(rbind, results)
    if(is.null(results) || nrow(results) == 0) stop("No TF cooperation detected.")

    # Add motif-cluster information
    motif_info <- load_motif_annotation(motif_anno_file)
    motif_cluster <- motif_info$motif_anno
    motif_cluster <- motif_cluster[!(motif_cluster$motif_id %in% remove_motifs), ]
    results$Cluster1 <- motif_cluster$cluster[match(results$TF1, motif_cluster$tf_name)]
    results$Cluster2 <- motif_cluster$cluster[match(results$TF2, motif_cluster$tf_name)]

    # Add prior information
    results$prior <- mapply(function(a, b) prior_mat[a, b], results$TF2, results$TF1)

    # Retain the largest absolute coefficient for each unordered pair
    results$pairs <- paste(pmin(results$TF1, results$TF2), pmax(results$TF1, results$TF2), sep = "_")
    df_sorted <- results[order(results$pairs, -abs(results$beta)), ]
    results <- df_sorted[!duplicated(df_sorted$pairs), ]
    results <- results[order(results$TF1, -abs(results$beta)), ]

    if(return.pos) results <- results[!is.na(results$beta) & results$beta > 0, , drop = FALSE]
    rownames(results) <- NULL
    results
}


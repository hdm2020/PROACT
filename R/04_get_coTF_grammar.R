# 1. Convert peak regions to GRanges and optionally resize around peak centers
get_center_peaks <- function(features, region_length = NULL, genome = NULL){
    parts <- do.call(rbind, strsplit(features, "-"))
    peaks_gr <- GRanges(seqnames = parts[, 1], ranges = IRanges(start = as.numeric(parts[, 2]), end = as.numeric(parts[, 3])))
    peaks_gr$peak_ID <- features
    #peaks_gr <- peaks_gr[!seqnames(peaks_gr) %in% c("chrY", "chrM")]
    peaks_gr <- peaks_gr[!as.character(seqnames(peaks_gr)) %in% c("chrY", "chrM")]
    message("Dealing with ", length(peaks_gr), " features!")

    if(!is.null(region_length)){
        if(is.null(genome)) stop("genome must be provided when region_length is not NULL.")
        chrom_sizes <- seqlengths(genome)
        chrom_sizes <- chrom_sizes[grep("^(chr[0-9X]+)$", names(chrom_sizes))]
        new_seqinfo <- Seqinfo(seqnames = names(chrom_sizes), seqlengths = chrom_sizes)
        seqlevels(peaks_gr) <- seqlevels(new_seqinfo)
        seqinfo(peaks_gr) <- new_seqinfo
        peaks_gr <- trim(resize(peaks_gr, width = region_length, fix = "center"))
    }

    peaks_gr
}

# 2. motif positions
prepare_motif_positions <- function(motif_gr_list, name_map, remove_motifs){
    keep_motifs <- setdiff(names(name_map), remove_motifs)
    motif_gr_sub <- motif_gr_list[intersect(names(motif_gr_list), keep_motifs)]
    if(length(motif_gr_sub) == 0) stop("No annotated motif positions remain after motif filtering.")

    # Add TF and motif columns
    motif_gr <- mapply(function(gr, motif_id){
        gr$TF <- name_map[[motif_id]]
        gr$motif <- motif_id
        gr
    }, motif_gr_sub, names(motif_gr_sub), SIMPLIFY = FALSE)

    unlist(GRangesList(motif_gr), use.names = FALSE)
}


# 3. Public API: identify TF-pair spatial grammar
#' Identify spatial grammar for TF pairs
#'
#' Extract motif instances for candidate TF pairs within selected chromatin
#' regions and summarize their relative spacing, strand geometry, and order.
#'
#' @param object Optional Signac/Seurat object whose selected assay contains motif
#'   positions. Used for single-cell chromatin accessibility data.
#' @param features Character vector of genomic regions formatted as
#'   `chr-start-end`. Required when `object` is supplied.
#' @param pairs Character vector of TF pairs formatted as `TF1_TF2`.
#' @param assay Name of the assay containing motif positions. Required when
#'   `object` is supplied.
#' @param motif_gr_list Optional named GRangesList containing motif positions
#'   for bulk data. List names should be motif IDs, and each motif instance
#'   must contain `peak_id` and `score` metadata columns.
#' @param strategy Pair-selection strategy within each peak: `"best"` selects
#'   the highest-scoring motif instance for each TF; `"nearest"` selects the
#'   motif pair with the smallest center-to-center distance.
#' @param motif_anno_file Optional path to the motif annotation file with columns
#'   `motif_id`, `cluster`, and `tf_name`. If `NULL`, the bundled
#'   `extdata/metadata.tsv` file is used.
#' @param remove_motifs Character vector of motif IDs to exclude. If `NULL`,
#'   `DEFAULT_REMOVE_MOTIFS` is used. Use `character(0)` to disable removal.
#' @param region_length Optional width in base pairs for resizing each feature
#'   around its center before motif overlap analysis.
#' @param genome Genome object providing chromosome lengths. Required when
#'   `region_length` is not `NULL`.
#'
#' @return A named list with one element per requested TF pair. Each non-empty
#'   element is a data frame describing motif spacing, strand combination,
#'   geometry, arrangement, order, and motif scores for representative instances
#'   across peaks.
#'
#' @details `edge_distance` is defined as the downstream motif start minus the
#'   upstream motif end. `center_distance` is the distance between motif centers.
#'   Geometry is classified as `HH`, `TT`, or `HT`; `HT` cases are further split
#'   into `A,B_HT` and `B,A_HT` according to motif order.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' grammar <- find_grammar(
#'     object = atac,
#'     features = da_peaks,
#'     pairs = c("JUN_FOS", "CEBPB_ATF4"),
#'     assay = "peaks",
#'     genome = genome
#' )
#' }
find_grammar <- function(
    object = NULL,
    features = NULL,
    pairs,
    assay = NULL,
    motif_gr_list = NULL,
    strategy = c("best", "nearest"),
    motif_anno_file = NULL,
    remove_motifs = NULL,
    region_length = NULL,
    genome = NULL
){
    strategy <- match.arg(strategy)
    if(is.null(remove_motifs)) remove_motifs <- DEFAULT_REMOVE_MOTIFS

    # Motif positions
    motif_info <- load_motif_annotation(motif_anno_file)
    name_map <- motif_info$name_map

    if(!is.null(object) && !is.null(motif_gr_list)){
        stop("Provide either 'object' or 'motif_gr_list', not both.")
    }

    if(is.null(object) && is.null(motif_gr_list)){
        stop("Either 'object' or 'motif_gr_list' must be provided.")
    }

    if(!is.null(object)){
        if(is.null(assay)) stop("'assay' must be provided when using a Signac object.")
        if(is.null(features)) stop("'features' must be provided when using a Signac object.")

        motif_gr <- object[[assay]]@motifs@positions
        motif_gr <- prepare_motif_positions(motif_gr, name_map, remove_motifs)

	# Peaks and motif-peak overlap
        peak_gr <- get_center_peaks(features, region_length, genome)
        hits <- findOverlaps(motif_gr, peak_gr, type = "within")
        qh <- queryHits(hits)
        sh <- subjectHits(hits)

        motif_df <- data.frame(
            peak_id = peak_gr$peak_ID[sh],
            chr = as.character(seqnames(motif_gr))[qh],
            start = start(motif_gr)[qh],
            end = end(motif_gr)[qh],
            strand = as.character(strand(motif_gr))[qh],
            TF = motif_gr$TF[qh],
            score = motif_gr$score[qh]
        )
    }

    if(!is.null(motif_gr_list)){
        motif_gr <- prepare_motif_positions(motif_gr_list, name_map, remove_motifs)

        if(!"peak_id" %in% colnames(mcols(motif_gr))){
            stop("motif_gr_list must contain a 'peak_id' metadata column.")
        }
        if(!"score" %in% colnames(mcols(motif_gr))){
            stop("motif_gr_list must contain a 'score' metadata column.")
        }    

        motif_df <- data.frame(
            peak_id = motif_gr$peak_id,
            chr = as.character(seqnames(motif_gr)),
            start = start(motif_gr),
            end = end(motif_gr),
            strand = as.character(strand(motif_gr)),
            TF = motif_gr$TF,
            score = motif_gr$score
        )
    }
    
    motif_df$center <- (motif_df$start + motif_df$end) / 2

    # Calculate the spatial relation of one motif pair
    compute_pair <- function(m1, m2, tf1, tf2){
        if(m1$center <= m2$center){
            up <- m1
            down <- m2
            motif_order <- paste0(tf1, "->", tf2)
            tf_up <- tf1
        } else {
            up <- m2
            down <- m1
            motif_order <- paste0(tf2, "->", tf1)
            tf_up <- tf2
        }

        edge_dist <- down$start - up$end
        center_dist <- down$center - up$center
        strand_pair <- paste0(as.character(up$strand), as.character(down$strand))
        geometry <- if(strand_pair == "+-") "HH" else if(strand_pair == "-+") "TT" else if(strand_pair %in% c("++", "--")) "HT" else NA_character_
        arrangement <- geometry
        if(!is.na(geometry) && geometry == "HT") arrangement <- if(tf_up == tf1) "A,B_HT" else "B,A_HT"

        data.frame(
            chr = up$chr,
            start1 = up$start, end1 = up$end, start2 = down$start, end2 = down$end,
            center1 = up$center, center2 = down$center,
            edge_distance = edge_dist, center_distance = center_dist,
            strand_pair = strand_pair, geometry = geometry, arrangement = arrangement, order = motif_order,
            score1 = up$score, score2 = down$score, pair_score = up$score + down$score
        )
    }

    # Select one representative motif pair per peak
    pair_by_strategy <- function(df, tf1, tf2){
        m1 <- df[df$TF == tf1, , drop = FALSE]
        m2 <- df[df$TF == tf2, , drop = FALSE]
        if(nrow(m1) == 0 || nrow(m2) == 0) return(NULL)

        if(strategy == "nearest"){
            comb <- expand.grid(i = seq_len(nrow(m1)), j = seq_len(nrow(m2)))
            comb$distance <- abs(m1$center[comb$i] - m2$center[comb$j])
            k <- which.min(comb$distance)
            return(compute_pair(m1[comb$i[k], ], m2[comb$j[k], ], tf1, tf2))
        }

        i <- which.max(m1$score)
        j <- which.max(m2$score)
        compute_pair(m1[i, ], m2[j, ], tf1, tf2)
    }

    peak_list <- split(motif_df, motif_df$peak_id)
    result_list <- lapply(pairs, function(p){
        tf <- strsplit(p, "_", fixed = TRUE)[[1]]
        if(length(tf) != 2) stop("Each pair must be formatted as 'TF1_TF2': ", p)
        tf1 <- tf[1]
        tf2 <- tf[2]
        res <- lapply(peak_list, pair_by_strategy, tf1 = tf1, tf2 = tf2)
        res <- Filter(Negate(is.null), res)
        if(length(res) == 0) return(NULL)
        res <- do.call(rbind, res)
        rownames(res) <- NULL
        res$pair <- p
        res
    })
    names(result_list) <- pairs
    result_list
}

# 3. Public API: test spacing, geometry, and arrangement preferences
#' Test spacing and orientation preferences of TF pairs
#'
#' Summarize the output of [find_grammar()] and test whether observed spacing,
#' geometry, and arrangement counts deviate from a uniform distribution. When
#' expected cell counts are small, a Monte Carlo chi-squared test is used.
#'
#' @param result_list Output from [find_grammar()].
#' @param dist_bin Bin width in base pairs used to discretize spacing.
#' @param min_count Minimum number of motif-pair observations required for a TF
#'   pair to be analyzed.
#' @param B Number of Monte Carlo replicates used when expected chi-squared cell
#'   counts are below 5.
#'
#' @return A named list containing counts, proportions, and p-values for edge
#'   spacing, center spacing, geometry, and arrangement for each retained TF pair.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' grammar_stats <- analyze_tf_grammar(grammar, dist_bin = 5, min_count = 50)
#' }
analyze_tf_grammar <- function(result_list, dist_bin = 5, min_count = 5, B = 5000){
    if(dist_bin <= 0) stop("dist_bin must be positive.")
    if(min_count < 1) stop("min_count must be at least 1.")
    if(B < 1) stop("B must be at least 1.")

    smart_chisq <- function(tab, n_sim = B){
        if(length(tab) <= 1) return(NA_real_)
        test0 <- suppressWarnings(chisq.test(tab))
        if(any(test0$expected < 5)) chisq.test(tab, simulate.p.value = TRUE, B = n_sim)$p.value else test0$p.value
    }

    res_summary <- list()
    for(pair in names(result_list)){
        df <- result_list[[pair]]
        if(is.null(df) || nrow(df) < min_count) next

        df$edge_bin <- round(df$edge_distance / dist_bin) * dist_bin
        df$edge_bin[df$edge_bin < 0] <- 0
        edge_table <- table(df$edge_bin)
        edge_p <- tryCatch(smart_chisq(edge_table), error = function(e) NA_real_)

        df$center_bin <- round(df$center_distance / dist_bin) * dist_bin
        center_table <- table(df$center_bin)
        center_p <- tryCatch(smart_chisq(center_table), error = function(e) NA_real_)

        geometry_table <- table(df$geometry)
        geometry_p <- tryCatch(smart_chisq(geometry_table), error = function(e) NA_real_)

        arrangement_table <- table(df$arrangement)
        arrangement_p <- tryCatch(smart_chisq(arrangement_table), error = function(e) NA_real_)

        res_summary[[pair]] <- list(
            n = nrow(df),
            edge_table = edge_table, edge_pvalue = edge_p,
            center_table = center_table, center_pvalue = center_p,
            geometry_table = geometry_table, geometry_prop = prop.table(geometry_table), geometry_pvalue = geometry_p,
            arrangement_table = arrangement_table, arrangement_prop = prop.table(arrangement_table), arrangement_pvalue = arrangement_p
        )
    }
    res_summary
}

# 4. Public API: plot one pair
#' Plot spatial grammar for one TF pair
#'
#' Plot edge spacing, center spacing, geometry, and arrangement summaries for a
#' single TF pair returned by [analyze_tf_grammar()].
#'
#' @param res Output from [analyze_tf_grammar()].
#' @param pair Name of a TF pair contained in `res`.
#'
#' @return A patchwork object containing four grammar panels.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_tf_grammar(grammar_stats, "JUN_FOS")
#' }
plot_tf_grammar <- function(res, pair){
    if(!pair %in% names(res)) stop("Pair not found in result list!")
    x <- res[[pair]]

    df_edge <- data.frame(dist = as.numeric(names(x$edge_table)), count = as.numeric(x$edge_table))
    peak_edge <- df_edge$dist[which.max(df_edge$count)]
    p_edge <- ggplot(df_edge, aes(dist, count)) +
        geom_col(fill = "steelblue") +
        geom_vline(xintercept = peak_edge, linetype = "dashed", color = "red") +
        theme_classic() +
        labs(title = paste0(pair, " (edge, p=", signif(x$edge_pvalue, 2), ")"), x = "Edge distance", y = "Count")

    df_center <- data.frame(dist = as.numeric(names(x$center_table)), count = as.numeric(x$center_table))
    peak_center <- df_center$dist[which.max(df_center$count)]
    p_center <- ggplot(df_center, aes(dist, count)) +
        geom_col(fill = "orange") +
        geom_vline(xintercept = peak_center, linetype = "dashed", color = "red") +
        theme_classic() +
        labs(title = paste0(pair, " (center, p=", signif(x$center_pvalue, 2), ")"), x = "Center distance", y = "Count")

    df_geometry <- data.frame(geometry = names(x$geometry_prop), prop = as.numeric(x$geometry_prop))
    df_geometry$label <- paste0(round(df_geometry$prop * 100, 1), "%")
    p_geometry <- ggplot(df_geometry, aes(geometry, prop)) +
        geom_col(fill = "tomato") +
        geom_text(aes(label = label), vjust = -0.5, size = 4) +
        theme_classic() +
        ylim(0, max(df_geometry$prop) * 1.2) +
        labs(title = paste0(pair, " (geometry, p=", signif(x$geometry_pvalue, 2), ")"), x = "Geometry", y = "Proportion")

    df_arrangement <- data.frame(arrangement = names(x$arrangement_prop), prop = as.numeric(x$arrangement_prop))
    df_arrangement$label <- paste0(round(df_arrangement$prop * 100, 1), "%")
    p_arrangement <- ggplot(df_arrangement, aes(arrangement, prop)) +
        geom_col(fill = "darkgreen") +
        geom_text(aes(label = label), vjust = -0.5, size = 4) +
        theme_classic() +
        ylim(0, max(df_arrangement$prop) * 1.2) +
        labs(title = paste0(pair, " (arrangement, p=", signif(x$arrangement_pvalue, 2), ")"), x = "Arrangement", y = "Proportion")

    wrap_plots(p_edge, p_center, p_geometry, p_arrangement, ncol = 2)
}

# 5. Public API: summarize pair-level grammar features
#' Summarize TF-pair grammar features
#'
#' Convert grammar test results into a pair-level table containing normalized
#' spacing and arrangement entropy, dominant spacing and arrangement, associated
#' proportions, raw p-values, and Benjamini-Hochberg adjusted FDR values.
#'
#' @param res_summary Output from [analyze_tf_grammar()].
#' @param spacing_metric Spacing definition to summarize: `"edge"` or
#'   `"center"`.
#' @param min_count Minimum number of observations required for a TF pair.
#'
#' @return A data frame with one row per retained TF pair and pair-level grammar
#'   summary statistics.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' grammar_df <- summarize_tf_grammar(grammar_stats, spacing_metric = "edge", min_count = 50)
#' }
summarize_tf_grammar <- function(res_summary, spacing_metric = c("edge", "center"), min_count = 50){
    spacing_metric <- match.arg(spacing_metric)

    all_spacing_bins <- unique(unlist(lapply(res_summary, function(x){
        if(spacing_metric == "edge") as.numeric(names(x$edge_table)) else as.numeric(names(x$center_table))
    })))
    spacing_global_bins <- length(all_spacing_bins)
    arrangement_global_bins <- 4

    calc_entropy <- function(tab, global_bins){
        p <- prop.table(tab)
        p <- p[p > 0]
        if(length(p) <= 1 || global_bins <= 1) return(0)
        -sum(p * log2(p)) / log2(global_bins)
    }

    grammar_list <- lapply(names(res_summary), function(pair){
        x <- res_summary[[pair]]
        if(is.null(x$n) || x$n < min_count) return(NULL)

        spacing_table <- if(spacing_metric == "edge") x$edge_table else x$center_table
        spacing_prop <- prop.table(spacing_table)
        arrangement_prop <- prop.table(x$arrangement_table)

        data.frame(
            pair = pair, n = x$n, spacing_metric = spacing_metric,
            spacing_entropy = calc_entropy(spacing_table, spacing_global_bins),
            dominant_spacing = as.numeric(names(which.max(spacing_prop))),
            dominant_spacing_prop = max(spacing_prop),
            spacing_pvalue = if(spacing_metric == "edge") x$edge_pvalue else x$center_pvalue,
            arrangement_entropy = calc_entropy(x$arrangement_table, arrangement_global_bins),
            dominant_arrangement = names(which.max(arrangement_prop)),
            dominant_arrangement_prop = max(arrangement_prop),
            arrangement_pvalue = x$arrangement_pvalue,
            stringsAsFactors = FALSE
        )
    })

    grammar_list <- Filter(Negate(is.null), grammar_list)
    if(length(grammar_list) == 0) return(data.frame())

    grammar_df <- do.call(rbind, grammar_list)
    rownames(grammar_df) <- NULL
    grammar_df$spacing_FDR <- p.adjust(grammar_df$spacing_pvalue, method = "BH")
    grammar_df$arrangement_FDR <- p.adjust(grammar_df$arrangement_pvalue, method = "BH")
    grammar_df <- grammar_df[order(-grammar_df$n), , drop = FALSE]
    rownames(grammar_df) <- NULL
    grammar_df
}


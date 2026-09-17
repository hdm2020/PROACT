#' PROACT: Prior-Regularized Optimization for Analysis of Co-occurring TFs
#'
#' PROACT provides functions for constructing TF prior matrices, preparing binary
#' peak-by-TF matrices, inferring candidate cooperative TF associations with
#' priority LASSO, characterizing spatial motif grammar, and visualizing TF-TF networks.
#'
#' The package bundles its motif annotation as `inst/extdata/metadata.tsv`.
#' PPI evidence files are not bundled and must be supplied by the user.
#'
#' @section Main functions:
#' \itemize{
#'   \item [build_prior()] constructs the TF prior matrix.
#'   \item [prepare_peak_tf_matrix()] prepares the binary peak-by-TF matrix.
#'   \item [run_PROACT()] performs prior-regularized cooperative TF inference.
#'   \item [find_grammar()], [analyze_tf_grammar()], [summarize_tf_grammar()],
#'     and [plot_tf_grammar()] characterize TF-pair spatial grammar.
#'   \item [plot_PROACT_network()] visualizes inferred TF association networks.
#' }
#'
#' @keywords internal
"_PACKAGE"

# Centralized package imports
#' @importFrom utils read.delim read.table read.csv globalVariables
#' @importFrom stats setNames na.omit chisq.test p.adjust
#' @importFrom prioritylasso prioritylasso
#' @importFrom BiocParallel SerialParam SnowParam bplapply
#' @importFrom GenomicRanges GRanges GRangesList findOverlaps resize trim strand
#' @importFrom IRanges IRanges start end
#' @importFrom GenomeInfoDb Seqinfo seqnames seqlevels "seqlevels<-" seqlengths seqinfo "seqinfo<-"
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom ggplot2 ggplot aes geom_col geom_vline geom_text theme_classic labs ylim scale_size_continuous guide_legend scale_fill_manual theme element_text
#' @importFrom patchwork wrap_plots
#' @importFrom dplyr select filter bind_rows transmute group_by summarise slice_max ungroup
#' @importFrom igraph graph_from_data_frame degree V "V<-" cluster_louvain membership
#' @importFrom ggraph ggraph geom_edge_link scale_edge_width scale_edge_color_manual geom_node_point geom_node_text theme_graph
#' @importFrom scales hue_pal
NULL

# Avoid R CMD check notes from non-standard evaluation in ggplot2/dplyr/ggraph.
globalVariables(c(
    "TF1", "TF2", "beta", "Category", "node", "anchor_TF", "max_beta",
    "degree", "anchor_group", "is_anchor", "show_label", "module", "name",
    "dist", "count", "geometry", "prop", "label", "arrangement"
))

# Default redundant motifs. Use remove_motifs = character(0) to disable removal.
DEFAULT_REMOVE_MOTIFS <- c(
    "MA1929.2", "MA1930.2", "MA1951.2", "MA1527.2", "MA0492.2", "MA1631.2", "MA1470.2", "MA1635.2",
    "MA1636.2", "MA1475.2", "MA1494.2", "MA0656.2", "MA1140.3", "MA1639.2", "MA1640.2", "MA1642.2",
    "MA1528.2", "MA1536.2", "MA1537.2", "MA1538.1", "MA1546.2", "MA1549.2", "MA0730.1", "MA0072.2",
    "MA1555.1", "MA1556.1", "MA0829.3", "MA0828.3", "MA0810.2", "MA0872.1", "MA0812.2", "MA0813.1",
    "MA0814.3", "MA0815.1", "MA1570.1", "MA1575.2", "MA1576.2"
)


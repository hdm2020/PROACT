#' Plot a PROACT TF-TF network
#'
#' Visualize PROACT candidate TF pairs as an undirected network. When anchor TFs
#' are supplied, partner TFs are assigned to the anchor with the largest pair
#' coefficient. Without anchors, Louvain community detection defines network
#' modules.
#'
#' @param pair_candidates A data frame containing `TF1`, `TF2`, and `beta`.
#'   Output from [run_PROACT()] can be supplied directly because its `prior`
#'   column is converted to edge categories automatically. Alternatively, a
#'   pre-existing `Category` column may be supplied.
#' @param anchor_tfs Optional character vector of anchor TFs to highlight.
#'
#' @return A ggraph/ggplot object representing the TF association network.
#'
#' @details When `Category` is absent, `prior > 0` is labeled
#'   `"Prior-supported"` and `prior == 0` is labeled `"No prior"`. Existing
#'   `"Known PPI"` and `"Unknown PPI"` categories remain supported for backward
#'   compatibility.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' p <- plot_PROACT_network(pair_candidates, anchor_tfs = c("JUN", "CEBPB"))
#' p
#' }
plot_PROACT_network <- function(pair_candidates, anchor_tfs = NULL){
    # 1. Check input and edge annotation
    required_cols <- c("TF1", "TF2", "beta")
    missing_cols <- setdiff(required_cols, colnames(pair_candidates))
    if(length(missing_cols) > 0) stop("pair_candidates is missing columns: ", paste(missing_cols, collapse = ", "))

    if(!"Category" %in% colnames(pair_candidates)){
        if(!"prior" %in% colnames(pair_candidates)){
            stop("pair_candidates must contain either a 'prior' column from run_PROACT() or a 'Category' column.")
        }
        pair_candidates$Category <- ifelse(!is.na(pair_candidates$prior) & pair_candidates$prior > 0, "Prior-supported", "No prior")
    }

    edges <- pair_candidates |>
        select(TF1, TF2, beta, Category) |>
        filter(!is.na(TF1), !is.na(TF2), !is.na(beta), !is.na(Category))

    if(nrow(edges) == 0) stop("No valid TF pairs found in pair_candidates.")

    edge_colors <- c(
        "Prior-supported" = "#E64B35", "No prior" = "#B0B8B4",
        "Known PPI" = "#E64B35", "Unknown PPI" = "#B0B8B4"
    )
    extra_categories <- setdiff(unique(edges$Category), names(edge_colors))
    if(length(extra_categories) > 0){
        extra_colors <- hue_pal()(length(extra_categories))
        names(extra_colors) <- extra_categories
        edge_colors <- c(edge_colors, extra_colors)
    }
    edge_colors <- edge_colors[unique(edges$Category)]

    # 2. Construct graph
    g <- graph_from_data_frame(d = edges, directed = FALSE)
    V(g)$degree <- degree(g)

    # 3. Determine whether anchor mode should be used
    if(!is.null(anchor_tfs)){
        anchor_tfs <- unique(anchor_tfs[!is.na(anchor_tfs) & anchor_tfs != ""])
        anchor_in_network <- anchor_tfs[anchor_tfs %in% V(g)$name]
    } else {
        anchor_in_network <- character(0)
    }
    use_anchor <- length(anchor_in_network) > 0

    # 4. Anchor mode
    if(use_anchor){
        message("Anchor mode: ", length(anchor_in_network), " anchor TFs detected in network.")

        anchor_partner_links <- bind_rows(
            edges |>
                filter(TF1 %in% anchor_in_network, !TF2 %in% anchor_in_network) |>
                transmute(node = TF2, anchor_TF = TF1, beta = beta),
            edges |>
                filter(TF2 %in% anchor_in_network, !TF1 %in% anchor_in_network) |>
                transmute(node = TF1, anchor_TF = TF2, beta = beta)
        )

        anchor_partner_links2 <- anchor_partner_links |>
            group_by(node, anchor_TF) |>
            summarise(max_beta = max(beta, na.rm = TRUE), .groups = "drop")

        partner_assignment <- anchor_partner_links2 |>
            group_by(node) |>
            slice_max(order_by = max_beta, n = 1, with_ties = FALSE) |>
            ungroup()
        partner_to_anchor <- setNames(partner_assignment$anchor_TF, partner_assignment$node)

        V(g)$is_anchor <- V(g)$name %in% anchor_in_network
        V(g)$anchor_group <- NA_character_
        V(g)$anchor_group[V(g)$is_anchor] <- V(g)$name[V(g)$is_anchor]

        non_anchor_idx <- which(!V(g)$is_anchor)
        V(g)$anchor_group[non_anchor_idx] <- partner_to_anchor[V(g)$name[non_anchor_idx]]
        V(g)$anchor_group[is.na(V(g)$anchor_group)] <- "Unassigned"

        anchor_colors <- hue_pal(h = c(15, 375), c = 100, l = 60)(length(anchor_in_network))
        names(anchor_colors) <- anchor_in_network
        node_colors <- c(anchor_colors, "Unassigned" = "grey80")
        V(g)$anchor_group <- factor(V(g)$anchor_group, levels = c(anchor_in_network, "Unassigned"))

        prior_tfs <- edges |>
            filter(Category %in% c("Prior-supported", "Known PPI")) |>
            select(TF1, TF2) |>
            unlist(use.names = FALSE) |>
            unique()
        V(g)$has_prior <- V(g)$name %in% prior_tfs
        V(g)$show_label <- V(g)$is_anchor | V(g)$has_prior

        set.seed(42)
        p_net <- ggraph(g, layout = "fr") +
            geom_edge_link(aes(edge_width = beta, edge_color = Category), alpha = 0.8) +
            scale_edge_width(range = c(0.3, 1.3), name = "Lasso Beta") +
            scale_edge_color_manual(values = edge_colors, name = "Prior Annotation") +
            geom_node_point(aes(size = degree, fill = anchor_group), shape = 21, color = "#e0e0e0", stroke = 0.5, alpha = 1) +
            geom_node_point(aes(filter = is_anchor, size = degree), shape = 21, fill = NA, color = "black", stroke = 0.8, alpha = 1, show.legend = FALSE) +
            scale_size_continuous(
                range = c(5, 15), name = "Node Degree",
                guide = guide_legend(override.aes = list(shape = 21, fill = "grey70", color = "grey30", stroke = 0.8, alpha = 1))
            ) +
            scale_fill_manual(values = node_colors, breaks = anchor_in_network, labels = anchor_in_network, name = "Anchor TF") +
            geom_node_text(aes(filter = show_label, label = name), repel = TRUE, size = 4.2, fontface = "bold", color = "black", max.overlaps = Inf) +
            theme_graph(base_family = "sans") +
            labs(
                title = "TF Network",
                subtitle = paste0("Anchor TFs and their associated partners share colors; ", "node size indicates degree centrality.")
            ) +
            theme(
                plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30"),
                legend.position = "right"
            )

    # 5. No-anchor mode: Louvain modules
    } else {
        message("No valid anchor TF supplied; using Louvain community mode.")

        communities <- cluster_louvain(g)
        V(g)$module <- as.character(membership(communities))
        module_levels <- sort(unique(V(g)$module))
        V(g)$module <- factor(V(g)$module, levels = module_levels)

        module_colors <- hue_pal()(length(module_levels))
        names(module_colors) <- module_levels

        set.seed(42)
        p_net <- ggraph(g, layout = "fr") +
            geom_edge_link(aes(edge_width = beta, edge_color = Category), alpha = 0.7) +
            scale_edge_width(range = c(0.8, 2.5), name = "Lasso Beta") +
            scale_edge_color_manual(values = edge_colors, name = "Prior Annotation") +
            geom_node_point(aes(size = degree, fill = module), shape = 21, color = "white", stroke = 1.2, alpha = 0.9) +
            scale_size_continuous(
                range = c(5, 15), name = "Node Degree",
                guide = guide_legend(override.aes = list(shape = 21, fill = "grey70", color = "grey30", stroke = 0.8, alpha = 1))
            ) +
            scale_fill_manual(values = module_colors, name = "Algorithm Module") +
            geom_node_text(aes(label = name), repel = TRUE, size = 4.5, fontface = "bold", color = "black", max.overlaps = Inf) +
            theme_graph(base_family = "sans") +
            labs(title = "TF Network") +
            theme(
                plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
                legend.position = "right"
            )
    }

    p_net
}


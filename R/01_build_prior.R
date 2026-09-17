# 1. Pair helper
make_pair <- function(tf1, tf2){
    paste(sort(c(tf1, tf2)), collapse = "_")
}

# 2. Batch pair generation
create_ppi_pairs <- function(df, col1, col2){
    if(nrow(df) == 0) return(character(0))
    apply(df[, c(col1, col2), drop = FALSE], 1, function(x) make_pair(x[1], x[2]))
}

# 3. Load PPI/prior evidence from four databases
load_ppi_database <- function(ppi_dir, jaspar_tfs){
    message("Loading PPI databases ...")

    # BioGRID
    biogrid <- read.delim(file.path(ppi_dir, "BioGRID.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    #biogrid <- biogrid[biogrid$`Organism Name Interactor A` == "Homo sapiens" & biogrid$`Organism Name Interactor B` == "Homo sapiens", ]
    #biogrid <- biogrid[biogrid$`Experimental System Type` == "physical", c("#BioGRID Interaction ID", "Official Symbol Interactor A", "Official Symbol Interactor B")]
    colnames(biogrid) <- c("BioGRID_ID", "Protein1", "Protein2")
    biogrid[biogrid$Protein1 == "T", "Protein1"] <- "TBXT"
    biogrid[biogrid$Protein2 == "T", "Protein2"] <- "TBXT"
    biogrid <- biogrid[biogrid$Protein1 %in% jaspar_tfs & biogrid$Protein2 %in% jaspar_tfs, ]
    biogrid_pairs <- create_ppi_pairs(biogrid, "Protein1", "Protein2")

    # STRING
    stringppi <- read.delim(file.path(ppi_dir, "STRING.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    stringppi <- stringppi[stringppi$experimentally_determined_interaction >= 0.15, ]
    stringppi <- stringppi[stringppi$node1 %in% jaspar_tfs & stringppi$node2 %in% jaspar_tfs, ]
    string_pairs <- create_ppi_pairs(stringppi, "node1", "node2")

    # RF2PPI
    rf2ppi <- read.delim(file.path(ppi_dir, "RF2PPI.tsv"), header = TRUE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
    rf2ppi <- rf2ppi[rf2ppi$Name1 %in% jaspar_tfs & rf2ppi$Name2 %in% jaspar_tfs, ]
    rf2_pairs <- create_ppi_pairs(rf2ppi, "Name1", "Name2")

    # CAP-SELEX
    capselex <- read.delim(file.path(ppi_dir, "CAPSELEX.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    capselex <- capselex[capselex$TF1 %in% jaspar_tfs & capselex$TF2 %in% jaspar_tfs, ]
    cap_pairs <- create_ppi_pairs(capselex, "TF1", "TF2")

    list(BioGRID = unique(biogrid_pairs), STRING = unique(string_pairs), RF2PPI = unique(rf2_pairs), CAPSELEX = unique(cap_pairs))
}

# 4. Assign prior levels
# Pairs supported by >=2 databases: priority 1; remaining supported pairs: priority 2.
classify_ppi_priority <- function(ppi_list){
    all_pairs <- unique(unlist(ppi_list))
    pair_count <- table(unlist(ppi_list))
    block1 <- names(pair_count[pair_count >= 2])
    block2 <- setdiff(all_pairs, block1)
    list(block1 = block1, block2 = block2)
}

# 5. Build prior matrix
build_prior_matrix <- function(tfs, block1, block2){
    prior <- matrix(0, nrow = length(tfs), ncol = length(tfs), dimnames = list(tfs, tfs))

    add_pairs <- function(prior, pair_vec, value){
        for(pair in pair_vec){
            tf_pair <- strsplit(pair, "_")[[1]]
            tf1 <- tf_pair[1]
            tf2 <- tf_pair[2]
            if(tf1 %in% tfs && tf2 %in% tfs){
                prior[tf1, tf2] <- value
                prior[tf2, tf1] <- value
            }
        }
        prior
    }

    prior <- add_pairs(prior, block2, 2)
    prior <- add_pairs(prior, block1, 1)
    diag(prior) <- 0
    prior
}

# 6. Public API
#' Build a TF prior matrix
#'
#' Construct a symmetric transcription factor (TF) prior matrix from BioGRID,
#' STRING, RF2PPI, and CAP-SELEX evidence. TF pairs supported by at least two
#' sources are assigned priority 1, pairs supported by one source are assigned
#' priority 2, and all remaining pairs are assigned 0.
#'
#' @param jaspar_tfs Character vector of TF names to include in the prior matrix.
#' @param ppi_dir Directory containing `BioGRID.tsv`, `STRING.tsv`,
#'   `RF2PPI.tsv`, and `CAPSELEX.tsv`.
#'
#' @return A symmetric numeric matrix with TF names as row and column names.
#'   Values are 1 for high-priority prior pairs, 2 for lower-priority prior
#'   pairs, and 0 when no prior evidence is assigned.
#'
#' @details CAP-SELEX provides DNA-guided TF-pair evidence and is treated as
#'   prior evidence rather than direct proof of a physical interaction in vivo.
#'   The returned matrix is used to define predictor blocks in [run_PROACT()].
#'
#' @export
#'
#' @examples
#' \dontrun{
#' tfs <- c("JUN", "FOS", "CEBPB")
#' prior_mat <- build_prior(tfs, ppi_dir = "PPI_data/human")
#' }
build_prior <- function(jaspar_tfs, ppi_dir){
    ppi_list <- load_ppi_database(ppi_dir = ppi_dir, jaspar_tfs = jaspar_tfs)
    blocks <- classify_ppi_priority(ppi_list)
    prior_mat <- build_prior_matrix(tfs = jaspar_tfs, block1 = blocks$block1, block2 = blocks$block2)
    message("Block1 pairs: ", length(blocks$block1))
    message("Block2 pairs: ", length(blocks$block2))
    prior_mat
}

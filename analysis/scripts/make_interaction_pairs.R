#!/usr/bin/env Rscript

# ======================== #
# make_interaction_pairs.R #
# ======================== #

# ----------- #
# Description #
# ----------- # ----------------------------------------------------------- #
# Script takes input file of interaction outcome matrix
# and generates a table listing results of all pair contests.
# The table is annotated with clade as well as phylogenetic distance,
# which is brought in from another file pre-computed.
#
# Assets produced as output:
# 1. analysis/data/interaction_pairs.csv
# ------------------------------------------------------------------------ #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
source("analysis/scripts/comp_utils.R")


# ------------- #
# function defs #
# ------------- #


parse_outcome <- function(ij, ji, str_i, str_j) {
    if (any(is.na(ij), is.na(ji))) {
        return("NA")
    }
    if (ij > ji) {
        return(str_j)
    }
    if (ij < ji) {
        return(str_i)
    }
    if (sum(ij, ji) == 0) {
        return("RNI")
    }
    if (ij == ji) {
        return("RI")
    }
}


pair_compare <- function(cmat, strains) {
    res <- list()
    for (i in seq_along(strains)) {
        for (j in seq_along(strains)) {
            i_str <- strains[i]
            j_str <- strains[j]
            outcome <- parse_outcome(cmat[i, j], cmat[j, i], i_str, j_str)
            str_names <- sort(c(i_str, j_str))
            res <- rbind(res, c("i" = str_names[1], "j" = str_names[2], "outcome" = outcome))
        }
    }
    results <- data.frame(unique(res))
    results$i <- sapply(results$i, function(x) gsub("^X", "", x))
    results$j <- sapply(results$j, function(x) gsub("^X", "", x))
    results$outcome <- sapply(results$outcome, function(x) gsub("^X", "", x))
    return(results)
}


main <- function(args) {
    cmat <- read_cmat(args$cfile)
    stopifnot(all(row.names(cmat) == names(cmat)))
    strains <- row.names(cmat)
    outcomes <- pair_compare(cmat, strains)
    phylodist <- readr::read_delim(args$phylodist, delim = "\t", col_types = readr::cols())
    phylodist %>% 
        dplyr::rename("i" = strain) %>%
        tidyr::pivot_longer(cols = names(phylodist)[!names(phylodist) %in% "strain"],
                            values_to = "pdist",
                            names_to = "j") ->
        phylodist_long
    strain_meta <- readr::read_csv(args$strain_meta, col_types = readr::cols()) %>%
        dplyr::select(-phylo_pos)
    outcomes %>%
        dplyr::left_join(phylodist_long, by = c("i", "j")) %>%
        dplyr::left_join(strain_meta, by = c("i" = "strain_id")) %>%
        dplyr::rename("i_clade" = clade) %>%
        dplyr::left_join(strain_meta, by = c("j" = "strain_id")) %>%
        dplyr::rename("j_clade" = clade) %>%
        dplyr::mutate(WB = ifelse(i_clade == j_clade, "W", "B"),
                      RI = ifelse(outcome == "RI", 1, 0),
                      RNI = ifelse(outcome == "RNI", 1, 0),
                      ASYM = ifelse(!outcome %in% c("RI", "RNI"), 1, 0)) %>%
        dplyr::filter(i != j,
                      outcome != "NA") %>%
        dplyr::arrange(i, j) ->
        outcomes_annot
    outcomes_annot %>%
        readr::write_csv(args$outfile)
}


# ==== #
# main #
# ==== #

"make_interaction_pairs.R

Usage:
    make_interaction_pairs.R [--help]
    make_interaction_pairs.R <cfile> <strain_meta> <phylodist> <outfile>

Arguments:
    cfile           Competitive outcomes matrix (txt)
    strain_meta     Inhibition outcomes matrix (txt)
    phylodist       Traits file for compiled traits (csv)
    outfile         Path to output file (csv)
" -> doc

args <- list()
args$cfile <- "analysis/data/c_matrix.txt"
args$strain_meta <- "analysis/data/strain_metadata.csv"
args$phylodist <- "analysis/data/phylogenetic_distance.txt"
args$outfile <- "analysis/data/interaction_pairs2.txt"

debug_status <- FALSE
arguments <- run_args_parse(args, debug_status, doc)

main(arguments)
#!/usr/bin/env Rscript

# ================== #
# calc_comp_traits.R #
# ================== #

# ----------- #
# Description #
# ----------- # ----------------------------------------------------------- #
# Script takes input files of competitive outcomes
# and calculates competitive trait values
# for exploitative (C_o, C_d, C_w)
# and interference (C_r, C_t, I_w) competition.
#
# Output element is .csv with all computed traits for all strains.
# ------------------------------------------------------------------------- #

# ------ #
# header #
# ------ #

suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
source("analysis/scripts/comp_utils.R")

# ------------- #
# function defs #
# ------------- #

main <- function(arguments) {

  # load comp matrix
  cmat <- read_cmat(arguments$cfile)

  # load interference matrix
  imat <- read_imat(arguments$ifile)

  # calculate indexes of competitiveness
  cmat_res <- c_calc(cmat)
  imat_res <- i_calc(imat)
  res_full <- dplyr::full_join(cmat_res,
                               imat_res,
                               by = "strain_id")

  write.table(res_full,
      file = file.path(arguments$outfile),
      sep = ",",
      col.names = T,
      row.names = F,
      quote = F)
}

# ==== #
# main #
# ==== #

"calc_comp_traits.R

Usage:
    calc_comp_traits.R [--help]
    calc_comp_traits.R <cfile> <ifile> <outfile>

Arguments:
    cfile        Competitive outcomes matrix (txt)
    ifile        Inhibition outcomes matrix (txt)
    outfile      Full outfile path for compiled competitive traits data (csv)
" -> doc

args <- list()
args$cfile   <- "analysis/data/c_matrix.txt"
args$ifile   <- "analysis/data/i_matrix.txt"
args$outfile <- "analysis/data/comp_traits.csv"

debug_status <- FALSE
arguments <- run_args_parse(args, debug_status, doc)

main(arguments)

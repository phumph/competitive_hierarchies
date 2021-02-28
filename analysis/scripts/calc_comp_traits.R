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

run_args_parse <- function(debug_status) {
  if (debug_status == TRUE) {
    arguments <- list()
    arguments$cfile   <- "analysis/data/c_matrix.txt"
    arguments$ifile   <- "analysis/data/i_matrix.txt"
    arguments$outfile <- "analysis/data/comp_traits.csv"
  } else if (debug_status == FALSE) {
    args <- commandArgs(trailingOnly = FALSE)
    arguments$cfile   <- args[1]
    arguments$ifile   <- args[2]
    arguments$outfile <- args[3]
  }
  return(arguments)
}


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

arguments <- run_args_parse(debug_status = TRUE)
res <- main(arguments)

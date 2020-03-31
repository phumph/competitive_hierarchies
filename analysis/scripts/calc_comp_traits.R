#!/usr/bin/env Rscript

# ------ #
# header #
# ------ #

library(dplyr)
library(tidyr)
source(file.path("scripts/utils.r"))

# ------------- #
# function defs #
# ------------- #

run_args_parse <- function(debug_status) {
  if (debug_status == TRUE) {
    arguments <- list()
    arguments$cfile   <- "data/c_matrix.txt"
    arguments$ifile   <- "data/i_matrix.txt"
    arguments$outfile <- "data/comp_traits.csv"
  } else if (debug_status == FALSE) {
    args <- commandArgs(trailingOnly = FALSE)
    arguments$cfile   <- args[1]
    arguments$ifile   <- args[2]
    arguments$outfile <- args[3]
  }
  return(arguments)
}

c_calc <- function(cmat) {
  cres <- data.frame(strain_id = row.names(cmat),
                      n_o = NA, n_d = NA, c_o = NA, c_d = NA, c_w = NA)
  # calculate the number of interactions in the defensive and offense positions:
  for (i in seq_len(nrow(cmat))) {
    # offense
    cres$n_o[i] <- sum(!is.na(cmat[, i])) - 1
    # defense
    cres$n_d[i] <- sum(!is.na(cmat[i, ])) - 1
  }
  cres$c_d_raw <- cres$n_d - rowSums(cmat, na.rm = TRUE)
  cres$c_o_raw <- colSums(cmat, na.rm = TRUE)
  cres$c_o     <- round(cres$c_o_raw / cres$n_o, 4)
  cres$c_d     <- round(cres$c_d_raw / cres$n_d, 4)
  cres$c_w     <- round(cres$c_o - (1 - cres$c_d), 4)
  return(cres)
}

i_calc <- function(imat) {
  ires <- data.frame(strain_id = row.names(imat),
                      n_t = NA, n_r = NA, c_r = NA, c_t = NA, i_w = NA)
  # subset to only those cols which have colSum > 0
  ires$c_t_raw <- colSums(imat, na.rm = TRUE)
  # sum num. strains each is killed by;
  # subtract from number of cT>0 strains to get cR);
  # define submatrix with strains with cT>0 as cols
  imat2 <- imat[, names(imat) %in% dplyr::filter(ires, c_t_raw > 0)[, "strain_id"]]
  ires$n_t <- c(rep(1, (nrow(imat2) - ncol(imat2))), colSums(!is.na(imat2)) - 1)
  ires$n_r <- rowSums(!is.na(imat2))
  ires$n_r[ires$strain_id %in% ires$strain_id[ires$c_t_raw > 0]] <-
    ires$n_r[ires$strain_id %in% ires$strain_id[ires$c_t_raw > 0]] - 1
  # tabulate raw interference score counts
  ires$c_r_raw <- ires$n_r - rowSums(imat2, na.rm = TRUE)
  ires$c_t_raw <- c(rep(0, (nrow(imat2) - ncol(imat2))),
                    colSums(imat2, na.rm = TRUE))
  # divide to recover proportions
  ires$c_r <- round(ires$c_r_raw / ires$n_r, 4)
  ires$c_t <- round(ires$c_t_raw / ires$n_t, 4)
  ires$i_w <- round((ires$c_t - 1) + ires$c_r, 4)
  ires$i_w <- round(ires$c_t - (1 - ires$c_r), 4)
  return(ires)
}

main <- function(arguments) {

  # setup comp matrix
  cmat <- read.table(file = file.path(arguments$cfile),
  row.names = 1, header = T, sep = "\t")
  row.names(cmat) <- names(cmat)
  diag(cmat) <- 0
  cmat <- cmat / 2

  # setup interference matrix
  imat <- read.table(file = file.path(arguments$ifile),
  row.names = 1, header = T, sep = "\t")
  row.names(imat) <- names(imat)
  diag(imat) <- 0

  # calculate indexes of competitiveness
  cmat_res <- c_calc(cmat)
  imat_res <- i_calc(imat)
  res_full <- dplyr::full_join(cmat_res, imat_res, by = "strain_id")

  return(res_full)
}

arguments <- run_args_parse(debug_status = TRUE)
res <- main(arguments)

write.table(res,
  file = file.path(arguments$outfile),
  sep = ",",
  col.names = T,
  row.names = F,
  quote = F)

#!/usr/bin/env Rscript

# ==================== #
# make_growth_curves.R #
# ==================== #

# ----------- #
# Description #
# ----------- # ----------------------------------------------------------- #
# Script takes input files of in vitro growth trajectories
# and plots strain-level growth curves.
#
# Output element is .png with all strains' growth trajectories plotted.
# ------------------------------------------------------------------------- #


# ------------- #
# function defs #
# ------------- #

run_args_parse <- function(debug_status) {
  if (debug_status == TRUE) {
    arguments <- list()
    arguments$g1_file   <- "data/MM_gcurvedata_Pflu.txt"
    arguments$g2_file   <- "data/MM_gcurvedata_Psyr.txt"
    arguments$plot_file <- "figs/growth_curves.pdf"
  } else if (debug_status == FALSE) {
    args <- commandArgs(trailingOnly = FALSE)
    arguments <- list()
    arguments$g1_file   <- args[1]
    arguments$g2_file   <- args[2]
    arguments$plot_file <- args[3]
  }
  return(arguments)
}

make_plot <- function(g1, g2) {

  # row 1
  par(mfrow = c(5, 3), mai = c(0.3, 0.4, 0.3, 0.3))
  plot(g1[, 1], g1[, 15], type = "l", lty = 1,
       xlim = c(0, 3600),
       ylim = c(0, 0.6),
       las = 2,
       cex.axis = 0.6)
  points(g1[, 1], g1[, 14], type = "l", lty = 2)
  points(g1[, 1], g1[, 13], type = "l", lty = 3)
  legend("topleft",
         c(dimnames(g1)[[2]][15],
           dimnames(g1)[[2]][14],
           dimnames(g1)[[2]][13]),
         cex = 0.66, lty = c(1, 2, 3), bty = "n")

  # row 2
  plot(g1[, 1], g1[, 12], type = "l", lty = 1,
       xlim = c(0, 3600),
       ylim = c(0, 0.6),
       las = 2,
       cex.axis = 0.66)
  points(g1[, 1], g1[, 11], type = "l", lty = 2)
  legend("topleft",
         c(dimnames(g1)[[2]][12],
           dimnames(g1)[[2]][11]),
         cex = 0.66, lty = c(1, 2), bty = "n")

  # row 3
  plot(g1[, 1], g1[, 10], type = "l", lty = 1,
       xlim = c(0, 3600),
       ylim = c(0, 0.6),
       las = 2,
       cex.axis = 0.66)
  points(g1[, 1], g1[, 9], type = "l", lty = 2)
  legend("topleft",
         c(dimnames(g1)[[2]][10],
           dimnames(g1)[[2]][9]),
         cex = 0.66, lty = c(1, 2), bty = "n")

  # row 4
  plot(g1[, 1], g1[, 7], type = "l", lty = 1,
       xlim = c(0, 3600),
       ylim = c(0, 0.6),
       las = 2,
       cex.axis = 0.66)
  points(g1[, 1], g1[, 6], type = "l", lty = 2)
  points(g1[, 1], g1[, 5], type = "l", lty = 3)
  points(g1[, 1], g1[, 4], type = "l", lty = 4)
  legend("topleft",
         c(dimnames(g1)[[2]][7],
           dimnames(g1)[[2]][6],
           dimnames(g1)[[2]][5],
           dimnames(g1)[[2]][4]),
         cex = 0.66, lty = c(1, 2, 3, 4), bty = "n")

  # row 5
  plot(g1[, 1], g1[, 3], type = "l", lty = 1,
       xlim = c(0, 3600),
       ylim = c(0, 0.6),
       las = 2,
       cex.axis = 0.66)
  points(g1[, 1], g1[, 2], type = "l", lty = 2)
  legend("topleft",
         c(dimnames(g1)[[2]][3],
           dimnames(g1)[[2]][2]),
         cex = 0.66, lty = c(1, 2), bty = "n")

  # row 6
  plot(g1[, 1], g1[, 8], type = "l", lty = 1,
       xlim = c(0, 3600),
       ylim = c(0, 0.6),
       las = 2,
       cex.axis = 0.66)
  legend("topleft",
         c(dimnames(g1)[[2]][8]),
         cex = 0.66, lty = c(1), bty = "n")

  # Now for P.syringae plots
  # first plot: 10A,  21B,  20A,  17A
  plot(g2[, 1], g2[, 9], type = "l", lty = 1,
       xlim = c(0, 3600),
       ylim = c(0, 0.8),
       las = 2,
       cex.axis = 0.6)
  points(g2[, 1], g2[, 11], type = "l", lty = 2)
  points(g2[, 1], g2[, 12], type = "l", lty = 3)
  points(g2[, 1], g2[, 10], type = "l", lty = 4)
  legend("topleft",
         c(dimnames(g2)[[2]][9],
           dimnames(g2)[[2]][11],
           dimnames(g2)[[2]][12],
           dimnames(g2)[[2]][10]),
         cex = 0.66, lty = c(1, 2, 3, 4), bty = "n")

  # second plot: 24A,  39C,  14B,  16A
  plot(g2[, 1], g2[, 13], type = "l", lty = 1,
       xlim = c(0, 3600),
       ylim = c(0, 0.8),
       las = 2,
       cex.axis = 0.6)
  points(g2[, 1], g2[, 14], type = "l", lty = 2)
  points(g2[, 1], g2[, 15], type = "l", lty = 3)
  points(g2[, 1], g2[, 16], type = "l", lty = 4)
  legend("topleft",
         c(dimnames(g2)[[2]][13],
           dimnames(g2)[[2]][14],
           dimnames(g2)[[2]][15],
           dimnames(g2)[[2]][16]),
         cex = 0.66, lty = c(1, 2, 3, 4), bty = "n")

  # thid plot: 46B 19A
  plot(g2[, 1], g2[, 25], type = "l", lty = 1,
       xlim = c(0, 3600),
       ylim = c(0, 0.8),
       las = 2,
       cex.axis = 0.6)
  points(g2[, 1], g2[, 26], type = "l", lty = 2)
  legend("topleft",
         c(dimnames(g2)[[2]][25],
           dimnames(g2)[[2]][26]),
         cex = 0.66, lty = c(1, 2), bty = "n")

  # fourth plot: 06C,  02A,  08A,  Psm4326
  plot(g2[, 1], g2[, 24], type = "l", lty = 1,
       xlim = c(0, 3600),
       ylim = c(0, 0.8),
       las = 2,
       cex.axis = 0.6)
  points(g2[, 1], g2[, 22], type = "l", lty = 2)
  points(g2[, 1], g2[, 23], type = "l", lty = 3)
  points(g2[, 1], g2[, 21], type = "l", lty = 4)
  legend("topleft",
         c(dimnames(g2)[[2]][24],
           dimnames(g2)[[2]][22],
           dimnames(g2)[[2]][23],
           dimnames(g2)[[2]][21]),
         cex = 0.66, lty = c(1, 2, 3, 4), bty = "n")

  # fifth plot: 22D,  07A,  22B
  plot(g2[, 1], g2[, 19], type = "l", lty = 1,
       xlim = c(0, 3600),
       ylim = c(0, 0.8),
       las = 2,
       cex.axis = 0.6)
  points(g2[, 1], g2[, 20], type = "l", lty = 2)
  points(g2[, 1], g2[, 18], type = "l", lty = 3)
  legend("topleft",
         c(dimnames(g2)[[2]][19],
           dimnames(g2)[[2]][20],
           dimnames(g2)[[2]][18]),
         cex = 0.66, lty = c(1, 2, 3), bty = "n")

  # sixth plot: 26B,  39F
  plot(g2[, 1], g2[, 17], type = "l", lty = 1,
       xlim = c(0, 3600),
       ylim = c(0, 0.8),
       las = 2,
       cex.axis = 0.6)
  points(g2[, 1], g2[, 2], type = "l", lty = 2)
  legend("topleft",
         c(dimnames(g2)[[2]][17],
           dimnames(g2)[[2]][2]),
         cex = 0.66, lty = c(1, 2), bty = "n")

  # seventh plot: 08C,  23A
  plot(g2[, 1], g2[, 7], type = "l", lty = 1,
       xlim = c(0, 3600),
       ylim = c(0, 0.8),
       las = 2,
       cex.axis = 0.6)
  points(g2[, 1], g2[, 8], type = "l", lty = 2)
  legend("topleft",
         c(dimnames(g2)[[2]][7],
           dimnames(g2)[[2]][8]),
         cex = 0.66, lty = c(1, 2), bty = "n")

  # eigth plot: 21A,  20B
  plot(g2[, 1], g2[, 4], type = "l", lty = 1,
       xlim = c(0, 3600),
       ylim = c(0, 0.8),
       las = 2,
       cex.axis = 0.6)
  points(g2[, 1], g2[, 5], type = "l", lty = 2)
  legend("topleft",
         c(dimnames(g2)[[2]][4],
           dimnames(g2)[[2]][5]),
         cex = 0.66, lty = c(1, 2), bty = "n")

  # ninth and last plot: 08B
  plot(g2[, 1], g2[, 6], type = "l", lty = 1,
       xlim = c(0, 3600),
       ylim = c(0, 0.8),
       las = 2,
       cex.axis = 0.6)
  legend("topleft",
         c(dimnames(g2)[[2]][6]),
         cex = 0.66, lty = c(1), bty = "n")
}


# -------- #
# main def #
# -------- #

main <- function(arguments) {

  g1 <- read.table(arguments$g1_file, header = T, sep = "\t")
  g2 <- read.table(arguments$g2_file, header = T, sep = "\t")

  pdf(file = arguments$plot_file, width = 8, height = 11)
    make_plot(g1, g2)
  dev.off()
}

# ==== #
# main #
# ==== #

arguments <- run_args_parse(debug_status = TRUE)
main(aguments)

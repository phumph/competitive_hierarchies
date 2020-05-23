#!/usr/bin/env Rscript

# ==================== #
# calc_facil_effects.R #
# ==================== #

# ----------- #
# Description #
# ----------- # ----------------------------------------------------------- #
# Script takes input files of competitive outcomes
# and calcuates competitive trait values
# for exploitative competition (C_o, C_d, C_w)
# in the context of each toxin-secreting third party.

# Output element is .csv with all delta C_w values for each strain
# in the context of each toxin producer, as well as the overall effect.
# ------------------------------------------------------------------------- #

# ------ #
# header #
# ------ #

suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(pals)))
source("analysis/scripts/comp_utils.R")

# ------------- #
# function defs #
# ------------- #

run_args_parse <- function(debug_status) {
  if (debug_status == TRUE) {
    arguments <- list()
    arguments$cfile    <- "analysis/data/c_matrix.txt"
    arguments$ifile    <- "analysis/data/i_matrix.txt"
    arguments$tfile    <- "analysis/data/all_traits.txt"
    arguments$outfile  <- "analysis/data/comp_facil_effects.csv"
    arguments$figs_dir <- "analysis/figs/"
  } else if (debug_status == FALSE) {
    args <- commandArgs(trailingOnly = FALSE)
    arguments$cfile    <- args[1]
    arguments$ifile    <- args[2]
    arguments$tfile    <- args[3]
    arguments$outfile  <- args[4]
    arguments$figs_dir <- args[5]
  }
  return(arguments)
}


calc_cw_deltas <- function(cmat, imat, focal_res_cols = c("strain_id", "c_w")) {

  toxin_cols <- colSums(imat, na.rm = TRUE)
  toxin_strains <- names(toxin_cols[toxin_cols > 0])
  imat_focal <- imat[, names(imat) %in% toxin_strains]

  # initialize results data.frame
  res_full <- c_calc(cmat)
  res_full <- res_full[, names(res_full) %in% focal_res_cols]

  for (strain in seq_along(toxin_strains)) {
    i_strain <- toxin_strains[strain]
    i_vec <- imat_focal[, strain]
    cmat_i <- knockout_by_ivec(cmat, i_vec)
    c_calc_i <- c_calc(cmat_i)
    res_tmp <- c_calc_i[, names(c_calc_i) %in% focal_res_cols]
    names(res_tmp)[2] <- i_strain
    res_full <- merge(res_full, res_tmp, by = "strain_id")
    res_full[i_vec == 1, i_strain] <- -1
  }

  deltas <- res_full[, !names(res_full) %in% focal_res_cols] - res_full$c_w
  res_full[, !names(res_full) %in% focal_res_cols] <- deltas
  ncols <- ncol(res_full[, !names(res_full) %in% focal_res_cols])
  res_full$net_delta <- rowSums(deltas)
  res_full$avg_delta <- res_full$net_delta / ncols

  return(res_full)
}


add_clade_info <- function(res_full, infile) {
  tfile <- read.table(file.path(arguments$tfile),
    sep = "\t",
    header = T,
    stringsAsFactors = FALSE
  )
  cols_to_join <- c("strain_id", "clade", "c_r")
  res_full$strain_id <- paste0(
    sapply(res_full$strain_id, function(x) gsub("X", "", x))
  )
  res_full <- merge(res_full, tfile[, cols_to_join], by = "strain_id")
  return(res_full)
}


calculate_rank_diffs <- function(df,
                                 clade = "all",
                                 deltas = TRUE,
                                 base_col = "c_w",
                                 avg_col = "avg_delta") {

  if (clade != "all") {
    stopifnot("clade" %in% names(df))
    stopifnot(clade %in% paste0(unique(df$clade)))
    df <- df[df$clade %in% clade, ]
  }

  xnames <- grep("^X", names(df), value = T)
  rank_res <- data.frame(strain_id = df$strain_id,
                         base_rank = rank(jitter(df[, base_col])))
  col_ranks <- list()

  for (icol in seq_along(xnames)) {
    if (deltas == TRUE) {
      col_ranks[[icol]] <- rank(jitter(df[, base_col])) -
        rank(jitter(df[, base_col] + df[, xnames[icol]]))
    } else {
      col_ranks[[icol]] <- rank(jitter(df[, base_col])) -
        rank(jitter(df[, xnames[icol]]))
    }
  }

  rank_df <- data.frame(do.call(cbind, col_ranks))
  names(rank_df) <- xnames
  tot_rank_impact <- sort(colSums(abs(rank_df)), decreasing = TRUE)
  rank_df <- rank_df[, names(tot_rank_impact)]
  rank_res <- cbind(rank_res, rank_df)
  rank_res$final_rank <- rank(jitter(df[, base_col] + df[, avg_col]))
  rank_res$final_rank_diff <- rank_res$final_rank - rank_res$base_rank
  rank_res <- rank_res[order(rank_res$base_rank), ]

  return(rank_res)
}


plot_facil_summaries <- function(res_full) {
  
  plot1 <- ggplot(res_full, aes(x = c_w, y = avg_delta, fill = clade)) +
    geom_point(pch = 21, col = "black", size = 2) +
    # facet_wrap(~ clade) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.x = element_text(
        size = 7, angle = 90,
        hjust = 1, vjust = 0.5
      ),
      axis.text.y = element_text(size = 7),
      axis.line = element_line(),
      axis.ticks = element_line()
    ) +
    geom_hline(yintercept = 0, lty = 3, col = "gray40") +
    scale_fill_manual(values = c("white", "black")) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    scale_y_continuous(
      limits = c(-0.65, 0.65),
      breaks = seq(-6, 0.6, 0.2)
    )

  plot2 <- ggplot(res_full, aes(x = c_w, y = c_w + avg_delta, fill = clade)) +
    geom_abline(slope = 1, intercept = 0, lty = 3, col = "gray60") +
    geom_point(pch = 21, col = "black", size = 2) +
    # facet_wrap(~ clade) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 7),
      axis.line = element_line(),
      axis.ticks = element_line()
    ) +
    geom_hline(yintercept = 0, lty = 3, col = "gray40") +
    scale_fill_manual(values = c("white", "black")) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5))

  return(
    ggpubr::ggarrange(
      plotlist = list(plot1, plot2),
      ncol = 1, nrow = 2,
      align = "hv"
    )
  )
}


plot_rank_diffs <- function(df) {
  
  factor_order <- grep("^X", names(df), value = TRUE) %>% gsub("^X", "", .)
  strain_id_order <- df$strain_id
  
  df %>%
    tidyr::gather(key = "i_strain", value = "rank_diff",
                  -strain_id,
                  -base_rank,
                  -final_rank,
                  -final_rank_diff) ->
    df2
  
  df2$i_strain <- sapply(df2$i_strain, function(x) gsub("^X", "", paste0(x)))
  df2$i_strain  <- factor(df2$i_strain, levels = factor_order)
  df2$strain_id <- factor(df2$strain_id, levels = rev(strain_id_order))
  
  df2 %>%
    ggplot(aes(x = i_strain, y = strain_id, fill = rank_diff)) +
    geom_tile() +
    #scale_fill_gradientn(colors = RColorBrewer::brewer.pal(11,"BrBG")) +
    scale_fill_gradientn(colors = pals::coolwarm()) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 7,
                                     angle = 90,
                                     hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 7)) ->
    rank_plot_1
  
  df %>%
    ggplot(aes(x = base_rank, y = final_rank, col = final_rank_diff)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, lty = 3) +
    theme_bw() +
    scale_color_gradientn(colors = pals::coolwarm()) ->
    rank_plot_2
  
  return(list(rank_plot_1, rank_plot_2))
}


# ======== #
# main def #
# ======== #

main <- function(arguments) {

  cmat <- read_cmat(arguments$cfile)
  cmat <- cmat[complete.cases(cmat), ]
  cmat <- cmat[, names(cmat) %in% row.names(cmat)]
  cmat <- cmat[order(row.names(cmat)), order(row.names(cmat))]

  # setup interference matrix
  imat <- read_imat(arguments$ifile)
  imat <- imat[complete.cases(imat), ]
  imat <- imat[order(row.names(imat)), order(row.names(imat))]

  stopifnot(all(row.names(imat) %in% row.names(cmat)))
  stopifnot(row.names(imat) == row.names(cmat))

  # perform calculations
  res_full <- calc_cw_deltas(cmat, imat)

  # add clade information
  res_full <- add_clade_info(res_full, arguments$tfile)

  # plot avg delta by clade
  res_full %>%
    plot_facil_summaries() ->
    facil_plots

  facil_plots %>%
    ggsave(filename = paste0(arguments$figs_dir,
                             "facil_plot_summaries.pdf"
                             ),
           width = 3,
           height = 5,
           #dpi = 300,
           device = "pdf"
         )

  # now calculate and plot rank changes
  res_full_rank_psyr <- calculate_rank_diffs(res_full, clade = "Psyr")
  res_full_rank_pflu <- calculate_rank_diffs(res_full, clade = "Pflu")
  
  res_full_rank_psyr %>%
    plot_rank_diffs() ->
    rank_plot_psyr
  
  res_full_rank_pflu %>%
    plot_rank_diffs() ->
    rank_plot_pflu
}


# ==== #
# main #
# ==== #

debug_status <- TRUE
arguments <- run_args_parse(debug_status)
main(arguments)

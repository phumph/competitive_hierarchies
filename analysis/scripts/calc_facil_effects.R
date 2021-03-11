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

  df %>%
    dplyr::arrange(desc(!!as.symbol(base_col))) %>%
    dplyr::mutate(base_rank = c(1:nrow(df))) %>%
    dplyr::select(strain_id, base_rank) ->
    rank_res

  for (icol in seq_along(xnames)) {
    df_tmp <-
      df %>%
      dplyr::select(strain_id, !!as.symbol(base_col), !!as.symbol(xnames[icol]))
    if (deltas == TRUE) {
      df_tmp %>%
        dplyr::mutate(c_w_i = !!as.symbol(base_col) + jitter(!!as.symbol(xnames[icol]))) %>%
        dplyr::arrange(desc(c_w_i)) %>%
        dplyr::mutate(rank = c(1:nrow(df_tmp))) %>%
        dplyr::select(strain_id, rank) ->
        df_tmp
    } else {
      df_tmp %>%
        dplyr::mutate(c_w_i = jitter(!!as.symbol(xnames[icol]))) %>%
        dplyr::arrange(desc(c_w_i)) %>%
        dplyr::mutate(rank = c(1:nrow(df_tmp))) %>%
        dplyr::select(strain_id, rank) ->
        df_tmp
    }

    suppressMessages(suppressWarnings(
        rank_res %>%
          dplyr::left_join(df_tmp, by = "strain_id") %>%
          dplyr::mutate(rank_diff = base_rank - rank) %>%
          dplyr::select(-rank) ->
          rank_res
    ))
    names(rank_res)[names(rank_res) == "rank_diff"] <- xnames[icol]
  }

  tot_rank_impact <- sort(colSums(abs(rank_res[, -c(1,2)])), decreasing = TRUE)
  rank_res <- rank_res[, c("strain_id", "base_rank", names(tot_rank_impact))]

  # add final_rank
  df %>%
    dplyr::select(strain_id,
                  !!as.symbol(base_col),
                  !!as.symbol(avg_col)) %>%
    dplyr::mutate(final_c_w = !!as.symbol(base_col) + !!as.symbol(avg_col)) %>%
    dplyr::arrange(desc(final_c_w)) %>%
    dplyr::mutate(final_rank = c(1:nrow(df))) %>%
    dplyr::select(-!!as.symbol(avg_col)) ->
    df_tmp_2

  rank_res %>%
    dplyr::left_join(df_tmp_2, by = "strain_id") %>%
    dplyr::mutate(final_rank_diff = base_rank - final_rank) ->
    rank_res

  return(rank_res)
}


plot_facil_summaries <- function(res_full) {

  plot1 <- ggplot(res_full, aes(x = c_w, y = avg_delta, fill = clade)) +
    geom_point(pch = 21, col = "black", size = 2) +
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

  factor_order <- grep("^X", names(df), value = TRUE) %>%
    gsub("^X", "", .)
  stopifnot("final_rank_diff" %in% names(df))
  names(df)[names(df) == "final_rank_diff"] <- "final"
  factor_order <- c(factor_order, "final")
  strain_id_order <- df$strain_id

  df %>%
    tidyr::gather(key = "i_strain", value = "rank_diff",
                  -strain_id,
                  -base_rank,
                  -final_rank,
                  -c_w,
                  -final_c_w) ->
    df2

  df2$i_strain <- sapply(df2$i_strain, function(x) gsub("^X", "", paste0(x)))
  df2$i_strain  <- factor(df2$i_strain, levels = factor_order)
  df2$strain_id <- factor(df2$strain_id, levels = rev(strain_id_order))

  value_range <- c(-max(abs(range(df2$rank_diff))),max(abs(range(df2$rank_diff))))
  palette_breaks <- seq(value_range[1], value_range[2], 1)
  color_palette  <- colorRampPalette(
    c("midnightblue", "white", "darkorange2"))(length(palette_breaks) - 1
    )

  df2 %>%
    ggplot(aes(x = i_strain, y = strain_id, fill = rank_diff)) +
    geom_tile(col = "white", lwd = 0.33) +
    scale_fill_gradientn(colors = color_palette, limits = value_range, name = "") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(axis.text.x = element_text(size = 7,
                                     angle = 90,
                                     hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 7)) ->
    rank_plot_1

  df %>%
    ggplot(aes(x = -base_rank, y = -final_rank, fill = final)) +
    geom_point(pch = 21, col = "black", size = 2.25) +
    geom_abline(intercept = 0, slope = 1, lty = 3) +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA),
          axis.text = element_text(size = 7),
          panel.grid = element_blank()) +
    scale_fill_gradientn(colors = color_palette, limits = value_range, name = "") ->
    rank_plot_2

  ggpubr::ggarrange(plotlist = list(rank_plot_1, rank_plot_2),
                    common.legend = T,
                    ncol = 2, widths = c(0.66, 1),
                    align = 'hv',
                    legend = "top",
                    labels = c("a", "b"),
                    vjust = 0.75) ->
    joint_plot

  return(joint_plot)
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
    ggsave(filename = file.path(arguments$figs_dir,
                                "facil_plot_summaries.pdf"),
           width = 2.5,
           height = 4.5,
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

  # save plots:
  rank_plot_psyr %>%
    ggsave(filename = file.path(arguments$figs_dir,
                                "facil_plot_Psyr_ranks.pdf"),
           device = "pdf",
           width = 6,
           height = 4)

  rank_plot_pflu %>%
    ggsave(filename = file.path(arguments$figs_dir,
                                "facil_plot_Pflu_ranks.pdf"),
           device = "pdf",
           width = 6,
           height = 4)
}


# ==== #
# main #
# ==== #


"calc_facil_effects.R

Usage:
    calc_facil_effects.R [--help]
    calc_facil_effects.R <cfile> <ifile> <traits_infile> <outfile> <figs_dir>

Arguments:
    cfile            Competitive outcomes matrix (txt)
    ifile            Inhibition outcomes matrix (txt)
    traits_infile    Traits file for compiled traits (csv)
    outfile          Full outfile path for compiled competitive traits data (csv)
" -> doc

args <- list()
args$cfile <- "analysis/data/c_matrix.txt"
args$ifile <- "analysis/data/i_matrix.txt"
args$traits_infile <- "analysis/data/all_traits.txt"
args$outfile <- "analysis/data/comp_facil_effects.csv"
args$figs_dir <- "analysis/figs"

debug_status <- FALSE
arguments <- run_args_parse(args, debug_status, doc)

main(arguments)


##

run_args_parse <- function(debug_status) {
  if (debug_status == TRUE) {
    arguments <- list()
    arguments$cfile    <- "analysis/data/c_matrix.txt"
    arguments$ifile    <- "analysis/data/i_matrix.txt"
    arguments$tfile    <- "analysis/data/all_traits.txt"
    arguments$outfile  <- "analysis/data/comp_facil_effects.csv"
    arguments$figs_dir <- "analysis/figs"
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

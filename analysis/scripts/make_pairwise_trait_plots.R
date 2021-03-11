#!/usr/bin/env Rscript

# =========================== #
# make_pairwise_trait_plots.R #
# =========================== #

# ----------- #
# Description #
# ----------- # ----------------------------------------------------------- #
# Script takes input files of compiled growth and competition traits
# and generates pairwise plots of all traits,
# broken down by comparison type.
#
# Output element is three .pdfs of each set of comparisons,
# which were manuaully edited for inclusion as figures S2-S4.
# ------------------------------------------------------------------------- #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggpubr)))
source("analysis/scripts/comp_utils.R")


# ============= #
# function defs #
# ============= #

trait_biplot <- function(x, fill_name = "clade") {
  clade_cols <- as.character(c("white", "black"))
  y_name <- names(x)[2]
  x_name <- names(x)[3]
  ggplot(x, aes_string(x = x_name, y = y_name, fill = fill_name)) +
    geom_point(size = 1.5, pch = 21, colour = "black", alpha = 0.66) +
    stat_smooth(method = lm, se = FALSE, size = 0.5, aes(lty = clade)) +
    theme(legend.position = "none") +
    scale_shape_manual(values = c(21, 21)) +
    scale_fill_manual(values = clade_cols) +
    theme_bw() +
    theme(legend.position = "none",
          legend.justification = "none",
          panel.background = element_blank(),
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 7, face = "italic"))
}


main <- function(arguments) {
  file_extn <- ".pdf"
  traits_all <- readr::read_csv(arguments$traits_infile, col_types = readr::cols())
  trait_vec <- names(traits_all)
  trait_vec <- trait_vec[!trait_vec %in% c("strain_id", "clade", "phylo_pos")]

  ## assemble sets of comparisons
  # first, growth comparisons
  growth_comparisons <- data.frame(x = c("r", "r", "L"),
                                   y = c("L", "K", "K"))
  growth_comp_plots <- list()
  for (i in 1:dim(growth_comparisons)[1]) {
    traits_all %>%
      dplyr::select(clade,
                    growth_comparisons[i, 1],
                    growth_comparisons[i, 2]) %>%
      trait_biplot() ->
      growth_comp_plots[[i]]
  }
  growth_plots <- suppressWarnings(
    suppressMessages(
      ggpubr::ggarrange(plotlist = growth_comp_plots,
                        nrow = 1, ncol = 3,
                        align = "hv")
    )
  )
  growth_plots %>%
    ggsave(filename = paste0(arguments$plot_file_base,
                             "_growth",
                             file_extn),
           width = 6.5,
           height = 2,
           #dpi = 300,
           device = "pdf"
         )

  # now competition comparisons
  comp_comparisons <- data.frame(expand.grid(trait_vec[-c(1:3)],
                                             trait_vec[-c(1:3)]))
  comp_comparisons[comp_comparisons[, 1] != comp_comparisons[, 2], ] ->
    comp_comparisons

  competition_comp_plots <- list()
  for (i in 1:dim(comp_comparisons)[1]) {
    traits_all %>%
      dplyr::select(clade,
                    comp_comparisons[i, 2],
                    comp_comparisons[i, 1]) %>%
      trait_biplot() ->
      competition_comp_plots[[i]]
  }
  competition_plots <- suppressWarnings(
    suppressMessages(
      ggpubr::ggarrange(plotlist = competition_comp_plots,
                        nrow = 6, ncol = 5,
                        align = "hv")
    )
  )

  competition_plots %>%
    ggsave(filename = paste0(arguments$plot_file_base,
                             "_comp",
                             file_extn),
           width = 10,
           height = 12,
           #dpi = 300,
           device = "pdf"
         )

  # finally, growth-v-competition comparisons
  growth_comp_comparisons <- expand.grid(c("r", "L", "K"),
                                         c(trait_vec[-c(1:3)]))
  growth_comp_comp_plots <- list()
  for (i in 1:dim(growth_comp_comparisons)[1]) {
    traits_all %>%
      dplyr::select(clade,
                    growth_comp_comparisons[i, 2],
                    growth_comp_comparisons[i, 1]) %>%
      trait_biplot() ->
      growth_comp_comp_plots[[i]]
  }
  growth_comp_plots <- suppressWarnings(
    suppressMessages(
      ggpubr::ggarrange(plotlist = growth_comp_comp_plots,
                        nrow = 6, ncol = 3,
                        align = "hv")
    )
  )

  growth_comp_plots %>%
    ggsave(filename = paste0(arguments$plot_file_base,
                             "_growth_comp",
                             file_extn),
           width = 5,
           height = 9,
           #dpi = 300,
           device = "pdf"
         )
}

# ============ #
# main routine #
# ============ #

"make_pairwise_trait_plots.R

Usage:
    make_pairwise_trait_plots.R [--help]
    make_pairwise_trait_plots.R <traits_infile> <plot_file_base>

Arguments:
    traits_infile         Input file of all traits (csv); output from make_traits_analysis.R
    plot_file_base      Path to figure outfile basename (excluding extension; will be appended to).
" -> doc

args <- list()
args$traits_infile <- "analysis/data/all_traits.txt"
args$plot_file_base <- "analysis/figs/pairwise_traits_biplots"

debug_status <- FALSE
arguments <- run_args_parse(args, debug_status, doc)

main(arguments)

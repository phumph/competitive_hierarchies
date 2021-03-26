#!/usr/bin/env Rscript

# ====================== #
# make_traits_analysis.R #
# ====================== #

# ----------- #
# Description #
# ----------- # ----------------------------------------------------------- #
# Script takes growth and competitive traits as input
# and generates (a) heatmap, (b) individual trait distributions, and
# (c) PCA + trait correlation plot for Psyr and Pfluo strains.
#
# Plot elements are manually assembled and edited as vector graphics
# prior to output as final versions.
#
# Table of summary statistics output as Table S1.
#
# Figures are not tracked in Makefile recipes.
# ------------------------------------------------------------------------ #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(gplots)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(psych)))

# To install ggbiplot: library(devtools); install_github("vqv/ggbiplot")
suppressWarnings(suppressMessages(library(ggbiplot)))
suppressWarnings(suppressMessages(library(ggrepel)))

source("analysis/scripts/comp_utils.R")

# ------------- #
# function defs #
# ------------- #

plot_heatmap <- function(traits_all, zero_to_na = FALSE) {

  if (zero_to_na == TRUE) {
    traits_all$c_t[traits_all$c_t == 0] <- NA
  }

  dat_mat <-
    traits_all %>%
    dplyr::select(-strain_id, -clade, -phylo_pos) %>%
    as.matrix()

  row.names(dat_mat) <- traits_all$strain_id
  the_cols <- colorRampPalette(c("midnightblue", "white", "darkorange2"))

  h1 <- gplots::heatmap.2(dat_mat,
                          dendrogram = c("none"),
                          trace = c("none"),
                          density.info = c("none"),
                          Rowv = FALSE,
                          Colv = FALSE,
                          na.color = "gray",
                          col = the_cols,
                          notecex = 0.6,
                          notecol = "black",
                          symbreaks = TRUE,
                          cexRow = 0.6,
                          cexCol = 0.6,
                          scale = "column",
                          colsep = c(1:9),
                          rowsep = c(1:40),
                          sepcolor = "white",
                          sepwidth = c(0.0025, 0.0025),
                          key.title = NA,
                          key.xlab = NA)

  return(h1)
}


plot_heatmap_clade <- function(traits_all, the_clade) {
  
  suppressWarnings(
    traits_all %>%
      dplyr::filter(clade == the_clade) %>%
      dplyr::select(-strain_id, -clade, -phylo_pos) %>%
      psych::corr.test(method = "pearson", adjust = "fdr") ->
      the_corrs
  )
  
  palette_breaks <- seq(-1, 1, 0.125)
  color_palette  <- colorRampPalette(
    c("midnightblue", "white", "darkorange2"))(length(palette_breaks) - 1
    )
  
  # capture only upper triangle
  the_corrs$r[lower.tri(the_corrs$r, diag = TRUE)] <- NA
  
  # generate labels: only cells with FDR p < 0.05 get text
  the_labels <- round(the_corrs$r, 2)
  the_labels[the_corrs$p >= 0.1] <- NA
  
  # produce heatmap
  h1 <- gplots::heatmap.2(the_corrs$r,
                  dendrogram = "none",
                  scale = "none",
                  trace = "none",
                  key = TRUE,
                  keysize = 1.5,
                  density.info = c("none"),
                  labCol = c(NA, colnames(the_corrs$r)[-1]),
                  labRow = c(row.names(the_corrs$r)[1:(length(row.names(the_corrs$r)) - 1)],
                             NA),
                  Rowv = FALSE,
                  Colv = FALSE,
                  col = color_palette,
                  breaks = palette_breaks,
                  symbreaks = TRUE,
                  cellnote = the_labels,
                  notecex = 0.9,
                  notecol = "black",
                  colsep,
                  rowsep,
                  sepcolor = "white",
                  sepwidth = c(0.1, 0.1)
  )
  return(h1)
}


plot_trait_distn <- function(dat, the_trait, the_lims, mult = NA) {
  
  if (!is.na(mult)) {
    dat[, the_trait] <- dat[, the_trait] * mult
  }
  
  tp <- ggplot(dat, aes_string(x = the_trait)) +
    geom_density(aes(linetype = clade)) +
    scale_linetype_manual(values = c(1, 2)) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()) +
    scale_x_continuous(limits = the_lims)
  
  return(tp)
}


do_pca <- function(traits_all) {

  traits_all %>%
    dplyr::filter(rowSums(is.na(.)) == 0) ->
    traits_all_complete
  
  traits_all_complete %>%
    dplyr::select(-strain_id, -clade, -phylo_pos) %>%
    as.matrix() %>%
    prcomp(center = TRUE, scale. = TRUE) ->
    pca_res

  pca_dat <- data.frame(pca_res$x,
                        clade = traits_all_complete$clade,
                        strain_id = traits_all_complete$strain_id)

  readr::write_csv(pca_dat, file.path(arguments$pca_outfile))

  ggbiplot::ggbiplot(pca_res,
                     obs.scale = 1,
                     var.scale = 1,
                     choices = 1:2,
                     alpha = 1,
                     size = 2,
                     groups = traits_all_complete$clade,
                     ellipse = TRUE,
                     ellipse.prob = 0.95,
                     circle = FALSE) +
    geom_point(aes(color = traits_all_complete$clade), size = 3) +
    theme_minimal() +
    theme(legend.position = c(1, 1),
          legend.justification = c(1, 1),
          panel.border = element_rect(fill = NA),
          panel.grid = element_blank()) +
    scale_color_manual(values = c("gray", "black")) +
    geom_text_repel(aes(label = pca_dat$strain_id), col = "gray40") +
    geom_hline(yintercept = 0, col = "gray40", lty = 3) +
    geom_vline(xintercept = 0, col = "gray40", lty = 3) +
    scale_x_continuous(limits = c(-5, 5)) +
    scale_y_continuous(limits = c(-5, 5)) ->
    pca_plot

  return(pca_plot)
}


produce_summary_stats <- function(traits_all) {
  
  traits_all %>%
    dplyr::mutate(r = r * 1000) %>%
    dplyr::select(-phylo_pos) %>%
    tidyr::gather(key = "trait", value = "trait_value", -strain_id, -clade) %>%
    dplyr::group_by(clade, trait) %>%
    dplyr::summarise(trait_mean = mean(trait_value, na.rm = T),
                     trait_sigma = sd(trait_value, na.rm = T),
                     trait_min = min(trait_value, na.rm = T),
                     trait_max = max(trait_value, na.rm = T)) %>%
    dplyr::mutate(mu = paste0(round(trait_mean, 2),
                              " (",
                              round(trait_sigma, 2), ")"),
                  range = paste0("[", round(trait_min, 2),
                                 "; ",
                                 round(trait_max, 2), "]")) %>%
    dplyr::select(clade, trait, mu, range) %>%
    tidyr::gather(key = "stat", value = "value", -clade, -trait) %>%
    tidyr::spread(key = "clade", value = "value") %>%
    dplyr::arrange(stat, desc(trait)) ->
    trait_table
  
  return(trait_table)
}


write_summary_stats_table <- function(summary_table) {
  
  rows_per_section <- nrow(summary_table) / 2
  the_footnote <- c("\\$r\\$ displayed as x 1000")
  
  summary_table %>%
    dplyr::select(-stat) %>%
    kableExtra::kable("latex",
                      booktabs = TRUE,
                      caption = "Summary statistics of life history and competitive trait distributions.") %>%
    kableExtra::kable_styling(
      position = "center",
      font_size = 8,
      bootstrap_options = c("striped", "condensed")
    ) %>%
    kableExtra::pack_rows("\\$\\mu (\\sigma)\\$",
                          1,
                          rows_per_section) %>%
    kableExtra::pack_rows("[min; max]",
                          rows_per_section + 1,
                          rows_per_section * 2) %>%
    kableExtra::add_footnote(the_footnote, notation = "alphabet") ->
    table_object
  
  file_conn <- file(file.path(arguments$table_s1), "w")
  cat(table_object, file = file_conn)
  close(file_conn)
}


main <- function(arguments) {

  comp      <- readr::read_csv(arguments$comp_traits, col_types = readr::cols())
  growth    <- readr::read_csv(arguments$growth_traits, col_types = readr::cols())
  meta_data <- readr::read_csv(arguments$strain_data, col_types = readr::cols())

  if ("clade" %in% names(growth)) {
    growth %>%
      dplyr::select(-clade) ->
      growth
  }

    if ("z4326" %in% growth$strain_id) {
    growth$strain_id[growth$strain_id == "z4326"] <- "4326"
  }

  comp$strain_id <- sapply(comp$strain_id, function(x) gsub("X", "", x))

  # combine data.frames and prepare for plot
  suppressWarnings(
    meta_data %>%
      dplyr::left_join(growth, by = "strain_id") %>%
      dplyr::left_join(comp, by = "strain_id") ->
      traits_all
  )

  traits_all %>%
    dplyr::select(strain_id, clade, phylo_pos,
                  r = r_max,
                  L = lambda,
                  K = K_max,
                  c_o, c_d, c_w, c_r, c_t, i_w) ->
    traits_all

  # save traits file as output for subsequent scripts
  readr::write_csv(traits_all, arguments$traits_outfile)

  # generate summary for table output
  suppressMessages(
    summary_table <- produce_summary_stats(traits_all)
  )
  write_summary_stats_table(summary_table)

  # plot and save heatmap
  pdf(file = "analysis/figs/heatmap_plot.pdf", width = 3.25, height = 6)
  plot_heatmap(traits_all)
  dev.off()

  # plot and save heatmap again to capture legend
  pdf(file = "analysis/figs/heatmap_plot_legend.pdf", width = 8, height = 6)
  plot_heatmap(traits_all)
  dev.off()

  # plot individual trait distributions
  growth_plotlist <- list(
    pdn_1 = plot_trait_distn(dat = traits_all,
                             the_trait = "r", the_lims = c(0, 15), mult = 10^4),
    pdn_2 = plot_trait_distn(dat = traits_all,
                             the_trait = "L", the_lims = c(0, 30), mult = 0.01),
    pdn_3 = plot_trait_distn(dat = traits_all,
                             the_trait = "K", the_lims = c(0, 1)),
    pdn_4 = plot_trait_distn(dat = traits_all,
                             the_trait = "c_o", the_lims = c(0, 1)),
    pdn_5 = plot_trait_distn(dat = traits_all,
                             the_trait = "c_d", the_lims = c(0, 1)),
    pdn_6 = plot_trait_distn(dat = traits_all,
                             the_trait = "c_w", the_lims = c(-1, 1)),
    pdn_7 = plot_trait_distn(dat = traits_all,
                             the_trait = "c_r", the_lims = c(0, 1)),
    pdn_8 = plot_trait_distn(dat = traits_all,
                             the_trait = "c_t", the_lims = c(0, 1)),
    pdn_9 = plot_trait_distn(dat = traits_all,
                             the_trait = "i_w", the_lims = c(-1, 1))
  )

  suppressWarnings(
    gps <- ggpubr::ggarrange(plotlist = growth_plotlist,
                             nrow = length(growth_plotlist),
                             align = "hv",
                             common.legend = TRUE)
  )

  ggsave(gps, file = file.path("analysis/figs/trait_distributions.pdf"), device = "pdf",
         width = 3, height = 8)

  # now generate PCA and produce output for panel (c)
  pdf(file = "analysis/figs/corr_plot_Psyr.pdf", width = 6, height = 6)
  plot_heatmap_clade(traits_all, the_clade = "Psyr")
  dev.off()

  pdf(file = "analysis/figs/corr_plot_Pfluo.pdf", width = 6, height = 6)
  plot_heatmap_clade(traits_all, the_clade = "Pflu")
  dev.off()

  # now for PCA
  ggsave(do_pca(traits_all),
         file = "analysis/figs/pca_plot.pdf",
         device = "pdf",
         width = 6, height = 6
  )
}

# ============ #
# main routine #
# ============ #

"make_traits_analysis.R

Usage:
    make_traits_analysis.R [--help]
    make_traits_analysis.R <comp_traits> <growth_traits> <strain_data> <traits_outfile> <pca_outfile> <table_s1>

Arguments:
    comp_traits         Input file of competitive traits (csv)
    growth_traits       Input file of growth traits (csv)
    strain_data         Strain meta-data file (csv)
    traits_outfile      Full path to full traits file
    pca_outfile         Full path to pca results file
    table_s1            Full path to Table S1 results (tex)
" -> doc

args <- list()
args$comp_traits    <- "analysis/data/comp_traits.csv"
args$growth_traits  <- "analysis/data/growth_traits_fitted.csv"
args$strain_data    <- "analysis/data/strain_metadata.csv"
args$traits_outfile <- "analysis/data/all_traits.txt"
args$pca_outfile    <- "analysis/data/pca_traits.txt"
args$table_s1       <- "analysis/tables/trait_summary_stats.tex"

debug_status <- FALSE
arguments <- run_args_parse(args, debug_status, doc)

main(arguments)

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
# ------------------------------------------------------------------------ #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(gplots)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(psych))) # for corr.test()
suppressWarnings(suppressMessages(library(ggbiplot)))
## To instal ggbiplot:
# library(devtools)
# install_github("vqv/ggbiplot")
suppressWarnings(suppressMessages(library(ggrepel)))

# ------------- #
# function defs #
# ------------- #

run_args_parse <- function(debug_status) {
  if (debug_status == TRUE) {
    arguments <- list()
    arguments$comp_traits    <- "analysis/data/comp_traits.csv"
    arguments$growth_traits  <- "analysis/data/growth_traits.csv"
    arguments$strain_data    <- "analysis/data/strain_metadata.csv"
    arguments$traits_outfile <- "analysis/data/all_traits.txt"
    arguments$pca_outfile    <- "analysis/data/pca_traits.txt"
  } else if (debug_status == FALSE) {
    args <- commandArgs(trailingOnly = FALSE)
    arguments <- list()
    arguments$comp_traits    <- args[1]
    arguments$growth_traits  <- args[2]
    arguments$strain_data    <- args[3]
    arguments$traits_outfile <- args[4]
    arguments$pca_outfile    <- args[5]
  }
  return(arguments)
}

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

  h1 <- heatmap.2(dat_mat,
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
  h1 <- heatmap.2(the_corrs$r,
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

  write.table(pca_dat,
              file = file.path(arguments$pca_outfile),
              row.names = F, col.names = T, sep = "\t")

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

  file_conn <- file(file.path("analysis/tables/trait_summary_stats.tex"), "w")
  cat(table_object, file = file_conn)
  close(file_conn)
}

main <- function(arguments) {

  cat("\n*** make_figure1.R ***\n\n")
  cat("Loading data...")

  comp      <- read.table(arguments$comp_traits, T, ",")
  growth    <- read.table(arguments$growth_traits, T, ",")
  meta_data <- read.table(arguments$strain_data, T, ",")

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
                  r = mu_mm,
                  L = lambda_mm,
                  K = A_mm,
                  c_o, c_d, c_w, c_r, c_t, i_w) ->
                  traits_all

  # save traits file as output for subsequent scripts
  write.table(traits_all,
              file = file.path(arguments$traits_outfile),
              sep = "\t",
              row.names = F,
              col.names = T)

  # generate summary for table output
  summary_table <- produce_summary_stats(traits_all)
  write_summary_stats_table(summary_table)

  cat("Done!\n")
  cat("Plotting heatmap...")

  # plot and save heatmap
  pdf(file = "analysis/figs/heatmap_plot.pdf", width = 3.25, height = 6)
    plot_heatmap(traits_all)
  dev.off()

  # plot and save heatmap again to capture legend
  pdf(file = "analysis/figs/heatmap_plot_legend.pdf", width = 8, height = 6)
    plot_heatmap(traits_all)
  dev.off()

  cat("Done!\n")
  cat("Generating trait plots...")
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
    nrow = length(growth_plotlist), align = "hv", common.legend = TRUE)
  )

  ggsave(gps, file = file.path("analysis/figs/trait_distributions.pdf"), device = "pdf",
         width = 3, height = 8)
  cat("Done!\n")
  cat("Generating pairwise trait correlations plot...")
  # now generate PCA and produce output for panel (c)
  pdf(file = "analysis/figs/corr_plot_Psyr.pdf", width = 6, height = 6)
    plot_heatmap_clade(traits_all, the_clade = "Psyr")
  dev.off()

  pdf(file = "analysis/figs/corr_plot_Pfluo.pdf", width = 6, height = 6)
    plot_heatmap_clade(traits_all, the_clade = "Pflu")
  dev.off()
  cat("Done!\n")
  cat("Performing principle components analysis...")
  # now for PCA
  pdf(file = "analysis/figs/pca_plot.pdf", width = 6, height = 6)
    do_pca(traits_all)
  dev.off()
  cat("Done!\n\n")
  cat("*** Script finished ***\n")
}

# ============ #
# main routine #
# ============ #

arguments <- run_args_parse(debug_status = TRUE)
main(arguments)

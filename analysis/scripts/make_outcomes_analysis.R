#!/usr/bin/env Rscript

# ======================== #
# make_outcomes_analysis.R #
# ======================== #

# ----------- #
# Description #
# ----------- # ----------------------------------------------------------- #
# Script takes input file of all traits
# and generates (a) Euclidean distances between all strains
# in order to predict interaction outcome between all pairs of strains
# as function of trait and/or genetic distance.
#
# Plot elements generated as output were manually assembled
# and edited as vector graphics prior to inclusion in manuscript.
# ------------------------------------------------------------------------ #

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(gplots)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(nnet)))

# ------------- #
# function defs #
# ------------- #

run_args_parse <- function(debug_status) {
  if (debug_status == TRUE) {
    arguments <- list()
    arguments$all_traits <- "data/all_traits.txt"
    arguments$pairs_file <- "data/interaction_pairs.txt"
    arguments$pca_file   <- "data/pca_traits.txt"
  } else if (debug_status == FALSE) {
    args <- commandArgs(trailingOnly = FALSE)
    arguments <- list()
    arguments$all_traits <- args[1]
    arguments$pairs_file <- args[2]
    arguments$pca_file   <- args[3]
  }
  return(arguments)
}

calc_trait_dists <- function(pairs, all_traits, focal_traits, as_vec = FALSE) {
  # take pairs, join traits for i and j; calc dist_fun
  # Euclidean distance
  dist_fun <- function(i, j) {
    sqrt(sum((i - j)^2, na.rm = T))
  }

  ij <- dplyr::select(pairs, i, j)

  do_dist_calc <- function(x, ij) {

    suppressWarnings(
    ij %>%
      dplyr::left_join(all_traits[, c("strain_id", x)],
                       by = c("i" = "strain_id")) %>%
      dplyr::left_join(all_traits[, c("strain_id", x)],
                       by = c("j" = "strain_id")) ->
      ij
    )

    ij$dist <- dist_fun(i = ij[, grep(".x", names(ij))],
                        j = ij[, grep(".y", names(ij))])
    ij %>%
      dplyr::select(dist) ->
      the_dist

    names(the_dist) <- paste0(paste0(x, collapse = "_"), "_dist")

    return(the_dist)
  }

  if (as_vec == FALSE) {
    dist_res  <- lapply(focal_traits, function(x) do_dist_calc(x, ij))
    all_dists <- cbind(pairs, do.call(cbind, dist_res))
  } else if (as_vec == TRUE) {
    dist_res <- c(NA)
    for (i in seq_len(nrow(ij))) {
      dist_res[i]  <- do_dist_calc(focal_traits, ij[i, ]) %>%
        do.call(rbind, .)
    }
    all_dists <- cbind(pairs, dist_res)
    names(all_dists)[names(all_dists) == "dist_res"] <-
      paste0(paste0(focal_traits, collapse = "_"), "_dist")
  }
  return(all_dists)
}

run_glms <- function(pairs_full) {

  # goal is to predict interaction type as function of genetic distance
  glm1 <- glm(RNI ~ PC1_dist + PC2_dist + PC3_dist,
              data = pairs_full,
              family = "binomial")
  broom::tidy(glm1)

  pairs_full %>%
    dplyr::filter(WB == "W", i_clade == "Psyr") %>%
    glm(ASYM ~ pdist + PC1_dist + PC2_dist + PC3_dist,
        family = "binomial", data = .) %>%
    broom::tidy() ->
    psyr_glm_coefs

  pairs_full %>%
    dplyr::filter(WB == "W", i_clade == "Pflu") %>%
    glm(ASYM ~ pdist + PC1_dist + PC2_dist + PC3_dist,
        family = "binomial", data = .) %>%
    broom::tidy() ->
    pflu_glm_coefs

  pairs_full %>%
    #dplyr::filter(WB == "W", i_clade == "Pflu") %>%
    glm(ASYM ~ i_clade : j_clade, family = "binomial", data = .) %>%
    broom::tidy() ->
    pflu_glm_coefs
    pairs_full %>%
      dplyr::filter(WB == "W", i_clade == "Psyr") %>%
      glm(ASYM ~ pdist + r_dist + K_dist + L_dist,
          family = "binomial", data = .) %>%
      broom::tidy() ->
      pflu_glm_coefs
}

plot_glm_res <- function(pairs_full) {

  ggplot(pairs_full, aes(x = PC1_dist, y = ASYM)) +
    geom_jitter(position = position_jitter(height = 0.02), alpha = 0.5) +
    stat_smooth(method = "glm", se = T) +
    theme_bw()

  gx2 <- ggplot(int_glm_s, aes(x = pdist, y = RI)) +
    geom_jitter(position = position_jitter(height = .025), alpha=0.5) +
    scale_y_continuous(limits = c(-0.025, 1.025)) +
    #geom_point() +
    stat_smooth( aes(y = RI),  method="glm", family="binomial", se=T) +
    theme_bw()

  gx3 <- ggplot(int_glm_s, aes(x = pdist, y = RNI)) +
    geom_jitter(position = position_jitter(height = .025), alpha=0.5) +
    scale_y_continuous(limits = c(-0.025, 1.025)) +
    #geom_point() +
    stat_smooth( aes(y = RNI), method="glm", family = "binomial", se=T) +
    theme_bw()

  ggpubr::ggarrange(plotlist = list(gx1, gx2, gx3), nrow = 3)

}

compare_mv_disp <- function(x,
                            dist_col,
                            clade_col = "i_clade",
                            write_out = TRUE) {

  # compute Welch's t-t.test
  x %>%
    dplyr::filter(WB == "W",
                  !!as.symbol(dist_col) > 0,
                  !!as.symbol(clade_col) == "Psyr") %>%
    dplyr::select(!!as.symbol(dist_col)) ->
    psyr_vec

  x %>%
    dplyr::filter(WB == "W",
                  !!as.symbol(dist_col) > 0,
                  !!as.symbol(clade_col) == "Pflu") %>%
    dplyr::select(!!as.symbol(dist_col)) ->
    pflu_vec

  t.test(psyr_vec, pflu_vec, var.equal = FALSE) %>%
    broom::tidy() -> t_res

  if (write_out == TRUE) {
    # export as markdown table
    kableExtra::kable(t_res) %>%
      kableExtra::kable_styling(
          bootstrap_options = c("striped", "hover", "condensed")
          ) ->
      table_object

    file_conn <- file(file.path("tables/mv_dist_res.md"), "w")
    cat("#### Multivariate pairwise trait distances by clade.\n",
        file = file_conn)
    cat("Results form Welch's two-sample unequal variance *t*-test:\n",
        file = file_conn, append = TRUE)
    cat(table_object, file = file_conn, append = TRUE)
    cat("\n<br>\n<center><img src='../figs/mv_trait_dists.png' width = 200></center>",
        file = file_conn, append = TRUE)
    close(file_conn)
  }

  # generate plot
  x %>%
    dplyr::filter(WB == "W", !!as.symbol(dist_col) > 0) %>%
    ggplot(aes_string(x = clade_col, y = dist_col)) +
    geom_jitter(width = 0.1, col = "gray40") +
    geom_boxplot(alpha = 0.5, col = "gray40", width = 0.25) +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA),
          panel.grid = element_blank()) +
    xlab("clade") +
    annotate("text", label = paste0("t = ",
                                    round(t_res$statistic, 2), "\np = ",
                                    round(t_res$p.value)), x = 1, y = 7) ->
    dist_plot

    ggplot2::ggsave(dist_plot, file = "figs/mv_trait_dists.png",
      device = "png",
      dpi = 300,
      width = 2.5,
      height = 4)

    # return objects ?

}

run_mn_outcomes <- function(pairs_full, outcome_col, predictor_col) {

  # produce table of outcomes
  input_dat <-
    pairs_full %>%
    dplyr::select(predictor_col = !!as.symbol(predictor_col),
                  outcome_col = !!as.symbol(outcome_col))

  count_table <-
    input_dat %>%
      dplyr::group_by(predictor_col,
                      outcome_col) %>%
      dplyr::summarise(outcome_count = n())

  # produce spread table object for output
  count_spread <-
    count_table %>%
      tidyr::spread(key = outcome_col, value = "outcome_count")

  freqs_spread <-
    count_spread %>%
    dplyr::mutate(tot = sum(ASYM, RNI, RI),
                  ASYM = round(ASYM / tot, 3),
                  RI = round(RI / tot, 3),
                  RNI = round(RNI / tot, 3)) %>%
    dplyr::select(-tot)

  # now run multinomial model
  # first set reference factor level:
  input_dat$outcome_col <- input_dat$outcome_col

  # run multinomial model
  mn <- nnet::multinom(outcome_col ~ predictor_col, data = input_dat)

  # generate fitted probabilities
  df_fitted <-
    data.frame(predictor_col = input_dat$predictor_col,
               round(fitted(mn), 3)) %>%
      dplyr::arrange(predictor_col) %>%
      unique()

  # generate p-values from 2-tailed z-test for all coefficients
  z <- summary(mn)$coefficients / summary(mn)$standard.errors
  p <- (1 - pnorm(abs(z), 0, 1)) * 2

  # write output for supplemental table
  coef_table <- broom::tidy(mn)
  names(coef_table) <- c("outcome", "term", "coefficient", "std.err", "z", "p")
  coef_table$term[grep("Intercept", coef_table$term)] <- "intercept"
  coef_table$term <- sapply(coef_table$term,
                            function(x) gsub("predictor_colW_", "", x))
  coef_table <-
    coef_table %>%
    dplyr::mutate(coefficient = round(coefficient, 3),
                  std.err = round(std.err, 3),
                  z = round(z, 3),
                  p = round(p, 4))

  # return observed counts, proportions, fitted, and coefs
  return(list(int_counts_obs = count_spread,
              int_freqs_fit  = df_fitted,
              model_coefs    = coef_table))
}

output_mn_res <- function(mn_res) {

  # barplot figure (Fig. 3a)
  mn_res$int_freqs_fit %>%
    tidyr::gather(key = "int", value = "int_freq", -predictor_col) %>%
    ggplot(aes(x = predictor_col, y = int_freq, fill = int)) +
    geom_bar(position = "stack", stat = "identity", col = "black") +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_manual(values = c("white", "gray40", "black"), name = "") +
    xlab("interaction type") +
    ylab("outcome frequency") ->
    barplot1

  ggsave(barplot1, file = "figs/interaction_barplot.pdf",
         device = "pdf",
         width = 3,
         height = 3
  )

  # tables (S1, S2)
  the_footnotes <- c("RNI = Reciprocal non-invasion",
                     "RI = Reciprocal invasion",
                     "ASYM = Asymmetric dominance",
                     "B = between-clade comparison",
                     "W = within-clade comparison")
  suppressWarnings(
      dplyr::bind_rows(mn_res$int_counts_obs, mn_res$int_freqs_fit) %>%
      dplyr::mutate()
    ) %>%
    knitr::kable() %>%
    kableExtra::kable_styling(bootstrap_options = c("condensed", "striped"),
                              full_width = FALSE,
                              position = "center") %>%
    kableExtra::pack_rows("Counts", 1, 3) %>%
    kableExtra::pack_rows("Frequencies", 4, 6) %>%
    kableExtra::add_footnote(the_footnotes, notation = "alphabet") ->
    outcome_dist_table

  file_conn <- file(file.path("tables/mn_outcomes_res.md"), "w")
  cat("#### Distribution of outcomes by pairing type.\n",
      file = file_conn)
  cat(outcome_dist_table, file = file_conn, append = TRUE)
  close(file_conn)

  # now output coefficients table
  the_footnotes_2 <- c("RNI = Reciprocal non-invasion",
                       "RI = Reciprocal invasion",
                       "coefficients = odds")
  mn_res$model_coefs %>%
    knitr::kable() %>%
    kableExtra::kable_styling(bootstrap_options = c("condensed", "striped"),
                              full_width = FALSE,
                              position = "center") %>%
    kableExtra::add_footnote(the_footnotes_2, notation = "alphabet") ->
    model_coef_output

  file_conn <- file(file.path("tables/mn_coef_res.md"), "w")
  cat("#### Multinomial model coefficients.\n",
      file = file_conn)
  cat(model_coef_output, file = file_conn, append = TRUE)
  close(file_conn)
}

# -------- #
# main def #
# -------- #

main <- function(arguments) {

  # load traits file
  all_traits <- read.table(arguments$all_traits,
                           header = T, sep = "\t", stringsAsFactors = F)

  # load pairs file
  pairs <- read.table(arguments$pairs_file,
                      header = T, sep = "\t", stringsAsFactors = F)
  # load pca results
  pca_res <- read.table(arguments$pca_file,
                        header = T, sep = "\t", stringsAsFactors = F)

  # add PCs to all_traits
  suppressWarnings(
    all_traits %>%
    dplyr::left_join(
      dplyr::select(pca_res, strain_id, PC1, PC2, PC3),
      by = "strain_id") ->
    all_traits
  )

  # calculate distances
  focal_traits <- c("r", "L", "K",
                    "c_o", "c_d", "c_w",
                    "c_r", "c_t", "i_w",
                    "PC1", "PC2", "PC3")

  pairs_full <- calc_trait_dists(pairs, all_traits, focal_traits)

  pairs_multivar <- calc_trait_dists(pairs, all_traits,
                                     focal_traits = c("PC1", "PC2", "PC3"),
                                     as_vec = TRUE)

  # generate multivariate dispersion comparison
  pairs_multivar %>%
    compare_mv_disp(dist_col = grep("_dist", names(pairs_multivar), value = T))

  # model interaction outcome distribution within versus between clades
  dist_col <- grep("_dist", names(pairs_multivar), value = T)

  pairs_full %>%
    dplyr::left_join(
      dplyr::select(pairs_multivar, i, j, !!as.symbol(dist_col)),
    by = c("i", "j")) ->
    pairs_full

  pairs_full$outcome_mn <- pairs_full$outcome
  pairs_full$outcome_mn[!pairs_full$outcome_mn %in% c("RI", "RNI")] <- "ASYM"
  pairs_full$pair_type <- paste0(pairs_full$WB, "_", pairs_full$i_clade)
  pairs_full$pair_type[grep("B", pairs_full$pair_type)] <- "B"

  # generate multinomial model on outcome proportion as function of pair type
  mn_res <-
    pairs_full %>%
      run_mn_outcomes(outcome_col = "outcome_mn", predictor_col = "pair_type")

  # outputs figures and tables from multinomial model
  output_mn_results(mn_res)
}

arguments <- run_args_parse(debug_status = TRUE)
main(arguments)
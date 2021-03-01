#!/usr/bin/env Rscript

# =================== #
# fit_growth_traits.R #
# =================== #

# ----------- #
# Description #
# ----------- # ----------------------------------------------------------- #
# Script fits growth trajectory data to find
# r_max, lag duration, and max yield.
# Max yield is taken as the maximum OD (minus background),
# rather than the end-point OD.

# Output is .csv with all growth traits fit for each strain.
# ------------------------------------------------------------------------- #

# ------ #
# header #
# ------ #

library(dplyr)
library(ggplot2)

base_dir <- "analysis/figs/growth_curves"

# ------------- #
# function defs #
# ------------- #

fit_spline <- function(time, y, label, lambda_fit = 1e-5, tsplit = 2000) {

    # step 1 ============================================
    # find r_max, lambda from early slice of trajectory =
    # ===================================================

    # remove zeros by adding pseudocounts.
    if (min(y) == 0) {
        y <- y + (sort(y)[sort(y) > 0][1] / 5)
    }

    # define initial interval
    tr <-  which(time == tsplit)

    # initial spline fit
    fit_spline <- smooth.spline(time[1:tr], y[1:tr], lambda = lambda_fit)

    # find higher resolution fit by decreasing x interval
    diff_interval <- 1
    xnew <- seq(min(time[1:tr]), max(time[1:tr]), diff_interval)
    dydt_fine <- predict(fit_spline, x = xnew, deriv = 1)
    max_index <- which.max(dydt_fine$y)
    t_max <- dydt_fine$x[max_index]

    optmax <- function(x) predict(fit_spline, x = x, deriv = 1)$y
    interval <- t_max + c(-2, 2) * diff_interval
    t_max2 <- optimize(f = optmax, interval = interval, maximum = TRUE)$maximum

    r_max <- predict(fit_spline, x = t_max2, deriv = 1)$y
    y_max <- predict(fit_spline, x = t_max2)$y

    # find y intercept
    b_max <- y_max - (t_max2 * r_max)

    # find x intercept
    lambda <- -b_max / r_max

    # step 2 ==================================
    # find K_max across the entire trajectory =
    # =========================================

    fit_for_k <- smooth.spline(time, y, lambda = lambda_fit)
    x_new_full <- seq(min(time), max(time), diff_interval)
    k_pred <- predict(fit_for_k, x = x_new_full)
    k_max_index <- which.max(k_pred$y)
    k_max_t <- k_pred$x[k_max_index]
    k_max <- max(k_pred$y)

    # produce and save plot of fit
    outfile <-  file.path(base_dir, paste0(label, "_", "curve.png"))
    png(file = outfile, width = 4.25, height = 4, res = 300, units = "in")
    plot(time, y, pch = 19, col = "gray30", cex = 0.5,
         las = 1, cex.axis = 0.66, xlab = "time",
         ylab = expression(OD[600]), bty = "n")
    points(c(0, lambda), c(0, 0), type = "l", col = "dodgerblue")
    points(fit_for_k, type = "l", col = "white", lwd = 1)
    abline(a = b_max, b = r_max, col = "dodgerblue")
    abline(h = k_max, col = "dodgerblue")
    points(t_max2, y_max, pch = 19, col = "darkorange2")
    points(lambda, 0, pch = 19, col = "darkorange2")
    points(k_max_t, k_max, pch = 19, col = "darkorange2")
    text(x = (t_max2 - 200), y = y_max, labels = expression(r[m]), cex = 0.75)
    text(x = lambda + 200, y = 0, labels = expression(lambda), cex = 0.75)
    text(x = k_max_t, y = k_max - (k_max / 10),
         labels = expression(K),
         cex = 0.75)
    text(x = 250, y = k_max - (k_max / 10), labels = label, cex = 0.66)
    dev.off()

    return(data.frame(
        "r_max" = r_max,
        "lag" = lambda,
        "K_max" = k_max
    ))
}


run_args_parse <- function(debug_status) {
    if (debug_status == TRUE) {
        arguments <- list()
        arguments$infile_pf <- "analysis/data/MM_gcurvedata_Pflu.txt"
        arguments$infile_ps <- "analysis/data/MM_gcurvedata_Psyr.txt"
        arguments$outfile <- "analysis/data/growth_traits_fitted.csv"
        arguments$tsplits <- "analysis/data/tsplits.csv"
    } else if (debug_status == FALSE) {
        args <- commandArgs(trailingOnly = FALSE)
        arguments$infile_pf <- args[1]
        arguments$infile_ps <- args[2]
        arguments$outfile <- args[3]
        arguments$tsplits <- args[4]
    }
    return(arguments)
}

df_to_long <- function(df, t_col = "min") {
    strain_cols <-  names(df)[!grepl(t_col, names(df))]
    df %>%
        tidyr::pivot_longer(cols = all_of(strain_cols),
                            names_to = "strain",
                            values_to = "OD") ->
        df_long
    return(df_long)
}


main <- function(arguments) {
    # load growth trajectories
    pf_traj <- readr::read_delim(
        arguments$infile_pf,
        delim = "\t"
    )
    ps_traj <- readr::read_delim(
        arguments$infile_ps,
        delim = "\t"
    )

    # load time splits file
    tsplits <- readr::read_csv(arguments$tsplits)

    # transform format for estimation
    pf_traj %>%
        df_to_long() %>%
        dplyr::mutate(clade = "Pf") -> pf_traj_long
    ps_traj %>%
        df_to_long() %>%
        dplyr::mutate(clade = "Ps") -> ps_traj_long
    pf_traj_long %>%
        dplyr::bind_rows(ps_traj_long) %>%
        dplyr::left_join(tsplits, by = c("clade", "strain")) %>%
        dplyr::mutate(strain_id = paste0(clade, "_", strain)) ->
        traj_long

    traj_split <- split(traj_long, traj_long$strain_id)
    # perform growth estimates
    lapply(traj_split, function(x) {
        fit_spline(
            x[["min"]], x[["OD"]],
            label = x[["strain_id"]][1],
            tsplit = x[["tsplit"]][1]
        )
    }) %>%
        do.call(rbind, .) -> curve_fits

    # gather strain and clade names for output df
    lapply(rownames(curve_fits),
           function(x) {
               s <- strsplit(x, "_") %>%
                   unlist()
               return(c(clade = s[1], strain_id = s[2]))
           }
    ) %>%
        do.call(rbind, .) %>%
        data.frame() -> names_df
    curve_fits %>%
        dplyr::bind_cols(names_df) %>%
        dplyr::mutate(r_max = round(r_max, 6),
                      lag = round(lag, 2),
                      K_max = round(K_max, 4)) %>%
        dplyr::rename(lambda = lag) ->
        curve_fits_out

    # write out
    curve_fits_out %>%
        readr::write_csv(file = arguments$outfile)
}

# ==== #
# main #
# ==== #

debug_status <- TRUE
arguments <- run_args_parse(debug_status)
main(arguments)

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

BASE_DIR <- "analysis/figs/growth_curves"

# ------------- #
# function defs #
# ------------- #

fit_spline <- function (time, y, label, lambda_fit = 1e-5, tsplit = 2000) {
    
    # step 1 ============================================
    # find r_max, lambda from early slice of trajectory =
    # ===================================================
    
    # remove zeros by adding pseudocounts.
    if (min(y) == 0) {
        y <- y + (sort(y)[sort(y) > 0][1] / 5)
    }
    
    # define initial interval
    tr = which(time == tsplit)
    
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
    lambda = -b_max / r_max
    
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
    outfile = file.path(BASE_DIR, paste0(label, "_", "curve.png"))
    png(file = outfile, width = 4.25, height = 4, res = 300, units = "in")
    plot(time, y, pch = 19, col = "gray30", cex = 0.5,
         las = 1, cex.axis = 0.66, xlab="time", ylab=expression(OD[600]), bty = "n")
    points(c(0, lambda), c(0, 0), type = "l", col = "dodgerblue")
    points(fit_for_k, type = "l", col = "white", lwd = 1)
    abline(a = b_max, b = r_max, col = "dodgerblue")
    abline(h = k_max, col = "dodgerblue")
    points(t_max2, y_max, pch=19, col = "darkorange2")
    points(lambda, 0, pch=19, col = "darkorange2")
    points(k_max_t, k_max, pch=19, col = "darkorange2")
    text(x = (t_max2 - 200), y = y_max, labels = expression(r[m]), cex = 0.75)
    text(x = lambda + 200, y = 0, labels = expression(lambda), cex = 0.75)
    text(x = k_max_t, y = k_max - (k_max / 10), labels = expression(K), cex = 0.75)
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
    } else if (debug_status == FALSE) {
        args <- commandArgs(trailingOnly = FALSE)
        arguments$infile_pf <- args[1]
        arguments$infile_ps <- args[2]
        arguments$outfile <- args[3]
    }
    return(arguments)
}


df_to_long <- function(df, t_col = "min") {
    strain_cols = names(df)[!grepl(t_col, names(df))]
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
    pf_traj %>% df_to_long() -> pf_traj_long
    ps_traj %>% df_to_long() -> ps_traj_long
    pf_split <- split(pf_traj_long, pf_traj_long$strain)
    ps_split <- split(ps_traj_long, ps_traj_long$strain)
    
    readr::write_csv(path = "analysis/data/tsplits.csv",
                     rbind(data.frame(clade = "Pf", strain = unique(pf_traj_long$strain)),
                           data.frame(clade = "Ps", strain = unique(ps_traj_long$strain))))
    
    tsplits <- data.frame(strain = strains, tsplit = tsplit)

    lapply(pf_split, function(x) {
        fit_spline(x[["min"]], x[["OD"]], label = paste0("Pf_", x[["strain"]][1]), tsplit = 2000)
    }) %>%
    do.call(rbind, .) -> pf_curve_fits
    lapply(ps_split, function(x) {
        fit_spline(x[["min"]], x[["OD"]], label = paste0("Ps_", x[["strain"]][1]), tsplit = 2000)
    }) %>%
    do.call(rbind, .) -> ps_curve_fits

    ps_curve_fits %>%
        dplyr::bind_rows(pf_curve_fits) %>%
        dplyr::mutate(strain_id = rownames(.)) -> df_joined   

    # read in known values to compare
    growth_traits <- readr::read_csv(file = file.path("analysis/data/growth_traits.csv"))
    growth_traits %>%
    dplyr::left_join(df_joined) -> master_df
    
    master_df %>%
        ggplot(aes(x = A_mm, y = K_max)) + 
        geom_point() +
        geom_abline(slope = 1, intercept = 0) +
        theme_bw()
}

# ==== #
# main #
# ==== #

debug_status <- TRUE
arguments <- run_args_parse(debug_status)
main(arguments)



# TO-DO
# 1. Fix label in plot script to include clade
# 2. Re-run growth curves fit
# 3. Inspect each plot; determine tsplit; add in tsplit
# 4. Re-run; re-inspect curves
# 5. Compare to previous values; adjust as needed
fit_spline <- function (time, y, label, lambda_fit = 1e-5) {
    
    # step 1 ============================================
    # find r_max, lambda from early slice of trajectory =
    # ===================================================

    # remove zeros by adding pseudocounts.
    if (min(y) == 0) {
        y <- y + (sort(y)[sort(y) > 0][1] / 5)
    }
    
    # define initial interval
    tr <- 200
    
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
    png(file = outfile, width = 4, height = 3.5, res = 300, units = "in")
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
        text(x = 10, y = k_max - (k_max / 10), labels = label, cex = 0.9)
    dev.off()
    
    return(data.frame(
        "r_max" = r_max,
        "lag" = lambda,
        "K_max" = k_max
    ))
}

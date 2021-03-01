# script to compare prior growth trait estimates to current
# PTH 2021-FEB-28

library(dplyr)
library(ggplot2)
library(ggrepel)

# read in new estimates
growth_traits_fitted <- readr::read_csv(file = file.path("analysis/data/growth_traits_fitted.csv"))
growth_traits_fitted$strain_id[growth_traits_fitted$strain_id == "z4326"] <- "4326"

# read in previous values to compare
growth_traits_old <- readr::read_csv(file = file.path("analysis/data/growth_traits_old.csv"))
growth_traits_old %>%
    dplyr::left_join(growth_traits_fitted, by = "strain_id") %>%
    dplyr::filter(!is.na(r_max)) -> master_df

master_df %>%
    ggplot(aes(x = A_mm, y = K_max, fill = clade)) + 
    geom_abline(slope = 1, intercept = 0, col = "gray40") +
    geom_point(pch = 21, col = "black") +
    theme_bw() +
    xlab("initial K estimate") +
    ylab("confirmed K estimate") +
    ggtitle("Confirmation of K estimates") +
    scale_fill_manual(values = c("darkorange2", "dodgerblue")) ->
    K_plot

ggsave(K_plot, file = "revision_scripts/k_confirm_plot.png",
       device = "png", width = 4, height = 3.5, units = "in")

master_df %>%
    ggplot(aes(x = lambda_mm, y = lambda, fill = clade, label = strain_id)) + 
    geom_abline(slope = 1, intercept = 0, col = "gray40") +
    geom_point(pch = 21, col = "black") +
    theme_bw() +
    xlab("initial lambda estiate") +
    ylab("confirmed lambda estimate") +
    ggtitle("Confirmation of lambda estimates") +
    scale_fill_manual(values = c("darkorange2", "dodgerblue")) +
    geom_text_repel(size = 2) ->
    lambda_plot

ggsave(lambda_plot, file = "revision_scripts/lambda_confirm_plot.png",
       device = "png", width = 4.5, height = 3.5, units = "in")

master_df %>%
    ggplot(aes(x = mu_mm, y = r_max, fill = clade, label = strain_id)) + 
    geom_abline(slope = 1, intercept = 0, col = "gray40") +
    geom_point(pch = 21, col = "black") +
    theme_bw() +
    xlab("initial r estimate") +
    ylab("confirmed r estimate") +
    ggtitle("Confirmation of r estimates") +
    scale_fill_manual(values = c("darkorange2", "dodgerblue")) +
    geom_text_repel(size = 2) ->
    r_plot

ggsave(r_plot, file = "revision_scripts/r_confirm_plot.png",
       device = "png", width = 4.5, height = 3.5, units = "in")

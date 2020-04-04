# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## load packages
library("dplyr")
library("ggplot2")
library("scales")

## control parameters
ns <- c(50, 100, 200)                  # number of observations
Ds <- c(5, 10, 20)                     # number of compositional parts
dims <- expand.grid(n = ns, D = Ds)    # all combinations

## load results
file <- "simulations/results/results_main_%d_%d.RData"
results_list <- mapply(function(n, D) {
  file_nD <- sprintf(file, n, D)
  if (file.exists(file_nD)) {
    load(file_nD)
    results_nD
  }
}, n = dims$n, D = dims$D, SIMPLIFY = FALSE, USE.NAMES = FALSE)
results <- do.call(rbind, results_list)

## nice labels for plotting
results <- results %>%
  mutate(n = factor(sprintf("n = %d", n), levels = sprintf("n = %d", ns)),
         k = sprintf("k = %d", k))

## loop over number of compositions and plot results
for (current_D in Ds) {

  # prepare MSE of coefficients
  core_methods <- c("OLS", "MM", "IF-MI", "BF-MI")
  df_MSE <- results %>%
    filter(D == current_D, Method %in% core_methods) %>%
    mutate(Method = factor(Method, levels = core_methods)) %>%
    group_by(n, Zeta, k, Method) %>%
    summarize(MSE = mean(MSE), Runs = n())
  # select colors and linetypes
  colors <- hue_pal()(6)[c(1, 6, 5, 4)]
  linetypes <- c("solid", "solid", "solid", "solid")
  # extract contamination levels to ensure gridlines and tick marks there
  zetas <- unique(df_MSE$Zeta)
  # create plot
  p_MSE <- ggplot() +
    geom_line(aes(x = Zeta, y = MSE, color = Method, linetype = Method),
              data = df_MSE) +
    facet_grid(n ~ k) +
    scale_color_manual("", values = colors) +
    scale_linetype_manual("", values = linetypes) +
    scale_x_continuous(breaks = zetas) +
    labs(x = "Contamination level", y = "Mean squared error") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          legend.direction = "horizontal",
          panel.grid.minor.x = element_blank(),
          panel.spacing.x = unit(0.5, "line"),
          panel.spacing.y = unit(0.5, "line"),
          plot.margin = margin(0, 0.4, 0.4, 0.4, unit = "line"))
  # plot to file
  file <- "simulations/figures/figure_MSE_%d.pdf"
  pdf(file = sprintf(file, current_D), width = 6.5, height = 6.75)
  print(p_MSE)
  dev.off()

  # prepare prediction results
  all_methods <- c("OLS", "MM", "3S", "shS-bi", "IF-MI", "BF-MI")
  df_MSEP <- results %>%
    filter(D == current_D, Method %in% all_methods) %>%
    mutate(Method = recode(factor(Method, levels = all_methods),
                           "shS-bi" = "ShS")) %>%
    group_by(n, Zeta, k, Method) %>%
    summarize(MSEP = mean(MSEP), Runs = n())
  # determine axis limits without pure cellwise methods (can be unstable)
  df_core <- df_MSEP %>% filter(Method %in% core_methods)
  ylim <- range(df_core$MSEP, na.rm = TRUE)
  # select colors and linetypes
  colors <- hue_pal()(6)[c(1, 6, 3, 2, 5, 4)]
  linetypes <- c("solid", "solid", "solid", "solid", "solid", "solid")
  # create plot
  p_MSEP <- ggplot() +
    geom_line(aes(x = Zeta, y = MSEP, color = Method, linetype = Method),
              data = df_MSEP) +
    facet_grid(n ~ k) +
    scale_color_manual("", values = colors, drop = FALSE) +
    scale_linetype_manual("", values = linetypes, drop = FALSE) +
    guides(color = guide_legend(nrow = 1),
           linetype = guide_legend(nrow = 1)) +
    coord_cartesian(ylim = ylim) +
    scale_x_continuous(breaks = zetas) +
    labs(x = "Contamination level", y = "Mean squared error of prediction") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          legend.direction = "horizontal",
          panel.grid.minor.x = element_blank(),
          panel.spacing.x = unit(0.5, "line"),
          panel.spacing.y = unit(0.5, "line"),
          plot.margin = margin(0, 0.4, 0.4, 0.4, unit = "line"))
  # plot to file
  file <- "simulations/figures/figure_MSEP_%d.pdf"
  pdf(file = sprintf(file, current_D), width = 6.5, height = 6.75)
  print(p_MSEP)
  dev.off()

}

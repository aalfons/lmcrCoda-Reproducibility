# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## load packages
library("dplyr")
library("ggplot2")
library("scales")
library("tidyr")

## load results
file <- "simulations/results/results_cpu_time.RData"
load(file)

## average CPU times
methods <- c("OLS", "MM", "3S", "shS-bi", "BF-MI")
df_CPU <- results %>%
  filter(Method %in% methods) %>%
  mutate(Method = recode(factor(Method, levels = methods),
                         "shS-bi" = "ShS")) %>%
  group_by(n, D, Method) %>%
  summarize(CPU = mean(CPU, na.rm = TRUE))

## compute relative CPU time (compared with our method)
dims <- unique(df_CPU[, c("n", "D")])
df_list <- mapply(function(n, D) {
  keep <- df_CPU$n == n & df_CPU$D == D
  df <- as.data.frame(df_CPU[keep, ])
  cpu_BF <- df[df$Method == "BF-MI", "CPU"]
  df$relCPU <- cpu_BF / df$CPU
  df
}, n = dims$n, D = dims$D, SIMPLIFY = FALSE, USE.NAMES = FALSE)
df_CPU <- do.call(rbind, df_list)

## prepare data for plotting
ns <- unique(df_CPU$n)
Ds <- unique(df_CPU$D)
df_plot <- df_CPU %>%
  rename(seconds = CPU, relative = relCPU) %>%
  gather(key = "Type", value = "CPU", -(1:3), factor_key = TRUE) %>%
  mutate(n = factor(sprintf("n = %d", n), levels = sprintf("n = %d", ns)))

# select colors and linetypes
colors <- hue_pal()(6)[c(1, 6, 3, 2, 4)]
linetypes <- c("solid", "solid", "solid", "solid", "solid")

## create plot
p_CPU <- ggplot() +
  geom_line(aes(x = D, y = CPU, color = Method), data = df_plot) +
  facet_grid(Type ~ n, scales = "free_y") +
  scale_color_manual("", values = colors) +
  scale_linetype_manual("", values = linetypes) +
  guides(color = guide_legend(nrow = 1),
         linetype = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = Ds) +
  labs(y = "CPU time") +
  theme_bw() +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        panel.grid.minor.x = element_blank(),
        panel.spacing.x = unit(0.5, "line"),
        panel.spacing.y = unit(0.5, "line"),
        plot.margin = margin(0, 0.4, 0.4, 0.4, unit = "line"))
# plot to file
pdf(file = "simulations/figures/figure_cpu_time.pdf",
    width = 6.5, height = 4)
print(p_CPU)
dev.off()

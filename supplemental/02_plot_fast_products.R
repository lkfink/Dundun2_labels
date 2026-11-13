# FAST product configuration: Experiments 1 and 2
# Run after 01_run_fast.R.
# Lauren Fink & chatGPT-6 at request of reviewer
# This script generates Figure S2

library(this.path)
library(tidyverse)
library(ggrepel)
library(patchwork)

out_dir <- file.path(this.path::this.dir(), "fast_outputs")

read_coords <- function(prefix) {
  x <- read_csv(file.path(out_dir, paste0(prefix, "_product_coordinates.csv")),
                show_col_types = FALSE, name_repair = "minimal")
  names(x)[1] <- "stimulus_id"
  stopifnot(all(c("stimulus_id", "Dim 1", "Dim 2") %in% names(x)))
  x %>% mutate(stimulus_class = factor(str_extract(stimulus_id, "^[SM]"),
                                        levels = c("S", "M")))
}

# FAST exports eigenvalue, percentage explained, and cumulative percentage.
read_variance <- function(prefix) {
  eig <- read.csv(file.path(out_dir, paste0(prefix, "_eig.csv")),
                  row.names = 1, check.names = FALSE)
  if (nrow(eig) < 2 || ncol(eig) < 2) stop("Malformed eigenvalue file: ", prefix)
  as.numeric(eig[1:2, 2])
}

colours <- c(S = "#CC6677", M = "#88CCEE")
shapes <- c(S = 16, M = 17)
class_labels <- c(S = "Speech", M = "Music")

plot_products <- function(coords, title, variance, seed) {
  ggplot(coords, aes(x = `Dim 1`, y = `Dim 2`)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3,
               colour = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
               colour = "grey70") +
    geom_point(aes(colour = stimulus_class, shape = stimulus_class),
               size = 3, alpha = 0.95) +
    geom_text_repel(aes(label = stimulus_id), size = 3.4, seed = seed,
                    box.padding = 0.35, point.padding = 0.2,
                    min.segment.length = 0, max.overlaps = Inf,
                    segment.colour = "grey60", segment.size = 0.3,
                    show.legend = FALSE) +
    scale_colour_manual(values = colours, labels = class_labels, name = "") +
    scale_shape_manual(values = shapes, labels = class_labels, name = "") +
    scale_x_continuous(expand = expansion(mult = 0.16)) +
    scale_y_continuous(expand = expansion(mult = 0.16)) +
    coord_equal(clip = "off") +
    labs(title = title,
         x = sprintf("FAST dimension 1 (%.2f%%)", variance[1]),
         y = sprintf("FAST dimension 2 (%.2f%%)", variance[2])) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(colour = "black"),
          legend.position = "top", plot.margin = margin(12, 18, 12, 12))
}

v1 <- read_variance("fast_exp1_wl")
v2 <- read_variance("fast_exp2_nl")
# Published values (rounded to two decimals): Exp 1 = 53.83, 6.39;
# Exp 2 = 58.33, 15.08. Flag changes without silently using stale labels.
if (any(abs(c(v1, v2) - c(53.83, 6.39, 58.33, 15.08)) > 0.05)) {
  warning("FAST variance percentages differ from the original manuscript plots; check inputs and package versions.")
}

p1 <- plot_products(read_coords("fast_exp1_wl"), "Experiment 1", v1, 101)
p2 <- plot_products(read_coords("fast_exp2_nl"), "Experiment 2", v2, 202)
p <- (p1 + p2 + plot_layout(guides = "collect")) &
  theme(legend.position = "top", legend.justification = "right",
        legend.box.just = "right")
p <- p + plot_annotation(tag_levels = "A",
                         theme = theme(plot.tag = element_text(face = "bold", size = 14)))

for (ext in c("png", "pdf")) {
  ggsave(file.path(out_dir, paste0("fast_product_configurations_exp1_exp2.", ext)),
         p, width = 12, height = 6.5, dpi = 600, bg = "white")
}
message("FAST product figure saved: ", out_dir)

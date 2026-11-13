# Experiment 2: participant-generated labels projected onto FAST products
# Run after 01_run_fast.R.
# Lauren Fink & chatGPT-6 at request of reviewer
# This script generates Figure 5 in the main manuscript

library(tidyverse)
library(ggrepel)
library(this.path)

out_dir <- file.path(this.path::this.dir(), "fast_outputs")
n_top_terms <- 10

# Retain the two-stage cleaning and phrase-level coding choices.
clean_label_phrase <- function(x) {
  replacements <- c(afro = "african", drums = "drum", dunun = "dundun",
                    percurssion = "percussion", base = "bass", load = "loud",
                    rythmic = "rhythmic", rhytmic = "rhythmic",
                    rithmic = "rhythmic", arythmic = "arhythmic",
                    arrhythmic = "arhythmic", arrythmic = "arhythmic",
                    arhythmic = "arhythmic")
  map_chr(x, function(s) {
    if (is.na(s) || str_squish(s) == "") return(NA_character_)
    words <- s %>% str_to_lower() %>%
      str_replace_all("[^a-z0-9 _-]", " ") %>%
      str_replace_all("[_-]+", " ") %>% str_squish() %>%
      str_split("\\s+") %>% pluck(1)
    change <- words %in% names(replacements)
    words[change] <- unname(replacements[words[change]])
    s <- paste(words, collapse = "_")
    if (s %in% c("", "1", "2")) return(NA_character_)
    dplyr::recode(s, non_rhythmic = "arhythmic", nonrhythmic = "arhythmic",
           bass = "bass_bassy", bassy = "bass_bassy", .default = s)
  })
}

read_fast_table <- function(filename, text_columns = FALSE) {
  # Explicitly read participant label/group matrices as text. Some columns may
  # contain only digits, but are still categorical labels rather than numbers.
  types <- if (text_columns) cols(.default = col_character()) else NULL
  x <- read_csv(file.path(out_dir, filename), show_col_types = FALSE,
                name_repair = "minimal", col_types = types)
  names(x)[1] <- "stimulus_id"
  x %>% select(-any_of("stim_acoustic_id"))
}

coords <- read_fast_table("fast_exp2_nl_product_coordinates.csv") %>%
  mutate(stimulus_class = factor(str_extract(stimulus_id, "^[SM]"),
                                 levels = c("S", "M")))
labels <- read_fast_table("fast_input_exp2_qualified_labels.csv", TRUE)
groups <- read_fast_table("fast_input_exp2_nl_group_ids.csv", TRUE)
stopifnot(all(c("Dim 1", "Dim 2") %in% names(coords)))
stopifnot(identical(names(labels), names(groups)))
stopifnot(setequal(coords$stimulus_id, labels$stimulus_id),
          setequal(coords$stimulus_id, groups$stimulus_id),
          !anyDuplicated(coords$stimulus_id),
          !anyDuplicated(labels$stimulus_id),
          !anyDuplicated(groups$stimulus_id),
          !anyNA(coords$stimulus_class))

groups_long <- groups %>%
  pivot_longer(-stimulus_id, names_to = "participant_id", values_to = "group_id") %>%
  mutate(group_id = as.character(group_id))
labels_long <- labels %>%
  pivot_longer(-stimulus_id, names_to = "participant_id", values_to = "label_raw") %>%
  mutate(label_raw = str_to_lower(label_raw)) %>%
  separate_rows(label_raw, sep = ";") %>%
  mutate(label = clean_label_phrase(label_raw)) %>%
  filter(!is.na(label), label != "") %>%
  left_join(groups_long, by = c("stimulus_id", "participant_id")) %>%
  distinct(participant_id, group_id, stimulus_id, label, .keep_all = TRUE)
if (anyNA(labels_long$group_id)) stop("Missing participant-box assignments.")
write_csv(labels_long, file.path(out_dir, "fast_exp2_labels_cleaned_long.csv"))

# Each participant-defined box contributes at most one use of a phrase.
term_frequency <- labels_long %>%
  distinct(participant_id, group_id, label) %>%
  count(label, name = "n_uses")

# Centroids are weighted by the number of stimulus-phrase associations.
label_centroids <- labels_long %>%
  count(label, stimulus_id, name = "n_stim_label") %>%
  left_join(coords %>% select(stimulus_id, `Dim 1`, `Dim 2`),
            by = "stimulus_id") %>%
  group_by(label) %>%
  summarise(dim1 = weighted.mean(`Dim 1`, n_stim_label),
            dim2 = weighted.mean(`Dim 2`, n_stim_label), .groups = "drop") %>%
  left_join(term_frequency, by = "label") %>%
  arrange(desc(n_uses), label)
if (anyNA(label_centroids$dim1) || anyNA(label_centroids$dim2)) {
  stop("Some labels cannot be matched to FAST stimulus coordinates.")
}
write_csv(label_centroids, file.path(out_dir, "fast_exp2_label_term_centroids.csv"))

top_terms <- label_centroids %>%
  slice_head(n = n_top_terms) %>%
  mutate(label_display = if_else(label == "bass_bassy", "bass / bassy",
                                 str_replace_all(label, "_", " ")))
write_csv(top_terms, file.path(out_dir, "fast_exp2_label_top10.csv"))
print(top_terms %>% select(label_display, n_uses, dim1, dim2))

# Single ggrepel layer prevents overlap between stimulus and phrase labels.
plot_text <- bind_rows(
  coords %>% transmute(x = `Dim 1`, y = `Dim 2`,
                       text = stimulus_id, fontface = "plain"),
  top_terms %>% transmute(x = dim1, y = dim2,
                          text = label_display, fontface = "bold")
)

# Use FAST eigenvalues
eig <- read.csv(file.path(out_dir, "fast_exp2_nl_eig.csv"),
                row.names = 1, check.names = FALSE)
variance <- as.numeric(eig[1:2, 2])


p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3,
             colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
             colour = "grey70") +
  geom_point(data = coords,
             aes(x = `Dim 1`, y = `Dim 2`, colour = stimulus_class,
                 shape = stimulus_class), size = 2.8, alpha = 0.9) +
  geom_point(data = top_terms,
             aes(x = dim1, y = dim2, size = n_uses),
             shape = 21, fill = "white", colour = "black", stroke = 0.6) +
  geom_text_repel(data = plot_text,
                  aes(x = x, y = y, label = text, fontface = fontface),
                  size = 3.2, seed = 456, box.padding = 0.4,
                  point.padding = 0.2, min.segment.length = 0,
                  max.overlaps = Inf, segment.colour = "grey65",
                  segment.size = 0.3, show.legend = FALSE) +
  scale_colour_manual(values = c(S = "#CC6677", M = "#88CCEE"),
                      labels = c(S = "Speech", M = "Music"), name = "") +
  scale_shape_manual(values = c(S = 16, M = 17),
                     labels = c(S = "Speech", M = "Music"), name = "") +
  scale_size_continuous(name = "Participant-box uses", range = c(2.5, 5)) +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 1),
         size = guide_legend(order = 2)) +
  scale_x_continuous(expand = expansion(mult = 0.14)) +
  scale_y_continuous(expand = expansion(mult = 0.14)) +
  coord_equal(clip = "off") +
  labs(x = sprintf("FAST dimension 1 (%.2f%%)", variance[1]),
       y = sprintf("FAST dimension 2 (%.2f%%)", variance[2])) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", legend.box = "vertical",
        legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        plot.margin = margin(15, 30, 15, 15))

for (ext in c("png", "pdf")) {
  ggsave(file.path(out_dir, paste0("fast_exp2_label_term_centroids_top10.", ext)),
         p, width = 7.5, height = 6.5, dpi = 600, bg = "white")
}
message("FAST label figure saved: ", out_dir)

# Analysis of participant-provided labels and reasons in Exp. 2
# Lauren Fink + chatGPT 5.6 Sol & 6
# This script reads in the anonymized labels and consensus categories 
# and creates supplemental figure S3


# Anchor paths to the project folder, even when run from supplemental/.
library(tidyverse)
library(this.path)

# Paths relative to this script
script_dir <- this.path::this.dir()

input_file <- file.path(
  script_dir, "..", "data",
  "exp2_consensus_categories_anonymized.csv"
)

out_dir <- file.path(script_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


categories <- readr::read_csv(input_file, show_col_types = FALSE)
expected <- c("instrument", "timing", "conceptual", "timbre",
              "pitch", "quality", "intensity")

if (nrow(readr::problems(categories)) > 0) stop("CSV parsing problems.")
if (!all(c("participant_id", "label_category") %in% names(categories))) {
  stop("Expected participant_id and label_category columns.")
}
if (nrow(categories) != 180 || anyNA(categories$participant_id) ||
    anyDuplicated(categories$participant_id) ||
    anyNA(categories$label_category) ||
    any(!categories$label_category %in% expected)) {
  stop("Unexpected participant count, IDs, or consensus categories.")
}

summary <- categories %>%
  count(label_category, name = "n") %>%
  mutate(percent = 100 * n / sum(n),
         category_display = stringr::str_to_sentence(label_category),
         value_display = sprintf("%d (%.1f%%)", n, percent)) %>%
  arrange(desc(n))

readr::write_csv(summary, file.path(out_dir, "exp2_consensus_category_summary.csv"))

p <- ggplot(summary, aes(x = forcats::fct_reorder(category_display, n), y = n)) +
  geom_col(width = 0.72, fill = "grey70", colour = "black", linewidth = 0.45) +
  geom_text(aes(label = value_display), hjust = -0.12, size = 3.8) +
  coord_flip(clip = "off") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6),
                     expand = expansion(mult = c(0, 0.22))) +
  labs(x = "Consensus category", y = "Number of participants") +
  theme_classic(base_size = 12) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 10, r = 45, b = 10, l = 10))

print(p)
for (ext in c("png", "pdf")) {
  ggsave(file.path(out_dir,
                   paste0("exp2_consensus_sorting_rationale_categories.", ext)),
         plot = p, width = 7, height = 4.8, units = "in",
         bg = "white", dpi = 600)
}
message("Saved figure and summary in: ", out_dir)

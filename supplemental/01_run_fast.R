# Dùndún Experiments 1 and 2: FAST analysis
# Lauren Fink & chatGPT 6
# Run this before 02_plot_fast_products.R and 03_plot_fast_labels.R.

# Set up
library(SensoMineR)
library(this.path)

script_dir <- this.path::this.dir()
data_dir <- file.path(script_dir, "..", "data")
out_dir <- file.path(script_dir, "fast_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
# Helper functions
############################################################
read_input <- function(name) {
  file <- file.path(data_dir, name)
  if (!file.exists(file)) stop("Missing input file: ", file)
  read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
}

# Stimulus columns are S1, M1, etc.; ignore any metadata columns, including IDs.
get_stimuli <- function(dat) {
  keys <- grep("^[SM][0-9]+$", names(dat), value = TRUE)
  if (!length(keys) || anyDuplicated(keys)) {
    stop("No unique S#/M# stimulus columns found; check input headers.")
  }
  pos <- dat[, keys, drop = FALSE]
  if (!all(vapply(pos, is.numeric, logical(1))) || anyNA(pos) ||
      !all(unlist(pos, use.names = FALSE) %in% c(1, 2))) {
    stop("Stimulus columns must contain only numeric group codes 1 or 2, without NAs.")
  }
  pos
}

# FAST uses a products x participants table with categorical group membership.
make_fast_input <- function(pos) {
  x <- as.data.frame(t(as.matrix(pos)))
  names(x) <- sprintf("P%03d", seq_len(nrow(pos)))
  x[] <- lapply(x, as.integer)  
  x
}

# Text normalization 
clean_phrase <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub(",", "", x)
  for (pair in list(c("afro", "african"), c("drums", "drum"),
                    c("dunun", "dundun"), c("percurssion", "percussion"),
                    c("base", "bass"), c("rythmic", "rhythmic"),
                    c("melodies", "melodic"), c("load", "loud"))) {
    x <- gsub(pair[1], pair[2], x)
  }
  trimws(gsub("\\s+", " ", x))
}

stim_to_acoustic_name <- function(x) {
  ifelse(substr(x, 1, 1) == "M", paste0("D", substring(x, 2), "M"),
         paste0("D", substring(x, 2), "s"))
}

save_matrix <- function(x, filename) {
  write.csv(cbind(stim_acoustic_id = stim_to_acoustic_name(rownames(x)), x),
            file.path(out_dir, filename), row.names = TRUE)
}

save_fast_results <- function(result, prefix) {
  write.csv(result$eig, file.path(out_dir, paste0(prefix, "_eig.csv")))
  write.csv(result$ind$coord,
            file.path(out_dir, paste0(prefix, "_product_coordinates.csv")))
  write.csv(result$cooccur,
            file.path(out_dir, paste0(prefix, "_reordered_cooccur.csv")))
}
############################################################

# Read in data
wl_raw <- read_input("wl_positions.csv")
nl_raw <- read_input("nl_positions.csv")
labels <- read_input("nl_labels.csv")
wl_pos <- get_stimuli(wl_raw)
nl_pos <- get_stimuli(nl_raw)

if (!setequal(names(wl_pos), names(nl_pos))) {
  stop("The two experiments must contain the same stimulus columns.")
}
nl_pos <- nl_pos[, names(wl_pos), drop = FALSE]
if (nrow(nl_pos) != nrow(labels) ||
    !all(c("left_box", "right_box") %in% names(labels))) {
  stop("nl_labels.csv must match nl_positions.csv by row and contain left_box/right_box.")
}
# If IDs exist in both original files, verify row alignment before discarding them.
shared_ids <- intersect(c("Subject", "subject", "participant_id", "ParticipantID"),
                        intersect(names(nl_raw), names(labels)))
for (id in shared_ids) {
  if (!identical(as.character(nl_raw[[id]]), as.character(labels[[id]]))) {
    stop("The Experiment 2 files have different participant ordering: ", id)
  }
}

wl_fast <- make_fast_input(wl_pos)
nl_fast <- make_fast_input(nl_pos)
message("Input participant counts: Experiment 1 = ", nrow(wl_pos),
        "; Experiment 2 = ", nrow(nl_pos))
if (nrow(wl_pos) != 108 || nrow(nl_pos) != 180) {
  warning("Input sample sizes differ from the expected 108/180; verify the final analysis samples.")
}

# Remove participants who never separated the stimuli into two groups.
keep_informative <- function(x, experiment) {
  keep <- vapply(x, function(p) length(unique(p)) >= 2, logical(1))
  if (!all(keep)) message(experiment, ": excluding ", sum(!keep),
                          " non-informative participant(s): ",
                          paste(names(x)[!keep], collapse = ", "))
  if (!any(keep)) stop(experiment, ": no informative participants remain.")
  x[, keep, drop = FALSE]
}
wl_fast <- keep_informative(wl_fast, "Experiment 1")
nl_fast <- keep_informative(nl_fast, "Experiment 2")

# Each stimulus receives its participant's label for the left or right box.
# Keep only participants included in the Experiment 2 FAST analysis.
nl_qualified_all <- matrix(NA_character_, nrow = nrow(nl_pos),
                           ncol = ncol(nl_pos),
                           dimnames = list(sprintf("P%03d", seq_len(nrow(nl_pos))),
                                           names(nl_pos)))
left_labels <- clean_phrase(labels$left_box)
right_labels <- clean_phrase(labels$right_box)
for (i in seq_len(nrow(nl_pos))) {
  nl_qualified_all[i, ] <- ifelse(unlist(nl_pos[i, ], use.names = FALSE) == 1,
                                left_labels[i], right_labels[i])
}
nl_qualified <- as.data.frame(t(nl_qualified_all), stringsAsFactors = FALSE)
nl_qualified <- nl_qualified[, names(nl_fast), drop = FALSE]

# Save data
save_matrix(wl_fast, "fast_input_exp1_wl_speech_music.csv")
save_matrix(nl_fast, "fast_input_exp2_nl_group_ids.csv")
save_matrix(nl_qualified, "fast_input_exp2_qualified_labels.csv")

set.seed(1234)
res_fast_wl <- SensoMineR::fast(wl_fast, graph = FALSE, axes = c(1, 2),
                              ncp = 5, B = 200)
res_fast_nl <- SensoMineR::fast(nl_fast, graph = FALSE, axes = c(1, 2),
                              ncp = 5, B = 200)

save_fast_results(res_fast_wl, "fast_exp1_wl")
save_fast_results(res_fast_nl, "fast_exp2_nl")
saveRDS(list(wl_fast = wl_fast, nl_fast = nl_fast,
             nl_qualified = nl_qualified,
             res_fast_wl = res_fast_wl, res_fast_nl = res_fast_nl),
        file.path(out_dir, "dundun_fast_results.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_session_info.txt"))
message("FAST analyses complete: ", out_dir)

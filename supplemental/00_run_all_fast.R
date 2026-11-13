# Run all FAST analyses and figures in order from any working directory.
# Lauren Fink & chatGPT 6, at request of reviewer
# packages required for each analysis are in each individual script

script_dir <- this.path::this.dir()
for (script in c("01_run_fast.R", "02_plot_fast_products.R",
                 "03_plot_fast_labels.R")) {
  message("Running: ", script)
  source(file.path(script_dir, script), local = TRUE)
}
message("All FAST analyses and figures completed.")

# run_all.R — reproduce every result and figure in order.
# Usage (from the repository folder):  Rscript run_all.R        or in R: source("run_all.R")
# Steps that download data skip themselves when their output files already exist.
# First full run downloads ~10 GB (Census slices, KPMP h5ad files, GTEx) and takes many hours.
# To make only the figures from saved results: source("R/00_setup.R"); source("R/16_figures.R")

repo <- normalizePath(getwd())
run <- function(f) {
  message("\n=========== ", f, " ===========")
  source(file.path(repo, "R", f), echo = FALSE, max.deparse.length = Inf)
}
steps <- c("00_setup.R",
           "01_ts_cohort_download.R", "02_ts_scores_subtypes.R", "03_ts_heart_relabel.R",
           "04_ts_main_sex_analyses.R", "05_ts_exploratory.R",
           "06_hca_validation.R", "07_hca_oxidative_stress.R",
           "09_hca_genome_wide.R",          # before 08: kidney models use Hs, cam2 and auto from 09
           "08_kpmp_kidney.R",
           "10_gtex_heart_arteries.R", "11_independent_lv_replication.R",
           "12_p2_sexchrom_receptor_panel.R", "13_p2_xy_dosage.R", "14_p2_loss_of_y.R",
           "15_p2_hormone_receptors.R",
           "16_figures.R")
for (f in steps) run(f)
message("Done. Figures are in ", file.path(getwd(), "figures"))

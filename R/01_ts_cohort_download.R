# =============================================================================
# TABULA SAPIENS: COHORT AND DOWNLOAD
# Tabula Sapiens 2.0: cohort, endothelial cells, download by donor
# Needs: 00_setup.R
# Makes: ts_obs_all.rds, ts_ec_obs.rds, ec_by_donor/<donor>.rds
# =============================================================================

# ---- 1. Find Tabula Sapiens 2.0 (use only "All Cells" so each cell counts once) ----
ds <- as.data.frame(census$get("census_info")$get("datasets")$read()$concat())
ts <- ds |> filter(str_detect(collection_name, regex("Tabula Sapiens", ignore_case = TRUE)))
ts |> select(dataset_title, dataset_total_cell_count) |> arrange(desc(dataset_total_cell_count))
all_id <- ts$dataset_id[ts$dataset_title == "Tabula Sapiens - All Cells"]   # 1,136,218 cells

# ---- 2. Cell metadata ----
if (!file.exists("ts_obs_all.rds")) {
  obs <- census$get("census_data")$get("homo_sapiens")$obs$read(
    value_filter = sprintf("dataset_id == '%s'", all_id),
    column_names = c("soma_joinid", "donor_id", "sex", "development_stage",
                     "tissue", "tissue_general", "cell_type", "assay", "raw_sum")
  )$concat() |> as.data.frame() |> as_tibble() |>
    mutate(age = parse_age(development_stage))
  saveRDS(obs, "ts_obs_all.rds")
}
obs <- readRDS("ts_obs_all.rds") |> as_tibble()      # older saved copies are plain data frames

obs |> distinct(donor_id, sex, age) |> arrange(sex, age) |> print(n = Inf)   # 24 donors: 13 F, 11 M
obs |> count(assay)

# ---- 3. Endothelial cells (10x only) ----
ec <- obs |>
  filter(str_detect(cell_type, "endotheli"), str_detect(assay, "10x")) |>
  mutate(age_grp = if_else(age >= 50, "50+", "<50"))
saveRDS(ec, "ts_ec_obs.rds")
nrow(ec)                                                                        # 69,271

ec |> count(cell_type, sex) |> pivot_wider(names_from = sex, values_from = n, values_fill = 0) |>
  arrange(desc(female + male)) |> print(n = Inf)
ec |> count(donor_id, sex, age) |> arrange(sex, age) |> print(n = Inf)
ec |> count(donor_id, sex, tissue_general) |> filter(n >= 50) |> count(tissue_general, sex) |>
  pivot_wider(names_from = sex, values_from = n, values_fill = 0) |>
  arrange(desc(pmin(female, male))) |> print(n = Inf)

# ---- 4. Donors with >= 100 endothelial cells (18 donors) ----
keep_donors <- ec |> count(donor_id) |> filter(n >= 100) |> pull(donor_id) |> as.character()
length(keep_donors)
saveRDS(keep_donors, "ts_keep_donors.rds")

# ---- 5. Download raw counts, one donor at a time (skips donors already saved) ----
dir.create("ec_by_donor", showWarnings = FALSE)
fetch_donor <- function(d) {
  out <- file.path("ec_by_donor", paste0(d, ".rds"))
  if (file.exists(out)) return(invisible("already done"))
  seu <- get_seurat(census, organism = "Homo sapiens",
                    obs_coords = ec$soma_joinid[ec$donor_id == d],
                    obs_column_names = c("soma_joinid", "donor_id", "sex", "development_stage",
                                         "tissue", "tissue_general", "cell_type", "assay", "raw_sum"))
  saveRDS(seu, out); rm(seu); gc(verbose = FALSE)
  invisible(out)
}
for (d in keep_donors) {
  message(format(Sys.time(), "%H:%M"), "  ", d)
  res <- tryCatch(fetch_donor(d), error = function(e) paste("ERROR:", conditionMessage(e)))
  if (startsWith(as.character(res), "ERROR")) message("   ", res)
}
length(list.files("ec_by_donor"))                                               # should be 18

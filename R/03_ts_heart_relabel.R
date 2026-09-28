# =============================================================================
# TABULA SAPIENS: FIXING HEART LABELS WITH THE HEART CELL ATLAS
# Fix heart vessel labels in Tabula Sapiens using the Heart Cell Atlas
# Why: TS expert heart labels were inconsistent between donors and confounded with sex
#      (men's "capillary" vs a woman's "venous" for similar cells).
# Needs: 00_setup.R, 02 outputs
# Makes: heart_dataset_ids.rds, hca_heart_ec_obs.rds, hca_heart_reference.rds,
#        hca_heart_lodo_accuracy.rds, ts_heart_labels.rds, ts_heart_markers.rds
# Methods facts: HCA whole-cell reference, leave-one-donor-out accuracy
#   capillary 87.1, venous 80.4, arterial 75.9, lymphatic 98.6 %
# =============================================================================
cells <- readRDS("ec_cells_scored.rds")

# ---- 1. The TS heart problem ----
heart_lab <- cells |> filter(tissue_general == "heart", !is.na(label_subtype))
heart_lab |> count(donor_id, sex, label_subtype) |>
  pivot_wider(names_from = label_subtype, values_from = n, values_fill = 0)

# ---- 2. Heart Cell Atlas datasets ----
ds <- as.data.frame(census$get("census_info")$get("datasets")$read()$concat())
heart_ids <- ds |> filter(str_detect(dataset_title, "^Vascular .* Cells of the adult human heart") |
                            str_detect(dataset_title, "Heart Global")) |> select(dataset_id, dataset_title)
saveRDS(heart_ids, "heart_dataset_ids.rds")
lit_id <- heart_ids$dataset_id[str_detect(heart_ids$dataset_title, "^Vascular")]   # Litvinukova 2020

# Expert-typed whole-cell endothelial cells (7 donors: 4 F, 3 M)
if (!file.exists("hca_heart_ec_obs.rds")) {
  lit <- census$get("census_data")$get("homo_sapiens")$obs$read(
    value_filter = sprintf("dataset_id == '%s'", lit_id),
    column_names = c("soma_joinid", "donor_id", "sex", "development_stage", "cell_type",
                     "assay", "suspension_type", "tissue"))$concat() |> as.data.frame() |> as_tibble() |>
    filter(suspension_type == "cell", str_detect(cell_type, "endotheli"), cell_type != "endothelial cell") |>
    mutate(subtype = case_when(str_detect(cell_type, "lymph") ~ "lymphatic",
                               str_detect(cell_type, "arter") ~ "arterial",
                               str_detect(cell_type, "vein")  ~ "venous",
                               str_detect(cell_type, "capillar") ~ "capillary"))
  saveRDS(lit, "hca_heart_ec_obs.rds")
}
lit <- readRDS("hca_heart_ec_obs.rds")
lit |> count(donor_id, sex, subtype) |> pivot_wider(names_from = subtype, values_from = n, values_fill = 0)

# ---- 3. Pseudobulk reference (donor x subtype) ----
if (!file.exists("hca_heart_reference.rds")) {
  href_list <- list()
  for (d in unique(as.character(lit$donor_id))) {
    message(format(Sys.time(), "%H:%M"), "  ", d)
    seu <- get_seurat(census, organism = "Homo sapiens", obs_coords = lit$soma_joinid[lit$donor_id == d],
                      obs_column_names = "soma_joinid")
    counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts"))
    sub <- lit$subtype[match(seu$soma_joinid, lit$soma_joinid)]; rm(seu)
    for (s in unique(sub)) {
      idx <- which(sub == s)
      if (length(idx) >= 20) href_list[[paste(d, s, sep = "|")]] <- Matrix::rowSums(counts[, idx, drop = FALSE])
    }
    rm(counts); gc(verbose = FALSE)
  }
  href <- list(counts = do.call(cbind, href_list))
  href$info <- tibble(id = colnames(href$counts)) |> separate(id, c("donor", "subtype"), sep = "\\|", remove = FALSE)
  saveRDS(href, "hca_heart_reference.rds")
}
href <- readRDS("hca_heart_reference.rds")
href_log <- log2(sweep(href$counts, 2, colSums(href$counts), "/") * 1e6 + 1)

# ---- 4. Leave-one-donor-out accuracy inside the Heart Cell Atlas ----
if (!file.exists("hca_heart_lodo_accuracy.rds")) {
  set.seed(1)
  test_ids <- lit |> group_by(donor_id, subtype) |> slice_sample(n = 300) |> ungroup()
  acc_list <- list()
  for (d in unique(as.character(test_ids$donor_id))) {
    tt <- test_ids |> filter(donor_id == d)
    seu <- get_seurat(census, organism = "Homo sapiens", obs_coords = tt$soma_joinid, obs_column_names = "soma_joinid")
    counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts"))
    truth <- tt$subtype[match(seu$soma_joinid, tt$soma_joinid)]; rm(seu)
    use <- href$info$donor != d
    pred <- SingleR(test = lognorm(counts), ref = href_log[, use], labels = href$info$subtype[use])
    acc_list[[d]] <- tibble(donor = d, truth = truth, pred = pred$labels)
  }
  saveRDS(bind_rows(acc_list), "hca_heart_lodo_accuracy.rds")
}
acc <- readRDS("hca_heart_lodo_accuracy.rds")
acc |> group_by(truth) |> summarise(cells = n(), pct_correct = round(100 * mean(pred == truth), 1))

# ---- 5. Relabel all TS heart endothelial cells + endocardium/venule markers ----
heart_cells <- cells |> filter(tissue_general == "heart")
if (!file.exists("ts_heart_labels.rds")) {
  hres <- list(); hm <- list()
  for (d in unique(as.character(heart_cells$donor_id))) {
    seu <- readRDS(file.path("ec_by_donor", paste0(d, ".rds")))
    seu <- seu[, seu$soma_joinid %in% heart_cells$soma_joinid]
    counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts"))
    pred <- SingleR(test = lognorm(counts), ref = href_log, labels = href$info$subtype)
    hres[[d]] <- tibble(soma_joinid = seu$soma_joinid, heart_subtype = pred$labels)
    hm[[d]]   <- tibble(soma_joinid = seu$soma_joinid,
                        NPR3 = as.numeric(counts["NPR3", ] > 0), ACKR1 = as.numeric(counts["ACKR1", ] > 0),
                        CA4 = as.numeric(counts["CA4", ] > 0), RGCC = as.numeric(counts["RGCC", ] > 0))
    rm(seu, counts, pred); gc(verbose = FALSE)
  }
  saveRDS(bind_rows(hres), "ts_heart_labels.rds")
  saveRDS(bind_rows(hm), "ts_heart_markers.rds")
}
hres <- readRDS("ts_heart_labels.rds"); hm <- readRDS("ts_heart_markers.rds")

heart_cells |> select(soma_joinid, donor_id, sex, tissue) |> left_join(hres, by = "soma_joinid") |>
  count(donor_id, sex, tissue, heart_subtype) |>
  pivot_wider(names_from = heart_subtype, values_from = n, values_fill = 0) |> print(n = Inf)
# Findings: TS "heart" includes coronary artery samples (TSP14 coronary only);
#           NPR3+ ACKR1- "venous" cells = endocardium (handled in 04).

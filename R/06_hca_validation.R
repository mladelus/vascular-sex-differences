# =============================================================================
# VALIDATION: HEART CELL ATLAS NUCLEI
# Independent validation in Heart Cell Atlas nuclei
# Litvinukova et al. 2020 vascular nuclei: 14 donors (7 F, 7 M), expert vessel labels,
# working myocardium only; doublets removed with scDblFinder (endothelial + mural nuclei).
# Needs: 00_setup.R, 03 output heart_dataset_ids.rds
# Makes: hca_vascular_nuclei_obs.rds, hca_vascular_nuclei_typed.rds, hg_nuclei/<donor>.rds,
#        hg_ec_nuclei_meta.rds, hg_donor_region_summaries.rds, hg_validation_results.rds
# Estimates are women minus men (positive = higher in women).
# =============================================================================
heart_ids <- readRDS("heart_dataset_ids.rds")
lit_id  <- heart_ids$dataset_id[str_detect(heart_ids$dataset_title, "^Vascular")]
glob_id <- heart_ids$dataset_id[str_detect(heart_ids$dataset_title, "Heart Global")]

# ---- Heart Global cohort and why we did not classify it ourselves ----
# Heart Global independent donors are nuclei only and unlabelled; a nucleus SingleR reference
# (hca_heart_nucleus_reference.rds) classified only 68.1% of expert capillary nuclei correctly
# (leave-one-donor-out), so expert labels from the vascular dataset were used instead.
if (file.exists("hca_heart_nucleus_lodo_accuracy.rds")) {
  readRDS("hca_heart_nucleus_lodo_accuracy.rds") |> group_by(truth) |>
    summarise(nuclei = n(), pct_correct = round(100 * mean(pred == truth), 1)) |> print()
}

# ---- 1. Expert-labelled vascular nuclei ----
if (!file.exists("hca_vascular_nuclei_obs.rds")) {
  vobs <- census$get("census_data")$get("homo_sapiens")$obs$read(
    value_filter = sprintf("dataset_id == '%s'", lit_id),
    column_names = c("soma_joinid", "donor_id", "sex", "development_stage", "tissue",
                     "cell_type", "assay", "suspension_type", "raw_sum"))$concat() |>
    as.data.frame() |> as_tibble() |> filter(suspension_type == "nucleus")
  saveRDS(vobs, "hca_vascular_nuclei_obs.rds")
}
vobs <- readRDS("hca_vascular_nuclei_obs.rds")
vt <- vobs |> mutate(type = case_when(str_detect(cell_type, "capillar") ~ "capillary",
                                      str_detect(cell_type, "artery")   ~ "arterial",
                                      str_detect(cell_type, "vein")     ~ "venous",
                                      str_detect(cell_type, "pericyte|smooth") ~ "mural",
                                      TRUE ~ NA_character_)) |> filter(!is.na(type))
saveRDS(vt, "hca_vascular_nuclei_typed.rds")
vt |> count(donor_id, sex, development_stage, type) |>
  pivot_wider(names_from = type, values_from = n, values_fill = 0) |> arrange(sex, development_stage) |> print(n = Inf)

# ---- 2. Per donor: doublets (scDblFinder), scores, pseudobulk ----
set.seed(1)
sel <- bind_rows(vt |> filter(type != "mural"),
                 vt |> filter(type == "mural") |> group_by(donor_id) |> slice_sample(n = 2000) |> ungroup())
dir.create("hg_nuclei", showWarnings = FALSE)
for (d in unique(as.character(sel$donor_id))) {
  out <- file.path("hg_nuclei", paste0(d, ".rds")); if (file.exists(out)) next
  message(format(Sys.time(), "%H:%M"), "  ", d)
  seu <- get_seurat(census, organism = "Homo sapiens", obs_coords = sel$soma_joinid[sel$donor_id == d],
                    obs_column_names = "soma_joinid")
  counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts"))
  meta <- sel[match(seu$soma_joinid, sel$soma_joinid), ]; rm(seu); gc(verbose = FALSE)

  sce <- scDblFinder(SingleCellExperiment(list(counts = counts)), clusters = meta$type)
  meta$dbl_class <- sce$scDblFinder.class; meta$dbl_score <- sce$scDblFinder.score; rm(sce)

  ec <- meta$type != "mural"; cE <- counts[, ec, drop = FALSE]; mE <- meta[ec, ]
  sc <- ScoreSignatures_UCell(cE, features = hg_sets, ncores = 1)
  cp <- function(g) as.numeric(Matrix::colSums(cE[intersect(g, rownames(cE)), , drop = FALSE]))
  g1 <- function(g) if (g %in% rownames(cE)) as.numeric(cE[g, ]) else 0
  mE <- mE |> mutate(peri_umi = cp(pan$peri), sm_umi = cp(pan$sm), str_umi = cp(pan$str),
                     CDKN2A = g1("CDKN2A"), CDKN1A = g1("CDKN1A"), MKI67 = g1("MKI67")) |>
    bind_cols(as_tibble(sc))

  sing <- mE$dbl_class == "singlet"; grp <- paste(mE$donor_id, mE$tissue, mE$type, sep = "|"); pbl <- list()
  for (g in unique(grp[sing])) { idx <- which(sing & grp == g)
    if (length(idx) >= 20) pbl[[g]] <- Matrix::rowSums(cE[, idx, drop = FALSE]) }
  saveRDS(list(meta = mE, pb = pbl), out)
  rm(counts, cE, mE, sc, pbl); gc(verbose = FALSE)
}
hg  <- lapply(list.files("hg_nuclei", full.names = TRUE), readRDS)
hgm <- bind_rows(lapply(hg, `[[`, "meta")); saveRDS(hgm, "hg_ec_nuclei_meta.rds")

# Doublet rates by donor and vessel type
hgm |> count(donor_id, sex, type, dbl_class) |>
  pivot_wider(names_from = dbl_class, values_from = n, values_fill = 0) |>
  mutate(pct_doublet = round(100 * doublet / (doublet + singlet), 1)) |> print(n = Inf)

# ---- 3. Donor x region x chemistry x vessel-type summaries (doublet-free) ----
hgs <- hgm |> filter(dbl_class == "singlet") |>
  mutate(age = parse_age(development_stage), age50 = age >= 50,
         sex = relevel(factor(as.character(sex)), ref = "male"),
         p16 = CDKN2A > 0 & MKI67 == 0, cdkn1a_cp10k = CDKN1A / raw_sum * 1e4) |>
  group_by(type) |> mutate(p21 = cdkn1a_cp10k > quantile(cdkn1a_cp10k, 0.9) & CDKN1A > 0) |> ungroup()
dr <- hgs |> group_by(donor_id, sex, age, age50, tissue, assay, type) |>
  summarise(n = n(),
            pericyte = mean(log1p(peri_umi / raw_sum * 1e4)), pericyte_ge3 = mean(peri_umi >= 3),
            smooth_muscle = mean(log1p(sm_umi / raw_sum * 1e4)), striated = mean(log1p(str_umi / raw_sum * 1e4)),
            inflammatory = mean(Inflammatory_UCell), tnfa = mean(TNFA_UCell), senmayo = mean(SenMayo_UCell),
            oxphos = mean(OXPHOS_UCell), stress = mean(Stress_UCell), p16 = mean(p16), p21 = mean(p21),
            .groups = "drop") |> filter(n >= 20)
saveRDS(dr, "hg_donor_region_summaries.rds")

# ---- 4. Women vs men (all ages: 7 vs 7; age 50+: 5 vs 7) ----
fitsex <- function(y, cls, older_only = FALSE) {
  d <- dr |> filter(type %in% cls); if (older_only) d <- filter(d, age50)
  f <- as.formula(paste(y, "~ sex + age + tissue + assay + stress",
                        if (length(cls) > 1) "+ type" else "", "+ (1 | donor_id)"))
  co <- summary(lmerTest::lmer(f, data = d))$coefficients["sexfemale", ]
  tibble(outcome = y, cells = paste(cls, collapse = "+"), ages = if (older_only) "50+" else "all",
         women = n_distinct(d$donor_id[d$sex == "female"]), men = n_distinct(d$donor_id[d$sex == "male"]),
         estimate = signif(co[1], 3), conf.low = signif(co[1] - 1.96 * co[2], 3),
         conf.high = signif(co[1] + 1.96 * co[2], 3), p = signif(co["Pr(>|t|)"], 2))
}
outs <- c("pericyte", "pericyte_ge3", "smooth_muscle", "striated",
          "inflammatory", "tnfa", "senmayo", "p16", "p21", "oxphos")
res <- bind_rows(
  lapply(outs, fitsex, cls = "capillary"),
  lapply(outs, fitsex, cls = "capillary", older_only = TRUE),
  fitsex("oxphos", c("capillary", "arterial", "venous")),
  fitsex("oxphos", c("capillary", "arterial", "venous"), older_only = TRUE),
  fitsex("oxphos", "arterial"), fitsex("oxphos", "venous"))
saveRDS(res, "hg_validation_results.rds")
print(res, n = Inf, width = Inf)
# "boundary (singular) fit" messages = donor variance estimated ~0 after covariates; estimates remain valid.

# Per-donor view: pericyte transcripts in capillary nuclei
dr |> filter(type == "capillary") |> group_by(donor_id, sex, age) |>
  summarise(pericyte = round(weighted.mean(pericyte, n), 3),
            pct_pericyte_ge3 = round(100 * weighted.mean(pericyte_ge3, n), 1), nuclei = sum(n), .groups = "drop") |>
  arrange(sex, age)

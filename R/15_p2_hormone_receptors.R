# =============================================================================
# HORMONE RECEPTORS: AGE, KIDNEY ESR1, KPMP SINGLE-CELL REPLICATION (Figure 8)
# Result: AR highest in kidney EC; ESR1 and PGR moderate; GPER1 low; CYP19A1 ~absent.
# ESR2 overlaps SYNE2 (not interpreted). No receptor changed with age in women (FDR >= 0.81).
# Kidney nuclei: ESR1 higher in women's glomerular EC (+2.15, p 0.0015 after tubular adjustment).
# Single cells: ESR1 higher in women in all 4 comparisons (+0.61 to +1.19); p = 0.043 for EHD3- cells
# (all donors), 0.07-0.42 otherwise; not glomerular-specific.
# =============================================================================

# ---- 1. Mean receptor expression by cohort x vessel type ----
long |> filter(gene %in% genes_hr) |> group_by(cohort, type, gene) |>
  summarise(mean_log2cpm = round(mean(logcpm), 1), .groups = "drop") |>
  pivot_wider(names_from = gene, values_from = mean_log2cpm) |> print(width = Inf)

# ---- 2. Checks: ESR2 vs overlapping SYNE2; kidney receptors vs tubular ambient RNA ----
wide_l <- long |> select(cohort, type, id, gene, logcpm) |> pivot_wider(names_from = gene, values_from = logcpm)
cat("ESR2 vs SYNE2 r =", round(cor(wide_l$ESR2, wide_l$SYNE2, use = "complete.obs"), 2), "\n")
wk <- wide_l |> filter(cohort == "KPMP kidney") |>
  mutate(tub = rowMeans(across(any_of(tub))))
print(sapply(c("AR", "ESR1", "PGR", "GPER1"), function(g) round(cor(wk[[g]], wk$tub), 2)))

# ---- 3. Age: women 50+ vs < 50, and sex x age50 (all FDR >= 0.81) ----
res_age <- long |> filter(gene %in% genes_hr, gene != "SYNE2") |> group_by(cohort, type, gene) |>
  group_modify(function(d, k) {
    w <- filter(d, sex == "female"); kid <- k$cohort == "KPMP kidney"
    a <- if (n_distinct(w$age50) == 2) summary(lm(if (kid) logcpm ~ age50 + disease2 else logcpm ~ age50,
                                                     data = w))$coefficients["age50TRUE", ] else rep(NA, 4)
    i <- summary(lm(if (kid) logcpm ~ sex * age50 + disease2 else logcpm ~ sex * age50, data = d))$coefficients
    ii <- if ("sexfemale:age50TRUE" %in% rownames(i)) i["sexfemale:age50TRUE", ] else rep(NA, 4)
    tibble(women_50plus_vs_under50 = a[1], p_age = a[4], sex_x_age50 = ii[1], p_int = ii[4])
  }) |> ungroup() |> mutate(FDR_age = p.adjust(p_age, "BH"), FDR_int = p.adjust(p_int, "BH"))
saveRDS(res_age, "p2_receptors_by_age.rds")
res_age |> arrange(p_age) |> mutate(across(where(is.numeric), ~ signif(.x, 2))) |> print(n = 20, width = Inf)

# ---- 4. Kidney nuclei: ESR1 by sex, adjusted for tubular ambient RNA ----
esr1_nuc <- wk |> select(id, ESR1, tub) |>
  left_join(s3 |> select(id, donor_id, type, sex, age, disease2), by = "id") |>
  mutate(sex = relevel(factor(sex), ref = "male"))
res_esr1_nuc <- esr1_nuc |> group_by(type) |> group_modify(function(d, k) {
  co <- summary(lm(ESR1 ~ sex + age + disease2 + tub, data = d))$coefficients["sexfemale", ]
  tibble(data = "nuclei", donors = "all", estimate = co[1], se = co[2], p = co[4])
}) |> ungroup()
print(res_esr1_nuc)

# ---- 5. KPMP single-cell replication (download ~ once) ----
# The single-cell file labels only subclass.l1 ("EC"), so glomerular-like cells = EHD3 > 0.
ds <- as.data.frame(census$get("census_info")$get("datasets")$read()$concat())
sc_id <- ds$dataset_id[ds$dataset_title == "Single-cell RNA-seq of the Adult Human Kidney (Version 1.5)"]
nuc_donors <- unique(km$donor_id)

fsc <- "kpmp/kpmp_scRNA_v1.5.h5ad"
if (!file.exists(fsc)) {
  u_sc <- sub("^s3://([^/]+)/", "https://\\1.s3.us-west-2.amazonaws.com/", get_source_h5ad_uri(sc_id)$uri)
  options(timeout = 3600); download.file(u_sc, fsc, mode = "wb")
}
library(rhdf5); filter <- dplyr::filter; select <- dplyr::select
rd2 <- function(col) {                        # same reader as script 08, for the single-cell file
  x <- h5read(fsc, paste0("/obs/", col))
  if (is.list(x) && all(c("categories", "codes") %in% names(x))) {
    i <- as.integer(x$codes) + 1; i[i < 1] <- NA; return(as.vector(x$categories[i])) }
  cats <- tryCatch(h5read(fsc, paste0("/obs/__categories/", col)), error = function(e) NULL)
  if (!is.null(cats)) { i <- as.integer(x) + 1; i[i < 1] <- NA; return(as.vector(cats[i])) }
  as.vector(x)
}
# Age from Census development_stage: "45-year-old human stage", "60-69 year-old ..." (midpoint) or decade stages
parse_age_sc <- function(x) {
  x <- as.character(x)
  one <- suppressWarnings(as.numeric(str_match(x, "^(\\d+)-year")[, 2]))
  rng <- str_match(x, "^(\\d+)\\s*-\\s*(\\d+)\\s*year")
  mid <- (suppressWarnings(as.numeric(rng[, 2])) + suppressWarnings(as.numeric(rng[, 3]))) / 2
  coalesce(one, mid, parse_age(x))
}
need <- c("soma_joinid", "donor_id", "sex", "age", "disease", "raw_sum")
sc_ok <- function(x) is.data.frame(x) && all(need %in% names(x)) && nrow(x) > 0
if (file.exists("kpmp_sc_ec_meta.rds") && !sc_ok(readRDS("kpmp_sc_ec_meta.rds"))) {
  message("kpmp_sc_ec_meta.rds is old or empty -> renamed to kpmp_sc_ec_meta_OLD.rds and rebuilt")
  file.rename("kpmp_sc_ec_meta.rds", "kpmp_sc_ec_meta_OLD.rds")
}
if (!file.exists("kpmp_sc_ec_meta.rds")) {
  stopifnot("Single-cell KPMP dataset not found in the Census" = length(sc_id) == 1)
  kp_sc2 <- tibble(observation_joinid = rd2("observation_joinid"), subclass.l1 = rd2("subclass.l1"))
  kc_sc <- census$get("census_data")$get("homo_sapiens")$obs$read(
    value_filter = sprintf("dataset_id == '%s'", sc_id),
    column_names = c("soma_joinid", "observation_joinid", "donor_id", "sex", "development_stage",
                     "disease", "raw_sum"))$concat() |> as.data.frame() |> as_tibble() |>
    mutate(across(c(observation_joinid, donor_id, sex, development_stage, disease), as.character))
  j <- kc_sc |> inner_join(kp_sc2 |> mutate(observation_joinid = as.character(observation_joinid)),
                           by = "observation_joinid")
  cat("Cells in Census:", nrow(kc_sc), "| matched to h5ad labels:", nrow(j), "\n")
  print(count(j, subclass.l1, sort = TRUE))
  sc_ec <- j |> filter(subclass.l1 == "EC") |> mutate(age = parse_age_sc(development_stage)) |>
    filter(!is.na(age), age >= 18, disease != "acute kidney failure")
  if (nrow(sc_ec) == 0) {
    print(count(j |> filter(subclass.l1 == "EC"), development_stage, disease))
    stop("No single-cell EC left after filtering: check the two tables above (age and disease labels may have changed in the Census).")
  }
  saveRDS(sc_ec, "kpmp_sc_ec_meta.rds")
}
sc_ec <- readRDS("kpmp_sc_ec_meta.rds")
cat("Single-cell kidney EC:", nrow(sc_ec), "cells from", n_distinct(sc_ec$donor_id), "donors\n")
print(sc_ec |> distinct(donor_id, sex, disease) |> count(disease, sex))

g_sc <- c("EHD3", "ESR1", "ESR2", "PGR", "AR", "GPER1", "PLVAP", "LRP2", "CUBN", "SLC34A1", "UMOD",
          "NPHS1", "NPHS2", "PECAM1")
get_genes <- function(ids) {                  # selected genes for many cells, 20,000 at a time
  ids <- sort(as.numeric(ids))
  bind_rows(lapply(split(ids, ceiling(seq_along(ids) / 2e4)), function(ch) {
    seu <- get_seurat(census, organism = "Homo sapiens", obs_coords = ch, obs_column_names = "soma_joinid",
                      var_index = "feature_name",
                      var_value_filter = sprintf("feature_name %%in%% c(%s)", paste0("'", g_sc, "'", collapse = ",")))
    cnt <- as.matrix(LayerData(seu, assay = "RNA", layer = "counts"))
    out <- as_tibble(t(cnt)) |> mutate(soma_joinid = seu$soma_joinid); rm(seu, cnt); gc(verbose = FALSE); out
  }))
}

# (a) Does EHD3 mark glomerular EC? Check on labelled kidney NUCLEI (EC-GC vs EC-PTC vs EC-AEA)
if (!file.exists("kpmp_nuc_ehd3.rds")) saveRDS(get_genes(km$soma_joinid), "kpmp_nuc_ehd3.rds")
km |> select(soma_joinid, subclass.l2) |> inner_join(readRDS("kpmp_nuc_ehd3.rds"), by = "soma_joinid") |>
  group_by(subclass.l2) |> summarise(nuclei = n(), pct_EHD3_pos = round(100 * mean(EHD3 > 0), 1)) |> print()
# Result: EHD3+ in 52.8% of EC-GC, 3.1% of EC-PTC, 2.1% of EC-AEA nuclei (in Methods)

# (b) Single cells: genes per cell, glomerular-like = EHD3 > 0; pseudobulk per donor x group
gfile <- "kpmp_sc_ec_genes.rds"
if (file.exists(gfile)) {                       # rebuild if it does not cover the current cells
  gg <- readRDS(gfile)
  if (!all(c("soma_joinid", g_sc) %in% names(gg)) || mean(sc_ec$soma_joinid %in% gg$soma_joinid) < 0.95) {
    message("kpmp_sc_ec_genes.rds does not match the current cells -> rebuilding (a few minutes)")
    file.rename(gfile, "kpmp_sc_ec_genes_OLD.rds")
  }
}
if (!file.exists(gfile)) saveRDS(get_genes(sc_ec$soma_joinid), gfile)
dsc <- sc_ec |> inner_join(readRDS(gfile) |> select(soma_joinid, all_of(g_sc)), by = "soma_joinid") |>
  mutate(gc_like = EHD3 > 0)
cat("Cells with gene counts:", nrow(dsc), "\n")
pbs <- dsc |> mutate(group = ifelse(gc_like, "glomerular-like (EHD3+)", "other EC (EHD3-)")) |>
  group_by(donor_id, sex, age, disease, group) |> filter(n() >= 20) |>
  summarise(n = n(), lib = sum(raw_sum), across(all_of(g_sc), sum), .groups = "drop") |>
  mutate(ESR1 = log2(ESR1 / lib * 1e6 + 1),
         tubular = log2((LRP2 + CUBN + SLC34A1 + UMOD) / lib * 1e6 + 1),
         podo = log2((NPHS1 + NPHS2) / lib * 1e6 + 1),
         sex = relevel(factor(as.character(sex)), ref = "male"),
         disease2 = relevel(factor(ifelse(disease == "normal", "healthy", "CKD")), ref = "healthy"),
         independent = !donor_id %in% nuc_donors)
print(pbs |> distinct(donor_id, sex, independent) |> count(independent, sex))
res_esr1_sc <- bind_rows(lapply(c("all", "independent"), function(w) {
  pbs |> filter(w == "all" | independent) |> group_by(type = group) |> group_modify(function(d, k) {
    co <- summary(lm(ESR1 ~ sex + age + disease2 + tubular + podo, data = d))$coefficients["sexfemale", ]
    tibble(data = "single cells", donors = w, estimate = co[1], se = co[2], p = co[4])
  }) |> ungroup()
}))
print(res_esr1_sc)
esr1_all <- bind_rows(res_esr1_nuc, res_esr1_sc)
saveRDS(esr1_all, "p2_esr1_kidney.rds")

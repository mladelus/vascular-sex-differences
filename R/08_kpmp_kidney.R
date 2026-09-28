# =============================================================================
# KIDNEY (KPMP): COHORT AND SOURCE FILE
# Kidney (cardiorenal): KPMP / Human Kidney Atlas v1.5 (Lake et al. 2023)
# Design: adults; healthy vs CKD; women vs men, esp. age 50+. AKI excluded (few women).
# Census labels are mostly generic "endothelial cell", so author subtypes (glomerular,
# peritubular, arteriolar, vasa recta) are read from the original .h5ad file.
# Needs: 00_setup.R
# Makes: kidney_ec_obs.rds, kpmp/kpmp_snRNA_v1.5.h5ad
# =============================================================================

# ---- 1. Candidate datasets and endothelial metadata ----
ds <- as.data.frame(census$get("census_info")$get("datasets")$read()$concat())
kid_ids <- ds |> filter(dataset_title %in% c(
  "Single-nucleus RNA-seq of the Adult Human Kidney (Version 1.5)",
  "Single-cell RNA-seq of the Adult Human Kidney (Version 1.5)",
  "Mature kidney dataset: full", "Living donor kidney")) |> select(dataset_id, dataset_title)

if (!file.exists("kidney_ec_obs.rds")) {
  kobs <- census$get("census_data")$get("homo_sapiens")$obs$read(
    value_filter = sprintf("dataset_id %%in%% c(%s)", paste0("'", kid_ids$dataset_id, "'", collapse = ",")),
    column_names = c("soma_joinid", "dataset_id", "donor_id", "sex", "development_stage", "disease",
                     "tissue", "cell_type", "assay", "suspension_type", "raw_sum"))$concat() |>
    as.data.frame() |> as_tibble() |> left_join(kid_ids, by = "dataset_id") |>
    filter(str_detect(cell_type, "endotheli"))
  saveRDS(kobs, "kidney_ec_obs.rds")
}
kobs <- readRDS("kidney_ec_obs.rds") |>
  mutate(age = parse_age(development_stage),                          # handles years and decades
         age_grp = case_when(is.na(age) ~ NA_character_, age < 18 ~ "child", age >= 50 ~ "50+", TRUE ~ "18-49"),
         data = if_else(str_detect(dataset_title, "nucleus"), "nuclei", "cells"))
v15 <- kobs |> filter(str_detect(dataset_title, "Version 1.5"), age_grp %in% c("18-49", "50+"))

# Adult donors by data type, disease, sex, age group
v15 |> distinct(data, disease, donor_id, sex, age_grp) |> count(data, disease, sex, age_grp) |>
  pivot_wider(names_from = c(sex, age_grp), values_from = n, values_fill = 0) |> print(n = Inf, width = Inf)
# Endothelial cells per donor
v15 |> count(data, disease, sex, donor_id) |> group_by(data, disease, sex) |>
  summarise(donors = n(), median_EC = median(n), min_EC = min(n), .groups = "drop") |> print(n = Inf)

# Do the cell and nucleus datasets share donors?
d_cells <- unique(as.character(v15$donor_id[v15$data == "cells"]))
d_nuc   <- unique(as.character(v15$donor_id[v15$data == "nuclei"]))
cat("cell donors:", length(d_cells), " nucleus donors:", length(d_nuc),
    " shared:", length(intersect(d_cells, d_nuc)), "\n")

# ---- 2. Download the original KPMP nucleus file (3.03 GB) ----
sn_id <- kid_ids$dataset_id[kid_ids$dataset_title == "Single-nucleus RNA-seq of the Adult Human Kidney (Version 1.5)"]
https_url <- sub("^s3://([^/]+)/", "https://\\1.s3.us-west-2.amazonaws.com/", get_source_h5ad_uri(sn_id)$uri)
options(timeout = 3600)
dir.create("kpmp", showWarnings = FALSE)
f <- "kpmp/kpmp_snRNA_v1.5.h5ad"
if (!file.exists(f) || file.size(f) < 3.0e9) download.file(https_url, f, mode = "wb")
round(file.size(f) / 1e9, 2)                                          # ~3.03

# ---- 3. List the label columns (reads names only, not data) ----
library(rhdf5)
filter <- dplyr::filter; select <- dplyr::select
h5ls(f, recursive = 2) |> filter(group == "/obs") |> pull(name)

# ---- Read author labels + clinical columns from the KPMP file ----
# Result: 304,989 nuclei; EC subtypes EC-PTC 7,567, EC-GC 4,029, EC-AEA 1,723, EC-DVR 4,611, EC-AVR 5,277.
# eGFR is stored as 10-unit bands (e.g. "50-59"); healthy donors mostly "unknown".
rd <- function(col) {                       # reads one obs column (either .h5ad category format)
  x <- h5read(f, paste0("/obs/", col))
  if (is.list(x) && all(c("categories", "codes") %in% names(x))) {
    i <- as.integer(x$codes) + 1; i[i < 1] <- NA; return(as.vector(x$categories[i]))
  }
  cats <- tryCatch(h5read(f, paste0("/obs/__categories/", col)), error = function(e) NULL)
  if (!is.null(cats)) { i <- as.integer(x) + 1; i[i < 1] <- NA; return(as.vector(cats[i])) }
  as.vector(x)
}
if (!file.exists("kpmp_obs_labels.rds")) {
  cols <- c("observation_joinid", "subclass.l2", "specimen", "eGFR", "hypertension", "diabetes_history")
  saveRDS(as_tibble(setNames(lapply(cols, rd), cols)), "kpmp_obs_labels.rds")
}
kp <- readRDS("kpmp_obs_labels.rds")
kp |> filter(str_detect(subclass.l2, "^(d|cyc)?EC")) |> count(subclass.l2) |> print(n = Inf)
kp |> count(hypertension); kp |> count(diabetes_history); kp |> count(eGFR) |> print(n = Inf)

# ---- Link labels to Census donor / sex / age / disease ----
# Result: all 304,989 nuclei matched by observation_joinid.
kc <- census$get("census_data")$get("homo_sapiens")$obs$read(
  value_filter = sprintf("dataset_id == '%s'", sn_id),
  column_names = c("soma_joinid", "observation_joinid", "donor_id", "sex", "development_stage",
                   "disease", "tissue", "assay", "raw_sum"))$concat() |> as.data.frame() |> as_tibble()
kk <- kc |> inner_join(kp, by = "observation_joinid")
c(census = nrow(kc), file = nrow(kp), matched = nrow(kk))

kec <- kk |> filter(subclass.l2 %in% c("EC-PTC", "EC-GC", "EC-AEA", "EC-DVR", "EC-AVR")) |>
  mutate(age = parse_age(development_stage), age50 = age >= 50,
         sex = relevel(factor(as.character(sex)), ref = "male"), disease = as.character(disease))
saveRDS(kec, "kpmp_ec_meta.rds")
kec |> distinct(donor_id, disease, sex, age50) |> count(disease, sex, age50) |> print(n = Inf)
kec |> count(donor_id, disease, sex, age, eGFR, hypertension, diabetes_history, subclass.l2) |>
  pivot_wider(names_from = subclass.l2, values_from = n, values_fill = 0) |>
  arrange(disease, sex, age) |> print(n = Inf, width = Inf)

# ---- Tissue sources and usable donors ----
# Adults only; acute kidney injury excluded (few women).
# Result (>= 20 nuclei, age 50+): EC-PTC  CKD 12 F / 14 M, healthy 5 F / 6 M
#                                 EC-GC   CKD 11 F / 8 M,  healthy 4 F / 4 M
#                                 EC-AEA  CKD 10 F / 6 M,  healthy 5 F / 3 M
kec <- kec |> filter(age >= 18, disease != "acute kidney failure")
kec |> distinct(donor_id, disease, sex, age, specimen, tissue) |>
  count(disease, specimen, tissue) |> print(n = Inf, width = Inf)          # tissue source by disease (procurement check)
kec |> filter(disease == "normal") |> distinct(donor_id, sex, age, specimen, tissue) |>
  arrange(specimen, sex, age) |> print(n = Inf, width = Inf)
elig <- kec |> filter(subclass.l2 %in% c("EC-PTC", "EC-GC", "EC-AEA")) |>
  count(donor_id, disease, sex, age50, subclass.l2) |> filter(n >= 20)
elig |> count(subclass.l2, disease, age50, sex) |>
  pivot_wider(names_from = sex, values_from = n, values_fill = 0) |> print(n = Inf)

# ---- Download + score kidney EC nuclei, one donor at a time ----
# Skips donors already saved: if it stops, just run this block again.
k_sets <- c(hg_sets, ox_sets[setdiff(names(ox_sets), names(hg_sets))])   # senescence, inflammation, oxidative, stress
tub <- c("LRP2", "CUBN", "SLC34A1", "UMOD", "SLC12A1", "SLC12A3", "AQP2")    # tubule genes: ambient-RNA check
ksel <- kec |> filter(subclass.l2 %in% c("EC-PTC", "EC-GC", "EC-AEA"))
dir.create("kp_ec", showWarnings = FALSE)
for (d in unique(as.character(ksel$donor_id))) {
  out <- file.path("kp_ec", paste0(gsub("[^A-Za-z0-9-]", "_", d), ".rds")); if (file.exists(out)) next
  message(format(Sys.time(), "%H:%M"), "  ", d)
  seu <- get_seurat(census, organism = "Homo sapiens", obs_coords = ksel$soma_joinid[ksel$donor_id == d],
                    obs_column_names = "soma_joinid")
  counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts"))
  meta <- ksel[match(seu$soma_joinid, ksel$soma_joinid), ]; rm(seu); gc(verbose = FALSE)
  cp <- function(g) as.numeric(Matrix::colSums(counts[intersect(g, rownames(counts)), , drop = FALSE]))
  g1 <- function(g) if (g %in% rownames(counts)) as.numeric(counts[g, ]) else 0
  sc <- ScoreSignatures_UCell(counts, features = k_sets, ncores = 1)
  meta <- meta |> mutate(tub_umi = cp(tub), CDKN2A = g1("CDKN2A"), CDKN1A = g1("CDKN1A"), MKI67 = g1("MKI67")) |>
    bind_cols(as_tibble(sc))
  pbl <- list()                                              # pseudobulk per donor x subtype (>= 20 nuclei)
  for (s in unique(meta$subclass.l2)) { idx <- which(meta$subclass.l2 == s)
    if (length(idx) >= 20) pbl[[paste(d, s, sep = "|")]] <- Matrix::rowSums(counts[, idx, drop = FALSE]) }
  saveRDS(list(meta = meta, pb = pbl), out)
  rm(counts, sc, pbl); gc(verbose = FALSE)
}
length(list.files("kp_ec"))

# ---- Kidney — does CKD change women's and men's capillaries differently? ----
# Needs Hs, cam2 and auto from script 09 (run_all.R runs 09 before 08).
# Result (EC-PTC: healthy 11 F / 9 M, CKD 16 / 16): sex x CKD interaction for interferon-g (D 7e-4) and
# TNF/NF-kB (D 0.012): women higher than men when healthy (IFN-g 0.003), lower in CKD (TNF 0.01).
# EC-GC: IFN-g interaction 6e-10 (women higher healthy 1e-4, lower CKD 3e-5). EC-AEA same direction (0.057).
# PTC age 50+: TNF interaction 1e-5, IFN-g 4e-4. OXPHOS lower in women with CKD (1e-5).
# CKD vs healthy in men: TNF and hypoxia LOWER in CKD -> healthy tissue source differs (procurement).
# Donor-level UCell scores (SenMayo, inflammatory, NRF2, p16): no sex or interaction effects.
# Adults, healthy vs CKD (AKI excluded). One pseudobulk per donor x EC subtype (>= 20 nuclei).
# Model: expression ~ sex * disease + age + tubular ambient RNA. Key term = sex x disease interaction.
# Primary: peritubular capillaries (EC-PTC); secondary: glomerular (EC-GC); non-capillary: arterioles (EC-AEA).
filter <- dplyr::filter; select <- dplyr::select; count <- dplyr::count; first <- dplyr::first
kp_all <- lapply(list.files("kp_ec", full.names = TRUE), readRDS)
km  <- bind_rows(lapply(kp_all, `[[`, "meta")) |>
  mutate(donor_id = as.character(donor_id), sex = relevel(factor(as.character(sex)), ref = "male"),
         disease2 = relevel(factor(ifelse(disease == "normal", "healthy", "CKD")), ref = "healthy"))
kpb <- unlist(lapply(kp_all, `[[`, "pb"), recursive = FALSE)
Kset <- c(Hs[c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
               "HALLMARK_HYPOXIA", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
               "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_ANGIOGENESIS")],
          list(SenMayo = hg_sets$SenMayo, NRF2 = ox_sets$NRF2, Antioxidant = ox_sets$Antioxidant))
camk <- function(st) { idx <- ids2indices(Kset, names(st)); idx <- idx[lengths(idx) >= 5]
  cameraPR(st, idx) |> tibble::rownames_to_column("pathway") |> as_tibble() }

run_kid <- function(sub, older = FALSE) {
  p <- kpb[str_detect(names(kpb), paste0("\\|", sub, "$"))]
  g <- Reduce(intersect, lapply(p, names)); m <- sapply(p, function(v) v[g])
  s <- tibble(id = colnames(m), donor_id = sub("\\|.*", "", colnames(m))) |>
    left_join(km |> filter(subclass.l2 == sub) |> group_by(donor_id) |>
                summarise(sex = first(sex), age = first(age), disease2 = first(disease2),
                          tub = mean(log1p(tub_umi / raw_sum * 1e4)), .groups = "drop"), by = "donor_id")
  if (older) { k <- s$age >= 50; m <- m[, k]; s <- s[k, ] }
  X  <- model.matrix(~ sex * disease2 + age + tub, data = s)
  yy <- DGEList(m[rownames(m) %in% auto, ])
  yy <- calcNormFactors(yy[filterByExpr(yy, X), , keep.lib.sizes = FALSE])
  f0 <- lmFit(voom(yy, X), X); ff <- eBayes(f0)
  cv <- setNames(rep(0, ncol(X)), colnames(X)); cv[c("sexfemale", "sexfemale:disease2CKD")] <- 1
  fc <- eBayes(contrasts.fit(f0, cv))
  cat("\n==", sub, if (older) "(age 50+)" else "(all adults)", "==\n"); print(table(s$disease2, s$sex))
  list(women_vs_men_healthy = ff$t[, "sexfemale"], women_vs_men_CKD = fc$t[, 1],
       CKD_effect_in_men = ff$t[, "disease2CKD"], sex_x_CKD = ff$t[, "sexfemale:disease2CKD"])
}
show_kid <- function(r) bind_rows(lapply(r, camk), .id = "test") |>
  mutate(val = paste(substr(Direction, 1, 1), signif(FDR, 2)), pathway = sub("HALLMARK_", "", pathway)) |>
  select(test, pathway, val) |> pivot_wider(names_from = test, values_from = val) |> print(width = Inf)
for (sub in c("EC-PTC", "EC-GC", "EC-AEA")) show_kid(run_kid(sub))
show_kid(run_kid("EC-PTC", older = TRUE))                              # sensitivity: age 50+ only

# Donor-level scores in peritubular capillaries (plain numbers + interaction test)
kd <- km |> filter(subclass.l2 == "EC-PTC") |> group_by(donor_id, sex, age, disease2) |>
  summarise(n = n(), across(c(SenMayo_UCell, Inflammatory_UCell, TNFA_UCell, NRF2_UCell, Antioxidant_UCell), mean),
            p16 = mean(CDKN2A > 0 & MKI67 == 0), tub = mean(log1p(tub_umi / raw_sum * 1e4)), .groups = "drop") |>
  filter(n >= 20)
kd |> group_by(disease2, sex) |> summarise(donors = n(), across(c(SenMayo_UCell, Inflammatory_UCell, NRF2_UCell, p16),
                                           ~ signif(mean(.x), 3)), .groups = "drop") |> print(width = Inf)
bind_rows(lapply(c("SenMayo_UCell", "Inflammatory_UCell", "TNFA_UCell", "NRF2_UCell", "Antioxidant_UCell", "p16"), function(y) {
  co <- summary(lm(as.formula(paste(y, "~ sex * disease2 + age + tub")), data = kd))$coefficients
  tibble(outcome = y, women_vs_men_healthy = signif(co["sexfemale", 1], 2), p_healthy = signif(co["sexfemale", 4], 2),
         sex_x_CKD = signif(co["sexfemale:disease2CKD", 1], 2), p_interaction = signif(co["sexfemale:disease2CKD", 4], 2))
})) |> print(width = Inf)

# ---- Last checks on the kidney interferon/TNF interaction ----
# Result: healthy tissue from several collection series, mix somewhat unbalanced by sex (one series 6 F / 1 M).
# Interaction unchanged with immune adjustment: IFN-g PTC 7e-4, GC 1e-8; TNF PTC 0.015.
# Drivers: IFI27, TXNIP, CASP8, IRF9, SERPING1, SELP, HLA-DRB1, TNFAIP3, HLA-DMA, IL15RA, IFI44L, CIITA.
# 1. Tissue source by sex (healthy tissue may come from different sources than CKD biopsies)
km |> distinct(donor_id, disease2, sex, specimen) |> count(disease2, specimen, sex) |>
  pivot_wider(names_from = sex, values_from = n, values_fill = 0) |> print(n = Inf, width = Inf)
# 2. Does the interaction survive adjusting for immune-cell (ambient) RNA? Which genes drive it?
run_kid2 <- function(sub, adj_imm = FALSE) {
  p <- kpb[str_detect(names(kpb), paste0("\\|", sub, "$"))]
  g <- Reduce(intersect, lapply(p, names)); m <- sapply(p, function(v) v[g])
  s <- tibble(id = colnames(m), donor_id = sub("\\|.*", "", colnames(m))) |>
    left_join(km |> filter(subclass.l2 == sub) |> group_by(donor_id) |>
                summarise(sex = first(sex), age = first(age), disease2 = first(disease2),
                          tub = mean(log1p(tub_umi / raw_sum * 1e4)), .groups = "drop"), by = "donor_id")
  lc <- edgeR::cpm(DGEList(m), log = TRUE)
  s$imm <- colMeans(lc[intersect(c("PTPRC", "CD3E", "CD68", "LYZ", "CD163"), rownames(lc)), ])
  X <- model.matrix(if (adj_imm) ~ sex * disease2 + age + tub + imm else ~ sex * disease2 + age + tub, data = s)
  yy <- DGEList(m[rownames(m) %in% auto, ])
  yy <- calcNormFactors(yy[filterByExpr(yy, X), , keep.lib.sizes = FALSE])
  ff <- eBayes(lmFit(voom(yy, X), X))
  list(t = ff$t[, "sexfemale:disease2CKD"], tt = topTable(ff, coef = "sexfemale:disease2CKD", n = Inf))
}
k_ptc <- run_kid2("EC-PTC"); k_ptc_imm <- run_kid2("EC-PTC", TRUE); k_gc_imm <- run_kid2("EC-GC", TRUE)
bind_rows(list(PTC = camk(k_ptc$t), PTC_plus_immune = camk(k_ptc_imm$t), GC_plus_immune = camk(k_gc_imm$t)), .id = "test") |>
  filter(pathway %in% c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")) |>
  mutate(val = paste(substr(Direction, 1, 1), signif(FDR, 2)), pathway = sub("HALLMARK_", "", pathway)) |>
  select(test, pathway, val) |> pivot_wider(names_from = test, values_from = val) |> print(width = Inf)
g_ifn <- intersect(Hs$HALLMARK_INTERFERON_GAMMA_RESPONSE, rownames(k_ptc$tt))
k_ptc$tt[g_ifn, ] |> tibble::rownames_to_column("gene") |> as_tibble() |> arrange(t) |>
  select(gene, logFC, P.Value) |> mutate(across(where(is.numeric), ~ signif(.x, 2))) |> head(12)

# ---- Stricter gene-set tests for the sex x CKD interaction ----
# Result: IFN-g interaction glomerular fry p = 0.007 (FDR 0.014), camera (estimated correlation) 0.0046;
# peritubular fry 0.06, camera 0.072.
# TNF-a/NF-kB not supported by fry (p >= 0.19). A simple per-donor mean score did not show it (p 0.38 / 0.87).
kid_check <- function(sub) {
  p <- kpb[str_detect(names(kpb), paste0("\\|", sub, "$"))]
  g <- Reduce(intersect, lapply(p, names)); m <- sapply(p, function(v) v[g])
  s <- tibble(donor_id = sub("\\|.*", "", colnames(m))) |>
    left_join(km |> filter(subclass.l2 == sub) |> group_by(donor_id) |>
                summarise(sex = first(sex), age = first(age), disease2 = first(disease2),
                          tub = mean(log1p(tub_umi / raw_sum * 1e4)), .groups = "drop"), by = "donor_id")
  X <- model.matrix(~ sex * disease2 + age + tub, data = s)
  yy <- DGEList(m); yy <- calcNormFactors(yy[filterByExpr(yy, X), , keep.lib.sizes = FALSE])
  v <- voom(yy, X); cf <- which(colnames(X) == "sexfemale:disease2CKD")
  idx <- ids2indices(Hs[c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")], rownames(v))
  list(fry = fry(v, idx, design = X, contrast = cf), camera_estimated_cor = camera(v, idx, design = X, contrast = cf, inter.gene.cor = NA))   # NA = estimate from data
}
kid_check("EC-PTC"); kid_check("EC-GC")

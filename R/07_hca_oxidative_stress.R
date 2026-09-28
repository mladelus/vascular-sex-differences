# =============================================================================
# OXIDATIVE STRESS: PERICYTES AND CAPILLARY ECs
# Pericytes and capillary ECs: oxidative-stress programs
# Hypothesis (pre-specified): postmenopausal women's capillary-pericyte units differ from
# age-matched men's in oxidative-stress defence/production genes. RNA measures response
# programs, not ROS or oxidative damage.
# Primary comparison: women 50+ vs men 50+ (5 vs 7), BH-corrected across programs.
# Needs: 00_setup.R, 06 outputs (hca_vascular_nuclei_typed.rds, hg_ec_nuclei_meta.rds)
# Makes: hg_ox/<donor>.rds, hg_oxidative_scores.rds
# =============================================================================
vt  <- readRDS("hca_vascular_nuclei_typed.rds")
hgm <- readRDS("hg_ec_nuclei_meta.rds")

set.seed(2)
peri_sel <- vt |> filter(cell_type == "pericyte") |> group_by(donor_id) |> slice_sample(n = 3000) |> ungroup()
cap_ids  <- hgm |> filter(type == "capillary", dbl_class == "singlet") |> pull(soma_joinid)
sel2 <- bind_rows(peri_sel |> mutate(kind = "pericyte"),
                  vt |> filter(soma_joinid %in% cap_ids) |> mutate(kind = "capillary EC"))

dir.create("hg_ox", showWarnings = FALSE)
for (d in unique(as.character(sel2$donor_id))) {
  out <- file.path("hg_ox", paste0(d, ".rds")); if (file.exists(out)) next
  message(format(Sys.time(), "%H:%M"), "  ", d)
  seu <- get_seurat(census, organism = "Homo sapiens", obs_coords = sel2$soma_joinid[sel2$donor_id == d],
                    obs_column_names = "soma_joinid")
  counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts"))
  meta <- sel2[match(seu$soma_joinid, sel2$soma_joinid), ]; rm(seu); gc(verbose = FALSE)
  ec_umi <- Matrix::colSums(counts[intersect(c("PECAM1", "CDH5", "VWF"), rownames(counts)), , drop = FALSE])
  keep <- !(meta$kind == "pericyte" & ec_umi >= 2)          # drop pericytes carrying endothelial RNA (likely doublets)
  sc <- ScoreSignatures_UCell(counts[, keep], features = ox_sets, ncores = 1)
  saveRDS(bind_cols(meta[keep, ], as_tibble(sc)), out)
  rm(counts, sc); gc(verbose = FALSE)
}
ox <- bind_rows(lapply(list.files("hg_ox", full.names = TRUE), readRDS)) |>
  mutate(age = parse_age(development_stage), age50 = age >= 50,
         sex = relevel(factor(as.character(sex)), ref = "male"))
saveRDS(ox, "hg_oxidative_scores.rds")

# Donor x region summaries; women vs men aged 50+
oxr <- ox |> group_by(kind, donor_id, sex, age, age50, tissue, assay) |>
  summarise(n = n(), across(ends_with("_UCell"), mean), .groups = "drop") |> filter(n >= 20)
progs <- c("NRF2_UCell", "Antioxidant_UCell", "ROS_sources_UCell", "ROS_hallmark_UCell", "Contractile_UCell", "eNOS_UCell")
res_ox <- bind_rows(lapply(c("pericyte", "capillary EC"), function(k) bind_rows(lapply(progs, function(y) {
  d <- filter(oxr, kind == k, age50)
  co <- summary(lmerTest::lmer(as.formula(paste(y, "~ sex + age + tissue + assay + Stress_UCell + (1 | donor_id)")),
                               data = d))$coefficients["sexfemale", ]
  tibble(cells = k, program = sub("_UCell", "", y),
         women = n_distinct(d$donor_id[d$sex == "female"]), men = n_distinct(d$donor_id[d$sex == "male"]),
         estimate = signif(co[1], 3), conf.low = signif(co[1] - 1.96 * co[2], 3),
         conf.high = signif(co[1] + 1.96 * co[2], 3), p = co["Pr(>|t|)"])
}))))
res_ox |> mutate(p_adj = signif(p.adjust(p, "BH"), 2), p = signif(p, 2)) |> print(n = Inf, width = Inf)

# Per-donor view
ox |> group_by(kind, donor_id, sex, age) |>
  summarise(across(c(NRF2_UCell, Antioxidant_UCell, ROS_sources_UCell), ~ round(mean(.x), 4)), n = n(),
            .groups = "drop") |> arrange(kind, sex, age) |> print(n = Inf)

# ---- Sensitivity — without outlier D6; adjusted for donor series (D vs H) ----
# Result: pericyte NRF2 estimate halved and antioxidant went to ~0 without D6; series changed little.
oxr <- oxr |> mutate(series = substr(as.character(donor_id), 1, 1))   # "D" or "H" (likely collection site)
sens <- function(k, y, data, extra = "") {
  d <- filter(data, kind == k, age50)
  co <- summary(lmerTest::lmer(as.formula(paste(y, "~ sex + age + tissue + assay + Stress_UCell", extra, "+ (1 | donor_id)")),
                               data = d))$coefficients["sexfemale", ]
  tibble(estimate = signif(co[1], 3), p = signif(co["Pr(>|t|)"], 2))
}
res_sens <- bind_rows(lapply(c("pericyte", "capillary EC"), function(k) bind_rows(lapply(progs, function(y) {
  bind_cols(tibble(cells = k, program = sub("_UCell", "", y)),
            sens(k, y, oxr) |> rename(est_main = estimate, p_main = p),
            sens(k, y, filter(oxr, donor_id != "D6")) |> rename(est_noD6 = estimate, p_noD6 = p),
            sens(k, y, oxr, "+ series") |> rename(est_series = estimate, p_series = p))
}))))
print(res_sens, n = Inf, width = Inf)
oxr |> filter(age50) |> distinct(donor_id, sex, series) |> dplyr::count(series, sex)   # D: 4 M / 4 F; H: 3 M / 1 F

# ---- Stress-response coupling (exploratory) ----
# Do cells under more stress switch on more defence genes, and is that weaker in women?
# Result: coupling present in both sexes; no sex difference (all p >= 0.33).
ov <- sapply(ox_sets, function(a) sapply(ox_sets, function(b) length(intersect(a, b))))
print(ov)                                                   # NRF2/Antioxidant share no genes with ROS_sources/Stress
cz <- function(x) if (sd(x) > 0) (x - mean(x)) / sd(x) else rep(NA_real_, length(x))
sl <- ox |> filter(age50) |> group_by(kind, donor_id, tissue, assay) |> filter(n() >= 50) |>
  mutate(across(ends_with("_UCell"), cz)) |> ungroup()
pairs <- expand.grid(defence = c("NRF2_UCell", "Antioxidant_UCell"),
                     input   = c("ROS_sources_UCell", "Stress_UCell"), stringsAsFactors = FALSE)
don <- sl |> group_by(kind, donor_id, sex, age) |>
  group_modify(function(d, key) bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
    x <- d[[pairs$input[i]]]; y <- d[[pairs$defence[i]]]; ok <- is.finite(x) & is.finite(y)
    tibble(defence = sub("_UCell", "", pairs$defence[i]), input = sub("_UCell", "", pairs$input[i]),
           n = sum(ok), slope = if (sum(ok) >= 50) unname(coef(lm(y[ok] ~ x[ok]))[2]) else NA_real_)
  }))) |> ungroup()
res_slope <- don |> filter(!is.na(slope)) |> group_by(kind, defence, input) |>
  group_modify(function(d, key) {
    co <- summary(lm(slope ~ sex + age, data = d))$coefficients["sexfemale", ]
    tibble(women = sum(d$sex == "female"), men = sum(d$sex == "male"),
           mean_women = round(mean(d$slope[d$sex == "female"]), 3),
           mean_men = round(mean(d$slope[d$sex == "male"]), 3),
           diff = signif(co[1], 3), p = signif(co[4], 2))
  }) |> ungroup() |> mutate(p_adj = signif(p.adjust(p, "BH"), 2))
print(res_slope, n = Inf, width = Inf)

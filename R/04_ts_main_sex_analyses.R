# =============================================================================
# TABULA SAPIENS: MAIN SEX ANALYSES (Table 2)
# Tabula Sapiens: main sex analyses (women = reference group)
# Needs: 00_setup.R, 02 + 03 outputs
# Makes: ec_cells_scored.rds (final labels), ts_analysis_dat.rds, donor_senescence_p16.rds
# Estimates are women minus men ("sexfemale"): OR > 1 or estimate > 0 = higher in women.
# In R the vessel classes are named "microvascular" (= capillary) and "large vessel"
# (= non-capillary: arterial + venous, incl. post-capillary venules), as worded in the manuscript.
# =============================================================================
cells <- readRDS("ec_cells_scored.rds")
ec    <- readRDS("ts_ec_obs.rds")
hres  <- readRDS("ts_heart_labels.rds"); hm <- readRDS("ts_heart_markers.rds")

# ---- 1. Final vessel classes ----
cells <- cells |>
  select(-any_of(c("age", "age_grp", "heart_subtype", "NPR3", "ACKR1", "CA4", "RGCC"))) |>
  left_join(select(ec, soma_joinid, age), by = "soma_joinid") |>
  left_join(hres, by = "soma_joinid") |>
  left_join(select(hm, soma_joinid, NPR3, ACKR1), by = "soma_joinid") |>
  mutate(
    age_grp      = if_else(age >= 50, "50+", "<50"),
    heart_region = case_when(tissue_general != "heart" ~ NA_character_,
                             tissue == "coronary artery" ~ "coronary artery", TRUE ~ "myocardium"),
    subtype_v2   = case_when(
      tissue_general == "heart" & heart_subtype == "venous" & NPR3 == 1 & ACKR1 == 0 ~ "endocardial",
      tissue_general == "heart" ~ heart_subtype,
      TRUE ~ subtype_final),
    source_v2    = case_when(tissue_general == "heart" ~ "heart atlas",
                             !is.na(label_subtype) ~ "expert label",
                             subtype_final != "uncertain" ~ "consensus", TRUE ~ "uncertain"),
    vessel_class = case_when(subtype_v2 == "capillary" ~ "microvascular",
                             subtype_v2 %in% c("arterial", "venous") ~ "large vessel",
                             TRUE ~ "excluded"))
saveRDS(cells, "ec_cells_scored.rds")
cells |> count(source_v2, vessel_class)

# Descriptive: per-donor % p16+ (all organs, before excluding reproductive organs)
donor_sen <- cells |> filter(vessel_class != "excluded") |>
  group_by(donor_id, sex, age, vessel_class) |>
  summarise(n = n(), pct_p16 = 100 * mean(sen_p16), .groups = "drop") |> filter(n >= 20)
saveRDS(donor_sen, "donor_senescence_p16.rds")
donor_sen |> group_by(vessel_class, sex) |>
  summarise(donors = n(), median_pct_p16 = round(median(pct_p16), 2), .groups = "drop")

# ---- 2. Analysis set: exclude reproductive organs; women as reference ----
dat <- cells |>
  filter(vessel_class != "excluded", !str_detect(as.character(tissue), REPRO)) |>
  mutate(sex          = relevel(factor(as.character(sex)), ref = "male"),
         vessel_class = factor(vessel_class, c("large vessel", "microvascular")),
         age10 = age / 10, logdepth = log10(raw_sum),
         organ = droplevels(factor(tissue_general)), chem = factor(assay),
         postmeno = factor(if_else(age >= 50, "50+", "<50"), c("<50", "50+"))) |>
  group_by(subtype_v2) |>
  mutate(sen_p21 = CDKN1A_cp10k > quantile(CDKN1A_cp10k, 0.9) & CDKN1A_cp10k > 0) |> ungroup() |>
  mutate(senescent_any = sen_p16 | sen_p21)
saveRDS(dat, "ts_analysis_dat.rds")
nrow(dat)                                                                     # 51,083
dat |> distinct(donor_id, sex, age) |> arrange(sex, age) |> print(n = Inf)   # 6 F, 8 M

nc  <- filter(dat, vessel_class == "large vessel")
cap <- filter(dat, vessel_class == "microvascular")

# ---- 3. Senescence: women vs men WITHIN each vessel class (Table 2) ----
for (vc in levels(dat$vessel_class)) {
  d <- filter(dat, vessel_class == vc)
  cat("\n######", vc, "######\n-- p16 --\n")
  print(show(glmer(sen_p16 ~ sex + age10 + logdepth + chem + (1 | donor_id) + (1 | organ),
                   family = binomial, data = d, control = glmer_ctrl)))
  cat("-- p21-high --\n")
  print(show(glmer(sen_p21 ~ sex + age10 + logdepth + chem + (1 | donor_id) + (1 | organ),
                   family = binomial, data = d, control = glmer_ctrl)))
  cat("-- SenMayo --\n")
  print(st(lmerTest::lmer(SenMayo_UCell ~ sex + age10 + logdepth + chem + (1 | donor_id) + (1 | organ), data = d),
           "sex|age10"))
}

# Does the sex difference differ between vessel classes? (random slopes for vessel class)
fit_p16_rs <- glmer(sen_p16 ~ sex * vessel_class + age10 + logdepth + chem +
                      (1 + vessel_class | donor_id) + (1 | organ),
                    family = binomial, data = dat, control = glmer_ctrl)
show(fit_p16_rs)

# Sensitivity: expert + heart-atlas labels only
show(glmer(sen_p16 ~ sex * vessel_class + age10 + logdepth + chem + (1 | donor_id) + (1 | organ),
           family = binomial, data = filter(dat, source_v2 != "consensus"), control = glmer_ctrl))

# Heart myocardial capillaries: descriptive only (3 F, 2 M)
cells |> filter(heart_region == "myocardium", subtype_v2 == "capillary") |>
  group_by(donor_id, sex, age) |>
  summarise(cells = n(), pct_p16 = round(100 * mean(sen_p16), 2),
            senmayo = round(mean(SenMayo_UCell), 4), .groups = "drop") |> arrange(sex, age)

# ---- 4. Inflammation: sex difference, stress adjustment, senescent cells removed (Results text) ----
run <- function(y, d, extra = "") {
  f <- as.formula(paste(y, "~ sex + age10 + logdepth + chem", extra, "+ (1 | donor_id) + (1 | organ)"))
  st(lmerTest::lmer(f, data = d), "^sexfemale$")
}
out <- list()
for (vc in levels(dat$vessel_class)) for (y in c("Inflammatory_UCell", "TNFA_UCell")) {
  d <- filter(dat, vessel_class == vc)
  out[[length(out) + 1]] <- bind_rows(
    cbind(model = "1 basic",             run(y, d)),
    cbind(model = "2 + stress",          run(y, d, "+ Stress_UCell")),
    cbind(model = "3 senescent removed", run(y, filter(d, !senescent_any), "+ Stress_UCell"))) |>
    mutate(vessel_class = vc, score = y, .before = 1)
}
bind_rows(out) |> select(-term) |> as_tibble() |> print(n = Inf, width = Inf)

# Robustness: one value per donor x organ (key check — the non-capillary difference did NOT hold)
dl <- nc |> group_by(donor_id, sex, age, postmeno, organ) |>
  summarise(infl = mean(Inflammatory_UCell), stress = mean(Stress_UCell), depth = mean(logdepth),
            n = n(), .groups = "drop") |> filter(n >= 20)
st(lmerTest::lmer(infl ~ sex + I(age / 10) + stress + depth + (1 | donor_id) + (1 | organ), data = dl))

# ---- 5. Menopause (exploratory; only 2 women and 3 men under 50) ----
for (nm in c("non-capillary", "capillary")) {
  d <- if (nm == "non-capillary") nc else cap
  cat("\n==", nm, ": sex x postmenopausal proxy ==\n")
  print(st(lmerTest::lmer(Inflammatory_UCell ~ sex * postmeno + logdepth + chem + Stress_UCell +
                            (1 | donor_id) + (1 | organ), data = d), "sex|postmeno"))
  cat("==", nm, ": sex x age per decade ==\n")
  print(st(lmerTest::lmer(Inflammatory_UCell ~ sex * age10 + logdepth + chem + Stress_UCell +
                            (1 | donor_id) + (1 | organ), data = d), "sex|age10"))
}
# Donor x organ level
st(lmerTest::lmer(infl ~ sex * postmeno + stress + depth + (1 | donor_id) + (1 | organ), data = dl), "sex|postmeno")
# Per-donor means (non-capillary)
nc |> group_by(donor_id, sex, age) |>
  summarise(inflammatory = round(mean(Inflammatory_UCell), 4), cells = n(), .groups = "drop") |> arrange(sex, age)
# Leave out each donor under 50
for (dd in unique(as.character(nc$donor_id[nc$age < 50]))) {
  r <- st(lmerTest::lmer(Inflammatory_UCell ~ sex * postmeno + logdepth + chem + Stress_UCell +
                           (1 | donor_id) + (1 | organ), data = filter(nc, donor_id != dd)), "sexfemale:postmeno")
  cat("without", dd, ": "); print(r[, -1], row.names = FALSE)
}

# ---- 6. Within women: change with age (6 women) ----
for (nm in c("non-capillary", "capillary")) {
  d <- filter(if (nm == "non-capillary") nc else cap, sex == "female")
  cat("\n== Women only,", nm, "(per decade) ==\n")
  for (y in c("Inflammatory_UCell", "SenMayo_UCell")) {
    r <- st(lmerTest::lmer(as.formula(paste(y, "~ age10 + logdepth + chem + Stress_UCell + (1 | donor_id) + (1 | organ)")),
                           data = d), "age10"); cat(y, ": "); print(r[, -1], row.names = FALSE)
  }
  print(show(glmer(sen_p16 ~ age10 + logdepth + chem + (1 | donor_id) + (1 | organ),
                   family = binomial, data = d, control = glmer_ctrl), "age10"))
}

# ---- 7. Female-specific tissues (descriptive) ----
cells |> filter(sex == "female", vessel_class != "excluded",
                str_detect(as.character(tissue), "ovary|uter|endometri|myometri|mammary|fallopian|vagin|cervix")) |>
  group_by(tissue, vessel_class) |>
  summarise(donors = n_distinct(donor_id), cells = n(), pct_p16 = round(100 * mean(sen_p16), 2),
            inflammatory = round(mean(Inflammatory_UCell), 4), senmayo = round(mean(SenMayo_UCell), 4),
            .groups = "drop") |> print(n = Inf)

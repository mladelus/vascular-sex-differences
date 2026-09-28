# =============================================================================
# SEX-CHROMOSOME + HORMONE-RECEPTOR PANEL (saved pseudobulks; no downloads)
# Uses the per-donor pseudobulks already on disk (they keep ALL genes, incl. X and Y):
#   hg_nuclei/*.rds (Heart Cell Atlas), rep_heart/*.rds (independent LV), kp_ec/*.rds (KPMP nuclei).
# One sample = donor x vessel type (>= 20 nuclei; >= 30 in the independent LV cohort).
# Expression = log2(CPM + 1), library size = all genes in the pseudobulk.
# Model per gene, cohort and vessel type: log2CPM ~ sex + age (+ CKD in kidney); BH within cohort.
# Makes: p2_panel_long.rds (long), p2_sex_panel_results.rds (res_p2)
# =============================================================================
filter <- dplyr::filter; select <- dplyr::select; count <- dplyr::count
first  <- dplyr::first;  rename <- dplyr::rename; desc <- dplyr::desc

genes_x  <- c("XIST", "JPX", "KDM6A", "KDM5C", "DDX3X", "EIF1AX", "ZFX", "USP9X", "RPS4X", "NLGN4X")
genes_y  <- c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY", "USP9Y", "NLGN4Y")
genes_hr <- c("AR", "ESR1", "ESR2", "GPER1", "PGR", "CYP19A1", "SYNE2")      # SYNE2 overlaps ESR2 (check)
tub      <- c("LRP2", "CUBN", "SLC34A1", "UMOD", "SLC12A1", "SLC12A3", "AQP2")
panel    <- unique(c(genes_x, genes_y, genes_hr, tub))

# ---- 1. Pseudobulk matrices: donor x vessel type ----
# (a) Heart Cell Atlas: sum the donor x region pseudobulks within donor x type
if (!exists("hg"))  hg  <- lapply(list.files("hg_nuclei", full.names = TRUE), readRDS)
if (!exists("hgm")) hgm <- readRDS("hg_ec_nuclei_meta.rds")
agg <- function(pbl, key_fun) {                     # sum a named list of count vectors by key
  keys <- key_fun(names(pbl)); g <- Reduce(intersect, lapply(pbl, names))
  sapply(split(names(pbl), keys), function(k) Reduce(`+`, lapply(pbl[k], function(v) v[g])))
}
p_hca <- unlist(lapply(hg, `[[`, "pb"), recursive = FALSE)
m1 <- agg(p_hca, function(n) sub("^([^|]+)\\|[^|]+\\|([^|]+)$", "\\1|\\2", n))   # donor|type
s1 <- tibble(id = colnames(m1), donor_id = sub("\\|.*", "", id), type = sub(".*\\|", "", id)) |>
  left_join(hgm |> mutate(donor_id = as.character(donor_id)) |> group_by(donor_id) |>
              summarise(sex = first(as.character(sex)), age = parse_age(first(development_stage)), .groups = "drop"),
            by = "donor_id") |> mutate(cohort = "Heart Cell Atlas", disease2 = "healthy")

# (b) Independent LV cohort
if (!exists("rep") || is.function(rep)) rep <-   # base R also has a function called rep
   lapply(list.files("rep_heart", full.names = TRUE), readRDS)
if (!exists("dinfo")) dinfo <- readRDS("rep_heart_lv_ec_obs.rds") |> distinct(donor_id, sex, age, assay) |>
  mutate(donor_id = as.character(donor_id))
p_rep <- list()
for (r in rep) for (tp in names(r$pb)) if (r$n[[tp]] >= 30) p_rep[[paste(r$donor, tp, sep = "|")]] <- r$pb[[tp]]
g2 <- Reduce(intersect, lapply(p_rep, names)); m2 <- sapply(p_rep, function(v) v[g2])
s2 <- tibble(id = colnames(m2), donor_id = sub("\\|.*", "", id), type = sub(".*\\|", "", id)) |>
  left_join(dinfo |> mutate(sex = as.character(sex)) |> select(donor_id, sex, age), by = "donor_id") |>
  mutate(cohort = "Independent LV", disease2 = "healthy")

# (c) KPMP nuclei (adults, healthy vs CKD; AKI already excluded in kp_ec)
if (!exists("km") || !exists("kpb")) {
  kp_all <- lapply(list.files("kp_ec", full.names = TRUE), readRDS)
  km  <- bind_rows(lapply(kp_all, `[[`, "meta")) |>
    mutate(donor_id = as.character(donor_id), sex = relevel(factor(as.character(sex)), ref = "male"),
           disease2 = relevel(factor(ifelse(disease == "normal", "healthy", "CKD")), ref = "healthy"))
  kpb <- unlist(lapply(kp_all, `[[`, "pb"), recursive = FALSE)
}
if (!"disease2" %in% names(km))                          # older sessions: add healthy/CKD label
  km <- km |> mutate(disease2 = relevel(factor(ifelse(disease == "normal", "healthy", "CKD")), ref = "healthy"))
g3 <- Reduce(intersect, lapply(kpb, names)); m3 <- sapply(kpb, function(v) v[g3])
s3 <- tibble(id = colnames(m3), donor_id = sub("\\|.*", "", id), type = sub(".*\\|", "", id)) |>
  left_join(km |> group_by(donor_id) |>
              summarise(sex = first(as.character(sex)), age = first(age),
                        disease2 = first(as.character(disease2)), .groups = "drop"), by = "donor_id") |>
  mutate(cohort = "KPMP kidney",
         type = recode(type, "EC-PTC" = "peritubular", "EC-GC" = "glomerular", "EC-AEA" = "arteriolar"))

# ---- 2. Long table: one row per sample x gene ----
to_long <- function(m, s) {
  cp <- t(t(m) / colSums(m)) * 1e6
  g  <- intersect(panel, rownames(cp))
  as_tibble(t(cp[g, , drop = FALSE])) |> mutate(id = colnames(m)) |>
    pivot_longer(-id, names_to = "gene", values_to = "cpm") |>
    left_join(s, by = "id") |> mutate(logcpm = log2(cpm + 1))
}
long <- bind_rows(to_long(m1, s1), to_long(m2, s2), to_long(m3, s3)) |>
  mutate(sex = relevel(factor(sex), ref = "male"), disease2 = relevel(factor(disease2), ref = "healthy"),
         age50 = age >= 50)
saveRDS(long, "p2_panel_long.rds")
long |> distinct(cohort, type, id, sex) |> count(cohort, type, sex) |>
  pivot_wider(names_from = sex, values_from = n) |> print(n = Inf)

# ---- 3. Women vs men for every panel gene ----
fit_sex <- function(d, kidney, y = "logcpm") {        # kidney = TRUE adds CKD as a covariate
  f <- as.formula(paste(y, if (kidney) "~ sex + age + disease2" else "~ sex + age"))
  co <- summary(lm(f, data = d))$coefficients
  tibble(estimate = co["sexfemale", 1], se = co["sexfemale", 2], p = co["sexfemale", 4])
}
res_p2 <- long |> group_by(cohort, type, gene) |> filter(n_distinct(sex) == 2) |>
  group_modify(~ fit_sex(.x, .y$cohort == "KPMP kidney")) |> ungroup() |>
  group_by(cohort) |> mutate(FDR = p.adjust(p, "BH")) |> ungroup()
saveRDS(res_p2, "p2_sex_panel_results.rds")
show_tab <- function(r, genes) r |> filter(gene %in% genes) |>
  mutate(val = paste0(sprintf("%+.2f", estimate), ifelse(FDR < 0.05, "*", ""))) |>
  select(cohort, type, gene, val) |> pivot_wider(names_from = gene, values_from = val) |> print(n = Inf, width = Inf)
show_tab(res_p2, c("XIST", genes_x[-1]))                 # X-escape genes (XIST = sex-label check)
show_tab(res_p2, genes_hr)                               # hormone receptors

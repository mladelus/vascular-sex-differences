# =============================================================================
# LOSS OF CHROMOSOME Y (Figure 7)
# Per-nucleus Y-gene counts (8 broadly expressed Y genes; TXLNGY removed: duplicated feature name).
# A nucleus can lack Y reads by chance when shallow, so the estimate uses only nuclei with
# >= 4.6 expected Y reads (P(zero | Y present) < 1%). Women = negative control.
# Result: kidney 8/493 deep nuclei Y-negative in 20 men (1.6%; upper 95% 3.2%); women 98.8% Y-negative.
#         Heart nuclei too shallow for this estimate. No association with CKD; weak age trend.
# Makes: p2_loy_counts.rds (allY), p2_loy_estimates.rds
# =============================================================================
ygenes <- c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY", "USP9Y", "NLGN4Y")
ctrl   <- c("XIST", "ACTB", "MALAT1")

get_y <- function(ids) {                         # ids = soma_joinids (numeric, sorted, non-empty)
  ids <- sort(as.numeric(ids))
  seu <- get_seurat(census, organism = "Homo sapiens", obs_coords = ids,
                    obs_column_names = "soma_joinid", var_index = "feature_name",
                    var_value_filter = sprintf("feature_name %%in%% c(%s)",
                                               paste0("'", c(ygenes, ctrl), "'", collapse = ",")))
  cnt <- LayerData(seu, assay = "RNA", layer = "counts")
  g1 <- function(g) if (g %in% rownames(cnt)) as.numeric(cnt[g, ]) else 0
  out <- tibble(soma_joinid = seu$soma_joinid,
                Y_umi = as.numeric(Matrix::colSums(cnt[intersect(ygenes, rownames(cnt)), , drop = FALSE])),
                XIST_umi = g1("XIST"), ACTB_umi = g1("ACTB"))
  rm(seu, cnt); gc(verbose = FALSE); out
}

if (!file.exists("p2_loy_counts.rds")) {
  lv <- readRDS("rep_heart_lv_ec_obs.rds")
  nuc <- bind_rows(
    hgm |> filter(dbl_class == "singlet", type != "mural") |>
      transmute(cohort = "Heart Cell Atlas", soma_joinid, donor_id = as.character(donor_id), sex = as.character(sex),
                age = parse_age(development_stage), type, disease2 = "healthy", raw_sum),
    lv |> transmute(cohort = "Independent LV", soma_joinid, donor_id = as.character(donor_id), sex = as.character(sex),
                    age, type, disease2 = "healthy", raw_sum),
    km |> transmute(cohort = "KPMP kidney", soma_joinid, donor_id, sex = as.character(sex), age,
                    type = subclass.l2, disease2 = as.character(disease2), raw_sum))
  dir.create("p2_loy", showWarnings = FALSE)
  for (d in unique(paste(nuc$cohort, nuc$donor_id, sep = "__"))) {
    out <- file.path("p2_loy", paste0(gsub("[^A-Za-z0-9_-]", "_", d), ".rds")); if (file.exists(out)) next
    message(format(Sys.time(), "%H:%M"), "  ", d)
    ids <- nuc$soma_joinid[paste(nuc$cohort, nuc$donor_id, sep = "__") == d]
    saveRDS(get_y(ids), out)
  }
  allY <- nuc |> inner_join(bind_rows(lapply(list.files("p2_loy", full.names = TRUE), readRDS)), by = "soma_joinid")
  saveRDS(allY, "p2_loy_counts.rds")
}
allY <- readRDS("p2_loy_counts.rds")
stopifnot(all(c("cohort", "donor_id", "sex", "raw_sum", "Y_umi") %in% names(allY)))   # columns used below
# Older saved files use short cohort names and lack age / disease2: standardise and add them from donor metadata
allY <- allY |> mutate(cohort = dplyr::recode(as.character(cohort), HCA = "Heart Cell Atlas",
                                              LV_rep = "Independent LV", Kidney = "KPMP kidney"),
                       donor_id = as.character(donor_id), sex = as.character(sex))
miss <- setdiff(c("age", "disease2"), names(allY))
if (length(miss) > 0) {
  if (!exists("hgm")) hgm <- readRDS("hg_ec_nuclei_meta.rds")
  if (!"disease2" %in% names(km))
    km <- km |> mutate(disease2 = relevel(factor(ifelse(disease == "normal", "healthy", "CKD")), ref = "healthy"))
  don <- bind_rows(
    hgm |> transmute(cohort = "Heart Cell Atlas", donor_id = as.character(donor_id),
                     age = as.numeric(parse_age(development_stage)), disease2 = "healthy"),
    readRDS("rep_heart_lv_ec_obs.rds") |> transmute(cohort = "Independent LV", donor_id = as.character(donor_id),
                                                   age = as.numeric(age), disease2 = "healthy"),
    km |> transmute(cohort = "KPMP kidney", donor_id = as.character(donor_id),
                    age = as.numeric(age), disease2 = as.character(disease2))) |>
    distinct(cohort, donor_id, .keep_all = TRUE)
  allY <- allY |> left_join(don |> select(cohort, donor_id, all_of(miss)), by = c("cohort", "donor_id"))
  message("Added ", paste(miss, collapse = ", "), " to saved loss-of-Y counts; donors unmatched: ",
          n_distinct(allY$donor_id[is.na(allY[[miss[1]]])]))
}

# ---- 1. Y-negative nuclei by sequencing depth (most "Y-negative" nuclei are just shallow) ----
allY <- allY |> mutate(Yneg = Y_umi == 0,
                       depth = cut(raw_sum, c(0, 1000, 2000, 3000, 5000, 10000, Inf), dig.lab = 6,
                                   labels = c("<1k", "1-2k", "2-3k", "3-5k", "5-10k", ">10k")))
loy_depth <- allY |> group_by(cohort, sex, depth) |>
  summarise(nuclei = n(), pct_Yneg = 100 * mean(Yneg), .groups = "drop")
print(loy_depth, n = Inf)
saveRDS(allY |> group_by(cohort, donor_id, sex, age, disease2) |>
          summarise(nuclei = n(), pct_Yneg = 100 * mean(Yneg), median_umi = median(raw_sum), .groups = "drop"),
        "p2_loy_by_donor.rds")

# ---- 2. Deep nuclei only: expected Y reads >= 4.6 (P(0) < 1% if Y present) ----
rate <- allY |> filter(sex == "male") |> group_by(cohort) |> summarise(Y_per_umi = sum(Y_umi) / sum(raw_sum))
deep <- allY |> left_join(rate, by = "cohort") |> mutate(expected_Y = Y_per_umi * raw_sum) |> filter(expected_Y >= 4.6)
loy_est <- deep |> group_by(cohort, sex) |>
  summarise(donors = n_distinct(donor_id), nuclei = n(), Yneg = sum(Yneg), .groups = "drop") |>
  rowwise() |> mutate(pct = 100 * Yneg / nuclei,
                      lo = 100 * binom.test(Yneg, nuclei)$conf.int[1],
                      hi = 100 * binom.test(Yneg, nuclei)$conf.int[2]) |> ungroup()
print(loy_est, width = Inf)
saveRDS(list(rate = rate, loy_depth = loy_depth, loy_est = loy_est), "p2_loy_estimates.rds")

# ---- 3. In men: Y-negative nuclei vs age and CKD (adjusted for depth) ----
suppressPackageStartupMessages(library(lme4))
if (!exists("glmer_ctrl")) glmer_ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
or_table <- function(fit, pat) {                 # odds ratios (Wald 95% CI) for terms matching pat
  co <- summary(fit)$coefficients; ci <- confint(fit, method = "Wald", parm = "beta_")
  k <- grepl(pat, rownames(co))
  data.frame(term = rownames(co)[k], OR = round(exp(co[k, 1]), 2),
             conf.low = round(exp(ci[rownames(co)[k], 1]), 2), conf.high = round(exp(ci[rownames(co)[k], 2]), 2),
             p.value = signif(co[k, 4], 2), row.names = NULL)
}
for (co in unique(allY$cohort)) {
  d <- allY |> filter(cohort == co, sex == "male") |> mutate(age10 = age / 10, ld = log10(raw_sum))
  f <- if (co == "KPMP kidney") Yneg ~ age10 + disease2 + ld + (1 | donor_id) else Yneg ~ age10 + ld + (1 | donor_id)
  cat("\n==", co, "(men) ==\n")
  print(or_table(glmer(f, family = binomial, data = d, control = glmer_ctrl), "age10|disease2"))
}

# =============================================================================
# INDEPENDENT HEART REPLICATION (single-nucleus, healthy left ventricle)
# Cohort: "Heart" dataset, collection "Automatic cell-type harmonization and integration across
# Human Cell Atlas datasets" (dataset 364bd0c7-...). 29 donors with LV (15 women, 14 men);
# donor IDs and exact ages differ from the Heart Cell Atlas -> independent. Chemistry 10x 5' v1
# (TWCM donors) or 10x 3' v2 (four-digit donors) = also marks study -> adjusted as 'assay'.
# =============================================================================

# ---- Healthy adult hearts in the Census (metadata only) ----
ds <- as.data.frame(census$get("census_info")$get("datasets")$read()$concat())
hobs <- census$get("census_data")$get("homo_sapiens")$obs$read(
  value_filter = "tissue_general == 'heart' & disease == 'normal' & is_primary_data == TRUE",
  column_names = c("soma_joinid", "dataset_id", "donor_id", "sex", "development_stage", "tissue",
                   "cell_type", "assay", "suspension_type", "raw_sum"))$concat() |> as.data.frame() |> as_tibble() |>
  left_join(ds |> select(dataset_id, dataset_title, collection_name), by = "dataset_id") |>
  mutate(age = parse_age(development_stage),
         is_ec = str_detect(cell_type, "endotheli"), is_cap = str_detect(cell_type, "capillar"))
hobs |> filter(is_ec, age >= 18) |>
  group_by(collection_name, dataset_title, suspension_type) |>
  summarise(donors = n_distinct(donor_id), women = n_distinct(donor_id[sex == "female"]),
            men = n_distinct(donor_id[sex == "male"]), ECs = n(), capillary_ECs = sum(is_cap), .groups = "drop") |>
  arrange(desc(donors)) |> print(n = Inf, width = Inf)

# ---- The candidate cohort, donor by donor ----
cand <- hobs |> filter(dataset_title == "Heart", suspension_type == "nucleus",
                       collection_name == "Automatic cell-type harmonization and integration across Human Cell Atlas datasets")
cand |> group_by(donor_id, sex, age) |>
  summarise(regions = paste(sort(unique(as.character(tissue))), collapse = "; "),
            assay = paste(unique(as.character(assay)), collapse = "; "), capillary = sum(is_cap), .groups = "drop") |>
  arrange(sex, age) |> print(n = Inf, width = Inf)

# ---- Download LV endothelial nuclei, pseudobulk per donor x vessel type ----
lv <- cand |>
  filter(tissue == "heart left ventricle",
         cell_type %in% c("capillary endothelial cell", "endothelial cell of artery", "vein endothelial cell")) |>
  mutate(type = case_when(cell_type == "capillary endothelial cell" ~ "capillary",
                          cell_type == "endothelial cell of artery" ~ "arterial", TRUE ~ "venous"))
saveRDS(lv, "rep_heart_lv_ec_obs.rds")
dir.create("rep_heart", showWarnings = FALSE)
for (d in unique(as.character(lv$donor_id))) {
  out <- file.path("rep_heart", paste0(d, ".rds")); if (file.exists(out)) next
  message(format(Sys.time(), "%H:%M"), "  ", d)
  seu <- get_seurat(census, organism = "Homo sapiens", obs_coords = lv$soma_joinid[lv$donor_id == d],
                    obs_column_names = "soma_joinid")
  counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts"))
  meta <- lv[match(seu$soma_joinid, lv$soma_joinid), ]; rm(seu); gc(verbose = FALSE)
  pb <- lapply(split(seq_len(ncol(counts)), meta$type), function(i) Matrix::rowSums(counts[, i, drop = FALSE]))
  saveRDS(list(pb = pb, n = table(meta$type), donor = d), out)
  rm(counts, pb); gc(verbose = FALSE)
}

# ---- Does the Heart Cell Atlas capillary result replicate? ----
# Result (29 donors; 15 F / 14 M): capillary hypoxia D 0.11, TNF U 0.83, TGF-beta U 0.39; EMT DOWN in women 0.002;
# UPR, mTORC1, MYC, glycolysis DOWN ~0.04-0.05. Same without donor 1221 (88 capillary nuclei, outlier).
# 55-gene check: 36/50 same direction (binomial p ~0.003); HCA-up set Up (FDR 0.018); effects much smaller;
# HIF genes not replicated; genome-wide correlation of sex effects 0.002.
lv  <- readRDS("rep_heart_lv_ec_obs.rds")
rep <- lapply(list.files("rep_heart", full.names = TRUE), readRDS)
dinfo <- lv |> distinct(donor_id, sex, age, assay) |> mutate(donor_id = as.character(donor_id))

run_rep <- function(tp, min_n = 30) {
  rr <- rep[sapply(rep, function(r) tp %in% names(r$n) && r$n[[tp]] >= min_n)]
  g <- Reduce(intersect, lapply(rr, function(r) names(r$pb[[tp]])))
  m <- sapply(rr, function(r) r$pb[[tp]][g]); colnames(m) <- sapply(rr, `[[`, "donor")
  s <- tibble(donor_id = colnames(m)) |> left_join(dinfo, by = "donor_id") |>
    mutate(sex = relevel(factor(as.character(sex)), ref = "male"), assay = factor(as.character(assay)))
  lcs <- edgeR::cpm(DGEList(m), log = TRUE)
  s$stress <- colMeans(lcs[intersect(STRESS, rownames(lcs)), ])          # same stress covariate idea as HCA
  yy <- DGEList(m[rownames(m) %in% auto, ]); de <- model.matrix(~ sex + age + assay + stress, data = s)
  yy <- calcNormFactors(yy[filterByExpr(yy, de), , keep.lib.sizes = FALSE])
  ff <- eBayes(lmFit(voom(yy, de), de))
  cat(tp, ":", sum(s$sex == "female"), "women,", sum(s$sex == "male"), "men;", nrow(yy), "genes tested\n")
  list(s = s, y = yy, t = ff$t[, "sexfemale"], tt = topTable(ff, coef = "sexfemale", n = Inf))
}
rc <- run_rep("capillary"); ra <- run_rep("arterial"); rv <- run_rep("venous")
hca <- readRDS(if (file.exists("hg_capillary_sex_DE_pathways.rds")) "hg_capillary_sex_DE_pathways.rds" else
                "hg_endothelium_sex_DE_pathways.rds")

# (1) Pathways: Heart Cell Atlas vs replication (U = higher in women)
top8 <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
          "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_MYC_TARGETS_V1",
          "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_GLYCOLYSIS")
bind_rows(list(HCA_capillary = hca$hallmark, rep_capillary = cam2(rc$t),
               rep_arterial = cam2(ra$t), rep_venous = cam2(rv$t)), .id = "set") |>
  filter(pathway %in% top8) |>
  mutate(val = paste(substr(Direction, 1, 1), signif(FDR, 2)), pathway = sub("HALLMARK_", "", pathway)) |>
  select(set, pathway, val) |> pivot_wider(names_from = set, values_from = val) |> print(width = Inf)
print(head(cam2(rc$t), 10), width = Inf)            # what the new cohort itself ranks highest

# (2) The 55 Heart Cell Atlas genes (FDR < 0.05): same direction in the new cohort?
h55 <- hca$tt |> tibble::rownames_to_column("gene") |> as_tibble() |>
  filter(adj.P.Val < 0.05, !gene %in% c("HLA-DRB5", "NPPA"))          # genotype / ambient genes excluded
cmp <- h55 |> select(gene, logFC_HCA = logFC) |>
  inner_join(rc$tt |> tibble::rownames_to_column("gene") |> as_tibble() |>
               select(gene, logFC_rep = logFC, P_rep = P.Value), by = "gene")
same <- sign(cmp$logFC_HCA) == sign(cmp$logFC_rep)
cat("testable:", nrow(cmp), "| same direction:", sum(same), "| same direction and p < 0.05:",
    sum(same & cmp$P_rep < 0.05), "| binomial p:", signif(binom.test(sum(same), nrow(cmp))$p.value, 2), "\n")
cmp |> mutate(across(where(is.numeric), ~ signif(.x, 2))) |> arrange(P_rep) |> print(n = Inf)
cameraPR(rc$t, ids2indices(list(HCA_up_in_women = h55$gene[h55$logFC > 0],
                                HCA_down_in_women = h55$gene[h55$logFC < 0]), names(rc$t)))
common <- intersect(names(rc$t), rownames(hca$tt))
cat("genome-wide correlation of sex effects (Spearman):", round(cor(rc$t[common], hca$tt[common, "t"], method = "spearman"), 3), "\n")

# (3) Per donor: hypoxia, TGF-beta, TNF scores in capillaries
zr <- t(scale(t(edgeR::cpm(rc$y, log = TRUE))))
scr <- sapply(c(hypoxia = "HALLMARK_HYPOXIA", tgfb = "HALLMARK_TGF_BETA_SIGNALING", tnf = "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
              function(p) colMeans(zr[intersect(Hs[[p]], rownames(zr)), ]))
bind_cols(rc$s, as_tibble(round(scr, 2))) |> select(donor_id, sex, age, assay, hypoxia, tgfb, tnf) |>
  arrange(sex, age) |> print(n = Inf)
saveRDS(list(capillary = rc$tt, arterial = ra$tt, venous = rv$tt, cmp = cmp), "rep_heart_sex_DE.rds")

# Sensitivity: without donor 1221 (low nuclei count, outlying scores)
rep_all <- rep; rep <- Filter(function(r) r$donor != "1221", rep_all); rc2 <- run_rep("capillary"); rep <- rep_all
print(head(cam2(rc2$t), 10), width = Inf)

# ---- Mural-cell composition, both heart cohorts ----
# Result: no replicated difference (pericytes per capillary combined fold 0.95, p 0.66;
# SMC per arterial EC combined 0.78, p 0.13). HCA ratios differ strongly by D/H series.
comp_rep <- cand |> filter(tissue == "heart left ventricle") |>
  mutate(grp = case_when(cell_type == "capillary endothelial cell" ~ "cap",
                         cell_type == "endothelial cell of artery" ~ "art_ec",
                         cell_type == "vein endothelial cell" ~ "vein_ec",
                         cell_type == "endothelial cell" ~ "ec_unlabelled",
                         cell_type == "pericyte" ~ "peri",
                         cell_type %in% c("smooth muscle cell", "vascular associated smooth muscle cell") ~ "smc",
                         str_detect(cell_type, "cardiac myocyte|cardiac muscle") ~ "cm",
                         cell_type == "fibroblast" ~ "fib", TRUE ~ "other")) |>
  count(donor_id, sex, age, assay, grp) |> pivot_wider(names_from = grp, values_from = n, values_fill = 0) |>
  mutate(sex = relevel(factor(as.character(sex)), ref = "male"),
         peri_per_cap = log((peri + 1) / (cap + 1)), smc_per_artEC = log((smc + 1) / (art_ec + 1)))
vt <- readRDS("hca_vascular_nuclei_typed.rds")
comp_hca <- vt |> filter(tissue %in% c("heart left ventricle", "apex of heart", "interventricular septum")) |>
  mutate(grp = case_when(str_detect(cell_type, "pericyte") ~ "peri", str_detect(cell_type, "smooth") ~ "smc", TRUE ~ type)) |>
  count(donor_id, sex, development_stage, grp) |> pivot_wider(names_from = grp, values_from = n, values_fill = 0) |>
  mutate(age = parse_age(development_stage), sex = relevel(factor(as.character(sex)), ref = "male"),
         peri_per_cap = log((peri + 1) / (capillary + 1)), smc_per_artEC = log((smc + 1) / (arterial + 1)))
fitc <- function(d, y, f) { co <- summary(lm(as.formula(paste(y, f)), data = d))$coefficients["sexfemale", ]; co[1:2] }
for (y in c("peri_per_cap", "smc_per_artEC")) {
  a <- fitc(comp_rep, y, "~ sex + age + assay"); b <- fitc(comp_hca, y, "~ sex + age")
  w <- 1 / c(a[2], b[2])^2; est <- sum(w * c(a[1], b[1])) / sum(w); se <- sqrt(1 / sum(w))
  cat(y, ": new cohort fold", round(exp(a[1]), 2), "| HCA fold", round(exp(b[1]), 2),
      "| combined fold", round(exp(est), 2), " p =", signif(2 * pnorm(-abs(est / se)), 2), "\n")
}

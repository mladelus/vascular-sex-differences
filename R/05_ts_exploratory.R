# =============================================================================
# TABULA SAPIENS: EXPLORATORY ANALYSES
# Tabula Sapiens exploratory analyses (supplementary / not headline)
# Needs: 00_setup.R, 04 outputs (ts_analysis_dat.rds)
# A. genome-wide sex differences + pathways   B. EndMT
# C. muscle-type genes in capillaries (by organ; heart myocardium; doublet check)
# D. postmenopausal women vs age-matched men (4 vs 4)
# logFC / estimates > 0 = higher in women.
# =============================================================================
dat <- readRDS("ts_analysis_dat.rds")
nc <- filter(dat, vessel_class == "large vessel"); cap <- filter(dat, vessel_class == "microvascular")

# ============================ A. Genome-wide ================================
if (!file.exists("pseudobulk_donor_organ_class.rds")) {
  keep_ids <- dat |> select(soma_joinid, donor_id, sex, age, organ, vessel_class)
  pb_list <- list(); pb_info <- list()
  for (d in unique(as.character(keep_ids$donor_id))) {
    seu <- readRDS(file.path("ec_by_donor", paste0(d, ".rds")))
    ids <- seu$soma_joinid; counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts")); rm(seu)
    m <- keep_ids[match(ids, keep_ids$soma_joinid), ]; ok <- !is.na(m$soma_joinid)
    grp <- paste(m$donor_id, m$organ, m$vessel_class, sep = "|")
    for (g in unique(grp[ok])) {
      idx <- which(ok & grp == g)
      if (length(idx) >= 20) {
        pb_list[[g]] <- Matrix::rowSums(counts[, idx, drop = FALSE])
        pb_info[[g]] <- m[idx[1], ] |> mutate(n_cells = length(idx), id = g)
      }
    }
    rm(counts); gc(verbose = FALSE)
  }
  saveRDS(list(counts = do.call(cbind, pb_list), info = bind_rows(pb_info)), "pseudobulk_donor_organ_class.rds")
}
pbd <- readRDS("pseudobulk_donor_organ_class.rds"); pb <- pbd$counts
info <- pbd$info |> mutate(sex = relevel(factor(as.character(sex)), ref = "male"))

# Mean dissociation-stress score per pseudobulk sample (for stress-adjusted DE)
info <- info |> select(-any_of("stress")) |>
  left_join(dat |> mutate(id = paste(donor_id, organ, vessel_class, sep = "|")) |>
              group_by(id) |> summarise(stress = mean(Stress_UCell), .groups = "drop"), by = "id")

run_de <- function(cls, adjust_stress = FALSE) {
  i <- info$vessel_class == cls
  inf <- info[i, ] |> mutate(age10 = age / 10, organ = droplevels(factor(organ)))
  design <- model.matrix(as.formula(paste("~ sex + age10 + organ", if (adjust_stress) "+ stress")), data = inf)
  design <- design[, !colnames(design) %in% nonEstimable(design), drop = FALSE]
  y <- DGEList(pb[, i]); keep <- filterByExpr(y, design)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  v <- voom(y, design); cor <- duplicateCorrelation(v, design, block = inf$donor_id)$consensus
  v <- voom(y, design, block = inf$donor_id, correlation = cor)
  fit <- eBayes(lmFit(v, design, block = inf$donor_id, correlation = cor))
  topTable(fit, coef = "sexfemale", number = Inf) |> tibble::rownames_to_column("gene")
}

# Chromosome of each gene (EnsDb v86)
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
filter <- dplyr::filter; select <- dplyr::select          # re-assert after loading ensembldb
gchr <- ensembldb::genes(EnsDb.Hsapiens.v86, return.type = "data.frame", columns = c("gene_id", "seq_name"))
chr_map <- gene_map |> left_join(gchr, by = c("feature_id" = "gene_id")) |>
  group_by(gene = feature_name) |> summarise(chr = first(na.omit(seq_name)), .groups = "drop")
annot <- function(t) t |> left_join(chr_map, by = "gene") |>
  mutate(chr_class = case_when(chr == "X" ~ "X", chr == "Y" ~ "Y", chr %in% c(1:22) ~ "autosomal", TRUE ~ "unknown"))

de  <- lapply(c(microvascular = "microvascular", `large vessel` = "large vessel"), \(x) annot(run_de(x)))
saveRDS(de, "sex_DE_by_vessel_class_annotated.rds")

for (cls in names(de)) {
  t <- de[[cls]]
  cat("\n######", cls, "######\n")
  print(t |> filter(adj.P.Val < 0.05) |> count(chr_class))
  print(t |> filter(gene %in% c("XIST", "RPS4Y1", "DDX3Y", "KDM5D")) |> select(gene, logFC, adj.P.Val))  # sanity check
  print(t |> filter(chr_class == "autosomal") |> head(25) |> select(gene, chr, logFC, AveExpr, adj.P.Val))
}

# Pathways (Hallmark), autosomal genes only; stress genes removed. Up = higher in women.
hdf <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
hl <- split(hdf$gene_symbol, hdf$gs_name)
pathways <- function(t) {
  t <- t |> filter(chr_class == "autosomal", !gene %in% STRESS)
  cameraPR(setNames(t$t, t$gene), ids2indices(hl, t$gene)) |> tibble::rownames_to_column("pathway") |>
    mutate(pathway = sub("HALLMARK_", "", pathway))
}
for (cls in names(de)) { cat("\n##", cls, "\n"); print(pathways(de[[cls]]) |> filter(FDR < 0.1), n = 50) }

# NOT YET RUN: stress-adjusted DE and pathways (tests whether the activation signal is dissociation stress)
# de_adj <- lapply(c(microvascular = "microvascular", `large vessel` = "large vessel"), \(x) annot(run_de(x, TRUE)))
# saveRDS(de_adj, "sex_DE_stress_adjusted.rds")
# for (cls in names(de_adj)) print(pathways(de_adj[[cls]]) |> filter(FDR < 0.1), n = 50)

# ================================ B. EndMT ==================================
if (!"EndMT" %in% names(dat)) {
  gs <- list(EndMT = c("ACTA2","TAGLN","CNN1","COL1A1","COL1A2","COL3A1","FN1",
                       "SNAI1","SNAI2","TWIST1","ZEB2","S100A4","POSTN"),
             EC_id = c("PECAM1","CDH5","VWF","TIE1","ERG"))
  res <- list()
  for (d in unique(as.character(dat$donor_id))) {
    seu <- readRDS(file.path("ec_by_donor", paste0(d, ".rds")))
    seu <- seu[, seu$soma_joinid %in% dat$soma_joinid]
    counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts"))
    sc <- ScoreSignatures_UCell(counts, features = gs, ncores = 1)
    res[[d]] <- tibble(soma_joinid = seu$soma_joinid, EndMT = sc[, "EndMT_UCell"], EC_id = sc[, "EC_id_UCell"],
                       ASPN_umi = as.numeric(counts["ASPN", ]))
    rm(seu, counts); gc(verbose = FALSE)
  }
  dat <- dat |> left_join(bind_rows(res), by = "soma_joinid"); saveRDS(dat, "ts_analysis_dat.rds")
  nc <- filter(dat, vessel_class == "large vessel"); cap <- filter(dat, vessel_class == "microvascular")
}
cor(dat$EndMT, dat$EC_id, method = "spearman")                     # ~ -0.06: no loss of EC identity
for (vc in levels(dat$vessel_class)) {
  d <- filter(dat, vessel_class == vc)
  print(st(lmerTest::lmer(EndMT ~ sex * age10 + logdepth + chem + Stress_UCell + (1 | donor_id) + (1 | organ),
                          data = d), "sex|age10"))
}

# ===================== C. Muscle-type genes in capillaries ====================
if (!file.exists("capillary_muscle_panels.rds")) {
  res3 <- list()
  for (d in unique(as.character(cap$donor_id))) {
    seu <- readRDS(file.path("ec_by_donor", paste0(d, ".rds")))
    seu <- seu[, seu$soma_joinid %in% cap$soma_joinid]; if (ncol(seu) == 0) next
    counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts"))
    cp <- function(g) Matrix::colSums(counts[intersect(g, rownames(counts)), , drop = FALSE])
    res3[[d]] <- tibble(soma_joinid = seu$soma_joinid, sm_umi = cp(pan$sm), peri_umi = cp(pan$peri), str_umi = cp(pan$str))
    rm(seu, counts); gc(verbose = FALSE)
  }
  capm <- cap |> left_join(bind_rows(res3), by = "soma_joinid") |>
    mutate(sm_cp10k = sm_umi / raw_sum * 1e4, peri_cp10k = peri_umi / raw_sum * 1e4, str_cp10k = str_umi / raw_sum * 1e4)
  saveRDS(capm, "capillary_muscle_panels.rds")
}
capm <- readRDS("capillary_muscle_panels.rds")

# Within each organ (organs with >= 2 donors of each sex): one row per donor, descriptive
don <- capm |> count(organ, donor_id, sex) |> filter(n >= 20)
keep_org <- don |> count(organ, sex) |> filter(n >= 2) |> count(organ) |> filter(n == 2) |> pull(organ) |> as.character()
org_don <- capm |> filter(as.character(organ) %in% keep_org) |>
  group_by(organ, donor_id, sex, age, postmeno) |>
  summarise(cells = n(), smooth_muscle = mean(sm_cp10k), pericyte = mean(peri_cp10k), striated = mean(str_cp10k),
            .groups = "drop") |> filter(cells >= 20)
print(org_don, n = Inf)
p <- org_don |> pivot_longer(c(smooth_muscle, pericyte, striated), names_to = "program", values_to = "cp10k") |>
  ggplot(aes(sex, cp10k, colour = sex, shape = postmeno)) +
  geom_point(size = 3, position = position_jitter(width = 0.1, height = 0)) +
  facet_grid(program ~ organ, scales = "free_y") +
  labs(y = "Mean counts per 10,000 (capillary cells)", x = NULL, shape = "Age") + theme_bw()
ggsave("capillary_muscle_by_organ.pdf", p, width = 3 + 2 * length(keep_org), height = 7)

# Heart myocardium only, with region and doublet check (total-UMI ratio)
hc <- capm |> filter(heart_region == "myocardium")
hc |> group_by(donor_id, sex, age, tissue) |>
  summarise(cells = n(), pericyte = round(mean(peri_cp10k), 2),
            pct_pericyte_ge3 = round(100 * mean(peri_umi >= 3), 1), .groups = "drop") |> arrange(sex, age)
hc |> mutate(peri_group = case_when(peri_umi == 0 ~ "0", peri_umi < 3 ~ "1-2", TRUE ~ ">=3")) |>
  group_by(donor_id, sex, age, peri_group) |>
  summarise(cells = n(), median_total_UMIs = median(raw_sum), EC_identity = round(mean(EC_id), 3), .groups = "drop") |>
  group_by(donor_id) |> mutate(ratio_vs_neg = round(median_total_UMIs / median_total_UMIs[peri_group == "0"], 2)) |>
  ungroup() |> arrange(sex, age, peri_group) |> print(n = Inf)
# (Not replicated in the Heart Cell Atlas nuclei — see 06.)

# ============== D. Postmenopausal women vs age-matched men (4 vs 4) ============
pm <- dat |> filter(donor_id %in% c("TSP27", "TSP1", "TSP2", "TSP7", "TSP14", "TSP17", "TSP25", "TSP6")) |> droplevels()
pm |> distinct(donor_id, sex, age) |> group_by(sex) |> summarise(n = n(), mean_age = mean(age))
for (vc in levels(pm$vessel_class)) {
  d <- filter(pm, vessel_class == vc); cat("\n##", vc, "\n")
  for (y in c("EndMT", "SenMayo_UCell", "Inflammatory_UCell", "TNFA_UCell")) {
    r <- st(lmerTest::lmer(as.formula(paste(y, "~ sex + logdepth + chem + Stress_UCell + (1 | donor_id) + (1 | organ)")),
                           data = d), "^sexfemale$"); cat(y, ": "); print(r[, -1], row.names = FALSE)
  }
}

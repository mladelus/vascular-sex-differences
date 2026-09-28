# =============================================================================
# GTEx: HYPOXIA / TGF-BETA BY SEX, ADJUSTED FOR ISCHEMIC TIME
# Bulk RNA-seq, open access (GTEx v8): heart left ventricle, atrial appendage, coronary artery.
# Question: are hypoxia, TGF-beta, TNF/NF-kB and matrix programs higher in women AFTER adjusting
# for age, ischemic time (SMTSISCH), manner of death (Hardy scale, DTHHRDY) and RNA quality (SMRIN)?
# Caveat: bulk tissue = mostly cardiomyocytes/smooth muscle, so this tests the whole tissue.
# =============================================================================

# ---- Download GTEx files (about 65 MB; skips files already saved) ----
# If a download fails: gtexportal.org -> Downloads -> Adult GTEx -> Bulk tissue expression (v8,
# "gene read counts by tissue") and Metadata (v8 SampleAttributesDS, SubjectPhenotypesDS);
# save them into ts2_project/gtex/ with the same file names.
library(data.table)
dir.create("gtex", showWarnings = FALSE); options(timeout = 3600)
base <- "https://storage.googleapis.com/adult-gtex"
gtex_urls <- c(
  attr  = paste0(base, "/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"),
  pheno = paste0(base, "/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"),
  lv    = paste0(base, "/bulk-gex/v8/rna-seq/counts-by-tissue/gene_reads_2017-06-05_v8_heart_left_ventricle.gct.gz"),
  aa    = paste0(base, "/bulk-gex/v8/rna-seq/counts-by-tissue/gene_reads_2017-06-05_v8_heart_atrial_appendage.gct.gz"),
  cor   = paste0(base, "/bulk-gex/v8/rna-seq/counts-by-tissue/gene_reads_2017-06-05_v8_artery_coronary.gct.gz"))
for (k in names(gtex_urls)) {
  f <- file.path("gtex", basename(gtex_urls[k]))
  if (!file.exists(f)) { message("downloading ", basename(f)); download.file(gtex_urls[k], f, mode = "wb") }
}
data.frame(file = list.files("gtex"), MB = round(file.size(list.files("gtex", full.names = TRUE)) / 1e6, 1))

# ---- Women vs men per tissue, adjusted for procurement ----
# Result: LV 137 F / 294 M; atrium 136 F / 293 M; coronary 94 F / 146 M.
# Women vs men: hypoxia not higher anywhere (FDR >= 0.55); HIF genes not higher (STC1 lower);
# TGF-beta lower in atrium (0.032); TNF/NF-kB LOWER in women: coronary 4e-8, atrium 0.033.
# Ischaemic time: strongly alters mTORC1, glycolysis, UPR, EMT (FDR down to 2e-5).
gattr <- fread("gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")[, .(SAMPID, SMTSISCH, SMRIN)]
gph   <- fread("gtex/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
gph[, `:=`(sex = ifelse(SEX == 2, "female", "male"), age = as.numeric(substr(AGE, 1, 2)) + 5)]   # decade midpoint

run_gtex <- function(fname) {
  g  <- fread(file.path("gtex", fname), skip = 2)
  m  <- as.matrix(g[, grep("^GTEX-", names(g)), with = FALSE])
  m  <- rowsum(m, g$Description)                                   # Ensembl IDs -> gene symbols
  m  <- m[rownames(m) %in% auto, ]                                  # autosomal only
  s  <- data.table(SAMPID = colnames(m))
  s[, SUBJID := sub("^(GTEX-[^-]+)-.*", "\\1", SAMPID)]
  s  <- merge(s, gattr, by = "SAMPID", sort = FALSE)
  s  <- merge(s, gph[, .(SUBJID, sex, age, DTHHRDY)], by = "SUBJID", sort = FALSE)
  s  <- s[!is.na(SMTSISCH) & !is.na(DTHHRDY) & !is.na(SMRIN)]
  s[, `:=`(sex = relevel(factor(sex), ref = "male"), DTHHRDY = factor(DTHHRDY), isch_h = SMTSISCH / 60)]
  m  <- m[, s$SAMPID]
  yy <- DGEList(m); de <- model.matrix(~ sex + age + isch_h + DTHHRDY + SMRIN, data = s)
  yy <- calcNormFactors(yy[filterByExpr(yy, de), , keep.lib.sizes = FALSE])
  ff <- eBayes(lmFit(voom(yy, de), de))
  cat(fname, ": ", sum(s$sex == "female"), " women, ", sum(s$sex == "male"), " men\n", sep = "")
  list(s = s, y = yy, t = ff$t[, "sexfemale"], tt = topTable(ff, coef = "sexfemale", n = Inf),
       isch = ff$t[, "isch_h"])
}
gt <- list(LV = run_gtex(basename(gtex_urls["lv"])), atrium = run_gtex(basename(gtex_urls["aa"])),
           coronary = run_gtex(basename(gtex_urls["cor"])))

# Donors by sex, age, death type
lapply(gt, function(x) x$s[, .N, by = .(sex, age50 = age >= 50)][order(sex, age50)])
lapply(gt, function(x) x$s[, .N, by = .(sex, DTHHRDY)][order(sex, DTHHRDY)])

# (1) Women vs men, key pathways, adjusted for ischemic time + death type (U = higher in women)
key <- c("HALLMARK_HYPOXIA", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
         "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_GLYCOLYSIS",
         "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
bind_rows(lapply(gt, function(x) cam2(x$t)), .id = "tissue") |> filter(pathway %in% key) |>
  mutate(val = paste(substr(Direction, 1, 1), signif(FDR, 2)), pathway = sub("HALLMARK_", "", pathway)) |>
  select(tissue, pathway, val) |> pivot_wider(names_from = tissue, values_from = val) |> print(width = Inf)

# (2) Does ischemic time itself raise these programs? (tells us if procurement matters)
bind_rows(lapply(gt, function(x) cam2(x$isch)), .id = "tissue") |> filter(pathway %in% key) |>
  mutate(val = paste(substr(Direction, 1, 1), signif(FDR, 2)), pathway = sub("HALLMARK_", "", pathway)) |>
  select(tissue, pathway, val) |> pivot_wider(names_from = tissue, values_from = val) |> print(width = Inf)

# (3) The heart-capillary genes, women vs men in GTEx (logFC; positive = higher in women)
genes_hc <- c("ANKRD37", "PDK1", "SLC16A3", "ERO1A", "STC1", "AKAP12", "BHLHE40", "WT1",
              "SMAD7", "SMAD6", "SKIL", "ID2", "BMPR2", "SOX17")
sapply(gt, function(x) round(x$tt[genes_hc, "logFC"], 2)) |> `rownames<-`(genes_hc)
saveRDS(lapply(gt, function(x) list(s = x$s, tt = x$tt)), "gtex_sex_results.rds")

# ---- Coronary artery — lower TNF/NF-kB in women: age-dependent? immune cells? ----
# Result: under 50 women much lower (TNF 2e-5, IFN-g 2e-12, inflammatory 7e-7, IL6 6e-5); 50+ only TNF 0.012;
# sex x age50 significant (IFN-g 2e-6, TNF 0.015). Narrowing = older men lower; women flat.
# + immune: TNF survives (2e-5); IFN-g/inflammatory/IL6 attenuate. Top genes IL1B, CCL4, CCL5, IL7R, CD69.
# Links to Tabula Sapiens: lower inflammatory activity in women's non-capillary endothelium, under 50.
co <- gt$coronary
s15 <- copy(co$s); s15[, age50 := age >= 50]
lc <- edgeR::cpm(co$y, log = TRUE)
imm_genes <- intersect(c("PTPRC", "CD68", "CD14", "CD163", "LYZ", "CD3E"), rownames(lc))
s15[, immune := colMeans(lc[imm_genes, SAMPID, drop = FALSE])]            # immune-cell content per sample
fit_t <- function(keep, form, coef = "sexfemale") {
  d <- droplevels(s15[keep]); X <- model.matrix(form, data = d)
  f <- eBayes(lmFit(voom(co$y[, keep], X), X)); f$t[, coef]
}
key2 <- c(key, "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING")
f1 <- ~ sex + age + isch_h + DTHHRDY + SMRIN
f2 <- ~ sex + age + isch_h + DTHHRDY + SMRIN + immune
all_s <- rep(TRUE, nrow(s15))
tests <- list(all = cam2(fit_t(all_s, f1)), all_plus_immune = cam2(fit_t(all_s, f2)),
              under50 = cam2(fit_t(!s15$age50, f1)), age50plus = cam2(fit_t(s15$age50, f1)),
              sex_x_age50 = cam2(fit_t(all_s, ~ sex * age50 + isch_h + DTHHRDY + SMRIN, "sexfemale:age50TRUE")))
bind_rows(tests, .id = "analysis") |> filter(pathway %in% key2) |>
  mutate(val = paste(substr(Direction, 1, 1), signif(FDR, 2)), pathway = sub("HALLMARK_", "", pathway)) |>
  select(analysis, pathway, val) |> pivot_wider(names_from = analysis, values_from = val) |> print(width = Inf)

# Plain view: TNF/NF-kB score (mean z) by sex and age group
zc <- t(scale(t(lc)))
s15[, tnf := colMeans(zc[intersect(Hs$HALLMARK_TNFA_SIGNALING_VIA_NFKB, rownames(zc)), SAMPID])]
s15[, .(n = .N, mean_TNF = round(mean(tnf), 2), mean_immune = round(mean(immune), 2)), by = .(age50, sex)][order(age50, sex)]

# Genes behind lower TNF/NF-kB in women (most negative first)
g <- intersect(Hs$HALLMARK_TNFA_SIGNALING_VIA_NFKB, rownames(co$tt))
co$tt[g, ] |> tibble::rownames_to_column("gene") |> as_tibble() |> arrange(t) |>
  select(gene, logFC, adj.P.Val) |> mutate(across(where(is.numeric), ~ signif(.x, 2))) |> head(15)

# ---- Replicate in aorta + tibial artery; one-death-type check; endothelial genes ----
# Result: aorta 152 F / 279 M; tibial 206 F / 445 M. Coronary: lower inflammation in women survives
# ventilator-only (TNF 2e-8, IFN-g 3e-5) and immune adjustment (TNF 2e-5). Aorta: no difference.
# Tibial: lower in women (TNF/IFN-g 0.001) but not after immune adjustment or in ventilator-only;
# difference at 50+ (TNF 2e-8), not < 50. Coronary endothelial genes all slightly lower in women
# (SELE -0.19, VCAM1 -0.23, CCL2 -0.40, IL6 -0.30 log2FC).
more <- c(aorta = "gene_reads_2017-06-05_v8_artery_aorta.gct.gz", tibial = "gene_reads_2017-06-05_v8_artery_tibial.gct.gz")
for (f in more) if (!file.exists(file.path("gtex", f)))
  download.file(paste0(base, "/bulk-gex/v8/rna-seq/counts-by-tissue/", f), file.path("gtex", f), mode = "wb")
gt$aorta <- run_gtex(more[["aorta"]]); gt$tibial <- run_gtex(more[["tibial"]])

art_tests <- function(x) {
  s <- copy(x$s); s[, age50 := age >= 50]
  lcx <- edgeR::cpm(x$y, log = TRUE)
  s[, immune := colMeans(lcx[intersect(c("PTPRC", "CD68", "CD14", "CD163", "LYZ", "CD3E"), rownames(lcx)), SAMPID, drop = FALSE])]
  ft <- function(keep, form, coef = "sexfemale") {
    d <- droplevels(s[keep]); X <- model.matrix(form, data = d)
    eBayes(lmFit(voom(x$y[, keep], X), X))$t[, coef] }
  a <- rep(TRUE, nrow(s))
  list(all             = cam2(ft(a, ~ sex + age + isch_h + DTHHRDY + SMRIN)),
       plus_immune     = cam2(ft(a, ~ sex + age + isch_h + DTHHRDY + SMRIN + immune)),
       ventilator_only = cam2(ft(s$DTHHRDY == "0", ~ sex + age + isch_h + SMRIN)),
       under50         = cam2(ft(!s$age50, ~ sex + age + isch_h + DTHHRDY + SMRIN)),
       age50plus       = cam2(ft(s$age50, ~ sex + age + isch_h + DTHHRDY + SMRIN)),
       sex_x_age50     = cam2(ft(a, ~ sex * age50 + isch_h + DTHHRDY + SMRIN, "sexfemale:age50TRUE")))
}
infl <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
          "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING")
for (tis in c("coronary", "aorta", "tibial")) {
  x <- gt[[tis]]
  cat("\n==", tis, ":", sum(x$s$sex == "female"), "women,", sum(x$s$sex == "male"), "men; under 50:",
      sum(x$s$sex == "female" & x$s$age < 50), "women,", sum(x$s$sex == "male" & x$s$age < 50), "men ==\n")
  bind_rows(art_tests(x), .id = "analysis") |> filter(pathway %in% infl) |>
    mutate(val = paste(substr(Direction, 1, 1), signif(FDR, 2)), pathway = sub("HALLMARK_", "", pathway)) |>
    select(analysis, pathway, val) |> pivot_wider(names_from = analysis, values_from = val) |> print(width = Inf)
}

# ---- Stricter gene-set tests (self-contained fry) for the artery results ----
# cameraPR with fixed inter-gene correlation 0.01 can overstate significance, so the manuscript reports fry.
# Result (fry p): coronary all ages TNF 0.053, IFN-g 0.073; under 50 IFN-g 0.008 (FDR 0.03), inflammatory 0.026,
# IL-6 0.023, TNF 0.07; 50+ >= 0.37; sex x age >= 0.11. Aorta all >= 0.18. Tibial 50+ TNF 0.006 (sex x age 0.03).
sets <- Hs[c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
             "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING")]
art_fry <- function(x) {
  s <- data.table::copy(x$s); s[, age50 := age >= 50]
  lcx <- edgeR::cpm(x$y, log = TRUE)
  s[, immune := colMeans(lcx[intersect(c("PTPRC","CD68","CD14","CD163","LYZ","CD3E"), rownames(lcx)), SAMPID, drop = FALSE])]
  ft <- function(keep, form, coef = "sexfemale") {
    d <- droplevels(s[keep]); X <- model.matrix(form, data = d); v <- voom(x$y[, keep], X)
    r <- fry(v, ids2indices(sets, rownames(v)), design = X, contrast = which(colnames(X) == coef))
    dplyr::as_tibble(tibble::rownames_to_column(r, "pathway"))
  }
  a <- rep(TRUE, nrow(s))
  dplyr::bind_rows(list(
    all             = ft(a, ~ sex + age + isch_h + DTHHRDY + SMRIN),
    plus_immune     = ft(a, ~ sex + age + isch_h + DTHHRDY + SMRIN + immune),
    ventilator_only = ft(s$DTHHRDY == "0", ~ sex + age + isch_h + SMRIN),
    under50         = ft(!s$age50, ~ sex + age + isch_h + DTHHRDY + SMRIN),
    age50plus       = ft(s$age50, ~ sex + age + isch_h + DTHHRDY + SMRIN),
    sex_x_age50     = ft(a, ~ sex * age50 + isch_h + DTHHRDY + SMRIN, "sexfemale:age50TRUE")), .id = "analysis")
}
fry_tab <- dplyr::bind_rows(lapply(c("coronary", "aorta", "tibial"),
                                   function(t) dplyr::mutate(art_fry(gt[[t]]), artery = t)))
fry_tab |>
  dplyr::mutate(val = paste(substr(Direction, 1, 1), signif(PValue, 2)), pathway = sub("HALLMARK_", "", pathway)) |>
  dplyr::select(artery, analysis, pathway, val) |>
  tidyr::pivot_wider(names_from = analysis, values_from = val) |> print(n = Inf, width = Inf)   # Table 4
saveRDS(fry_tab, "gtex_arteries_fry.rds")

# ---- Competitive camera test with the inter-gene correlation ESTIMATED from the data ----
# Result (camera p): coronary all ages TNF 0.081, IFN-g 0.10, inflammatory 0.20, IL-6 0.21; under 50 IFN-g 0.035
# (others 0.08-0.17); 50+ >= 0.33. Aorta all >= 0.49. Tibial 50+ TNF 0.056.
# Reported in the Results next to fry (coronary, all ages). inter.gene.cor = NA makes camera estimate the
# correlation; limma's default is a fixed 0.01. Same models as art_fry().
art_camera <- function(x) {
  s <- data.table::copy(x$s); s[, age50 := age >= 50]
  lcx <- edgeR::cpm(x$y, log = TRUE)
  s[, immune := colMeans(lcx[intersect(c("PTPRC","CD68","CD14","CD163","LYZ","CD3E"), rownames(lcx)), SAMPID, drop = FALSE])]
  ft <- function(keep, form, coef = "sexfemale") {
    d <- droplevels(s[keep]); X <- model.matrix(form, data = d); v <- voom(x$y[, keep], X)
    r <- camera(v, ids2indices(sets, rownames(v)), design = X, contrast = which(colnames(X) == coef),
                inter.gene.cor = NA)
    dplyr::as_tibble(tibble::rownames_to_column(r, "pathway"))
  }
  a <- rep(TRUE, nrow(s))
  dplyr::bind_rows(list(
    all         = ft(a, ~ sex + age + isch_h + DTHHRDY + SMRIN),
    plus_immune = ft(a, ~ sex + age + isch_h + DTHHRDY + SMRIN + immune),
    under50     = ft(!s$age50, ~ sex + age + isch_h + DTHHRDY + SMRIN),
    age50plus   = ft(s$age50, ~ sex + age + isch_h + DTHHRDY + SMRIN)), .id = "analysis")
}
cam_tab <- dplyr::bind_rows(lapply(c("coronary", "aorta", "tibial"),
                                   function(t) dplyr::mutate(art_camera(gt[[t]]), artery = t)))
cam_tab |>
  dplyr::mutate(val = paste(substr(Direction, 1, 1), signif(PValue, 2)), pathway = sub("HALLMARK_", "", pathway)) |>
  dplyr::select(artery, analysis, pathway, val) |>
  tidyr::pivot_wider(names_from = analysis, values_from = val) |> print(n = Inf, width = Inf)
saveRDS(cam_tab, "gtex_arteries_camera.rds")

# How did under-50 donors die, by sex? (Hardy: 0 ventilator, 1 violent fast, 2 natural fast, 3 intermediate, 4 slow)
lapply(gt[c("coronary", "aorta", "tibial")], function(x) dcast(x$s[age < 50, .N, by = .(sex, DTHHRDY)], sex ~ DTHHRDY, value.var = "N", fill = 0))

# Endothelial activation / protection genes, women vs men (logFC, all ages; positive = higher in women)
ecg <- c("SELE", "SELP", "VCAM1", "ICAM1", "CX3CL1", "CCL2", "IL6", "NFKBIA", "KLF2", "KLF4", "NOS3")
sapply(gt[c("coronary", "aorta", "tibial")], function(x) round(x$tt[ecg, "logFC"], 2)) |> `rownames<-`(ecg)
saveRDS(lapply(gt, function(x) list(s = x$s, tt = x$tt)), "gtex_sex_results.rds")

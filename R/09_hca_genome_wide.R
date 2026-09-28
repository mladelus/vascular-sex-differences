# =============================================================================
# HEART ENDOTHELIUM: GENOME-WIDE SEX COMPARISON
# Question: without pre-chosen gene sets, what differs between women's and men's
# myocardial capillary endothelium (Heart Cell Atlas, 7 women / 7 men, doublet-free)?
# Pseudobulk per donor x region (>= 20 nuclei), limma-voom with donor blocking,
# autosomal genes only (X/Y genes are analyzed in scripts 12-15), pathways by cameraPR (Hallmark, Reactome).
# =============================================================================

# ---- Every gene, then pathways ----
# Result: 77 samples (35 F / 42 M); 4,631 genes tested; 55 at FDR < 0.05.
# Hallmark UP in women: TNFA/NF-kB (FDR 1e-6), hypoxia (7e-5), UPR (9e-5), TGF-beta (8e-4),
# mTORC1 (8e-4), MYC (0.003), EMT (0.004), glycolysis (0.02).
# Top genes up in women: HIF targets ANKRD37, PDK1, SLC16A3, ERO1A, STC1; WT1; SOX17, BCL6B; AKAP12.
# Ignore as biology: HLA-DRB5 (genotype), NPPA (cardiomyocyte ambient RNA).
library(msigdbr); library(EnsDb.Hsapiens.v86)
filter <- dplyr::filter; select <- dplyr::select; count <- dplyr::count; rename <- dplyr::rename
first  <- dplyr::first;  desc <- dplyr::desc;     slice <- dplyr::slice     # EnsDb masks these too

hg  <- lapply(list.files("hg_nuclei", full.names = TRUE), readRDS)
hgm <- readRDS("hg_ec_nuclei_meta.rds")
chr  <- ensembldb::genes(EnsDb.Hsapiens.v86, return.type = "data.frame", columns = c("gene_name", "seq_name"))
auto <- unique(chr$gene_name[chr$seq_name %in% as.character(1:22)])
filter <- dplyr::filter; select <- dplyr::select          # ensembldb masks these again

# One function for any vessel type: pseudobulk -> limma-voom (donor blocking) -> women vs men
run_type <- function(tp) {
  p <- unlist(lapply(hg, `[[`, "pb"), recursive = FALSE); p <- p[str_detect(names(p), paste0("\\|", tp, "$"))]
  g <- Reduce(intersect, lapply(p, names)); m <- sapply(p, function(v) v[g])
  s <- tibble(id = colnames(m)) |> separate(id, c("donor_id", "tissue", "type"), sep = "\\|", remove = FALSE) |>
    left_join(hgm |> filter(type == tp, dbl_class == "singlet") |>
                mutate(donor_id = as.character(donor_id), tissue = as.character(tissue)) |>
                group_by(donor_id, tissue) |>
                summarise(sex = first(as.character(sex)), age = parse_age(first(development_stage)),
                          assay = names(which.max(table(assay))), stress = mean(Stress_UCell), .groups = "drop"),
              by = c("donor_id", "tissue")) |>
    mutate(sex = relevel(factor(sex), ref = "male"), series = substr(donor_id, 1, 1))
  yy <- DGEList(m[rownames(m) %in% auto, ]); de <- model.matrix(~ sex + age + tissue + assay + stress, data = s)
  yy <- calcNormFactors(yy[filterByExpr(yy, de), , keep.lib.sizes = FALSE])
  vv <- voom(yy, de); dd <- duplicateCorrelation(vv, de, block = s$donor_id)
  vv <- voom(yy, de, block = s$donor_id, correlation = dd$consensus)
  dd <- duplicateCorrelation(vv, de, block = s$donor_id)
  ff <- eBayes(lmFit(vv, de, block = s$donor_id, correlation = dd$consensus))
  list(m = m, s = s, y = yy, cor = dd$consensus, t = ff$t[, "sexfemale"],
       tt = topTable(ff, coef = "sexfemale", n = Inf))
}
capr <- run_type("capillary"); art <- run_type("arterial"); ven <- run_type("venous")
cap <- capr$m; sm <- capr$s; y <- capr$y; tt <- capr$tt; stat <- capr$t
count(sm, sex)

# Sanity check: X/Y genes separate the sexes
xy <- c("XIST", "RPS4Y1", "DDX3Y", "UTY", "KDM5D")
round(t(t(cap[intersect(xy, rownames(cap)), ]) / colSums(cap) * 1e6)) |> t() |>
  as.data.frame() |> tibble::rownames_to_column("sample") |> mutate(sex = sm$sex) |> print()

cat("genes tested:", nrow(tt), " FDR < 0.05:", sum(tt$adj.P.Val < 0.05), " FDR < 0.10:", sum(tt$adj.P.Val < 0.1), "\n")
tt |> tibble::rownames_to_column("gene") |> as_tibble() |> select(gene, logFC, AveExpr, P.Value, adj.P.Val) |>
  mutate(across(where(is.numeric), ~ signif(.x, 3))) |> print(n = 55)      # positive logFC = higher in women

# Pathways (Direction "Up" = higher in women)
getS <- function(...) { d <- msigdbr(species = "Homo sapiens", ...); split(d$gene_symbol, d$gs_name) }
Hs <- getS(category = "H"); Rs <- getS(category = "C2", subcategory = "CP:REACTOME")   # deprecation warning: harmless
cam2 <- function(st, sets = Hs, cor = 0.01) {
  idx <- ids2indices(sets, names(st)); idx <- idx[lengths(idx) >= 10 & lengths(idx) <= 500]
  cameraPR(st, idx, inter.gene.cor = cor) |> tibble::rownames_to_column("pathway") |> as_tibble() }
H <- cam2(stat); R <- cam2(stat, Rs)
print(head(H, 15), width = Inf); print(head(R, 25), width = Inf)
saveRDS(list(tt = tt, hallmark = H, reactome = R, arterial = art$tt, venous = ven$tt),
        "hg_endothelium_sex_DE_pathways.rds")

# ---- Is the signal capillary-specific? Same test in arterial + venous ECs ----
# Result: hypoxia UP in all three (FDR cap 7e-5, ven 0.026, art 0.076); mTORC1 in all three (<= 0.024);
# EMT cap 0.004 + arterial 0.012. HIF genes (STC1, NDRG1, P4HA1, VEGFA, ERO1A) up in women in all types.
# -> hypoxia looks tissue-wide (biology of women's myocardium OR procurement), not capillary-specific.
bind_rows(list(capillary = cam2(stat), arterial = cam2(art$t), venous = cam2(ven$t)), .id = "cells") |>
  filter(pathway %in% H$pathway[1:10]) |>
  mutate(val = paste(substr(Direction, 1, 1), signif(FDR, 2)), pathway = sub("HALLMARK_", "", pathway)) |>
  select(cells, pathway, val) |> pivot_wider(names_from = cells, values_from = val) |> print(width = Inf)
hif <- c("ANKRD37", "PDK1", "SLC16A3", "ERO1A", "STC1", "ADM", "VEGFA", "BNIP3", "EGLN3", "NDRG1", "P4HA1", "WT1")
sapply(list(capillary = tt, arterial = art$tt, venous = ven$tt), function(x) round(x[hif, "logFC"], 2)) |>
  `rownames<-`(hif)

# ---- Are women's capillaries more arterial-like? ----
# Result: SOX17, BCL6B up and EMCN down, but DLL4/EFNB2 not up and CA4/RGCC unchanged -> no arterial shift.
zon <- c("SOX17", "GJA5", "HEY1", "EFNB2", "DLL4", "BCL6B", "CA4", "RGCC", "EMCN", "ACKR1", "NR2F2", "SELP")
tt[intersect(zon, rownames(tt)), c("logFC", "P.Value", "adj.P.Val")] |> signif(2)

# ---- Contamination (should NOT be higher in women) ----
# Result: DCN, ACTA2 flat; cardiomyocyte genes slightly lower; LUM/COL1A1/2/PDGFRA/FAP too low to test.
chk <- c("DCN", "LUM", "COL1A1", "COL1A2", "PDGFRA", "FAP", "RGS5", "PDGFRB", "ACTA2", "TTN", "MYH7", "TNNT2")
tt[intersect(chk, rownames(tt)), c("logFC", "P.Value", "adj.P.Val")] |> signif(3)

# ---- Genes driving each top pathway ----
# Result: hypoxia = HIF targets (AKAP12, ERO1A, PDK1, STC1, BHLHE40); TNF = immediate-early/flow genes;
# TGF-beta = SMAD7, SMAD6, SKIL, SMURF2, KLF10, ID2, BMPR2; "EMT" = integrins/basement membrane, not EndMT.
for (p in c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA",
            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_TGF_BETA_SIGNALING")) {
  cat("\n", p, "\n")
  g <- intersect(Hs[[p]], rownames(tt))
  print(tt[g, ] |> tibble::rownames_to_column("gene") |> as_tibble() |> arrange(desc(t)) |>
          select(gene, logFC, t, P.Value) |> mutate(across(where(is.numeric), ~ signif(.x, 2))) |> head(10))
}

# ---- Sensitivity — no stress genes / + D-H series / stricter correlation ----
# Result: all hold without stress genes and with series (balanced 4/4, 3/3); stricter test: TNF 0.046,
# UPR 0.057, hypoxia + TGF-beta 0.072, EMT 0.21.
ieg <- unique(c(STRESS, grep("^HSP|^DNAJ|^FOS|^JUN|^EGR|^IER|^NR4A|^ATF3$|^KLF[246]$|^DUSP|^ZFP36|^PPP1R15A$|^GADD45",
                             names(stat), value = TRUE)))
sm |> distinct(donor_id, sex, series) |> count(series, sex)
d2 <- model.matrix(~ sex + age + tissue + assay + stress + series, data = sm)
v2 <- voom(y, d2, block = sm$donor_id, correlation = capr$cor)
fit2 <- eBayes(lmFit(v2, d2, block = sm$donor_id, correlation = capr$cor))
stat2 <- fit2$t[, "sexfemale"]
comp <- list(main = cam2(stat), no_stress_genes = cam2(stat[!names(stat) %in% ieg]),
             plus_series = cam2(stat2), stricter = cam2(stat, cor = 0.05))
bind_rows(comp, .id = "analysis") |> filter(pathway %in% H$pathway[1:10]) |>
  mutate(val = paste(substr(Direction, 1, 1), signif(FDR, 2)), pathway = sub("HALLMARK_", "", pathway)) |>
  select(analysis, pathway, val) |> pivot_wider(names_from = analysis, values_from = val) |> print(width = Inf)

# ---- Per-donor pathway scores (most women, or one or two donors?) ----
# Result: hypoxia/TNF high in D1, D11, D4 (women) and D6 (man); TGF-beta higher in 5 of 7 women.
zs <- t(scale(t(edgeR::cpm(y, log = TRUE))))      # edgeR:: because cpm is masked
pw <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_TGF_BETA_SIGNALING")
psc <- sapply(pw, function(p) colMeans(zs[intersect(Hs[[p]], rownames(zs)), , drop = FALSE]))
as_tibble(psc) |> mutate(donor_id = sm$donor_id, sex = sm$sex, age = sm$age) |>
  group_by(donor_id, sex, age) |> summarise(across(starts_with("HALLMARK"), ~ round(mean(.x), 2)), .groups = "drop") |>
  rename_with(~ sub("HALLMARK_", "", .x)) |> arrange(sex, age) |> print(n = Inf, width = Inf)

# =============================================================================
# FIGURES 1-8, SUPPLEMENTARY FIGURES S1-S3, SESSION INFO
# Run R/00_setup.R first. Then run the "figure setup" block, then any figure block, in any order.
# Each block lists the objects it needs. Objects already in your session are used as they are;
# anything missing is rebuilt from the saved .rds files, or the block tells you which script to run.
# Every figure is one multi-panel file (panels A, B) saved to fig_dir (set below) as
#   PDF (vector, for the journal) + PNG (300 dpi, for Word) + TIFF (600 dpi, LZW; journal upload).
# Colors are the same in every figure: women = purple, men = blue; cohorts have fixed colours.
# =============================================================================

# ---- Figure setup (run once per session) ----
for (pk in c("patchwork", "scales", "ggrepel")) if (!requireNamespace(pk, quietly = TRUE)) install.packages(pk)
library(ggplot2); library(patchwork)
filter <- dplyr::filter; select <- dplyr::select; count <- dplyr::count
first  <- dplyr::first;  rename <- dplyr::rename; desc <- dplyr::desc

col_sex    <- c(female = "#8E24AA", male = "#1E88E5")                 # women purple, men blue
lab_sex    <- c(female = "Women", male = "Men")
col_cohort <- c("Heart Cell Atlas" = "#FB8C00", "Independent LV" = "#43A047", "KPMP kidney" = "#E53935")
col_div    <- c(low = "#1E88E5", mid = "#F5F5F5", high = "#8E24AA")    # blue = lower in women, purple = higher

theme_pub <- function(base = 9) {
  theme_classic(base_size = base) +
    theme(plot.margin = margin(3, 4, 3, 3),
          plot.tag = element_text(face = "bold", size = base + 4),
          plot.title = element_text(face = "bold", size = base + 1, margin = margin(b = 2)),
          strip.background = element_rect(fill = "#EDE7F6", colour = NA),
          strip.text = element_text(face = "bold", size = base, margin = margin(2, 2, 2, 2)),
          legend.position = "top", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -4, 0),
          legend.key.size = unit(3.5, "mm"), legend.title = element_text(size = base - 1),
          legend.text = element_text(size = base - 1),
          axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3))
}
theme_set(theme_pub())

# Where figures are saved: <data folder>/figures (00_setup.R sets the working directory to the data folder).
fig_dir <- file.path(getwd(), "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

pdf_device <- function(filename, width, height, ...) {    # vector PDF that keeps Greek letters
  if (Sys.info()[["sysname"]] == "Darwin")                # Mac: Quartz (no XQuartz needed)
    grDevices::quartz(file = filename, type = "pdf", width = width, height = height, bg = "white")
  else grDevices::cairo_pdf(filename = filename, width = width, height = height)
}
save_fig <- function(p, name, w, h, tag = TRUE) { # w, h in inches; BMC full page 170 mm = 6.69 in, max height 225 mm incl. legend
  if (tag) p <- p + plot_annotation(tag_levels = "A")   # panel letters (not for single-panel supplementary figures)
  f <- function(ext) file.path(fig_dir, paste0(name, ext))
  ggsave(f(".pdf"),  p, width = w, height = h, device = pdf_device)
  ggsave(f(".png"),  p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(f(".tiff"), p, width = w, height = h, dpi = 600, bg = "white", compression = "lzw")
  message("saved ", f(""), " (.pdf, .png, .tiff)")
  invisible(p)
}
score_set <- function(y, genes) {                  # mean z-score of a gene set per sample
  z <- t(scale(t(edgeR::cpm(y, log = TRUE)))); colMeans(z[intersect(genes, rownames(z)), , drop = FALSE])
}
if (!exists("Hs")) {
  hdf <- msigdbr::msigdbr(species = "Homo sapiens", category = "H"); Hs <- split(hdf$gene_symbol, hdf$gs_name)
}
if (!exists("cam2")) cam2 <- function(st, sets = Hs, cor = 0.01) {
  idx <- ids2indices(sets, names(st)); idx <- idx[lengths(idx) >= 10 & lengths(idx) <= 500]
  cameraPR(st, idx, inter.gene.cor = cor) |> tibble::rownames_to_column("pathway") |> as_tibble() }
tstat  <- function(tt) setNames(tt$t, rownames(tt))                   # t-statistics from a saved topTable
pretty_hallmark <- function(x) str_to_sentence(gsub("_", " ", sub("HALLMARK_", "", x))) |>
  str_replace("Tnfa signaling via nfkb", "TNF-α/NF-κB") |> str_replace("Tgf beta signaling", "TGF-β") |>
  str_replace("Mtorc1 signaling", "mTORC1") |> str_replace("Myc targets v1", "MYC targets") |>
  str_replace("Interferon gamma response", "Interferon-γ response")

# =============================================================================
# FIGURE 1 — Study design
# =============================================================================
# Drawn outside R (figures/make_fig1.py in the repository). Use fig1_study_design.png / .tiff. Nothing to run.

# =============================================================================
# FIGURE 2 — Healthy myocardial capillaries do not differ by sex (Heart Cell Atlas)
# Needs: hg_validation_results.rds, hg_donor_region_summaries.rds (script 06)
# =============================================================================
f2res <- readRDS("hg_validation_results.rds") |> filter(cells == "capillary")
dr    <- readRDS("hg_donor_region_summaries.rds")
lab2  <- c(pericyte = "Pericyte transcripts", pericyte_ge3 = "Pericyte-positive nuclei",
           smooth_muscle = "Smooth-muscle transcripts", striated = "Striated-muscle transcripts",
           inflammatory = "Inflammatory response", tnfa = "TNF-α/NF-κB", senmayo = "SenMayo senescence",
           p16 = "p16-positive", p21 = "p21-high", oxphos = "Oxidative phosphorylation")
sds <- dr |> filter(type == "capillary") |> summarise(across(all_of(names(lab2)), ~ sd(.x, na.rm = TRUE))) |>
  pivot_longer(everything(), names_to = "outcome", values_to = "sd")
f2a <- f2res |> left_join(sds, by = "outcome") |>
  mutate(std = estimate / sd, lo = conf.low / sd, hi = conf.high / sd,
         outcome = factor(lab2[outcome], levels = rev(lab2)),
         ages = factor(ifelse(ages == "all", "All ages (7 W / 7 M)", "Age ≥ 50 (5 W / 7 M)"),
                       levels = c("All ages (7 W / 7 M)", "Age ≥ 50 (5 W / 7 M)")))
f2a <- f2a |> mutate(y = as.numeric(outcome) + ifelse(as.integer(ages) == 1, 0.17, -0.17))   # manual dodge
p2a <- ggplot(f2a, aes(std, y, colour = ages)) +
  annotate("rect", xmin = -0.5, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = "#F3E5F5", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.7) +
  geom_point(size = 2) +
  geom_text(aes(x = hi, label = paste0(" p=", signif(p, 2))), hjust = 0, size = 2.1, show.legend = FALSE) +
  scale_colour_manual(values = c("#5E35B1", "#EC407A"), name = NULL) +
  scale_y_continuous(breaks = seq_along(levels(f2a$outcome)), labels = levels(f2a$outcome),
                     expand = expansion(add = 0.5)) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  labs(x = "Women minus men (SD units, 95% CI)", y = NULL, title = "Capillary measures, women vs men")
f2b <- dr |> filter(type == "capillary") |>
  select(donor_id, sex, age50, inflammatory, tnfa, senmayo) |>
  pivot_longer(c(inflammatory, tnfa, senmayo), names_to = "score", values_to = "value") |>
  mutate(score = factor(lab2[score], levels = lab2[c("inflammatory", "tnfa", "senmayo")]),
         sex = factor(as.character(sex), levels = c("female", "male")),
         age = ifelse(age50, "≥ 50", "< 50"))
p2b <- ggplot(f2b, aes(sex, value, colour = sex, fill = sex)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.55, linewidth = 0.35, alpha = 0.25) +
  geom_point(aes(shape = age), position = position_jitter(width = 0.13, height = 0, seed = 1), size = 1.8, alpha = 0.85) +
  facet_wrap(~ score, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = col_sex, labels = lab_sex, guide = "none") +
  scale_fill_manual(values = col_sex, guide = "none") + scale_shape_manual(values = c(16, 17), name = "Age") +
  scale_x_discrete(labels = lab_sex) + labs(x = NULL, y = "Score (per donor × region)", title = "Per-donor scores")
save_fig(p2a / p2b + plot_layout(heights = c(1.35, 1)), "Figure2_capillary_null", 6.69, 7)

# =============================================================================
# FIGURE 3 — Single-cohort differences do not replicate
# Needs: hca (scripts 09/11), rc (script 11), gt$LV (script 10), cmp — or the saved .rds files
# =============================================================================
if (!exists("hca")) hca <- readRDS(if (file.exists("hg_capillary_sex_DE_pathways.rds"))
  "hg_capillary_sex_DE_pathways.rds" else "hg_endothelium_sex_DE_pathways.rds")
rep_t <- if (exists("rc")) rc$t else tstat(readRDS("rep_heart_sex_DE.rds")$capillary)
lv_t  <- if (exists("gt") && !is.null(gt$LV$t)) gt$LV$t else tstat(readRDS("gtex_sex_results.rds")$LV$tt)
cmp   <- if (exists("cmp")) cmp else readRDS("rep_heart_sex_DE.rds")$cmp
top8 <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
          "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_MYC_TARGETS_V1",
          "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_GLYCOLYSIS")
f3a <- bind_rows(list(`Heart Cell\nAtlas\n7 W / 7 M` = hca$hallmark, `Independent\nLV\n15 W / 14 M` = cam2(rep_t),
                      `GTEx\nLV\n137 W / 294 M` = cam2(lv_t)), .id = "cohort") |>
  filter(pathway %in% top8) |>
  mutate(signed = pmax(pmin(ifelse(Direction == "Up", 1, -1) * -log10(FDR), 6), -6),
         pathway = factor(pretty_hallmark(pathway), levels = rev(pretty_hallmark(top8))),
         cohort = factor(cohort, levels = unique(cohort)),
         lab = ifelse(FDR < 0.001, formatC(FDR, format = "e", digits = 0), formatC(FDR, format = "fg", digits = 2)))
p3a <- ggplot(f3a, aes(cohort, pathway, fill = signed)) +
  geom_tile(colour = "white", linewidth = 0.8) + geom_text(aes(label = lab), size = 2.4) +
  scale_fill_gradient2(low = col_div[["low"]], mid = col_div[["mid"]], high = col_div[["high"]], limits = c(-6, 6),
                       name = "Signed −log10 FDR (+ higher in women)") +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Hallmark programs") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", legend.key.width = unit(8, "mm"))
lim3 <- max(abs(c(cmp$logFC_HCA, cmp$logFC_rep)), na.rm = TRUE)
p3b <- ggplot(cmp, aes(logFC_HCA, logFC_rep)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "#F3E5F5") +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "#E3F2FD") +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_abline(slope = 1, linetype = 2, colour = "grey40", linewidth = 0.4) +
  geom_point(aes(fill = P_rep < 0.05), shape = 21, size = 2.2, colour = "grey20", stroke = 0.3) +
  ggrepel::geom_text_repel(data = ~ dplyr::filter(.x, P_rep < 0.05), aes(label = gene), size = 2.4,
                           fontface = "italic", nudge_x = -0.6, nudge_y = 0.6, min.segment.length = 0,
                           segment.size = 0.3, box.padding = 0.4) +
  scale_fill_manual(values = c(`TRUE` = "#FDD835", `FALSE` = "grey75"), labels = c("p ≥ 0.05", "p < 0.05"),
                    name = "Independent cohort") +
  coord_equal(xlim = c(-lim3, lim3), ylim = c(-lim3, lim3)) +
  labs(x = "Heart Cell Atlas (log2FC, women − men)", y = "Independent LV (log2FC, women − men)",
       title = "The 50 Heart Cell Atlas genes")
save_fig(p3a + p3b + plot_layout(widths = c(1, 1.05)), "Figure3_replication", 6.69, 3.9)

# =============================================================================
# FIGURE 4 — Lower inflammatory signaling in younger women's coronary arteries (GTEx)
# Needs: gt$coronary, gt$aorta, gt$tibial WITH $y (session, script 10). If missing, run script 10 (~5 min).
# =============================================================================
stopifnot("Run R/10_gtex_heart_arteries.R first (gt with $y needed)" = exists("gt") && all(c("coronary", "aorta", "tibial") %in% names(gt)))
arts <- c(coronary = "Coronary", aorta = "Aorta", tibial = "Tibial")
f4a <- bind_rows(lapply(names(arts), function(tis) {
  x <- gt[[tis]]
  tibble(artery = arts[[tis]], sex = factor(as.character(x$s$sex), levels = c("female", "male")),
         age = ifelse(x$s$age >= 50, "≥ 50", "< 50"),
         TNF = score_set(x$y, Hs$HALLMARK_INTERFERON_GAMMA_RESPONSE))   # column kept as TNF; now interferon-g score
})) |> mutate(age = factor(age, levels = c("< 50", "≥ 50")), artery = factor(artery, levels = arts))
f4n <- f4a |> count(artery, age, sex)
p4a <- ggplot(f4a, aes(age, TNF, fill = sex, colour = sex)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75, seed = 1), size = 0.5, alpha = 0.35) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.35, linewidth = 0.4, position = position_dodge(width = 0.75)) +
  geom_text(data = f4n, aes(y = -Inf, label = n), position = position_dodge(width = 0.75), vjust = -0.4,
            size = 2.1, show.legend = FALSE) +
  facet_wrap(~ artery, nrow = 1) +
  scale_fill_manual(values = col_sex, labels = lab_sex, name = NULL) +
  scale_colour_manual(values = col_sex, labels = lab_sex, name = NULL) +
  labs(x = "Age (years)", y = "Interferon-γ response score (mean z)", title = "Interferon-γ response by sex and age")
ecg <- c("SELE", "SELP", "VCAM1", "ICAM1", "CX3CL1", "CCL2", "IL6", "NFKBIA", "KLF2", "KLF4", "NOS3")
f4b <- bind_rows(lapply(names(arts), function(tis) {
  tt <- gt[[tis]]$tt; g <- intersect(ecg, rownames(tt))
  tibble(artery = arts[[tis]], gene = g, logFC = tt[g, "logFC"], FDR = tt[g, "adj.P.Val"])
})) |> mutate(artery = factor(artery, levels = arts), gene = factor(gene, levels = rev(ecg)),
              star = ifelse(FDR < 0.05, "*", ""))
p4b <- ggplot(f4b, aes(artery, gene, fill = logFC)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = paste0(sprintf("%+.2f", logFC), star)), size = 2.3) +
  scale_fill_gradient2(low = col_div[["low"]], mid = col_div[["mid"]], high = col_div[["high"]],
                       limits = c(-0.6, 0.6), oob = scales::squish, name = "log2FC\n(W − M)") +
  scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Endothelial genes") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "right", axis.text.x = element_text(face = "bold"))
save_fig(p4a + p4b + plot_layout(widths = c(2.4, 1)), "Figure4_coronary_inflammation", 6.69, 3.4)

# =============================================================================
# FIGURE 5 — Modest sex x CKD interaction in glomerular capillary interferon-g signaling (KPMP)
# Needs: kpb, km (scripts 08/12) — rebuilt from kp_ec/ if missing
# =============================================================================
if (!exists("km") || !exists("kpb")) {
  kp_all <- lapply(list.files("kp_ec", full.names = TRUE), readRDS)
  km  <- bind_rows(lapply(kp_all, `[[`, "meta")) |>
    mutate(donor_id = as.character(donor_id), sex = relevel(factor(as.character(sex)), ref = "male"),
           disease2 = relevel(factor(ifelse(disease == "normal", "healthy", "CKD")), ref = "healthy"))
  kpb <- unlist(lapply(kp_all, `[[`, "pb"), recursive = FALSE)
}
if (!"disease2" %in% names(km))                          # older sessions: add healthy/CKD label
  km <- km |> mutate(disease2 = relevel(factor(ifelse(disease == "normal", "healthy", "CKD")), ref = "healthy"))
pb_sub <- function(sub) {                               # count matrix + donor info for one EC subtype
  p <- kpb[str_detect(names(kpb), paste0("\\|", sub, "$"))]
  g <- Reduce(intersect, lapply(p, names)); m <- sapply(p, function(v) v[g])
  s <- tibble(donor_id = sub("\\|.*", "", colnames(m))) |>
    left_join(km |> filter(subclass.l2 == sub) |> group_by(donor_id) |>
                summarise(sex = first(sex), age = first(age), disease2 = first(disease2),
                          tub = mean(log1p(tub_umi / raw_sum * 1e4)), .groups = "drop"), by = "donor_id")
  list(m = m, s = s)
}
f5a <- bind_rows(lapply(c(`EC-PTC` = "Peritubular capillaries", `EC-GC` = "Glomerular capillaries"), function(lbl) {
  sub <- ifelse(lbl == "Peritubular capillaries", "EC-PTC", "EC-GC"); x <- pb_sub(sub)
  x$s |> mutate(IFNG = score_set(DGEList(x$m), Hs$HALLMARK_INTERFERON_GAMMA_RESPONSE), capillary = lbl)
})) |> mutate(sex = factor(as.character(sex), levels = c("female", "male")),
              disease2 = factor(as.character(disease2), levels = c("healthy", "CKD"), labels = c("Healthy", "CKD")),
              capillary = factor(capillary, levels = c("Peritubular capillaries", "Glomerular capillaries")))
# Check: does the sex x CKD interaction show up in the plotted score itself?
# Same covariates as the pathway model; positive interaction = CKD change larger in men.
chk5 <- f5a |> mutate(sex = relevel(sex, ref = "male")) |> group_by(capillary) |>
  group_modify(function(d, k) {
    co <- summary(lm(IFNG ~ sex * disease2 + age + tub, data = d))$coefficients
    tibble(term = rownames(co), estimate = round(co[, 1], 3), p = signif(co[, 4], 2))
  }) |> ungroup() |> filter(grepl("sex", term))
print(chk5)
# Unadjusted per-donor score = supplementary Figure S1
pS3 <- ggplot(f5a, aes(disease2, IFNG, colour = sex, fill = sex, group = sex)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.45, seed = 1), size = 1.6, alpha = 0.7) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.1, position = position_dodge(width = 0.45)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.12, linewidth = 0.6, position = position_dodge(width = 0.45)) +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 21, colour = "white", stroke = 0.8,
               position = position_dodge(width = 0.45)) +
  facet_wrap(~ capillary, nrow = 1) +
  scale_colour_manual(values = col_sex, labels = lab_sex, name = NULL) +
  scale_fill_manual(values = col_sex, labels = lab_sex, name = NULL) +
  labs(x = NULL, y = "Interferon-γ response score (mean z)")
save_fig(pS3, "FigureS1_kidney_IFNG_score", 6.69, 3, tag = FALSE)

# Panel A: gene-level sex x CKD interaction, interferon-g genes vs all other genes (+ self-contained fry p)
kid_fit <- function(sub) {
  x <- pb_sub(sub); X <- model.matrix(~ sex * disease2 + age + tub, data = x$s)
  yy <- DGEList(x$m); yy <- calcNormFactors(yy[filterByExpr(yy, X), , keep.lib.sizes = FALSE])
  v <- voom(yy, X); cf <- which(colnames(X) == "sexfemale:disease2CKD")
  ifn <- Hs$HALLMARK_INTERFERON_GAMMA_RESPONSE
  fr <- fry(v, ids2indices(list(IFNG = ifn), rownames(v)), design = X, contrast = cf)
  tt <- topTable(eBayes(lmFit(v, X)), coef = cf, n = Inf)
  list(tt = tt, genes = tibble(gene = rownames(tt), t = tt$t,
                               set = ifelse(rownames(tt) %in% ifn, "Interferon-γ genes", "All other genes")),
       p_fry = fr["IFNG", "PValue"])
}
kf <- list(`Peritubular capillaries` = kid_fit("EC-PTC"), `Glomerular capillaries` = kid_fit("EC-GC"))
f5A <- bind_rows(lapply(names(kf), function(n) mutate(kf[[n]]$genes, capillary = n))) |>
  mutate(capillary = factor(capillary, levels = names(kf)),
         set = factor(set, levels = c("All other genes", "Interferon-γ genes")))
f5p <- tibble(capillary = factor(names(kf), levels = names(kf)),
              lab = sapply(kf, function(k) paste0("fry p = ", signif(k$p_fry, 2))))
p5a <- ggplot(f5A, aes(set, t, fill = set)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey45", linewidth = 0.4) +
  geom_violin(colour = NA, alpha = 0.55, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.4, fill = "white") +
  geom_text(data = f5p, aes(x = 1.5, y = Inf, label = lab), inherit.aes = FALSE, vjust = 1.5, size = 2.6,
            fontface = "bold") +
  facet_wrap(~ capillary, nrow = 1) +
  scale_fill_manual(values = c("All other genes" = "grey70", "Interferon-γ genes" = "#E53935"), guide = "none") +
  coord_cartesian(ylim = c(-4.5, 4)) +
  labs(x = NULL, y = "Sex × CKD interaction (t, per gene)", title = "Sex × CKD: interferon-γ genes vs all genes")
tt5 <- kf$`Peritubular capillaries`$tt
f5b <- tt5[intersect(Hs$HALLMARK_INTERFERON_GAMMA_RESPONSE, rownames(tt5)), ] |>
  tibble::rownames_to_column("gene") |> as_tibble() |> arrange(t) |> slice(1:15) |>
  mutate(gene = factor(gene, levels = rev(gene)), star = ifelse(adj.P.Val < 0.05, "*", ""))
p5b <- ggplot(f5b, aes(logFC, gene, fill = logFC)) +
  geom_col(width = 0.75) + geom_text(aes(label = star, x = logFC), hjust = 1.3, size = 3) +
  scale_fill_gradient(low = "#1E88E5", high = "#90CAF9", guide = "none") +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(x = "Sex × CKD (log2FC)", y = NULL, title = "Driver genes")
save_fig(p5a + p5b + plot_layout(widths = c(2, 1)), "Figure5_kidney_CKD_interferon", 6.69, 3.5)

# =============================================================================
# FIGURE 6 — Sex-chromosome dosage of chromatin regulators (3 cohorts)
# Needs: dres, yshare (script 13) — or p2_xy_dosage.rds
# =============================================================================
ok_df <- function(x) exists(x, envir = globalenv()) && is.data.frame(get(x, envir = globalenv()))
if (!ok_df("dres") || !ok_df("yshare")) {
  xyd <- readRDS("p2_xy_dosage.rds")
  if (is.data.frame(xyd) || is.null(xyd$dres) || is.null(xyd$yshare))
    stop("p2_xy_dosage.rds is the OLD version. Run R/13_p2_xy_dosage.R first, then re-run this block.")
  dres <- xyd$dres; yshare <- xyd$yshare
}
pair_order <- c("KDM6A/UTY", "KDM5C/KDM5D", "USP9X/USP9Y", "DDX3X/DDX3Y", "EIF1AX/EIF1AY", "ZFX/ZFY",
                "RPS4X/RPS4Y1", "NLGN4X/NLGN4Y")
short_c <- c("Heart Cell Atlas" = "HCA", "Independent LV" = "LV", "KPMP kidney" = "Kidney")
f6a <- dres |> mutate(where = paste(short_c[cohort], type),
                      measure = factor(measure, c("X_only", "X_plus_Y"), c("X copy alone", "X + Y combined")),
                      pair = factor(pair, levels = rev(intersect(pair_order, unique(pair)))),
                      sig = ifelse(FDR < 0.05, "FDR < 0.05", "n.s."))
f6a$where <- factor(f6a$where, levels = unique(f6a$where[order(f6a$cohort, f6a$type)]))
p6a <- ggplot(f6a, aes(where, pair)) +
  geom_point(aes(size = pmin(-log10(FDR), 8), fill = estimate, shape = sig), colour = "grey25", stroke = 0.3) +
  scale_shape_manual(values = c("FDR < 0.05" = 21, "n.s." = 24), name = NULL) +
  scale_fill_gradient2(low = col_div[["low"]], mid = col_div[["mid"]], high = col_div[["high"]], limits = c(-1.5, 1.5),
                       oob = scales::squish, name = "W − M\n(log2 CPM)") +
  scale_size_continuous(range = c(1.2, 5), name = "−log10 FDR") +
  facet_wrap(~ measure, nrow = 1) + labs(x = NULL, y = NULL, title = "Women minus men, by cohort and vessel") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "right",
        panel.grid.major = element_line(colour = "grey94", linewidth = 0.3))
f6b <- yshare |> filter(pair %in% pair_order) |> mutate(pair = factor(pair, levels = pair_order))
p6b <- ggplot(f6b, aes(pair, 100 * Y_share)) +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey45", linewidth = 0.4) +
  geom_boxplot(outlier.shape = NA, fill = "#E3F2FD", colour = "#1565C0", width = 0.6, linewidth = 0.4) +
  geom_point(aes(colour = cohort), position = position_jitter(width = 0.12, seed = 1), size = 1.8) +
  scale_colour_manual(values = col_cohort, name = NULL) +
  labs(x = NULL, y = "Y copy share in men (%)", title = "Share supplied by the Y copy") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p6a / p6b + plot_layout(heights = c(1.45, 1)), "Figure6_XY_dosage", 6.69, 7)

# =============================================================================
# FIGURE 7 — Loss of chromosome Y is uncommon in men's kidney endothelium
# Needs: loy_depth, loy_est (script 14) — or p2_loy_estimates.rds
# =============================================================================
if (!exists("ok_df")) ok_df <- function(x) exists(x, envir = globalenv()) && is.data.frame(get(x, envir = globalenv()))
if (!ok_df("loy_depth") || !ok_df("loy_est")) {
  if (!file.exists("p2_loy_estimates.rds")) stop("p2_loy_estimates.rds not found. Run R/14_p2_loss_of_y.R first.")
  le <- readRDS("p2_loy_estimates.rds")
  if (is.data.frame(le) || is.null(le$loy_depth)) stop("p2_loy_estimates.rds is the OLD version. Run R/14_p2_loss_of_y.R first.")
  loy_depth <- le$loy_depth; loy_est <- le$loy_est
}
f7a <- loy_depth |> filter(nuclei >= 20) |> mutate(sex = factor(as.character(sex), levels = c("female", "male")))
p7a <- ggplot(f7a, aes(depth, pct_Yneg, colour = sex, group = sex)) +
  geom_line(linewidth = 0.9) + geom_point(aes(size = nuclei), alpha = 0.9) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_colour_manual(values = col_sex, labels = lab_sex, name = NULL) +
  scale_size_area(max_size = 4, name = "Nuclei", breaks = c(100, 1000, 10000), labels = scales::comma) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0.01, 0.03))) +
  labs(x = "UMIs per nucleus", y = "Nuclei without Y reads (%)", title = "Y-negative nuclei by sequencing depth") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
f7b <- loy_est |> mutate(sex = factor(as.character(sex), levels = c("female", "male"), labels = lab_sex),
                         label = paste0(Yneg, "/", nuclei, " (", donors, " donors)"))
p7b <- ggplot(f7b, aes(pct, cohort, colour = sex)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y", linewidth = 0.7) +
  geom_point(size = 2.8) +
  geom_text(aes(x = Inf, label = label), hjust = 1.05, vjust = -1.2, size = 2.3, show.legend = FALSE) +
  facet_wrap(~ sex, scales = "free_x") +
  scale_colour_manual(values = setNames(col_sex, lab_sex), guide = "none") +
  labs(x = "Y-negative among deep nuclei (%, exact 95% CI)", y = NULL,
       title = "Deep nuclei only (≥ 4.6 expected Y reads)")
save_fig(p7a / p7b + plot_layout(heights = c(1.3, 1)), "Figure7_loss_of_Y", 6.69, 5.6)

# =============================================================================
# FIGURE 8 — Sex-hormone receptors in vascular endothelium
# Needs: long (script 12) or p2_panel_long.rds; esr1_all (script 15) or p2_esr1_kidney.rds
# =============================================================================
if (!exists("ok_df")) ok_df <- function(x) exists(x, envir = globalenv()) && is.data.frame(get(x, envir = globalenv()))
if (!ok_df("long") || !"type" %in% names(long)) long <- readRDS("p2_panel_long.rds")
if (!"type" %in% names(long)) stop("p2_panel_long.rds is the OLD version. Run R/12_p2_sexchrom_receptor_panel.R first.")
if (!ok_df("esr1_all")) {
  if (!file.exists("p2_esr1_kidney.rds")) stop("p2_esr1_kidney.rds not found. Run R/15_p2_hormone_receptors.R first.")
  esr1_all <- readRDS("p2_esr1_kidney.rds")
}
recs <- c("AR", "ESR1", "ESR2", "PGR", "GPER1", "CYP19A1")
f8a <- long |> filter(gene %in% recs) |> group_by(cohort, type, gene) |>
  summarise(mean = mean(logcpm), .groups = "drop") |>
  mutate(where = paste(short_c[cohort], type), gene = factor(gene, levels = recs))
f8a$where <- factor(f8a$where, levels = rev(unique(f8a$where[order(f8a$cohort, f8a$type)])))
p8a <- ggplot(f8a, aes(gene, where, fill = mean)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f", mean), colour = mean < 5.5), size = 2.3, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_fill_viridis_c(option = "plasma", name = "Mean\nlog2 CPM") +
  scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Receptor expression") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, face = "bold"),
        legend.position = "right", legend.key.height = unit(5, "mm"))
f8b <- esr1_all |>
  mutate(sub = case_when(grepl("EHD3\\+", type) ~ "EHD3+", grepl("EHD3", type) ~ "EHD3\u2212",
                         TRUE ~ as.character(type)),
         lab = paste0(ifelse(data == "nuclei", "Nuclei, ", "Cells, "), sub,
                      ifelse(donors == "independent", " (new donors)", "")),
         lab = factor(lab, levels = rev(unique(lab))))
p8b <- ggplot(f8b, aes(estimate, lab, colour = data)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#E3F2FD", alpha = 0.6) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#F3E5F5", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45", linewidth = 0.4) +
  geom_errorbar(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se), width = 0.25, orientation = "y", linewidth = 0.7) +
  geom_point(size = 2.6) +
  geom_text(aes(label = paste0("p=", signif(p, 2))), vjust = -1, size = 2.2, show.legend = FALSE) +
  scale_colour_manual(values = c(nuclei = "#E53935", `single cells` = "#8E24AA"), name = NULL) +
  labs(x = "ESR1, women \u2212 men (log2 CPM, 95% CI)", y = NULL, title = "ESR1 in kidney endothelium") +
  theme(plot.title.position = "plot")
save_fig(p8a + p8b + plot_layout(widths = c(1, 1.1)), "Figure8_hormone_receptors", 6.69, 3.8)

# =============================================================================
# SUPPLEMENTARY FIGURES (Additional file 1). S1 (kidney interferon-g score) is saved in the Figure 5 block.
# =============================================================================
# S2 — Y-negative nuclei per man by age
lbd <- readRDS("p2_loy_by_donor.rds") |> filter(sex == "male")
pS1 <- ggplot(lbd, aes(age, pct_Yneg)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", fill = "grey85", linewidth = 0.6) +
  geom_point(aes(size = nuclei, colour = disease2), alpha = 0.85) +
  facet_wrap(~ cohort, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c(healthy = "#43A047", CKD = "#E53935"), name = NULL) +
  scale_size_area(max_size = 4, name = "Nuclei") +
  labs(x = "Age (years)", y = "Nuclei without Y reads (%, all depths)")
save_fig(pS1, "FigureS2_Ynegative_by_age_men", 6.69, 2.8, tag = FALSE)

# S3 — Receptors by age in women (no change; all FDR >= 0.81)
pS2 <- long |> filter(gene %in% c("AR", "ESR1", "PGR"), sex == "female") |>
  ggplot(aes(age, logcpm)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", fill = "grey85", linewidth = 0.5) +
  geom_point(aes(colour = cohort), size = 1.3, alpha = 0.8) +
  facet_grid(gene ~ cohort, scales = "free_y") +
  scale_colour_manual(values = col_cohort, guide = "none") +
  labs(x = "Age (years)", y = "log2 CPM (women)")
save_fig(pS2, "FigureS3_receptors_by_age_women", 6.69, 5, tag = FALSE)

# ---- Session info for the Methods ----
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
list.files(fig_dir)
message("Figures saved in ", normalizePath(fig_dir))


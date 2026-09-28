# =============================================================================
# OPTIONAL FUTURE WORK (not run; not in the manuscript) — clustering, PMVEC signature, pericyte coverage
if (FALSE) {   # optional: not run for the manuscript (downloads ~14 donors)
# ---- Download capillary nuclei counts for clustering ----
# Up to 4,000 doublet-free capillary nuclei per donor; one file per donor; skips saved donors.
set.seed(3)
capsel <- hgm |> filter(type == "capillary", dbl_class == "singlet") |>
  group_by(donor_id) |> slice_sample(n = 4000) |> ungroup()
count(capsel, sex, donor_id) |> print(n = Inf)
saveRDS(capsel, "hg_cap_selection.rds")
dir.create("hg_cap_counts", showWarnings = FALSE)
for (d in unique(as.character(capsel$donor_id))) {
  out <- file.path("hg_cap_counts", paste0(d, ".rds")); if (file.exists(out)) next
  message(format(Sys.time(), "%H:%M"), "  ", d)
  seu <- get_seurat(census, organism = "Homo sapiens", obs_coords = capsel$soma_joinid[capsel$donor_id == d],
                    obs_column_names = "soma_joinid")
  counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts"))
  colnames(counts) <- as.character(seu$soma_joinid)
  rm(seu); gc(verbose = FALSE)
  saveRDS(counts, out); rm(counts); gc(verbose = FALSE)
}
length(list.files("hg_cap_counts"))     # should be 14
}

# ---- Unsupervised clustering + which capillary states differ by sex ----
# Code will be added after 14d: normalise, autosomal variable genes, PCA, integrate by chemistry,
# cluster, UMAP, cluster markers (+ fibroblast contamination check), donor-level proportions by sex.

# ---- Your PMVEC signature in heart capillaries ----
# her_up <- c("FAP","LAMC2","COL1A2","SULF1","CREB3L1","HEY2","NFATC4","IL17D")
# his_up <- c("MYCN","PRKCZ","RARG","ZKSCAN3")
# Better: use the full GSE233547 DE table (female-up / male-up at padj < 0.05) as signatures.

# ---- Pericyte coverage / microvessel composition ----
# vt <- readRDS("hca_vascular_nuclei_typed.rds")  -> pericytes per capillary EC, SMC per arterial EC, by sex.

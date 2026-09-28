# =============================================================================
# TABULA SAPIENS: CELL SCORES AND VESSEL SUBTYPES
# Score each endothelial cell; assign vessel subtype
# Needs: 00_setup.R, 01_ts_data.R outputs
# Makes: cell_meta/<donor>.rds, ec_cells_scored.rds, subtype_reference.rds
# Methods facts produced here:
#   marker agreement with expert labels: capillary 82.4, venous 83.7, lymphatic 62.5, arterial 44.5 %
#   SingleR (leave-one-donor-out): venous 84.7, arterial 78.9, capillary 78.8, lymphatic 73.5 % (79.4 overall)
#   final: 24,393 expert-labelled, 35,166 consensus, 9,531 uncertain
# =============================================================================
keep_donors <- readRDS("ts_keep_donors.rds")

label_subtype <- function(ct) case_when(
  str_detect(ct, "lymph")                                    ~ "lymphatic",
  str_detect(ct, "arter")                                    ~ "arterial",
  str_detect(ct, "vein|venous|venule")                       ~ "venous",
  str_detect(ct, "capillar|microvascular|sinusoid|glomerul") ~ "capillary",
  TRUE                                                       ~ NA_character_)

# ---- 1. Per-donor scoring (UCell), marker-based subtype, senescence flags ----
process_donor <- function(d) {
  seu <- readRDS(file.path("ec_by_donor", paste0(d, ".rds")))
  meta <- seu[[]]
  counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts"))
  rm(seu); gc(verbose = FALSE)
  keep <- str_detect(meta$assay, "10x")
  meta <- meta[keep, ]; counts <- counts[, keep, drop = FALSE]

  sc <- ScoreSignatures_UCell(counts, features = ts_sets, ncores = 1, chunk.size = 1000)
  meta <- cbind(meta, as.data.frame(sc))

  mk <- as.matrix(meta[, c("Arterial_UCell", "Capillary_UCell", "Venous_UCell", "Lymphatic_UCell")])
  meta$marker_subtype <- c("arterial", "capillary", "venous", "lymphatic")[max.col(mk, ties.method = "first")]
  meta$marker_score   <- apply(mk, 1, max)
  meta$label_subtype  <- label_subtype(as.character(meta$cell_type))

  g <- function(x) if (x %in% rownames(counts)) as.numeric(counts[x, ]) else rep(0, ncol(counts))
  meta$CDKN2A_umi   <- g("CDKN2A")
  meta$MKI67_umi    <- g("MKI67")
  meta$CDKN1A_cp10k <- g("CDKN1A") / meta$raw_sum * 1e4
  meta$sen_p16      <- meta$CDKN2A_umi > 0 & meta$MKI67_umi == 0     # TS2.0 definition
  meta
}

dir.create("cell_meta", showWarnings = FALSE)
for (d in keep_donors) {
  out <- file.path("cell_meta", paste0(d, ".rds")); if (file.exists(out)) next
  message(format(Sys.time(), "%H:%M"), "  ", d)
  saveRDS(process_donor(d), out); gc(verbose = FALSE)
}
cells <- bind_rows(lapply(list.files("cell_meta", full.names = TRUE), readRDS))
nrow(cells)                                                                   # 69,090

# Agreement of marker calls with expert labels
lab <- cells |> filter(!is.na(label_subtype))
lab |> group_by(label_subtype) |>
  summarise(cells = n(), pct_agree = round(100 * mean(marker_subtype == label_subtype), 1))

# ---- 2. SingleR reference: expert-labelled cells, pseudobulk per donor x subtype ----
if (!file.exists("subtype_reference.rds")) {
  lab_ids <- lab |> select(soma_joinid, label_subtype)
  ref_list <- list()
  for (d in keep_donors) {
    seu <- readRDS(file.path("ec_by_donor", paste0(d, ".rds")))
    meta <- seu[[]]; counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts")); rm(seu)
    l <- lab_ids$label_subtype[match(meta$soma_joinid, lab_ids$soma_joinid)]
    keep <- !is.na(l) & str_detect(meta$assay, "10x")
    for (s in unique(l[keep])) {
      idx <- which(keep & l == s)
      if (length(idx) >= 20) ref_list[[paste(d, s, sep = "|")]] <- Matrix::rowSums(counts[, idx, drop = FALSE])
    }
    rm(counts); gc(verbose = FALSE)
  }
  ref <- list(counts = do.call(cbind, ref_list))
  ref$info <- tibble(id = colnames(ref$counts)) |> separate(id, c("donor", "subtype"), sep = "\\|", remove = FALSE)
  saveRDS(ref, "subtype_reference.rds")
}
ref <- readRDS("subtype_reference.rds")
ref$info |> count(subtype)                                                    # 6-7 donors per subtype
ref_log <- log2(sweep(ref$counts, 2, colSums(ref$counts), "/") * 1e6 + 1)

# ---- 3. SingleR classification, leaving each donor's own profiles out ----
if (file.exists("ec_cells_scored.rds")) {                     # reuse earlier results if present
  prev <- readRDS("ec_cells_scored.rds")
  if ("singler" %in% names(prev)) cells <- cells |> left_join(select(prev, soma_joinid, singler), by = "soma_joinid")
  rm(prev)
}
if (!"singler" %in% names(cells)) {
  res_list <- list()
  for (d in keep_donors) {
    message(format(Sys.time(), "%H:%M"), "  ", d)
    seu <- readRDS(file.path("ec_by_donor", paste0(d, ".rds")))
    meta <- seu[[]]; counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts")); rm(seu)
    keep <- str_detect(meta$assay, "10x"); meta <- meta[keep, ]; counts <- counts[, keep, drop = FALSE]
    use <- ref$info$donor != d
    pred <- SingleR(test = lognorm(counts), ref = ref_log[, use], labels = ref$info$subtype[use])
    res_list[[d]] <- tibble(soma_joinid = meta$soma_joinid, singler = pred$labels)
    rm(counts, pred); gc(verbose = FALSE)
  }
  cells <- cells |> left_join(bind_rows(res_list), by = "soma_joinid")
}
lab <- cells |> filter(!is.na(label_subtype))
lab |> group_by(label_subtype) |>
  summarise(cells = n(), pct_agree = round(100 * mean(singler == label_subtype), 1))

# ---- 4. Final subtype: expert label > agreement of both methods > uncertain ----
cells <- cells |> mutate(
  subtype_final  = case_when(!is.na(label_subtype) ~ label_subtype,
                             singler == marker_subtype ~ singler, TRUE ~ "uncertain"),
  subtype_source = case_when(!is.na(label_subtype) ~ "expert label",
                             subtype_final != "uncertain" ~ "consensus", TRUE ~ "uncertain"))
saveRDS(cells, "ec_cells_scored.rds")
cells |> count(subtype_source)

# Organ composition of expert labels (48% of expert capillaries are lung -> reference bias)
cells |> filter(!is.na(label_subtype)) |> count(tissue_general, label_subtype) |>
  pivot_wider(names_from = label_subtype, values_from = n, values_fill = 0) |> print(n = Inf)

# Note: an organ-specific reference (donor x organ x subtype) was also tested and NOT used:
# arterial agreement fell to 62%; heart capillary agreement stayed ~25% (see 03_ts_heart.R).

# =============================================================================
# X-Y DOSAGE OF CHROMATIN REGULATORS (Table 5, Figure 6)
# For each X-Y pair: X copy alone vs combined X + Y (log2 of summed CPM), women minus men.
# Y share = Y CPM / (X + Y CPM) in men.
# Result: X+Y lower in women (FDR < 0.05) for KDM6A/UTY in 9/9 cohort x vessel combinations,
# USP9X/USP9Y 8/9, KDM5C/KDM5D 7/9; Y copy supplies 52-62% (UTY), 60-74% (KDM5D), 37-56% (USP9Y), 32-41% (DDX3Y).
# Makes: p2_xy_dosage.rds (list dres, dos)
# =============================================================================
pairs <- tibble(X = c("KDM6A", "KDM5C", "USP9X", "DDX3X", "EIF1AX", "ZFX", "RPS4X", "NLGN4X"),
                Y = c("UTY", "KDM5D", "USP9Y", "DDX3Y", "EIF1AY", "ZFY", "RPS4Y1", "NLGN4Y")) |>
  mutate(pair = paste(X, Y, sep = "/"))
wide <- long |> select(cohort, type, id, donor_id, sex, age, disease2, gene, cpm) |>
  pivot_wider(names_from = gene, values_from = cpm)
dos <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
  wide |> transmute(cohort, type, id, donor_id, sex, age, disease2, pair = pairs$pair[i],
                    X_cpm = .data[[pairs$X[i]]], Y_cpm = .data[[pairs$Y[i]]])
})) |> mutate(X_only = log2(X_cpm + 1), X_plus_Y = log2(X_cpm + Y_cpm + 1), Y_share = Y_cpm / (X_cpm + Y_cpm))

dres <- dos |> pivot_longer(c(X_only, X_plus_Y), names_to = "measure", values_to = "value") |>
  group_by(cohort, type, pair, measure) |> filter(n_distinct(sex) == 2) |>
  group_modify(~ fit_sex(.x, .y$cohort == "KPMP kidney", "value")) |> ungroup() |>
  group_by(cohort, measure) |> mutate(FDR = p.adjust(p, "BH")) |> ungroup()
# Table 5: X alone / X + Y (* FDR < 0.05 for X + Y)
dres |> mutate(v = sprintf("%+.2f", estimate), star = ifelse(measure == "X_plus_Y" & FDR < 0.05, "*", "")) |>
  group_by(cohort, type, pair) |> summarise(val = paste0(v[measure == "X_only"], " / ", v[measure == "X_plus_Y"],
                                                        star[measure == "X_plus_Y"]), .groups = "drop") |>
  pivot_wider(names_from = c(cohort, type), values_from = val) |> print(width = Inf)
# How often is combined X + Y significantly lower in women?
dres |> filter(measure == "X_plus_Y") |> group_by(pair) |>
  summarise(lower_in_women = sum(estimate < 0 & FDR < 0.05), combos = n())
# Y share in men (median across donors)
yshare <- dos |> filter(sex == "male") |> group_by(cohort, type, pair) |>
  summarise(Y_share = median(Y_share, na.rm = TRUE), .groups = "drop")
yshare |> mutate(Y_share = round(100 * Y_share)) |> pivot_wider(names_from = pair, values_from = Y_share) |> print(width = Inf)
saveRDS(list(dres = dres, dos = dos, yshare = yshare), "p2_xy_dosage.rds")

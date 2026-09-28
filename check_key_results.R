# check_key_results.R — re-runs the steps behind the key numbers in the manuscript (kidney sex x CKD
# tests, GTEx artery tests / Table 4, X-Y dosage / Table 5, hormone receptors) and prints them.
# Needs the saved files from a full run (run_all.R); no large downloads. From the repository folder:
#   source("check_key_results.R")
repo <- normalizePath(getwd())
src  <- function(f) { message("\n----- running ", f, " -----"); source(file.path(repo, "R", f), echo = FALSE) }
hdr  <- function(x) cat("\n\n===================== ", x, " =====================\n")
options(width = 200)

src("00_setup.R")
src("09_hca_genome_wide.R")          # defines Hs, auto and cam2, which scripts 08 and 10 need
src("08_kpmp_kidney.R")

hdr("1. KIDNEY: fry and camera (estimated inter-gene correlation), sex x CKD interaction")
for (s in c("EC-PTC", "EC-GC")) { cat("\n--", s, "--\n"); print(kid_check(s)) }

src("10_gtex_heart_arteries.R")
hdr("2. GTEx ARTERIES: camera with estimated inter-gene correlation")
print(cam_tab |> dplyr::select(artery, analysis, pathway, NGenes, Correlation, Direction, PValue, FDR), n = Inf)
hdr("2b. GTEx ARTERIES: fry (Table 4)")
print(fry_tab |> dplyr::select(artery, analysis, pathway, Direction, PValue, FDR), n = Inf)

src("12_p2_sexchrom_receptor_panel.R")
src("13_p2_xy_dosage.R")
hdr("3. X-Y DOSAGE: Table 5, combos lower in women, Y share")
dres |> dplyr::mutate(v = sprintf("%+.2f", estimate), star = ifelse(measure == "X_plus_Y" & FDR < 0.05, "*", "")) |>
  dplyr::group_by(cohort, type, pair) |>
  dplyr::summarise(val = paste0(v[measure == "X_only"], " / ", v[measure == "X_plus_Y"], star[measure == "X_plus_Y"]),
                   .groups = "drop") |>
  tidyr::pivot_wider(names_from = c(cohort, type), values_from = val) |> print(width = Inf)
print(dres |> dplyr::filter(measure == "X_plus_Y") |> dplyr::group_by(pair) |>
        dplyr::summarise(lower_in_women = sum(estimate < 0 & FDR < 0.05), combos = dplyr::n()))
print(yshare |> dplyr::mutate(Y_share = round(100 * Y_share)) |>
        tidyr::pivot_wider(names_from = pair, values_from = Y_share), width = Inf)
print(res_p2 |> dplyr::filter(gene %in% c("KDM6A", "JPX", "DDX3X", "KDM5C", "EIF1AX", "ZFX", "USP9X")) |>
        dplyr::group_by(gene) |> dplyr::summarise(sig_higher_in_women = sum(estimate > 0 & FDR < 0.05),
                                                  combos = dplyr::n()))

src("15_p2_hormone_receptors.R")
hdr("4. RECEPTORS: age in women, kidney ESR1 (nuclei and single cells)")
print(res_age |> dplyr::summarise(min_FDR_age = min(FDR_age, na.rm = TRUE), min_FDR_int = min(FDR_int, na.rm = TRUE)))
print(esr1_all |> dplyr::mutate(dplyr::across(c(estimate, se), ~ round(.x, 2)), p = signif(p, 2)), width = Inf)
cat("\nDone. Copy the output under the ===== headings into the chat.\n")

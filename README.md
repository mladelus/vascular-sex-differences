# Sex differences in the human vasculature

Code for: **"Sex differences in human vascular endothelium are context-specific: consistent X–Y paralogue dosage differences but no reproducible baseline capillary programs"** (Adelus, 2026; bioRxiv DOI to be added).

**Author:** Maria Adelus, MD, PhD, MPH — Department of Internal Medicine, University of Michigan · ORCID [0000-0002-9676-9214](https://orcid.org/0000-0002-9676-9214) · [Google Scholar](https://scholar.google.com/citations?user=5ETXB4UAAAAJ) · [LinkedIn](https://www.linkedin.com/in/mariaadelus/) · [ResearchGate](https://www.researchgate.net/profile/Maria-Adelus-2)

We compared women's and men's vascular endothelium across public single-cell, single-nucleus and bulk RNA-seq datasets. Women are the reference group throughout: positive estimates mean higher in women. A sex difference was treated as established only if it replicated in an independent cohort.

## Main findings

1. **Healthy capillaries:** no reproducible sex difference in senescence, inflammatory, oxidative-stress or genome-wide programs (Tabula Sapiens 2.0, Heart Cell Atlas, an independent left-ventricle cohort and GTEx).
2. **Vascular bed:** before age 50, women's coronary arteries had lower interferon-γ and inflammatory signaling than men's (self-contained fry test, FDR 0.03). No difference was seen at 50 or older, or in the aorta (GTEx).
3. **Disease:** in glomerular capillaries, CKD was associated with modest, opposite shifts in interferon-γ signaling in women and men (sex × disease, fry p = 0.007; KPMP).
4. **Sex-chromosome signature:** KDM6A and JPX were consistently higher in women across three single-nucleus cohorts. In men, however, the Y paralogues supplied 52–74% of KDM6A/UTY and KDM5C/KDM5D expression, so combined X + Y dosage of these chromatin regulators was *lower* in women. This held in all three single-nucleus cohorts.
5. Loss of Y was rare in kidney endothelium (1.6% of deep nuclei; upper 95% limit 3.2%); heart nuclei were too shallow to estimate it. Endothelium expressed AR, ESR1 and PGR but not aromatase, and receptor levels were stable with age.

## Data (all public)

| Source | Access | Used for |
|---|---|---|
| Tabula Sapiens 2.0 | CZ CELLxGENE Census, release `2025-11-08` | Discovery, many organs |
| Heart Cell Atlas (Litviňuková 2020), vascular nuclei | Census | Heart validation, genome-wide |
| Independent healthy left ventricle ("Heart", HCA harmonization collection) | Census | Replication |
| KPMP / Human Kidney Atlas v1.5, single-nucleus + single-cell | Census + source `.h5ad` (author labels) | Kidney sex × CKD, ESR1 |
| GTEx v8 gene read counts: LV, atrial appendage, coronary, aorta, tibial | `storage.googleapis.com/adult-gtex` (open access) | Bulk replication, arteries |

The scripts download everything themselves. No data files are stored in this repository.

## Scripts (`R/`)

| Script | What it does | Manuscript |
|---|---|---|
| `00_setup.R` | Packages, Census, helpers, gene sets. **Run first every session.** | — |
| `01`–`05` | Tabula Sapiens: cohort, cell scores, vessel subtypes, heart relabeling, sex analyses | Table 2 |
| `06_hca_validation.R` | Heart Cell Atlas nuclei: capillary programs by sex | Fig 2 |
| `07_hca_oxidative_stress.R` | Pericytes and capillary ECs, oxidative stress | Results |
| `08_kpmp_kidney.R` | KPMP kidney: sex × CKD | Fig 5 |
| `09_hca_genome_wide.R` | Genome-wide capillary comparison (limma-voom with donor blocking; gene sets by cameraPR) | Fig 3 |
| `10_gtex_heart_arteries.R` | GTEx heart and arteries, adjusted for ischemic time and manner of death | Table 4, Fig 4 |
| `11_independent_lv_replication.R` | Independent LV replication, mural-cell ratios | Table 3, Fig 3 |
| `12_p2_sexchrom_receptor_panel.R` | X/Y genes and hormone receptors from saved pseudobulks | Results |
| `13_p2_xy_dosage.R` | X–Y dosage of chromatin regulators | Table 5, Fig 6 |
| `14_p2_loss_of_y.R` | Loss of chromosome Y (depth-aware) | Fig 7 |
| `15_p2_hormone_receptors.R` | Receptors by age, kidney ESR1, single-cell replication | Fig 8 |
| `16_figures.R` | All figures as multi-panel PDF, PNG and TIFF, and `sessionInfo.txt` (S1 kidney score, S2 loss of Y, S3 receptors) | Figs 2–8, S1–S3 |
| `99_optional_future_work.R` | Not run; not in the manuscript | — |

`05_ts_exploratory.R` contains one commented-out block (stress-adjusted genome-wide test) that was not run for the manuscript.

## How to run

```r
# everything, in order (long first run: downloads ~10 GB)
source("run_all.R")

# re-check the key manuscript numbers from saved results (after a full run)
source("check_key_results.R")

# figures only, from saved results
source("R/00_setup.R"); source("R/16_figures.R")
```

Data go to `~/Documents/ts2_project` by default. To change this, set the environment variable `VSD_DATA` or edit `DATA_DIR` in `R/00_setup.R`. Figures are written to `<data folder>/figures/`. Figure 1, the study design, is drawn by `figures/make_fig1.py` (Python, matplotlib) and supplied as `figures/fig1_study_design.png`.

## Requirements

- R ≥ 4.3.
- CRAN: `cellxgene.census` (r-universe), `Seurat`, `Matrix`, `lme4`, `lmerTest`, `msigdbr`, `dplyr`, `tidyr`, `stringr`, `ggplot2`, `patchwork`, `ggrepel`, `scales`, `data.table`.
- Bioconductor: `UCell`, `SingleR`, `scDblFinder`, `SingleCellExperiment`, `edgeR`, `limma`, `EnsDb.Hsapiens.v86`, `rhdf5`.

Installation commands are at the end of `R/00_setup.R`. A laptop with 16 GB RAM is enough, because downloads are done one donor at a time.

## Notes

- `dplyr` verbs are masked by several Bioconductor packages. The scripts re-assign `filter`, `select`, `count`, `first`, `rename` and `desc` to the `dplyr` versions, and call `edgeR::cpm` explicitly.
- `TXLNGY` is excluded from the Y-gene set because its feature name is duplicated in the Census.
- The KPMP single-cell file labels only `subclass.l1`, so glomerular-like endothelial cells are defined as EHD3+.

## How to cite

See `CITATION.cff` (GitHub shows a "Cite this repository" button).

## License

MIT (see `LICENSE`).

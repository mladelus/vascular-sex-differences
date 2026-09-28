# =============================================================================
# SETUP: run this first in every new R session
# Sex differences in human capillary endothelium (Tabula Sapiens 2.0 + HCA heart + KPMP kidney)
# Sets the folder, loads packages, opens CELLxGENE Census, defines helpers and gene sets.
# =============================================================================

# Data folder (downloads, .rds files). Change here or set the env variable VSD_DATA.
DATA_DIR <- Sys.getenv("VSD_DATA", "~/Documents/ts2_project")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(DATA_DIR)
CENSUS_VERSION <- "2025-11-08"          # pinned for reproducibility (report in Methods)

# ---- Packages (install once; see bottom of file) ----------------------------
suppressPackageStartupMessages({
  library(cellxgene.census)
  library(Seurat); library(Matrix)
  library(UCell); library(SingleR)
  library(scDblFinder); library(SingleCellExperiment)
  library(edgeR); library(limma)
  library(lme4); library(lmerTest)       # lmerTest::lmer gives small-sample (Satterthwaite) p-values
  library(ggplot2); library(tidyr); library(stringr); library(dplyr)
})

# Bioconductor packages mask several dplyr verbs; always use dplyr's versions
filter <- dplyr::filter; select <- dplyr::select; count <- dplyr::count
first  <- dplyr::first;  rename <- dplyr::rename; slice <- dplyr::slice; desc <- dplyr::desc

census <- open_soma(census_version = CENSUS_VERSION)

# ---- Gene ID -> symbol table (created once) ----------------------------------
if (!file.exists("gene_map.rds")) {
  genes <- census$get("census_data")$get("homo_sapiens")$ms$get("RNA")$var$read(
    column_names = c("feature_id", "feature_name"))$concat() |> as.data.frame()
  saveRDS(genes, "gene_map.rds")
}
gene_map <- readRDS("gene_map.rds")

# ---- Helpers -----------------------------------------------------------------

# Convert Ensembl-ID rows to gene symbols, summing IDs that share a symbol
to_symbols <- function(m) {
  sym <- gene_map$feature_name[match(rownames(m), gene_map$feature_id)]
  f <- factor(sym)
  out <- fac2sparse(f) %*% m
  rownames(out) <- levels(f)
  out
}

# log2(counts per 10,000 + 1) for a sparse cells x genes matrix (genes in rows)
lognorm <- function(m) {
  m <- t(t(m) / Matrix::colSums(m)) * 1e4
  m@x <- log2(m@x + 1)
  m
}

# MSigDB gene sets (loaded once per session)
.msig <- list()
get_set <- function(name, coll) {
  if (is.null(.msig[[coll]])) {
    .msig[[coll]] <<- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = coll),
                               error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = coll))
  }
  unique(.msig[[coll]]$gene_symbol[.msig[[coll]]$gs_name == name])
}

# Age from "56-year-old stage" or "sixth decade stage" (decade -> midpoint)
dec_mid <- c("third decade stage" = 25, "fourth decade stage" = 35, "fifth decade stage" = 45,
             "sixth decade stage" = 55, "seventh decade stage" = 65, "eighth decade stage" = 75,
             "ninth decade stage" = 85)
parse_age <- function(x) {
  x <- as.character(x)
  coalesce(as.numeric(str_extract(x, "^\\d+(?=-year)")), unname(dec_mid[x]))
}

# Odds ratios (Wald 95% CI) from a glmer fit, for terms matching `pat`
show <- function(fit, pat = "sex|age10|vessel|postmeno") {
  co <- summary(fit)$coefficients
  ci <- confint(fit, method = "Wald", parm = "beta_")
  k <- grepl(pat, rownames(co))
  data.frame(term = rownames(co)[k],
             OR = round(exp(co[k, 1]), 2),
             conf.low = round(exp(ci[rownames(co)[k], 1]), 2),
             conf.high = round(exp(ci[rownames(co)[k], 2]), 2),
             p.value = signif(co[k, 4], 2), row.names = NULL)
}

# Estimates with small-sample (Satterthwaite) p-values from a lmerTest::lmer fit
st <- function(m, pat = "sex") {
  co <- summary(m)$coefficients; k <- grepl(pat, rownames(co))
  data.frame(term = rownames(co)[k], estimate = signif(co[k, 1], 3),
             conf.low = signif(co[k, 1] - 1.96 * co[k, 2], 3),
             conf.high = signif(co[k, 1] + 1.96 * co[k, 2], 3),
             df = round(co[k, "df"], 1), p = signif(co[k, "Pr(>|t|)"], 2), row.names = NULL)
}

glmer_ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# ---- Gene sets (single source of truth) -------------------------------------
STRESS <- c("FOS","FOSB","JUN","JUNB","JUND","ATF3","EGR1","IER2","IER3","DUSP1",
            "ZFP36","HSPA1A","HSPA1B","HSPA8","HSPH1","DNAJB1","HSP90AA1","HSPE1","SOCS3","NR4A1")

# Tabula Sapiens cell scoring + vessel-type markers (final versions; VWF/SELP deliberately excluded)
ts_sets <- list(
  SenMayo      = get_set("SAUL_SEN_MAYO", "C2"),
  Inflammatory = get_set("HALLMARK_INFLAMMATORY_RESPONSE", "H"),
  TNFA         = get_set("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "H"),
  Stress       = STRESS,
  Arterial     = c("GJA5","HEY1","SEMA3G","GJA4","DKK2","FBLN5"),
  Capillary    = c("CA4","RGCC","BTNL9","CD36"),
  Venous       = c("ACKR1","NR2F2","CPE"),
  Lymphatic    = c("PROX1","LYVE1","CCL21","PDPN","FLT4","MMRN1"))

# Heart-atlas nucleus scoring
hg_sets <- list(
  SenMayo      = ts_sets$SenMayo,
  Inflammatory = ts_sets$Inflammatory,
  TNFA         = ts_sets$TNFA,
  OXPHOS       = get_set("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "H"),
  Stress       = STRESS)

# Muscle-type / mural panels
pan <- list(peri = c("RGS5","PDGFRB","KCNJ8","ABCC9","NOTCH3"),
            sm   = c("ACTA2","TAGLN","CNN1","MYH11","MYL9","DES","SMTN","LMOD1"),
            str  = c("MYL2","MB","ACTA1","COX6A2","TNNT2","MYH7","CKM"))

# Oxidative-stress programs (pericytes and capillary ECs)
ox_sets <- list(
  NRF2         = c("NQO1","HMOX1","GCLM","GCLC","TXNRD1","SRXN1","SLC7A11","G6PD","PGD","ME1","GSR","FTH1","FTL"),
  Antioxidant  = c("SOD1","SOD2","SOD3","CAT","GPX1","GPX3","GPX4","PRDX1","PRDX2","PRDX3","PRDX5","PRDX6","TXN","TXNRD1","GSR"),
  ROS_sources  = c("NOX4","CYBB","NOX1","XDH","CYBA","NCF1","NCF2"),
  ROS_hallmark = get_set("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "H"),
  Contractile  = c("ACTA2","MYH11","CNN1","TAGLN","MYL9","DES","KCNJ8","ABCC9"),
  eNOS         = c("NOS3","GCH1","SLC7A1","CAV1","HSP90AA1"),
  Stress       = STRESS)

# Reproductive organs (excluded from sex comparisons)
REPRO <- "ovary|uter|endometri|myometri|prostate|testis|mammary|fallopian|vagin|cervix"

message("Setup complete. Census ", CENSUS_VERSION, " open.")

# ---- One-time installation (run manually if a package is missing) -----------
# install.packages(c("BiocManager","dplyr","tidyr","stringr","ggplot2","Matrix","Seurat",
#                    "lme4","lmerTest","msigdbr","httr"))
# install.packages("cellxgene.census",
#                  repos = c("https://chanzuckerberg.r-universe.dev", "https://cloud.r-project.org"))
# options(timeout = 600)
# BiocManager::install(c("UCell","SingleR","scDblFinder","SingleCellExperiment","edgeR","limma",
#                        "rtracklayer","EnsDb.Hsapiens.v86","rhdf5"), update = FALSE, ask = FALSE)

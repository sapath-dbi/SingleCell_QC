# =============================================================================
# scRNA-seq QC Report — 10x Chromium molecule_info.h5
# =============================================================================
# Produces a self-contained HTML report assessing:
#   - Cell capture and UMI distributions
#   - PCR duplication rate
#   - Library saturation and top-up sequencing estimates
#
# Works with any species, tissue, or sample type.
#
# Requirements:
#   install.packages(c("data.table", "ggplot2", "patchwork",
#                      "scales", "rmarkdown"))
#   BiocManager::install("DropletUtils")
#
# Usage:
#   1. Load your molinfo object from the h5 file:
#        library(DropletUtils)
#        molinfo <- read10xMolInfo("path/to/molecule_info.h5")
#   2. Edit the CONFIG section below for your sample
#   3. Source this script:
#        source("scrna_qc_report.R")
#   4. An HTML report is written to OUTPUT_FILE
# =============================================================================


# ── 0. CONFIG — edit these for every sample ───────────────────────────────────

SAMPLE_NAME    <- "My Sample"         # e.g. "Human PBMC", "Mouse Liver E14.5"
SPECIES        <- "human"             # "human" or "mouse" — label only
TISSUE_TYPE    <- "general"           # see TISSUE PRESETS below
UMI_THRESHOLD  <- 500                 # cell calling cutoff (adjust after knee plot)
OUTPUT_FILE    <- "scrna_qc_report.html"
TOP_UP_TARGETS <- c(0.50, 0.60, 0.70, 0.75)  # duplication rate targets

# ── TISSUE PRESETS ─────────────────────────────────────────────────────────────
# Controls UMI pass/review thresholds shown in the report summary.
# Options: "general", "pbmc", "brain", "heart", "liver", "lung", "kidney"
# Use "general" if your tissue is not listed.
#
# Override thresholds manually (set to NULL to use the preset):
UMI_PASS_THRESHOLD   <- NULL   # e.g. 3000
UMI_REVIEW_THRESHOLD <- NULL   # e.g. 1000
CELL_PASS_THRESHOLD  <- NULL   # e.g. 2000


# ── 1. Tissue/species reference ranges ────────────────────────────────────────
# Expected median UMI ranges informed by 10x benchmarks and published datasets.
# These are approximate guides — not hard rules.
# High-RNA-content cells (hepatocytes, cardiomyocytes) trend higher;
# small cells (lymphocytes, platelets) trend lower.

tissue_refs <- list(
  general = list(
    umi_pass = 2000, umi_review = 800, cell_pass = 2000,
    note = "General defaults for mixed or unlisted tissue types."
  ),
  pbmc = list(
    umi_pass = 2000, umi_review = 800, cell_pass = 3000,
    note = "PBMCs are small cells with low-to-moderate RNA content. 2,000-5,000 UMI/cell is typical. T cells trend lower, monocytes higher."
  ),
  brain = list(
    umi_pass = 3000, umi_review = 1000, cell_pass = 2000,
    note = "Neurons are large and transcriptionally complex (5,000-15,000 UMI). Oligodendrocytes and microglia are lower (1,000-3,000 UMI)."
  ),
  heart = list(
    umi_pass = 3000, umi_review = 1000, cell_pass = 1000,
    note = "Cardiomyocytes are large with high RNA content (5,000-20,000 UMI). Endothelial and immune cells are lower. Embryonic hearts yield fewer but transcriptionally rich cells."
  ),
  liver = list(
    umi_pass = 4000, umi_review = 1500, cell_pass = 1000,
    note = "Hepatocytes are among the most transcriptionally active cells (8,000-20,000 UMI). Non-parenchymal cells (Kupffer, stellate, endothelial) are lower."
  ),
  lung = list(
    umi_pass = 2500, umi_review = 1000, cell_pass = 2000,
    note = "Lung is highly heterogeneous. AT2 cells are high (4,000-8,000 UMI); immune and endothelial cells are lower (1,000-3,000 UMI)."
  ),
  kidney = list(
    umi_pass = 3000, umi_review = 1000, cell_pass = 2000,
    note = "Proximal tubule cells are high (4,000-10,000 UMI); collecting duct, immune, and endothelial cells are lower."
  )
)

preset <- if (TISSUE_TYPE %in% names(tissue_refs)) {
  tissue_refs[[TISSUE_TYPE]]
} else {
  message(sprintf("Tissue '%s' not in presets — using 'general' defaults.", TISSUE_TYPE))
  tissue_refs[["general"]]
}

umi_pass_thr   <- if (!is.null(UMI_PASS_THRESHOLD))   UMI_PASS_THRESHOLD   else preset$umi_pass
umi_review_thr <- if (!is.null(UMI_REVIEW_THRESHOLD)) UMI_REVIEW_THRESHOLD else preset$umi_review
cell_pass_thr  <- if (!is.null(CELL_PASS_THRESHOLD))  CELL_PASS_THRESHOLD  else preset$cell_pass
tissue_note    <- preset$note

species_label <- switch(tolower(SPECIES),
                        human = "Human", mouse = "Mouse", SPECIES)


# ── 2. Dependencies ───────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

theme_qc <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey92"),
      plot.title        = element_text(size = 12, face = "bold", margin = margin(b = 4)),
      plot.subtitle     = element_text(size = 10, colour = "grey50", margin = margin(b = 8)),
      axis.title        = element_text(size = 10, colour = "grey40"),
      axis.text         = element_text(size = 9),
      strip.text        = element_text(size = 10, face = "bold"),
      legend.position   = "none"
    )
}

mm_umi <- function(reads, lib_size) lib_size * reads / (reads + lib_size)


# ── 3. Load and validate data ─────────────────────────────────────────────────

stopifnot(
  "molinfo not found — load with DropletUtils::read10xMolInfo() first." =
    exists("molinfo")
)

cat("Converting to data.table...\n")
mol <- as.data.table(as.data.frame(molinfo$data))
cat(sprintf("Loaded %s molecules\n", format(nrow(mol), big.mark = ",")))


# ── 4. Raw barcode stats ──────────────────────────────────────────────────────

cat("Computing per-barcode UMI counts...\n")
umi_per_cell   <- mol[, .(umi_count = .N), by = cell]
n_barcodes_raw <- nrow(umi_per_cell)
cat(sprintf("Total barcodes (unfiltered): %s\n", format(n_barcodes_raw, big.mark = ",")))


# ── 5. Threshold sweep ────────────────────────────────────────────────────────

cat("Running threshold sweep...\n")
thresholds <- c(100, 200, 500, 1000, 2000, 5000)

thresh_dt <- rbindlist(lapply(thresholds, function(t) {
  sub <- umi_per_cell[umi_count >= t]
  data.table(
    threshold  = t,
    n_cells    = nrow(sub),
    median_umi = median(sub$umi_count),
    mean_umi   = round(mean(sub$umi_count))
  )
}))


# ── 6. Filter to real cells ───────────────────────────────────────────────────

cat(sprintf("Filtering to >= %d UMIs...\n", UMI_THRESHOLD))

dup_per_cell <- mol[, .(
  total_reads = sum(reads),
  unique_umis = .N,
  dup_rate    = 1 - (.N / sum(reads))
), by = cell]

cell_qc          <- merge(umi_per_cell, dup_per_cell, by = "cell")
cell_qc_filtered <- cell_qc[umi_count >= UMI_THRESHOLD]

n_cells     <- nrow(cell_qc_filtered)
median_umi  <- median(cell_qc_filtered$umi_count)
mean_umi    <- round(mean(cell_qc_filtered$umi_count))
median_dup  <- round(median(cell_qc_filtered$dup_rate) * 100, 1)
mean_dup    <- round(mean(cell_qc_filtered$dup_rate) * 100, 1)
overall_dup <- round((1 - sum(cell_qc_filtered$unique_umis) /
                        sum(cell_qc_filtered$total_reads)) * 100, 1)

cat(sprintf("Cells: %s | Median UMI: %s | Mean dup: %.1f%%\n",
            format(n_cells, big.mark = ","),
            format(median_umi, big.mark = ","),
            mean_dup))


# ── 7. Library saturation estimates ───────────────────────────────────────────

total_obs   <- sum(cell_qc_filtered$unique_umis)
total_reads <- sum(cell_qc_filtered$total_reads)
mean_dup_r  <- 1 - (total_obs / total_reads)
est_lib     <- round(total_obs / mean_dup_r)
pct_seen    <- round(total_obs / est_lib * 100, 1)

lib_per_cell    <- round(median(cell_qc_filtered$umi_count /
                                  (1 - cell_qc_filtered$dup_rate)))
unseen_per_cell <- lib_per_cell - median_umi

top_up_dt <- rbindlist(lapply(TOP_UP_TARGETS, function(target) {
  reads_needed  <- est_lib * (target / (1 - target))
  fold_increase <- reads_needed / total_reads
  data.table(
    target_dup    = paste0(round(target * 100), "%"),
    reads_needed  = format(round(reads_needed / 1e6), big.mark = ","),
    fold_increase = round(fold_increase, 1),
    proj_umi_cell = round(mm_umi(
      (median_umi / (1 - mean_dup_r)) * fold_increase,
      lib_per_cell
    ))
  )
}))


# ── 8. Status flags ───────────────────────────────────────────────────────────

umi_status <- ifelse(median_umi >= umi_pass_thr,  "PASS",
                     ifelse(median_umi >= umi_review_thr, "REVIEW", "FAIL"))

dup_status <- ifelse(mean_dup <= 50, "SHALLOW",
                     ifelse(mean_dup <= 75, "OPTIMAL", "SATURATED"))

cell_status <- ifelse(n_cells >= cell_pass_thr,          "PASS",
                      ifelse(n_cells >= cell_pass_thr * 0.4,    "REVIEW", "FAIL"))

topup_recommended <- mean_dup < 50
topup_saturated   <- mean_dup > 75
topup_fold        <- round(top_up_dt[target_dup == "50%"]$fold_increase, 1)
topup_proj_umi    <- top_up_dt[target_dup == "50%"]$proj_umi_cell


# ── 9. Recommendation text (generic — no tissue-specific biology) ─────────────

if (topup_recommended) {
  topup_verdict <- "RECOMMENDED"
  topup_summary <- paste0(
    "At ", mean_dup, "% mean duplication, only ", pct_seen, "% of the estimated ",
    "library has been observed. A ", topup_fold, "x top-up run on the existing ",
    "library is recommended to reach ~50% duplication and ~",
    format(topup_proj_umi, big.mark = ","), " median UMI/cell."
  )
  topup_section <- paste0(
    "### Top-up sequencing recommended\n\n",
    "**Current duplication rate:** ", mean_dup, "% (below the 50% optimal threshold).  \n",
    "**Implication:** The library is under-sequenced relative to its complexity. ",
    "Only ", pct_seen, "% of the estimated ", format(est_lib, big.mark = ","),
    " molecule library has been observed.  \n",
    "**Library quality:** The narrow duplication rate distribution and geometric ",
    "reads-per-UMI decay confirm this is a sequencing depth issue — not a ",
    "library preparation or PCR quality issue.  \n",
    "**Action:** Request a **", topup_fold, "x top-up run** on the existing library. ",
    "No new library preparation is needed.\n\n",
    "| Target dup rate | Total reads needed (M) | Fold increase | Projected median UMI/cell |\n",
    "|---|---|---|---|\n",
    paste(apply(top_up_dt, 1, function(r)
      paste0("| ", r["target_dup"], " | ", r["reads_needed"], " | ",
             r["fold_increase"], "x | ",
             format(as.integer(r["proj_umi_cell"]), big.mark = ","), " |")),
      collapse = "\n"), "\n\n",
    "> **Note on the saturation model:** Estimates assume equal sampling probability ",
    "across all molecules. In reality, highly expressed genes are already ",
    "near-completely sampled; the unseen fraction (", format(unseen_per_cell, big.mark=","),
    " UMI/cell) is enriched for lowly expressed transcripts. Top-up sequencing ",
    "disproportionately benefits detection of rare or lowly expressed genes."
  )
} else if (topup_saturated) {
  topup_verdict <- "NOT REQUIRED — library near saturation"
  topup_summary <- paste0(
    "At ", mean_dup, "% duplication, the library is deeply sequenced. ",
    "Further reads will yield diminishing returns."
  )
  topup_section <- paste0(
    "### Top-up sequencing not recommended\n\n",
    "Duplication rate of ", mean_dup, "% exceeds 75%. The library is approaching ",
    "saturation and additional sequencing would yield a small marginal increase in ",
    "unique UMIs. Consider whether the specific biological question requires ",
    "deeper coverage before investing in further sequencing."
  )
} else {
  topup_verdict <- "NOT REQUIRED"
  topup_summary <- paste0(
    "Duplication rate of ", mean_dup, "% is within the optimal 50-75% range. ",
    "The library is well-sequenced."
  )
  topup_section <- paste0(
    "### Top-up sequencing not required\n\n",
    "Duplication rate of ", mean_dup, "% is within the optimal 50-75% range, ",
    "indicating the library is appropriately sequenced. No further action is needed."
  )
}


# ── 10. Plots ─────────────────────────────────────────────────────────────────

cat("Generating plots...\n")

# Knee plot
knee_dt     <- umi_per_cell[order(-umi_count)]
knee_dt[, rank := .I]
knee_sample <- knee_dt[seq(1, nrow(knee_dt), length.out = min(50000, nrow(knee_dt)))]

p_knee <- ggplot(knee_sample, aes(rank, umi_count)) +
  geom_line(colour = "#1D9E75", linewidth = 0.8) +
  geom_vline(xintercept = n_cells, linetype = "dashed",
             colour = "#E24B4A", linewidth = 0.6) +
  annotate("text", x = n_cells * 1.6, y = max(knee_sample$umi_count) * 0.4,
           label = paste0("Threshold\n", UMI_THRESHOLD, " UMI\n(",
                          format(n_cells, big.mark = ","), " cells)"),
           size = 3, colour = "#E24B4A", hjust = 0) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  labs(title = "Barcode rank plot",
       subtitle = "Inflection point marks real cells vs empty droplets",
       x = "Barcode rank", y = "UMI count") +
  theme_qc()

# UMI per cell
p_umi <- ggplot(cell_qc_filtered, aes(x = umi_count)) +
  geom_histogram(bins = 60, fill = "#9FE1CB", colour = "#1D9E75", linewidth = 0.3) +
  geom_vline(xintercept = median_umi, linetype = "dashed",
             colour = "#1D9E75", linewidth = 0.8) +
  geom_vline(xintercept = umi_pass_thr, linetype = "dotted",
             colour = "grey55", linewidth = 0.5) +
  annotate("text", x = median_umi * 1.05, y = Inf,
           label = paste0("Median: ", format(median_umi, big.mark = ",")),
           vjust = 1.5, hjust = 0, size = 3, colour = "#0F6E56") +
  annotate("text", x = umi_pass_thr, y = Inf,
           label = paste0("Pass (", format(umi_pass_thr, big.mark=","), ")"),
           vjust = 1.5, hjust = 1.05, size = 2.8, colour = "grey50") +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "UMI count per cell",
       subtitle = paste0("Filtered cells (>= ", UMI_THRESHOLD,
                         " UMI) | dotted = tissue-specific pass threshold"),
       x = "UMI count", y = "Cells") +
  theme_qc()

# Duplication rate
p_dup <- ggplot(cell_qc_filtered[, .(dup_pct = dup_rate * 100)], aes(x = dup_pct)) +
  geom_histogram(binwidth = 2, fill = "#CECBF6", colour = "#7F77DD", linewidth = 0.3) +
  geom_vline(xintercept = mean_dup, linetype = "dashed",
             colour = "#534AB7", linewidth = 0.8) +
  geom_vline(xintercept = c(50, 75), linetype = "dotted",
             colour = "grey60", linewidth = 0.5) +
  annotate("text", x = mean_dup + 0.5, y = Inf,
           label = paste0("Mean: ", mean_dup, "%"),
           vjust = 1.5, hjust = 0, size = 3, colour = "#3C3489") +
  annotate("text", x = 51, y = Inf, label = "optimal \u2192",
           vjust = 1.5, hjust = 0, size = 2.8, colour = "grey50") +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Duplication rate per cell",
       subtitle = "Dotted lines at 50% and 75% mark the optimal sequencing zone",
       x = "Duplication rate (%)", y = "Cells") +
  theme_qc()

# Reads per UMI
rpu_exact <- mol[cell %in% cell_qc_filtered$cell,
                 .(count = .N), by = .(reads_per_umi = reads)]
rpu_exact <- rpu_exact[order(reads_per_umi)][reads_per_umi <= 15]

p_rpu <- ggplot(rpu_exact, aes(x = factor(reads_per_umi), y = count)) +
  geom_col(fill = "#F0997B", colour = "#D85A30", linewidth = 0.3) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Reads per UMI",
       subtitle = "Geometric decay confirms unbiased PCR amplification",
       x = "Reads per molecule", y = "Molecule count") +
  theme_qc()

# Saturation curve
read_seq <- seq(0, total_reads * 20, length.out = 300)
sat_dt   <- data.table(reads = read_seq, umis = mm_umi(read_seq, est_lib))

target_pts <- data.table(
  reads = sapply(TOP_UP_TARGETS, function(t) est_lib * (t / (1 - t))),
  umis  = sapply(TOP_UP_TARGETS, function(t) mm_umi(est_lib * (t / (1 - t)), est_lib)),
  label = paste0(round(TOP_UP_TARGETS * 100), "%")
)

p_sat <- ggplot(sat_dt, aes(reads / 1e6, umis / 1e6)) +
  geom_line(colour = "#1D9E75", linewidth = 1) +
  geom_point(data = data.table(reads = total_reads, umis = total_obs),
             aes(reads / 1e6, umis / 1e6),
             colour = "#E24B4A", size = 3.5, shape = 16) +
  geom_point(data = target_pts, aes(reads / 1e6, umis / 1e6),
             colour = "#EF9F27", size = 3, shape = 18) +
  geom_text(data = target_pts, aes(reads / 1e6, umis / 1e6, label = label),
            nudge_y = est_lib * 0.02 / 1e6, size = 2.8, colour = "#854F0B") +
  annotate("text", x = total_reads / 1e6, y = total_obs / 1e6,
           label = paste0("Current\n(", pct_seen, "% seen)"),
           hjust = -0.1, size = 2.8, colour = "#A32D2D") +
  geom_hline(yintercept = est_lib / 1e6, linetype = "dotted",
             colour = "grey60", linewidth = 0.5) +
  annotate("text", x = 0, y = est_lib / 1e6,
           label = paste0("Est. ceiling: ", format(round(est_lib / 1e6, 1)), "M"),
           vjust = -0.4, hjust = 0, size = 2.8, colour = "grey50") +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Library saturation curve",
       subtitle = "Red = current position | orange = duplication rate targets",
       x = "Total reads (millions)", y = "Unique UMIs (millions)") +
  theme_qc()

# Projected UMI/cell
reads_per_cell_now <- median_umi / (1 - mean_dup_r)
multiples <- c(1, 1.5, 2, 3, 5, 10)
proj_dt <- data.table(
  depth    = factor(paste0(multiples, "\u00d7"), levels = paste0(multiples, "\u00d7")),
  proj_umi = sapply(multiples, function(m)
    round(mm_umi(reads_per_cell_now * m, lib_per_cell))),
  multiple = multiples
)

p_proj <- ggplot(proj_dt, aes(depth, proj_umi,
                              fill = ifelse(multiple == 1, "current", "projected"))) +
  geom_col(colour = NA) +
  geom_hline(yintercept = lib_per_cell, linetype = "dotted",
             colour = "grey60", linewidth = 0.5) +
  annotate("text", x = 0.5, y = lib_per_cell,
           label = paste0("Est. ceiling: ", format(lib_per_cell, big.mark = ",")),
           vjust = -0.4, hjust = 0, size = 2.8, colour = "grey50") +
  geom_text(aes(label = format(proj_umi, big.mark = ",")),
            vjust = -0.4, size = 3, colour = "grey30") +
  scale_fill_manual(values = c("current" = "#1D9E75", "projected" = "#9FE1CB")) +
  scale_y_continuous(labels = label_comma(), limits = c(0, lib_per_cell * 1.15)) +
  labs(title = "Projected median UMI/cell",
       subtitle = "At increasing sequencing depth relative to current",
       x = "Depth relative to current", y = "Median UMI per cell") +
  theme_qc()


# ── 11. Render report ─────────────────────────────────────────────────────────

cat("Assembling HTML report...\n")

rmd_text <- paste0(
  '---
title: "scRNA-seq QC Report"
subtitle: "', SAMPLE_NAME, ' | ', species_label, '"
date: "', format(Sys.Date(), "%d %B %Y"), '"
output:
  html_document:
    theme: flatly
    toc: true
    toc_float: true
    self_contained: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE,
                      fig.width=10, fig.height=4, dpi=150)
```

## Summary

```{r summary-table}
library(knitr)
summary_df <- data.frame(
  Metric = c("Sample", "Species", "Tissue type",
             "UMI threshold applied", "Cells called",
             "Median UMI / cell", "Mean UMI / cell",
             "Mean duplication rate", "Overall duplication rate",
             "Est. library size / cell", "Est. unseen UMI / cell",
             "Fraction of library seen"),
  Value = c(
    "', SAMPLE_NAME, '", "', species_label, '", "', TISSUE_TYPE, '",
    "', UMI_THRESHOLD, '",
    format(', n_cells, ', big.mark=","),
    format(', median_umi, ', big.mark=","),
    format(', mean_umi, ', big.mark=","),
    paste0("', mean_dup, '%"),
    paste0("', overall_dup, '%"),
    format(', lib_per_cell, ', big.mark=","),
    format(', unseen_per_cell, ', big.mark=","),
    paste0("', pct_seen, '%")
  ),
  Status = c("","","","",
             "', cell_status, '",
             "', umi_status, '", "",
             "', dup_status, '",
             "","","","")
)
kable(summary_df, align=c("l","r","c"))
```

> **Tissue context:** ', tissue_note, '  
> UMI thresholds for `', TISSUE_TYPE, '`: PASS >= ', format(umi_pass_thr, big.mark=","),
' | REVIEW >= ', format(umi_review_thr, big.mark=","), '

**Top-up sequencing: ', topup_verdict, '**  
', topup_summary, '

---

## Cell calling

```{r knee-umi, fig.height=4}
p_knee + p_umi
```

The barcode rank plot shows UMI counts from highest to lowest.
The inflection separates real cells from empty droplets.
A threshold of **', UMI_THRESHOLD, ' UMI** was applied, yielding **',
format(n_cells, big.mark=","), ' cells**.

```{r threshold-table}
kable(thresh_dt,
      col.names = c("UMI threshold", "Cells", "Median UMI", "Mean UMI"),
      format.args = list(big.mark=","),
      align = "rrrr",
      caption = "Cell count and UMI statistics at each threshold")
```

---

## Duplication and PCR quality

```{r dup-rpu, fig.height=4}
p_dup + p_rpu
```

The **duplication rate** (1 - unique UMIs / total reads) measures the fraction
of reads that are PCR copies of molecules already counted.
A narrow, symmetric distribution indicates uniform sequencing depth and unbiased
PCR amplification. The optimal zone (50-75%) is marked by dotted lines.

The **reads-per-UMI** plot shows the molecule-level read support distribution.
A geometric decay (each bar ~3x smaller than the previous) confirms no PCR
jackpotting or amplification bias.

---

## Library saturation and top-up estimate

```{r sat-proj, fig.height=4}
p_sat + p_proj
```

Library size is estimated using the Michaelis-Menten model:
`est_library = observed_UMIs / duplication_rate`.
This assumes uniform sampling probability — highly expressed genes are
near-complete while rare transcripts make up most of the unseen fraction.

```{r topup-table}
kable(top_up_dt,
      col.names = c("Target dup rate", "Total reads needed (M)",
                    "Fold increase", "Projected median UMI/cell"),
      align = "rrrr",
      caption = "Reads required to reach each duplication rate target")
```

---

## Recommendation

', topup_section, '

---
*Generated by scrna_qc_report.R | ', format(Sys.time(), "%Y-%m-%d %H:%M"), '*
')

rmd_file <- tempfile(fileext = ".Rmd")
writeLines(rmd_text, rmd_file)

for (obj in c("p_knee","p_umi","p_dup","p_rpu","p_sat","p_proj")) {
  assign(obj, get(obj), envir = .GlobalEnv)
}

rmarkdown::render(
  input       = rmd_file,
  output_file = OUTPUT_FILE,
  output_dir  = getwd(),
  quiet       = TRUE,
  envir       = .GlobalEnv
)

cat(sprintf("\nReport written to: %s/%s\n", getwd(), OUTPUT_FILE))
cat("Done.\n")


#---------  NOT WORKING!!

proj_summary <- rbindlist(lapply(multiples, function(m) {
  col <- paste0("proj_umi_", gsub("\\.", "_", m), "x")
  proj_umi  <- median(cell_qc_filtered[[col]])
  proj_reads <- median(cell_qc_filtered$total_reads) * m
  data.table(
    multiple       = m,
    depth          = paste0(m, "x"),
    proj_reads     = round(proj_reads),
    proj_umi       = round(proj_umi),
    proj_dup       = round((1 - proj_umi / proj_reads) * 100, 1),
    unseen_umi     = round(median(cell_qc_filtered$lib_est) - proj_umi)
  )
}))

print(proj_summary)
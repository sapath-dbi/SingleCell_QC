# SingleCell_QC
Process the molecule_info.h5 to assess Cell, UMI, Dduplication and Library Saturation

Needs setup instructions and R libraries required.

But just run the R script (and load required libraries) 


## scRNA-seq QC Report — 10x Chromium molecule_info.h5

 Produces a self-contained HTML report assessing:
   - Cell capture and UMI distributions
   - PCR duplication rate
   - Library saturation and top-up sequencing estimates

 Works with any species, tissue, or sample type.

 Requirements:
   install.packages(c("data.table", "ggplot2", "patchwork",
                      "scales", "rmarkdown"))
   BiocManager::install("DropletUtils")

 Usage:
   1. Load your molinfo object from the h5 file:
        library(DropletUtils)
        molinfo <- read10xMolInfo("path/to/molecule_info.h5")
   2. Edit the CONFIG section for your sample
   3. Source this script:
        source("scrna_qc_report.R")
   4. An HTML report is written to OUTPUT_FILE (config) default: scrna_qc_report.html
 =============================================================================

``` bash
For Research Use Only. Not for use in diagnostic procedures.
```

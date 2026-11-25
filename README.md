# neoadjuvant_chemoimmunotherapy

This repository contains the custom R (and Python via `reticulate`) code used for the
single-cell and translational analyses in:

> **Longitudinal single-cell profiling identifies histology-specific features associated
> with response to neoadjuvant chemoimmunotherapy in NSCLC.**

The study is an investigator-initiated, prospective, multicenter, single-arm clinical trial
(ClinicalTrials.gov [NCT05894889](https://clinicaltrials.gov/study/NCT05894889)) of
neoadjuvant pembrolizumab plus platinum-doublet chemotherapy in stage II–IIIB NSCLC.

---

## Contents

| Script | Purpose |
|---|---|
| `BEYOND_T_cell_type.R` | Converts sample × subcluster counts to proportions and builds the pre- to post-treatment alluvial flow of T-cell subpopulations; calls Scanpy/PHATE/Palantir through `reticulate`. |
| `find_co_change_cell_type_based_on_trajectory.R` | Identifies cell types that change together along the inferred trajectory and ranks their associations (with P values). |
| `forestplot.R` | Univariable logistic-regression forest plot of baseline clinical and PD-L1 variables against major pathological response (MPR). |
| `gene_search_by_similiar_trajectory.R` | Searches for genes whose expression trajectories parallel the change of a target cell type, using weighted and winsorized correlations. |
| `remodeling_index.R` | Computes the tumor-remodeling index against the continuous pathological response, with pre/post reshaping and `ComplexHeatmap` visualization. |
| `OAK_two_group.R` | External validation in the OAK cohort: CAF GSVA scoring, two-group comparisons, and survival (PFS/OS) analyses. |

---

## 1. System requirements

- **R version 4.1.2** (the code is expected to run on recent R ≥ 4.1 on Linux or macOS).
- A working **Python 3** environment used through `reticulate` (see Installation).
- No special hardware is required. A standard workstation with ≥16 GB RAM is recommended;
  loading the full single-cell expression matrices may benefit from additional memory.

### R packages

CRAN packages:

```
broom, circlize, dplyr, forcats, ggalluvial, ggplot2, ggpubr, ggrepel, ggsci, ggsurvfit,
gridExtra, matrixStats, mgcv, pheatmap, purrr, readxl, reticulate, scales, SHAPforxgboost,
slingshot, stringr, survival, survminer, tibble, tidyr, tidyverse, viridis, xgboost
```

Bioconductor packages:

```
ComplexHeatmap, GSVA, limma
```

Version-validated components in the manuscript: R **4.1.2**, Seurat **5.3.0**, Slingshot
**2.2.0**, randomForest **4.7.1.1**, survminer **0.5.0**, ggplot2 **3.5.2**, ggrepel
**0.9.6**, decoupler **1.6.0**. Other packages can be taken from a contemporary
CRAN/Bioconductor release; no non-standard hardware is required.

### Python components (via `reticulate`)

`scanpy` **1.9.1**, `anndata`, `pandas`, `numpy`, `phate`, and `palantir`. Raw FASTQ
processing used Cell Ranger **6.1.2** (10x Genomics), which is not part of this repository.

---

## 2. Installation guide

1. Install R 4.1.2 (or a compatible newer release) from https://cran.r-project.org.
2. Install the package dependencies:

```r
install.packages(c("broom","circlize","dplyr","forcats","ggalluvial","ggplot2","ggpubr",
  "ggrepel","ggsci","ggsurvfit","gridExtra","matrixStats","mgcv","pheatmap","purrr",
  "readxl","reticulate","scales","SHAPforxgboost","slingshot","stringr","survival",
  "survminer","tibble","tidyr","tidyverse","viridis","xgboost"))

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap","GSVA","limma"))
```

3. Configure the Python environment used by `reticulate` (a conda environment is recommended):

```r
library(reticulate)
# reticulate::py_install(c("scanpy==1.9.1","anndata","pandas","numpy"), pip=TRUE)
# reticulate::py_install(c("phate","palantir"), pip=TRUE)
use_condaenv("<your-env-name>")   # or point to your Python interpreter
```

**Typical install time:** approximately 15–30 minutes on a normal desktop, depending on
download speed and the number of packages already present.

---

## 3. Demo

The primary single-cell and TCR data are held under **controlled access** (GSA accession
[HRA019014](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA019014)) and therefore cannot be
redistributed here; a public demo dataset is consequently not bundled in this repository.
After obtaining access (or substituting your own equivalently formatted data):

- run any individual script with `Rscript <script>.R` (or source it interactively);
- expected outputs are the figures and summary tables reported in the manuscript, including
  the alluvial flow, forest plot, trajectory/remodeling heatmaps, and OAK survival curves.

**Expected run time:** from seconds to a few minutes per script on a normal desktop once the
input data and packages are available; the `reticulate` Python steps may take longer on the
first run.

---

## 4. Instructions for use

1. Obtain the controlled-access data (HRA019014) or prepare matched input files with the same
   column/object structure referenced in each script.
2. **Edit the input paths before running.** The scripts contain the original authors' local
   paths (for example `/home/zwj/...`); replace these with the locations of your downloaded
   data and annotation objects.
3. Source the scripts in an R session or run them with `Rscript`. The alluvial and trajectory
   steps additionally require the configured `reticulate` Python environment above.
4. Random seeds are fixed for the stochastic analyses; the statistical tests and
   multiple-testing settings (Benjamini–Hochberg FDR) are documented in the Methods section.

### Data availability

Single-cell and TCR data are deposited in the Genome Sequence Archive under accession
**HRA019014** (controlled access). TCGA data are from the [GDC](https://gdc.cancer.gov); the
OAK cohort is from [EGAS00001005013](https://ega-archive.org/studies/EGAS00001005013) under a
restricted license. Other data are available from the corresponding author upon reasonable
request.

### Code availability

All custom code is available in this repository at
https://github.com/zwj-tina/neoadjuvant_chemoimmunotherapy.

### License

<!-- TODO: add a LICENSE file; an OSI-approved license such as MIT is recommended. -->

### Citation

Please cite the manuscript above when using or adapting this code.

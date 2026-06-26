# LRGS-cervical-cancer-code

This repository contains the R scripts used for the bioinformatic analyses in the manuscript:

**Integrative single-cell analysis unveils a lactylation-based gene signature for prognosis and treatment prediction in cervical cancer**

The scripts correspond to the single-cell analysis, TCGA-CESC prognostic model construction, tumor immune microenvironment analysis, and IMvigor210 immunotherapy validation.

## Repository structure

```text
LRGS-cervical-cancer-code/
├── README.md
└── scripts/
    ├── 01_single_cell_LRG_screening.R
    ├── 02_TCGA_LRGS_model_nomogram.R
    ├── 03_TME_immune_infiltration.R
    └── 04_IMvigor210_immunotherapy_validation.R
```

## Script descriptions

| Script | Main analysis | Related manuscript figures |
|---|---|---|
| `01_single_cell_LRG_screening.R` | Single-cell RNA-seq processing of GSE279998, quality control, dimensional reduction, clustering, SingleR/manual cell annotation, lactylation activity scoring, high/low lactylation group comparison, LRG screening, GSVA, GO, and KEGG enrichment analyses. | Figure 3; Supplementary Figure 1 |
| `02_TCGA_LRGS_model_nomogram.R` | TCGA-CESC bulk RNA-seq processing, univariate Cox regression, LASSO regression, random forest feature selection, 5-gene LRGS model construction, risk score calculation, survival analysis, time-dependent ROC analysis, clinical Cox regression, nomogram construction, and functional annotation of the five model genes. | Figure 4; Figure 9; Supplementary Figure 2A |
| `03_TME_immune_infiltration.R` | Tumor immune microenvironment analysis using CIBERSORT, MCPcounter, EPIC, xCell, quanTIseq, ESTIMATE, TIMER, ssGSEA, and IPS-related analyses; comparison between high- and low-risk groups. | Figure 7; Figure 8B; Supplementary Figure 3 |
| `04_IMvigor210_immunotherapy_validation.R` | External immunotherapy validation using the IMvigor210 cohort, including LRGS risk score calculation, anti-PD-L1 treatment response analysis, survival analysis, and ROC analysis for treatment response prediction. | Figure 8D-F; Supplementary Figure 4 |

## Data sources

- Single-cell RNA-seq dataset: GEO GSE279998
- Bulk transcriptome and clinical information: TCGA-CESC cohort
- Immunotherapy validation cohort: IMvigor210

## Notes

Some file paths in the scripts reflect the local analysis environment used during manuscript preparation. Before rerunning the scripts, please modify the working directories and input file paths according to your local environment.

Some analyses described in the manuscript, such as experimental validation of core gene expression, mutation/TMB analysis, drug sensitivity prediction, TIDE analysis, and decision curve analysis, may require additional scripts, laboratory data processing, or web-based outputs that are not included in this minimal code package.

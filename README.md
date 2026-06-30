# **Integrative Multi-faceted Proteomics: Deconvoluting Expression, Stability, and Redox Dynamics in the Interferon Response**

Publication: TBA.

# Overview

This repository contains the R-based bioinformatics pipeline used to process, analyze, and visualize data generated via the PISA-REX mass spectrometry methodology. The pipeline integrates protein expression, thermal stability, and redox dynamics to characterize the interferon (IFN) response on a system-wide scale.

## Primary Data Analysis Scripts:

*[pisaRex_preprocessing.R](Data_processing_scripts/pisaRex_preprocessing.R):* Isolates and preprocesses the replicate data from the pair of interest for down-stream analysis.

*[pisaRex_pairwiseReport.R](Data_processing_scripts/pisaRex_pairwiseReport.R):* Provides a statistical/differential regulation report for the pair of interest. Calculated using the empirical Bayes method (limma).

*[GSEA_Ranked_List.R](Data_processing_scripts/GSEA_Ranked_List.R):* Outputs the pre-ranked lists used for GSEA analysis.

*[DIABLO.R](Data_processing_scripts/DIABLO.R):* Performs a multiblock (s)PLS-DA analysis using the mixOmics package.

## Figures:

All R scripts utilized to generate the primary manuscript visualizations are located in the Figures/ directory.

## Main Dependencies:

The analysis was performed in R, leveraging the tidyverse syntax alongside standard Bioconductor tools. Main packages utilized include:

```
Tidyverse 2.0.0
limma 3.66.0
mixOmics 6.32.0
ComplexHeatmap 2.26.1
dbscan 1.2.4
igraph 2.2.3
```

# **Integrative Multi-faceted Proteomics: Deconvoluting Expression, Stability, and Redox Dynamics in the Interferon Response**

R code to reproduce the main analysis and figures of the study.

Publications TBA.

## Primary Data Analysis Scripts:

*pisaRex_preprocessing.R:* Isolates and preprocesses the replicate data from the pair of interest for down-stream analysis.

*pisaRex_pairwiseReport.R:* Provides a statistical/differential regulation report for the pair of interest. Calculated using the empirical Bayes method (limma).

*DIABLO.R:* Performs a multiblock (s)PLS-DA analysis using the mixOmics package.

*GSEA_Ranked_List.R:* Outputs the pre-ranked lists used for GSEA analysis.

## Figures:

Code used to generate the main figures of the paper can be found in the "Figures" folder.

## Main R Packages used:

```
Tidyverse 2.0.0
limma 3.66.0
mixOmics 6.32.0
ComplexHeatmap 2.26.1
dbscan 1.2.4
igraph 2.2.3
```

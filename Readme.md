[![R-CMD-check](https://github.com/jrybarczyk/levi/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jrybarczyk/levi/actions/workflows/R-CMD-check.yaml)
![](https://bioconductor.org/shields/availability/release/levi.svg)
![](https://bioconductor.org/shields/downloads/release/levi.svg)
![](https://bioconductor.org/shields/years-in-bioc/levi.svg)
![](https://bioconductor.org/shields/build/release/bioc/levi.svg)
![](https://bioconductor.org/shields/dependencies/release/levi.svg)


# levi - Landscape Expression Visualization and Network Inference

**Authors:** José R. Pilan, Agnes A. S. Takeda, Jose L. Rybarczyk-Filho

**Maintainer:** José L. Rybarczyk-Filho

## Summary
The integration of biological data emerges as a powerful tool for extracting information capable
of describing the biological system in greater detail. The use of transcriptomic data has become
ubiquitous among researchers in addressing various biological questions. In the past decade, we
have witnessed the integration of biological network analysis to enrich and complement these
responses. This study proposes the synergistic integration of transcriptomic data and networks,
adopting an approach analogous to the Heatmap technique. The Landscape Expression
Visualization Interface (LEVI) is an open-source package developed for the R environment,
aiming to provide an enhanced visualization of gene expression projection onto a biological
network.

## Installation
`levi` is accessible on the bioconductor.org platform.

To install this package, start R and enter:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("levi")
```

### Alternative installation method using devtools
If you prefer, you can also install `levi` using `devtools`. This method might be useful if you want to install the development version directly from a GitHub repository or another source. To install using `devtools`, you'll first need to ensure that `devtools` is installed. If it's not, you can install it using the following command:

```r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("jrybarczyk/levi")
```

## Overview

**`levi`** (**Landscape Expression Visualization Interface**) is an R package developed to enable the visualization of gene expression projections on a biological network. It leverages two main modes of interaction: a Graphical User Interface (GUI) powered by the Shiny package for an accessible, user-friendly experience, and script-based usage for advanced analysis and automation. 

The GUI mode is designed to make the package approachable for users who may not be familiar with coding in R or prefer a visual approach to data analysis. It allows users to upload data, adjust visualization parameters, and interact with the results intuitively. This accessibility facilitates the exploration and interpretation of gene expression data within biological networks.

For users experienced in R or those requiring quicker, possibly automated analyses, `levi` can be operated through scripts. This method offers greater flexibility and customization, making it suitable for integrating `levi` into broader data analysis workflows. Scripting can significantly save time, especially when dealing with large datasets or conducting batch analyses.

`levi` requires two files for use: 
- A file containing the expression levels of the genes (microarray , RNA-seq, Single Cell data). [See example](https://github.com/jrybarczyk/levi/blob/devel/inst/extdata/expression.dat).
- A file containing the biological network. [See example](https://github.com/jrybarczyk/levi/blob/devel/inst/extdata/medusa.dat).

## Files

### Gene Expression Levels

This file should contain the genes of interest, previously normalized by the user. The expression file must have a column with gene identification (Gene Symbol, Entrez, etc.) and at least one column with gene expression levels (treatment, case, control, etc.) [see example](https://github.com/jrybarczyk/levi/blob/devel/inst/extdata/expression.dat). The user can compare expression levels between samples if there are more columns containing these data.

If the expression file does not have values for all genes in the network, a message will be displayed showing a log file path to a temporary directory with gene names. Missing measurements receive an input score of 0.5 and are recorded in `result$metadata$missing_genes`. This marks unavailable information; it does not establish absence of biological change. Smoothing also mixes neighbouring signals.

Datasets of gene expression can be obtained from online databases:
- [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/)
- [Array Express](https://www.ebi.ac.uk/arrayexpress/)
- [The Cancer Genome Atlas (TCGA)](https://cancergenome.nih.gov/)
- [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra)


### Biological Network

The user builds the biological network with a tool such as Cytoscape, RedeR or
Medusa, or passes an edge list directly with `leviFromEdges()` or a STRING
query with `leviFromSTRING()`. Interaction data can be obtained from online
repositories:
- [STRING database](https://string-db.org/)
- [STITCH database](http://stitch.embl.de/)
- [StarBase](http://starbase.sysu.edu.cn/)
- [miRBase](https://mirbase.org/)
- [lncRNAdb](http://www.lncrnadb.org/)

Supported file formats (`fileTypeInput`):

| Code | File type |
|------|-----------|
| dat  | Medusa (DAT), [example](https://github.com/jrybarczyk/levi/blob/devel/inst/extdata/medusa.dat) |
| dyn  | RedeR (DYN) |
| net  | Pajek (NET) |
| stg  | STRING / STITCH (coordinates file plus interactions file) |


## Viewing Modes

### Graphical User Interface (GUI)

The GUI is built with Shiny and can be launched inside RStudio or in the
default web browser. The heavy computation runs in C++ through Rcpp, so the
interface stays responsive while the landscape is built.

```r
library(levi)
LEVIui(browser=TRUE)  # Launch Levi to Browser.
LEVIui(browser=FALSE) # Launch Levi to R environment.
```

The **File** tab loads the network and expression files and selects the gene
identifier, test and control columns. The **Settings** tab controls contrast,
resolution, smoothing, zoom, the colour palette (multicolor or one of the
two-colour sets), contour lines, the signal mode and the permutation test,
which outlines and labels the significant regions as in script mode. The
landscape, its 3D surface and the tables of regions and selected genes can be
downloaded. A full walkthrough is in `vignette("levi")`.

### Script

The script mode offers enhanced flexibility for more detailed settings adjustments and comprehensive automation. This mode is particularly well-suited for advanced users who require precise control over their analysis parameters and workflows. Additionally, it is ideal for processing large datasets in bulk, enabling efficient data handling and analysis customization. This capability makes it an invaluable tool for researchers and data scientists looking to streamline their data analysis pipelines, ensuring both scalability and reproducibility in their work.

```r
library(levi)

template_network <- file.path(system.file(package="levi"),"extdata",
                                "medusa.dat", fsep = .Platform$file.sep)

template_expression <- file.path(system.file(package="levi"),
                                "extdata","expression.dat", 
                                fsep = .Platform$file.sep)

multicolor <- levi(networkCoordinatesInput = template_network,
                expressionInput = template_expression, fileTypeInput = "dat",
                geneSymbolInput = "ID", 
                readExpColumn=
                readExpColumn("TumorCurrentSmoker-NormalNeverSmoker"), 
                contrastValueInput = 50, resolutionValueInput  = 50, 
                zoomValueInput = 50, smoothValueInput = 50, contourLevi = TRUE)

twocolors <- levi(networkCoordinatesInput = template_network,
                expressionInput = template_expression, fileTypeInput = "dat",
                geneSymbolInput = "ID", 
                readExpColumn=
                readExpColumn("TumorCurrentSmoker-NormalNeverSmoker"),
                setcolor = "pink_green", contourLevi = FALSE)

```

The script mode allows the user to compare combinations between two 
experiments in the gene expression levels file. The readExpColumn function 
can be used to this task to inform the combination separating by dash (-) 
and to add more combinations separate by comma (,).

```r
library(levi)
base <- readExpColumn("TumorFormerSmoker-NormalFormerSmoker", 
                        "TumorNeverSmoker-NormalNeverSmoker")

template_network <- file.path(system.file(package="levi"),"extdata",
                                "medusa.dat", fsep = .Platform$file.sep)

template_expression <- file.path(system.file(package="levi"),
                                "extdata","expression.dat", 
                                fsep = .Platform$file.sep)

multicolor <- levi(networkCoordinatesInput = template_network,
                    expressionInput = template_expression, 
                    fileTypeInput = "dat",
                    geneSymbolInput = "ID", readExpColumn= base, 
                    contrastValueInput = 50, resolutionValueInput  = 50, 
                    zoomValueInput = 50, smoothValueInput = 50, 
                    contourLevi = FALSE)

twocolors <- levi(networkCoordinatesInput = template_network,
                expressionInput = template_expression, fileTypeInput = "dat",
                geneSymbolInput = "ID", 
                readExpColumn= base,
                setcolor = "pink_green", contourLevi = FALSE)

```
More examples are in `vignette("levi")`, `vignette("levi_inference")`,
`vignette("levi_validation")` and the scripts in `inst/scripts/`.


## Signal contracts and reproducibility

Signal interpretation: `ratio` preserves Test/(Test + Control), with no min-max
rescaling; equal nonzero inputs give 0.5. `expressionLog = TRUE` back-transforms
log2 inputs only in this mode. A single ratio column means abundance/(abundance + 1),
not a comparison with a control. `logfc` accepts two log-scale columns or one
already computed logFC, mapping zero to 0.5. `zscore` centres on the mean logFC
of measured network support points, not on biological absence of change.
Missing measurements are assigned 0.5 and listed in `result$metadata`.
Gaussian smoothing mixes neighbouring signals, so these baseline statements
apply to the input signals and to uniformly neutral networks.

Permutation inference is conditional on the fixed network and layout. Measured
gene values are shuffled as pairs and edge midpoints are recalculated; missing
positions stay fixed. Both tails over occupied cells form one multiple-testing
family per comparison. `p_adjust_method = "BY"` is the default; contours and
`result$pvalues` use adjusted values. `result$raw_pvalues` retains the unadjusted
values. `perm_side` selects the displayed tails, without changing that family.
This is not a test of differential expression between biological replicates;
increasing `n_perm` alone does not validate inferential use. Calibration across
networks, layouts and missingness patterns still requires simulation studies.

Call `set.seed()` before permutation runs. `result$metadata` records the RNG
state, network, coordinates, signal mode, grid settings and software versions.
`leviDiff()` rejects incompatible metadata or grid coordinates. For DESeq2,
edgeR and limma adapters, select the logFC column against itself with
`signal_mode = "logfc"`; abundance annotations are not control measurements.
For edgeR, select the contrast in `glmLRT()` or `glmQLFTest()` before calling
`leviFromEdgeR()` on the resulting test object. KEGG uses the supplied universe,
converting both selected genes and background to ENTREZID with the same OrgDb.

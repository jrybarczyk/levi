# levi usage examples

A series of self-contained scripts, from loading the data to interpreting the
results. They all use only the datasets shipped with the package, so they run
with no download and no preparation.

To run a script:

```r
scripts <- system.file("scripts", package = "levi")
source(file.path(scripts, "01-getting-started.R"), echo = TRUE)
```

Or, to open and read it:

```r
file.edit(file.path(system.file("scripts", package = "levi"),
                    "01-getting-started.R"))
```

## The scripts

| # | Script | What it covers |
|---|--------|----------------|
| 01 | `01-getting-started.R` | The two input files, the minimal call, and the anatomy of the object returned by `levi()`. |
| 02 | `02-reading-the-landscape.R` | The scale and the neutral point 0.5, the silhouette, and how to read `$scores` and `$peaks`. |
| 03 | `03-topology-and-landscape.R` | Five topologies under the same parameters: what comes from the expression and what comes from the shape of the network. |
| 04 | `04-signal-modes.R` | `signal_mode` (`ratio`, `logfc`, `zscore`), `logfc_k`, `expressionLog` and single-column input. |
| 05 | `05-batch-and-comparison.R` | Several comparisons in one call, `leviGrid()` and `leviDiff()`. |
| 06 | `06-significance.R` | Legacy per-cell permutation test (`inference_unit = "cell"`), p-value matrices, contours and the `1/(n_perm+1)` floor. The default regional test is covered in `vignette("levi_inference")`. |
| 07 | `07-bioconductor-data.R` | `SummarizedExperiment`, `ExpressionSet`, DE tables and STRING networks. |
| 08 | `08-visual-parameters.R` | Resolution, smoothing, contrast, zoom, palettes and exporting. |
| 09 | `09-reliability-validation.R` | Numerical regression diagnostics for the cell-wise test and the node-score contracts. |
| 10 | `10-calibration.R` | **Offline, about an hour.** Type I error, power and layout sensitivity by simulation; writes the summaries in `inst/extdata/validation/` used by `vignette("levi_validation")`. |
| 11 | `11-airway-preprocess.R` | **Offline, needs network access.** DESeq2 on `airway` and the STRING query; writes `inst/extdata/airway/` used by `vignette("levi_validation")`. |

The order is progressive, but each script works on its own.

## Integration tests on real data

`tests/testthat/test_integration_real_data.R` replays the analyses of the
supplement (airway RNA-seq, GSE10072 microarray with its smoking strata and
the KEGG focal-adhesion network, Kang 2018 single-cell monocytes with the
STRING and KEGG JAK-STAT networks) and checks that their conclusions still
hold: which region comes first, whether it is significant, exact p-values
where the test enumerates every arrangement, cluster sizes and lead genes.
They are skipped by default. To run them from the package root:

```sh
LEVI_INTEGRATION=1 LEVI_REAL_DATA=/path/to/real_data_tests \
    Rscript -e 'testthat::test_local(filter = "integration")'
```

`LEVI_REAL_DATA` must hold `GSE10072_series_matrix.txt.gz`, `GPL96.annot.gz`,
`kang18.rds` (muscData `Kang18_8vs8` as a `SingleCellExperiment`) and the
saved networks `gse10072*_{nodes,edges}.tsv`, `kang*_{nodes,edges}.tsv`. The
airway case needs only the files shipped in `inst/extdata/airway/`. The whole
file takes a few minutes; the Kang cases are the slow part.

## If you are short on time

Three ideas account for most of the correct use of the package:

**The scale is absolute.** `LandscapeScore` lives in `[0, 1]` and 0.5 means "no
change between test and control". This holds regardless of the smoothing and of
the network density, and it is what makes two landscapes comparable side by
side. The `flat` dataset, with no variation at all, gives exactly 0.5 — worth
using as a sanity check whenever a result looks odd. See script 02.

**The signal mode is the decision that matters most.** Data on a linear scale
(counts, TPM, FPKM) call for `ratio`; log2 data (microarray RMA, VST/rlog,
proteomics, any ready-made logFC) call for `logfc`. Using `ratio` on log2 data
is the most common mistake and compresses the real differences. See script 04.

**Neighbourhood is information.** An altered gene surrounded by genes altered in
the same direction produces a broad region; an isolated gene barely shows. That
is what the landscape adds to a ranked gene list, and it is why the network
layout is an analytical choice, not merely an aesthetic one. See script 03.

## Graphical interface

Everything here is also available in the GUI:

```r
LEVIui(browser = TRUE)   # in the browser
LEVIui(browser = FALSE)  # in the RStudio Viewer pane
```

## Datasets used

| Files | Description |
|-------|-------------|
| `hub_*` | Star of 9 nodes: over-expressed core, repressed periphery. |
| `gradient_*` | Linear chain of 6 nodes, monotonically increasing expression. |
| `bimodal_*` | Two opposing modules joined by a bridge. |
| `flat_*` | Null control: no variation between the conditions. |
| `sparse_*` | 15 nodes in the network, only 5 with a measured value. |
| `logfc_*` | Six nodes with values already on the log2 scale. |
| `hub_multicomp_expression.dat` | Three conditions over the `hub` network, for batch mode. |
| `medusa.dat` + `expression.dat` | Real network of 30 nodes and 325 interactions, with lung tumour expression data. |

Each dataset has a documentation page with the expected results:
`?hub_dataset`, `?gradient_dataset`, `?bimodal_dataset`, `?flat_dataset`,
`?sparse_dataset`, `?logfc_dataset`, `?multicomp_dataset`.

`09-reliability-validation.R` runs small null simulations and reports sensitivity
to smoothing, layout and missing measurements. Pass an output directory to keep
the CSV diagnostics. It is a numerical regression diagnostic, not a calibration
study establishing inferential validity.

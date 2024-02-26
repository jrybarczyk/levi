![](https://bioconductor.org/shields/availability/release/levi.svg)
![](https://bioconductor.org/shields/downloads/release/levi.svg)
![](https://bioconductor.org/shields/years-in-bioc/levi.svg)
![](https://bioconductor.org/shields/build/release/bioc/levi.svg)
![](https://bioconductor.org/shields/dependencies/release/levi.svg)


# levi - Landscape Expression Visualization Interface

**Authors:**  José R. Pilan, Isabelle M. da Silva, Agnes A. S. Takeda, Jose L. Rybarczyk-Filho  
**Maintainer:** José L. Rybarczyk-Filho

## Installation
`levi` is accessible on the bioconductor.org platform.

To install this package, start R and enter:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("levi")

```
### Alternative installation method using devtools
If you prefer, you can also install `levi` using `devtools`. This method might be useful if you want to install the development version directly from a **GitHub repository** or another source. To install using `devtools`, you'll first need to ensure that `devtools` is installed. If it's not, you can install it using the following command:

```r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("jrybarczyk/levi")

```


## Overview

**`levi`** **(Landscape Expression Visualization Interface)** is an R package developed to enable the visualization of gene expression projections on a biological network. It is based on two other software: Viacomplex (This tool is no longer supported)  and Galant (This tool is no longer supported), which corresponds to a plugin for Cytoscape 2.x software.

Two files are required to use `levi`: 
- The file containing the expression levels of the genes.
- A file containing the biological network.

## Files

### Gene Expression Levels

This file should contain the genes of interest, previously normalized by the user. The expression file must have a column with gene identification (Gene Symbol, Entrez, etc.) and at least one column with gene expression levels (treatment, case, control, etc.) [see example](https://github.com/jrybarczyk/levi/blob/devel/inst/extdata/expression.dat). The user can compare expression levels between samples if there are more columns containing these data.

If the expression file does not have values for all genes in the network, a message will be displayed showing a log file path to a temporary directory with gene names. In the **landscape** construction, genes with no expression value will be displayed with values close to 0.5, demonstrating that there were no changes (down-regulated or up-regulated genes).

Data sets of gene expression can be obtained from online databases:
- [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/)
- [Array Express](https://www.ebi.ac.uk/arrayexpress/)
- [The Cancer Genome Atlas (TCGA)](https://cancergenome.nih.gov/)
- [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra)

### Biological Network

The `levi` supports several extensions of biological network files (\*.net, \*.dyn, \*.txt, \*.dat). The user should build the biological network using specific tools such as Cytoscape, RedeR, Medusa, etc. It is recommended to obtain interaction data/biological associations from online repositories:
- [STRING database](https://string-db.org/)
- [StarBase](http://starbase.sysu.edu.cn/)
- [miRBase](http://www.mirbase.org/)
- [lncRNAdb](http://www.lncrnadb.org/)
- [HTRIdb](http://www.lbbc.ibb.unesp.br/htri)

| Extension | File Type | Link |
|-----------|-----------|-----|
| dat       | Medusa (DAT) | [Example DAT](https://github.com/jrybarczyk/levi/blob/devel/inst/extdata/medusa.dat) |
| dyn       | RedeR (DYN)  | [Example DYN](https://github.com/jrybarczyk/levi/blob/devel/inst/extdata/RedeR.dyn) |
| net       | Pajek (NET)  | [Example NET](https://github.com/jrybarczyk/levi/blob/devel/inst/extdata/Pajek.net) |
| stg       | STRING / STITCH | [Example coordinates](https://github.com/jrybarczyk/levi/blob/devel/inst/extdata/string_network_coordinates.txt) |
| stg       | STRING / STITCH | [Example interactions](https://github.com/jrybarczyk/levi/blob/devel/inst/extdata/string_interactions.tsv)|

## Viewing Modes

`levi` has two viewing modes: **Graphical User Interface (GUI)** and **script**.

### Graphical User Interface (GUI)

The GUI mode was developed using the Shiny package (Figure 1). This viewing mode can be used in the R environment (`browser=FALSE`) or in the user's default operating system browser (`browser=TRUE`).

```r
library(levi)
LEVIui(browser=TRUE)  # Launch Levi to Browser.
LEVIui(browser=FALSE) # Launch Levi to R environment.
```
<center><img src="https://github.com/jrybarczyk/levi/blob/master/vignettes/levi.png" width=150%></center>
<center> <b>Figure 2: </b> Graphical User Interface (GUI).</center>

* **1** - **Tabs**
    + **File** - General options for uploading the data files.
    + **Settings** - Options to build the **landscape**. 
* **2** - **Network input type:** - Selection of the biological network input 
format. The options are available on the table 1.
* **3** - **Upload the network file** - Input button to select the biological 
network file.
* **4** - **Upload the expression file** - Input button to select the gene 
expression levels file.
* **5** - **Expression values in log scale** - 
select when there is low variation between expression data
* **6** - **Selected fields** - Selecting the type of analysis. Comparison 
between two samples or a single sample.
* **7** - **Select gene symbol field** -  Selection of the gene ID contained 
in the gene expression levels file.  The gene ID in the biological network 
file must be the same in the gene expression levels file. **Select test field**
-  Select the case/test sample. **Select control field** - Select the 
control sample. If **levi** detects the case/test sample and control sample as
equal, then **levi** will apply the "single sample" analysis.
* **8** - **Chart Colors** - Color palette to build the **landscape**. There is 
two options the **Multicolor** has 20 color levels combined. The 
**Two colors** has two categories of color palettes, multicolor and two colors. 
For the selection of two colors the available options are: *purple_pink*, 
*green_blue*, *blue_yellow*, *pink_green*, *orange_purple*, *green_marine*. 
Both of them is targeted to show the biological network with higher expression. 
(Figure 3 and Figure 4).
* **9** - **Chart with contour** - Enable or disable the contour lines in 
the **landscape**.
* **10** -  **landscape** Image display area. The user can select specific 
areas of the image to inspect the gene name and position.
* **11** - **Expression area:** - Total expression value of the selected 
area in the image.
* **12** - **Download options**
    + **File format to landscape** - Select the output format of the
    **landscape**. The options available are: TIFF, BMP,  JPEG and PNG. We 
    recommend using this option in the browser.
    + **Download Plot** - Button to save the file.
* **13** -  Visualization of gene name and expression value of the selected 
area in the **landscape**. The user can save the table in a csv file 
(button **Download Data**)
* **14** - **Settings**
    + **Contrast** - Constrast value in the **landscape**. The variable 
    range is 1 to 100. The default value is 50.
    + **Resolution** -  Image  size of the **landscape**. The variable 
    range is 1 to 100. The default value is 50. If this parameter is higher, 
    then the total time required will be longer.
    + **Smoothing** - Smoothing of the **landscape**. The variable range 
    is 1 to 100. The default value is 50.I f this parameter is higher, then 
    the total time required will be longer.
    + **Zoom** - Zoom value for the **landscape**. The variable range is 
    1 to 100. The default value is 50.

## Script

The `levi` scripting mode also has the settings for **landscape** building 
(see example below).

```r
library(levi)

template_network <- file.path(system.file(package="levi"),"extdata",
                                "medusa.dat", fsep = .Platform$file.sep)

template_expression <- file.path(system.file(package="levi"),
                                "extdata","expression.dat", 
                                fsep = .Platform$file.sep)

multicolor <- levi(networkCoordinatesInput = template_network,
                expressionInput = template_expression, fileTypeInput = "dat",
                geneSymbolnput = "ID", 
                readExpColumn=
                readExpColumn("TumorCurrentSmoker-NormalNeverSmoker"), 
                contrastValueInput = 50, resolutionValueInput  = 50, 
                zoomValueInput = 50, smoothValueInput = 50, contourLevi = TRUE)

twocolors <- levi(networkCoordinatesInput = template_network,
                expressionInput = template_expression, fileTypeInput = "dat",
                geneSymbolnput = "ID", 
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
                        "TumorNeverSmoker-TumorNeverSmoker")

template_network <- file.path(system.file(package="levi"),"extdata",
                                "medusa.dat", fsep = .Platform$file.sep)

template_expression <- file.path(system.file(package="levi"),
                                "extdata","expression.dat", 
                                fsep = .Platform$file.sep)

multicolor <- levi(networkCoordinatesInput = template_network,
                    expressionInput = template_expression, 
                    fileTypeInput = "dat",
                    geneSymbolnput = "ID", readExpColumn= base, 
                    contrastValueInput = 50, resolutionValueInput  = 50, 
                    zoomValueInput = 50, smoothValueInput = 50, 
                    contourLevi = FALSE)

twocolors <- levi(networkCoordinatesInput = template_network,
                expressionInput = template_expression, fileTypeInput = "dat",
                geneSymbolnput = "ID", 
                readExpColumn= base,
                setcolor = "pink_green", contourLevi = FALSE)

```

[more examples](https://github.com/jrybarczyk/levi/blob/paper/Case/Case.pdf)

## Information about downloads of levi on Bioconductor.org
![](https://github.com/jrybarczyk/levi/blob/devel/miscellaneous/Rplot.png)

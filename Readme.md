# levi - Landscape Expression Visualization Interface

**Authors:** Isabelle M. da Silva, José R. Pilan, Agnes A. S. Takeda, Jose L. Rybarczyk-Filho  
**Package Version:** `r pkg_ver('levi')`  
**Date:** `r doc_date()`  
**Vignette:** [Using levi](#using-levi)  

## Overview

**levi (Landscape Expression Visualization Interface)** is an R package developed to enable the visualization of gene expression projections on a biological network. It is based on two other software: Viacomplex (@Castro2009) and Galant (@Camilo2013), which corresponds to a plugin for Cytoscape software.

Two files are required to use **levi**: 
- The file containing the expression levels of the genes.
- A file containing the biological network.

## Files

### Gene Expression Levels

This file should contain the genes of interest, previously normalized by the user. The expression file must have a column with gene identification (Gene Symbol, Entrez, etc.) and at least one column with gene expression levels (treatment, case, control, etc.). The user can compare expression levels between samples if there are more columns containing these data.

If the expression file does not have values for all genes in the network, a message will be displayed showing a log file path to a temporary directory with gene names. In the **landscape** construction, genes with no expression value will be displayed with values close to 0.5, demonstrating that there were no changes (down-regulated or up-regulated genes).

Data sets of gene expression can be obtained from online databases:
- [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/)
- [Array Express](https://www.ebi.ac.uk/arrayexpress/)
- [The Cancer Genome Atlas (TCGA)](https://cancergenome.nih.gov/)
- [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra)

### Biological Network

The **levi** supports several extensions of biological network files (\*.net, \*.dyn, \*.txt, \*.dat). The user should build the biological network using specific tools such as Cytoscape, RedeR, Medusa, etc. It is recommended to obtain interaction data/biological associations from online repositories:
- [STRING database](https://string-db.org/)
- [StarBase](http://starbase.sysu.edu.cn/)
- [miRBase](http://www.mirbase.org/)
- [lncRNAdb](http://www.lncrnadb.org/)
- [HTRIdb](http://www.lbbc.ibb.unesp.br/htri)

| Extension | File Type | Link |
|-----------|-----------|-----|
| dat       | Medusa (DAT) | [Example DAT](https://bit.ly/2sj6RJh) |
| dyn       | RedeR (DYN)  | [Example DYN](https://bit.ly/2IJGaZ7) |
| net       | Pajek (NET)  | [Example NET](https://bit.ly/2IKzaLv) |
| stg       | STRING / STITCH | [Example coordinates](https://bit.ly/2xjnaLU) |
| stg       | STRING / STITCH | [Example interactions](https://bit.ly/2J2OToH)|

## Viewing Modes

**levi** has two viewing modes: **Graphical User Interface (GUI)** and **script**.

### Graphical User Interface (GUI)

The GUI mode was developed using the Shiny package (Figure 1). This viewing mode can be used in the R environment (`browser=FALSE`) or in the user's default operating system browser (`browser=TRUE`).

```r
library(levi)
LEVIui(browser=TRUE)  # Launch Levi to Browser.
LEVIui(browser=FALSE) # Launch Levi to R environment.
``

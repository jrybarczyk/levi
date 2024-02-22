---
title: '**levi**: an R package for landscape expression visualization' 

tags: 
 - Biological Networks 
 - Gene Expression Levels 
 - R

authors: 
- name: José Rafael Pilan 
  affiliation: "1" 
  orcid: 0000-0002-9403-5030 
- name: Isabelle M. da Silva 
  affiliation: "1" 
  orcid: 0000-0003-1435-6545 
- name: Agnes A. Sekijima Takeda 
  affiliation: "1" 
  orcid: 0000-0001-5000-6652 
- name: José L. Rybarczyk-Filho 
  affiliation: "1, 2" 
  orcid: 0000-0002-0757-3608 

affiliations: 
 - name: Departament of Biophysics and Pharmacology, Institute of Bioscience of Botucatu, Universidade Estadual Paulista (UNESP) 
   index: 1 
 - name: For correspondence, contact [jose.luiz\@unesp.br](mailto:jose.luiz@unesp.br){.email}

date: "20 february 2024" 
bibliography: paper.bib
---

# Summary

The integration of biological data emerges as a powerful tool for extracting information capable of describing the biological system in greater detail. The use of transcriptomic data has become ubiquitous among researchers in addressing various biological questions. In the past decade, we have witnessed the integration of biological network analysis to enrich and complement these responses. This study proposes the synergistic integration of transcriptomic data and networks, adopting an approach analogous to the Heatmap technique. The Landscape Expression Visualization Interface (levi) is an open-source package developed for the R environment, aiming to provide an enhanced visualization of gene expression projection onto a biological network.

# Introduction

The high throughput technologies have been developed in the last decades for obtention of life sciences knowledge, whose result is the 'omics'-scale data [@gehlenborg2010]. Consequently, it is essential the development of methods and tools to allow the data integration and visualization for better understanding of a biological system, both qualitatively and quantitatively [@kannan2015]. The biological networks allows the visualization of the relationships among several molecules in a pathway. When additional information (e.g. expression data) is integrated to the network it can help to understand which molecules play central roles in complex systems [@charitou2016]. One approach to facilitate the visual identification of patterns in the biological data is its integration to biological networks, which will result in\`colored networks' as smoothed data maps [@castro2009,@camilo2013]. Here, we present the **L**andscape **E**xpression **V**isualization **I**nterface (levi), a package developed for the R environment that integrates data from biological networks with gene expression levels or others features (clustering coeficient, connectivity, stress, etc), and displays a heatmap with surface curves (landscape) to evidence the altered regions. It allows the user to perform the data processing and evaluation in a single environment.

# Implementation and main functions

Levi was developed in R language, using the Shiny package to create an graphical user interface (GUI). Normalized gene expression data can be adjusted to logarithmic scale for better visualization of the expression levels projected onto the network. The network upload requires files in formats Medusa (DAT), RedeR (DYN), Pajek (NET) or STRING / STITCH. The GUI mode allows the user the landscape visualization from a single gene expression data or a comparison of data expression. Moreover, it allows the interaction with the final image, selecting areas in the landscape to identify the genes in a chosen region. Levi settings contrast, resolution, smoothing, zoom and color palettes are also available. The script mode enables the batch execution to create landscapes for different comparisons between sample types for a single pathway/network. Both modes export the landscapes as images in TIFF, BMP, JPEG and PNG formats.


![Landscape analysis to Non-Small Cell Lung Cancer signaling pathway. (A) Biological network for Non-Small Cell Lung Cancer. (B) Landscape for normal lung and never smoked. (C) Landscape for tumor and current smoker. (D) Landscape for tumor and former smoker. (E) Landscape for tumor and never smoked.](figure1.png)

# Case Study

To illustrate the levi application, we built a Non-Small Cell Lung Cancer (NSCLC) network, which was obtained by prospection of the genes from the Non-Small Cell Lung Cancer signaling pathway (map05223) from the Kyoto Enciclopedia of Genes and Genomes(KEGG) database [@kanehisa2019new]. These genes were used as input for the STRING v11 [@szklarczyk2019string] to built the NSCLC network using the following parameters: active interaction sources - Experiments and database, confidence score 0.40 (considered medium confidence) Figure (A). The files "network-coordinates" and "as simple tabular text output" were selected and used as input for levi package. Then, we prepared the expression file. In our example, we selected a gene expression dataset of lung adenocarcinoma and lung tissue from current, former and never smokers (GSE10072) [@landi2008] available in Gene Expression Omnibus (GEO) [@barrett2012]. Before the upload in levi, the data was previously normalized with MAS5 and the probes were mapped with gene symbol using the R packages affy and hgu133a2.db, respectively.

The data processing in levi allowed the visualization of the expression data onto the biological newtork to evidence difference among samples (Figure 1 (B-E)). The comparison among images indicated a difference of colors in one area that corresponds to the  E2F2  and  CDK4 genes, which were more expressed in the samples of lung cancer patients than in the control samples. This result was corroborated by previous studies, which also suggested these genes as potential biomarkers for prognosis in detection of NSCLC [@chen2015e2f2, @wu2011elevated].

# Conclusion

The landscape visualization, provided by 'levi'— a tool available on the Bioconductor.org platform — results from the integration of expression data and the network. This allows for a more comprehensive qualitative analysis of the biological problem, due to the improved visualization of both upregulated and downregulated genes, as well as their interactions with the surrounding network.

# Funding

This study was supported by CNPq 473789/2013-2, 134469/2016-0, 134467/2016-7 and financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior - Brasil (CAPES) - Finance Code 001.

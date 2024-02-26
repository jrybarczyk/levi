---
title: 'LEVI: an R package for landscape expression visualization' 

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
 - name: For correspondence, contact jose.luiz@unesp.br

date: "20 february 2024" 
bibliography: paper.bib
---

# Summary

The integration of biological data emerges as a powerful tool for extracting information capable of describing the biological system in greater detail. The use of transcriptomic data has become ubiquitous among researchers in addressing various biological questions. In the past decade, we have witnessed the integration of biological network analysis to enrich and complement these responses. This study proposes the synergistic integration of transcriptomic data and networks, adopting an approach analogous to the Heatmap technique. The Landscape Expression Visualization Interface (LEVI) is an open-source package developed for the R environment, aiming to provide an enhanced visualization of gene expression projection onto a biological network.

# Introduction

The high throughput technologies have been developed over the last decades have revolutionized the acquisition of life sciences knowledge, resulting in 'omics'-scale data [@Gehlenborg2010]. Consequently, developing methods and tools  for data integration and visualization is essential to qualitatively and quantitatively enhance our understanding of biological systems [@Kannan2016]. Biological networks allow for the visualization of relationships among  molecules within a pathway. When additional information (e.g. expression data) is integrated to the network it can help to understand which molecules play central roles in complex systems [@Charitou2016]. An effective approach for visual identification of patterns in biological data involves its integration into biological networks, resulting in 'colored networks' that serve as smoothed data maps [@Castro2009; @Camilo2013]. We present the **L**andscape **E**xpression **V**isualization **I**nterface (LEVI), a package developed for the R environment that integrates gene expression levels or other features (e.g., clustering coefficient, connectivity, stress) with biological networks. LEVI displays a heatmap with surface curves (landscape) to highlight altered regions, facilitating data processing and evaluation within a single environment.

# Implementation and main functions

LEVI was developed in R  programming language, incorporating the Shiny package to create a graphical user interface (GUI). Normalized gene expression data can be adjusted to logarithmic scale for better visualization of the expression levels projected onto the network. The network upload requires files in formats Medusa (DAT), RedeR (DYN), Pajek (NET) or STRING / STITCH. The GUI mode allows the user the landscape visualization from a single gene expression data or a comparison of data expression. Moreover, it allows the interaction with the final image, selecting areas in the landscape to identify the genes in a chosen region. LEVI settings contrast, resolution, smoothing, zoom and color palettes are also available. The script mode enables the batch execution to create landscapes for different comparisons between sample types for a single pathway/network. Both modes export the landscapes as images in TIFF, BMP, JPEG and PNG formats.


![Landscape analysis to Non-Small Cell Lung Cancer signaling pathway. (A) Biological network for Non-Small Cell Lung Cancer. (B) Landscape for normal lung and never smoked. (C) Landscape for tumor and current smoker. (D) Landscape for tumor and former smoker. (E) Landscape for tumor and never smoked.](figure1.png)

# Case Study

To illustrate the LEVI application, we built a Non-Small Cell Lung Cancer (NSCLC) network, which was obtained by prospection of the genes from the Non-Small Cell Lung Cancer signaling pathway (map05223) from the Kyoto Enciclopedia of Genes and Genomes(KEGG) database [@Kanehisa2019]. These genes were used as input for the STRING v11 [@Szklarczyk2018] to built the NSCLC network using the following parameters: active interaction sources - Experiments and database, confidence score 0.40 (considered medium confidence) Figure (A). The files "network-coordinates" and "as simple tabular text output" were selected and used as input for LEVI package. Then, we prepared the expression file. In our example, we selected a gene expression dataset of lung adenocarcinoma and lung tissue from current, former and never smokers (GSE10072) [@Landi2008] available in Gene Expression Omnibus (GEO) [@Barrett2013]. Before the upload in LEVI, the data was previously normalized with MAS5 and the probes were mapped with gene symbol using the R packages affy and hgu133a2.db, respectively.

The data processing in LEVI allowed the visualization of the expression data onto the biological newtork to evidence difference among samples (Figure 1 (B-E)). The comparison among images indicated a difference of colors in one area that corresponds to the  E2F2  and  CDK4 genes, which were more expressed in the samples of lung cancer patients than in the control samples. This result was corroborated by previous studies, which also suggested these genes as potential biomarkers for prognosis in detection of NSCLC [@Du2022; @Wu2011].

# Conclusion

The landscape visualization, provided by LEVI — a tool available on the Bioconductor.org platform — results from the integration of expression data and the network. This allows for a more comprehensive qualitative analysis of the biological problem, due to the improved visualization of both upregulated and downregulated genes, as well as their interactions with the surrounding network.

# Funding

This study was supported by CNPq 473789/2013-2, 134469/2016-0, 134467/2016-7 and financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior - Brasil (CAPES) - Finance Code 001.


# References

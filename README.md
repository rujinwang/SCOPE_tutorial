# SCOPE
A normalization and copy number estimation method for single-cell DNA sequencing


## Authors
Rujin Wang, Danyu Lin, and Yuchao Jiang


## Maintainer
Rujin Wang <rujin@email.unc.edu>


## Installation
```r
install.packages('devtools')

devtools::install_github("rujinwang/SCOPE/package")
```

## Description
Whole genome single-cell DNA sequencing (scDNA-seq) enables characterization of copy number profiles at the cellular level. This circumvents the averaging effects associated with bulk-tissue sequencing and has increased resolution yet decreased ambiguity in deconvolving cancer subclones and elucidating cancer evolutionary history. ScDNA-seq data is, however, sparse, noisy, and highly variable even within a homogeneous cell population, due to the biases and artifacts that are introduced during the library preparation and sequencing procedure. Here, we propose SCOPE, a normalization and copy number estimation method for scDNA-seq data. The distinguishing features of SCOPE include: (i) utilization of cell-specific Gini coefficients for quality controls and for identification of normal/diploid cells, which are further used as negative control samples in a Poisson latent factor model for normalization; (ii) modeling of GC content bias using an expectation-maximization algorithm embedded in the Poisson generalized linear models, which accounts for the different copy number states along the genome; (iii) a cross-sample iterative segmentation procedure to identify breakpoints that are shared across cells from the same genetic background. We evaluate performance of SCOPE on real scDNA-seq data sets from cancer genomic studies. Compared to existing methods, SCOPE more accurately estimates subclonal copy number aberrations and is shown to have higher correlation with array-based copy number profiles of purified bulk samples from the same patient. We further demonstrate SCOPE on three recently released data sets using the 10X Genomics single-cell CNV pipeline and show that it can reliably recover 1% of the cancer cells from a background of normal.


## Manuscript
Rujin Wang, Danyu Lin, and Yuchao Jiang. SCOPE: A Normalization and Copy Number Estimation Method for Single-Cell DNA Sequencing. ***bioRxiv***, 2019. ([link](https://www.biorxiv.org/content/10.1101/594267v1))


## Pipeline overview
The flowchart below outlines the procedures for profiling single-cell CNV. The first step is bioinformatic pre-processing pipeline. Assembled BAM files are finally required. The second step is calculation of GC content, mappability and read depth using Rsamtools with QC measurements. An EM embedded normalization procedure is then applied to single cells to remove biases and artifacts along the whole genome. The cross-sample Poisson likelihood segmentation is performed to call CNVs, which can be further used to infer single-cell clusters or clones.

<p align="center">
  <img src='https://github.com/rujinwang/SCOPE_tutorial/blob/master/demo/Figure1_SCOPE_outline.jpg' width='250' height='400'>
</p>

## Running SCOPE
**Tutorial** is available [**here**](http://htmlpreview.github.io/?https://github.com/rujinwang/SCOPE_tutorial/blob/master/SCOPE_tutorial.html)


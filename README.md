
# 

## MDTToolset

### THIS R PACKAGE IS CURRENTLY UNDER DEVELOPMENT

### Overview

The [R](https://www.r-project.org) package **MDTToolset** is designed to
calculate the **M**ost **D**ominant **T**ranscript and their switches in
the RNA-Seq datasets. It enables users to:

- Remove redundant transcripts with identical protein sequences and
  aggregate their transcript counts.
- Calculate the most dominant transcripts in RNA-Seq data using flexible
  cutoffs, such as user-defined transcript expression and enrichment
  cutoffs.
- Identify MDTs present in a user-defined percentage of samples.
- Calculate MDT Switching events with flexible user-defined cutoffs.
- Integrate the isoform interaction network with the MDT switches to
  understand their functional impact on the protein interaction network.
- Visualize the expression of MDTs and disease-specific MDTs (dMDTs).
- Visualize the network of dMDTs.

### Background

Alternative splicing plays an essential role in development, tissue
specificity and essential cell functions. One gene can encode for many
transcripts, one of those transcripts might be expressed at a
significantly higher level than the other transcripts. These transcripts
are defined as most dominant transcripts (MDTs). In cancer, those
canonical MDTs might be switched to other transcripts which we call MDT
switching event. The *MDTToolset* aims to find the MDTs in RNA-seq data
and evaluated the MDT switching event in a flexible user-defined way.
The user can define cutoffs for the transcript expression, and MDT
enrichment.

### How to Install

You can install the development version of MDTToolset from
[GitHub](https://github.com/KarakulakTulay/MDTToolset) with:

``` r
# install.packages("devtools")
devtools::install_github("KarakulakTulay/MDTToolset")
```

### Vignettes

Please find a comprehensive workflow to MDTToolset in the
[vignette](https://karakulaktulay.github.io/MDTToolset/articles/MDTToolset-vignette.html).

### Documentation

The HTML documentation of the latest version is available at
[GitHub](https://github.com/KarakulakTulay/MDTToolset).

### Credits

The idea behind `MDTToolset` was developed by Abdullah Kahraman, and the
original perl codes for the calculation of MDTs and their switches can
be found [here](https://github.com/abxka/CanIsoNet).

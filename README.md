
# MDTToolset

The [R](https://www.r-project.org) package **MDTToolset** is designed to
calculate the **M**ost **D**ominant *T*ranscript and their switches in
the RNA-Seq datasets. It can calculate;

- Remove redundant transcripts that has the same protein sequences, and
  merges their transcript counts.
- Most dominant transcripts in RNA-Seq data with flexible cutoffs (e.g
  user-defined transcript expression cutoff and enrichment cutoff)
- MDTs found in user-defined percentage of samples
- MDT Switching events with flexible cutoffs (e.g. user-defined cutoffs)

## Installation

You can install the development version of MDTToolset from
[GitHub](https://github.com/KarakulakTulay/MDTToolset) with:

``` r
# install.packages("devtools")
devtools::install_github("KarakulakTulay/MDTToolset")
```

## Background

Alternative splicing plays an essential role in development, tissue
specificity and essential cell functions. One gene can encode for many
transcripts, one of those transcripts might be expressed at a
significantly higher level than the other transcripts. These transcripts
are defined as most dominant transcripts (MDTs). In cancer, those
canonical MDTs might be switched to other transcripts which we call MDT
switching event. This R package aims to find the MDTs in RNA-seq data
and evaluated the MDT switching event in a flexible user-defined way.
The user can define cutoffs for the transcript expression, and MDT
enrichment.

## Example

``` r
library(MDTToolset)
## basic example code
```

## Documentation

The HTML documentation of the latest version released to CRAN is
available at and [GitHub](https://github.com/KarakulakTulay/MDTToolset).

## Credits

The idea behind `MDTToolset` was developed by Abdullah Kahraman, and the
original perl codes for the calculation of MDTs and their switches can
be found [here](https://github.com/abxka/CanIsoNet).

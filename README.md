# r3Cseq
r3Cseq: an R package for the discovery of long-range genomic interactions from chromosome conformation capture and next-generation sequencing data.

# Installation

The r3Cseq package can be easily installed using the usual Bioconductor method:

```{r, eval = F}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("r3Cseq")
```
# Important update
The BioC3.9 in R3.6 will not anymore support the _RangedData_ which is mainly used by r3Cseq in previous versions. 
This version of r3Cseq (version 1.29.0) therefore has migrated codes to use the GRanges to support BioC3.9 in R3.6.  

# Documentation

The documentation for the version of this package can be found in https://github.com/supatt/r3Cseq/blob/master/vignettes folder. How to use to package is also described on http://r3cseq.genereg.net/Site/index.html.

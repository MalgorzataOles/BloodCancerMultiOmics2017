## Installation

To install *BloodCancerMultiOmics2017*, open R (>= 3.5.0) and run the `biocLite` installation script in order to resolve dependencies on the Bioconductor packages.

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("MalgorzataOles/BloodCancerMultiOmics2017") # You might need to first run `install.packages("devtools")`
```

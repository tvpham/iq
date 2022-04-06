# iq
An R package  for protein quantification in mass spectrometry-based proteomics

**Citation**

Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein abundances from ion quantification in DIA-MS-based proteomics, _Bioinformatics_ 2020 Apr 15;36(8):2611-2613.

## Installation

The package is hosted on CRAN. It is best to install from within R.

```
install.packages("iq")
```

If you have a small dataset and just want to try out the MaxLFQ algorithm, you can use the version in pure R without any installation. Just source the R code as follows and omit the ```iq::`` in your script.

```
source("https://github.com/tvpham/iq/releases/download/pureR/iq.R")
```


## Usage

The package can be loaded in the usual manner

```
library("iq")
```

See [a recent example](https://cran.r-project.org/web/packages/iq/vignettes/iq-fast.html).

Or [an older vignette](https://cran.r-project.org/web/packages/iq/vignettes/iq.html) with some visualization.
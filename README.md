# iq: an R package  for protein quantification

This R package provides an implementation of the MaxLFQ algorithm by Cox et al. (2014) in a comprehensive pipeline for DIA-MS (Pham et al. 2020). It also offers options for protein quantification using the N most intense fragment ions, using all fragment ions, and the Tukey's median polish algorithm. In general, the tool can be used to integrate multiple proportional observations into a single quantitative value.

**Citation**

Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein abundances from ion quantification in DIA-MS-based proteomics, _Bioinformatics_ 2020 Apr 15;36(8):2611-2613.
[https://doi.org/10.1093/bioinformatics/btz961](https://doi.org/10.1093/bioinformatics/btz961)

## Installation

The package is hosted on CRAN. It is best to install from within R.

```
install.packages("iq")
```

## Usage

See [a recent example](https://cran.r-project.org/web/packages/iq/vignettes/iq-fast.html) for processing a Spectronaut output. 

Or [an older vignette](https://cran.r-project.org/web/packages/iq/vignettes/iq.html) for processing output from Spectronaut, OpenSWATH and MaxQuant with some visualization.

The package can be loaded in the usual manner

```
library("iq")
```

***To process a DIA-NN output***

DIA-NN 1.9.2 uses the parquet data format for output. We can use the R package ```arrow``` to read the data. 

```
require("arrow")
# if the package "arrow" is not available, you can install it by 
# install.packages("arrow") 
```

The following is an ```iq``` function call to filter on the ```Q.Value```, ```PG.Q.Value```, ```Lib.Q.Value```, and ```Lib.PG.Q.Value``` for a match-between run (MBR) DIA-NN search as discussed [here](https://github.com/vdemichev/DiaNN/discussions/1172#discussioncomment-10680048).

```
process_long_format(arrow::read_parquet("report.parquet"), 
                    output_filename = "report-protein-group.txt", 
                    sample_id = "Run",
                    intensity_col = "Precursor.Normalised",
                    intensity_col_sep = NULL,
                    annotation_col = c("Protein.Ids","Protein.Names", "Genes"),
                    filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.05", 
                                           "Lib.Q.Value" = "0.01", "Lib.PG.Q.Value" = "0.01"))
```

Similarly, for a DIA-NN search without MBR

```
process_long_format(arrow::read_parquet("report.parquet"), 
                    output_filename = "report-protein-group.txt", 
                    sample_id = "Run",
                    intensity_col = "Precursor.Normalised",
                    intensity_col_sep = NULL,
                    annotation_col = c("Protein.Ids","Protein.Names", "Genes"),
                    filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.05", 
                                           "Global.Q.Value" = "0.01", "Global.PG.Q.Value" = "0.01"))
```

We can still process a tsv output file as follows, for example for a search with MBR

```
process_long_format("report.tsv", 
                    output_filename = "report-protein-group.txt", 
                    sample_id = "Run",
                    intensity_col = "Precursor.Normalised",
                    intensity_col_sep = NULL,
                    annotation_col = c("Protein.Ids","Protein.Names", "Genes"),
                    filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.05", 
                                           "Lib.Q.Value" = "0.01", "Lib.PG.Q.Value" = "0.01"))
```

Use the paramter `peptide_extractor` if you want to get the number of peptides per protein, again for example with MBR

```
process_long_format(arrow::read_parquet("report.parquet"), 
                    output_filename = "report-protein-group.txt", 
                    sample_id = "Run",
                    intensity_col = "Precursor.Normalised",
                    intensity_col_sep = NULL,
                    annotation_col = c("Protein.Ids","Protein.Names", "Genes"),
                    filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.05", 
                                           "Lib.Q.Value" = "0.01", "Lib.PG.Q.Value" = "0.01"),
                    peptide_extractor = function(x) gsub("[0-9].*$", "", x))
```

***To process a Spectronaut output***

Use this export schema [iq.rs](https://github.com/tvpham/iq/releases/download/v1.1/iq.rs) to make a long report, for example "Spectronaut_Report.xls".

```
process_long_format("Spectronaut_Report.xls",
                    output_filename = "iq-MaxLFQ.tsv", 
                    sample_id  = "R.FileName",
                    primary_id = "PG.ProteinGroups",
                    secondary_id = c("EG.Library", "FG.Id", "FG.Charge", "F.FrgIon", 
                                     "F.Charge", "F.FrgLossType"),
                    intensity_col = "F.PeakArea",
                    annotation_col = c("PG.Genes", "PG.ProteinNames", "PG.FastaFiles"),
                    filter_string_equal = c("F.ExcludedFromQuantification" = "False"),
                    filter_double_less = c("PG.Qvalue" = "0.01", "EG.Qvalue" = "0.01"),
                    log2_intensity_cutoff = 0)
```

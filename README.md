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

A [blog post](https://digitalbiologylab.github.io/posts/220301-parquet-tsv/) on converting a ```.parquet``` file to a ```.tsv``` file.

The package can be loaded in the usual manner

```
library("iq")
```

***To process a DIA-NN output***

For version of DIA-NN prior to 2.0, the following is an ```iq``` function call to filter on the ```Q.Value```, ```PG.Q.Value```, ```Lib.Q.Value```, and ```Lib.PG.Q.Value``` for a match-between run (MBR) DIA-NN search as discussed [here](https://github.com/vdemichev/DiaNN/discussions/1172#discussioncomment-10680048).

```
process_long_format("report.tsv", 
                    sample_id = "Run",
                    intensity_col = "Fragment.Quant.Raw",
                    output_filename = "report-pg-global.txt", 
                    annotation_col = c("Protein.Names", "Genes"),
                    filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.05", 
                                           "Lib.Q.Value" = "0.01", "Lib.PG.Q.Value" = "0.01"))
```

DIA-NN version 2.0 uses the parquet data format for output. We can use the R package ```arrow``` to read the data. However, the fragment intensities are not reported by the default setting. To perform quantification using MS/MS fragments, one must switch on the --export-quant option. Then we can use R to create a .tsv as an intermediate step as follows

```
require("arrow")
# if the package "arrow" is not available, you can install it by 
# install.packages("arrow") 

raw <- arrow::read_parquet("report.parquet")

# create a new column called "Intensities"
raw$Intensities = paste(raw$Fr.0.Quantity, raw$Fr.1.Quantity, raw$Fr.2.Quantity, 
                        raw$Fr.3.Quantity, raw$Fr.4.Quantity, raw$Fr.5.Quantity, 
                        raw$Fr.6.Quantity, raw$Fr.7.Quantity, raw$Fr.8.Quantity, 
                        raw$Fr.9.Quantity, raw$Fr.10.Quantity, raw$Fr.11.Quantity, sep = ";")

write.table(raw, "report.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# using the new column "Intensities"
iq::process_long_format("report.tsv", 
                        output_filename = "report-protein-group.txt", 
                        sample_id = "Run",
                        intensity_col = "Intensities",
                        annotation_col = c("Protein.Ids","Protein.Names", "Genes"),
                        filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.05", 
                                               "Lib.Q.Value" = "0.01", 
                                               "Lib.PG.Q.Value" = "0.01"))
```


Alternatively, we can use the aggregated intensities in ```Precursor.Normalised```as discussed [here](https://github.com/vdemichev/DiaNN/discussions/951#discussioncomment-8631014).

```
iq::process_long_format(arrow::read_parquet("report.parquet"), 
                        output_filename = "report-protein-group.txt", 
                        sample_id = "Run",
                        intensity_col = "Precursor.Normalised",
                        intensity_col_sep = NULL,
                        annotation_col = c("Protein.Ids","Protein.Names", "Genes"),
                        filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.05", 
                                               "Lib.Q.Value" = "0.01", 
                                               "Lib.PG.Q.Value" = "0.01"))
```

Similarly, for a DIA-NN search without MBR

```
iq::process_long_format(arrow::read_parquet("report.parquet"), 
                        output_filename = "report-protein-group.txt", 
                        sample_id = "Run",
                        intensity_col = "Precursor.Normalised",
                        intensity_col_sep = NULL,
                        annotation_col = c("Protein.Ids","Protein.Names", "Genes"),
                        filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.05", 
                                               "Global.Q.Value" = "0.01", 
                                               "Global.PG.Q.Value" = "0.01"))
```

Finally, use the parameter `peptide_extractor` if you want to get the number of peptides per protein, for example with MBR and the --export-quant option.

```
iq::process_long_format("report.tsv", 
                        output_filename = "report-protein-group.txt", 
                        sample_id = "Run",
                        intensity_col = "Intensities",
                        annotation_col = c("Protein.Ids","Protein.Names", "Genes"),
                        filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.05", 
                                               "Lib.Q.Value" = "0.01", 
                                               "Lib.PG.Q.Value" = "0.01"),
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

# iq: an R package  for protein quantification

This R package provides an implementation of the MaxLFQ algorithm by Cox et al. (2014) in a comprehensive pipeline for DIA-MS (Pham et al. 2020). It also offers options for protein quantification using the N most intense fragment ions, using all fragment ions, and the Tukey's median polish algorithm. In general, the tool can be used to integrate multiple proportional observations into a single quantitative value.

**Citation**

Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein abundances from ion quantification in DIA-MS-based proteomics, _Bioinformatics_ 2020 Apr 15;36(8):2611-2613.
[https://doi.org/10.1093/bioinformatics/btz961](https://doi.org/10.1093/bioinformatics/btz961)

For [version 2](#version2), please cite

Pham TV, Tran CT, Henneman AA, Pham LH, Le DG, Can AH, Bui PH, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. _bioRxiv_. 2025:2025.10.06.680769.
[https://doi.org/10.1101/2025.10.06.680769](https://doi.org/10.1101/2025.10.06.680769)

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

<a id="version2"></a>
## Version 2
***A simple example***

Quantify a data matrix for a single protein
```
X <- matrix(c(10, 12, 12, NA, NA, 10), nrow = 2, 
            dimnames = list(c("ion1", "ion2"), 
                            c("sample1", "sample2", "sample3")))

#      sample1 sample2 sample3
# ion1      10      12      NA
# ion2      12      NA      10

iq::process_matrix(X, method = "maxlfq_bit")

# [1] 11 13  9
```

The result by maxlfq-bit shows that sample2 is 4-fold higher than sample1 (difference of 2 in log2 space). This two samples share ion1 only whose quantitative difference is reflected in the result. Similarly, sample 1 is 4-fold higher than sample3.

***The new iq data format***

A new data format is introduced to support large-scale data processing. This example convert the usual tab-deliminated format file [Bruderer15-DIA-longformat-compact.txt](https://github.com/tvpham/iq/releases/download/v1.1/Bruderer15-DIA-longformat-compact.txt.gz) to the new iq format

```
primary_id <- "PG.ProteinGroups"
sample_id  <- "R.Condition" 
secondary_id <- c("EG.ModifiedSequence", "FG.Charge", "F.FrgIon", "F.Charge")
annotation_col <- c("PG.Genes", "PG.ProteinNames")

iq::long_format_to_iq_format("E:/devops/iq-dev/branches/_test/Bruderer15-DIA-longformat-compact.txt", "bruderer15.iq", 
                        primary_id = primary_id,
                        secondary_id = secondary_id,
                        sample_id = sample_id,
                        annotation_col = annotation_col,
                        intensity_col = "F.PeakArea",
                        intensity_col_sep = NULL,
                        normalization = "median")
```

The new data is contained in a folder **bruderer15.iq**. Then we can process the data using different method

```
iq::process_iq_format("bruderer15.iq", output_filename = "result-maxlfq_bit.tsv",
                      method = "maxlfq_bit")

iq::process_iq_format("bruderer15.iq", output_filename = "result-maxlfq.tsv",
                      method = "maxlfq")

a <- read.delim("result-maxlfq.tsv")
b <- read.delim("result-maxlfq_bit.tsv")

max(a[,4:27] - b[,4:27], na.rm = TRUE)
# [1] 1.012523e-13

min(a[,4:27] - b[,4:27], na.rm = TRUE)
# [1] -1.012523e-13
```
The new method _maxlfq-bit_ and the current implementation _maxlfq_ should give the same result within the tolerance of floating-point arithmetic.


***Cluster processing***

The iq data format makes it very easy to process subsets of proteins on different processing nodes. Here is an example of using the Dutch national grid for processing. You will need to adapt the shell script to your cluster processing and storage architecture.
```
#!/bin/bash
#SBATCH --array=0-119
#SBATCH --time=5:00
#SBATCH --mem=32G

RANGES_FILE=ranges.txt

LOCAL_DIR=/gpfs/work4/0/prjs0919/iq/cluster-processing

R_SCRIPT=$LOCAL_DIR/script.R

OUTPUT_DIR=$LOCAL_DIR/res_files
INPUT_DIR=$LOCAL_DIR/data.iq

RANGES_FILE=$LOCAL_DIR/ranges.txt

mapfile -t <$RANGES_FILE PROTEIN_RANGES

# This needs to be syncronized by hand with the
# --array length
# wc -l ranges.txt
## 120 ranges.txt

RANGE_VALUES=${PROTEIN_RANGES[$SLURM_ARRAY_TASK_ID]}

module load 2024
module load R/4.4.2-gfbf-2024a

OUTPUT_FILE_PREFIX=$OUTPUT_DIR/res
METHOD="maxlfq_bit"

srun $R_SCRIPT $INPUT_DIR $OUTPUT_FILE_PREFIX $METHOD $RANGE_VALUES

```

In this example, we process a dataset of 11987 proteins (**data.iq**) using the **maxlfq_bit** algorithm. We split the data into 120 partitions (file **ranges.txt**). An R script **script.R** is called in each node to process the data.

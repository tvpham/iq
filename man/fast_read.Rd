\name{fast_read}
\alias{fast_read}
\title{
Reading data from an input file
}
\description{
A highly efficient reading of a tab-separated text file for iq processing.
}
\usage{
fast_read(filename,
          sample_id = "R.Condition",
          primary_id = "PG.ProteinGroups",
          secondary_id = c("EG.ModifiedSequence", "FG.Charge", "F.FrgIon", "F.Charge"),
          intensity_col = "F.PeakArea",
          annotation_col = c("PG.Genes", "PG.ProteinNames"),
          filter_string_equal = c("F.ExcludedFromQuantification" = "False"),
          filter_double_less = c("PG.Qvalue" = "0.01", "EG.Qvalue" = "0.01"))
}
\arguments{
  \item{filename}{
A long-format tab-separated text file with a primary column of protein identification, secondary columns of fragment ions, a column of sample names, a column for quantitative intensities, and extra columns for annotation.
}
\item{primary_id}{Unique values in this column form the list of proteins to be quantified.}
\item{secondary_id}{A concatenation of these columns determines the fragment ions used for quantification.}
\item{sample_id}{Unique values in this column form the list of samples.}
\item{intensity_col}{The column for intensities.}
\item{annotation_col}{Annotation columns}
\item{filter_string_equal}{A named vector of strings. Only rows satisfying the filter are kept.}
\item{filter_double_less}{A named vector of strings. Only rows satisfying the filter are kept. Default PG.Qvalue < 0.01 and EG.Qvalue < 0.01.}
}
\value{
A list is returned with following components
\item{protein}{A table of proteins in the first column followed by annotation columns.}
\item{sample}{A vector of samples.}
\item{ion}{A vector of fragment ions to be used for quantification.}
\item{quant_table}{A list of four components: protein_list (index pointing to \code{protein})), sample_list (index pointing to \code{sample}), id (index pointing to \code{ion}), and quant (intensities).}
}
\references{
Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}

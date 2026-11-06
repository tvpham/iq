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
          filter_string_not_equal = NULL,
          filter_double_less = c("PG.Qvalue" = "0.01", "EG.Qvalue" = "0.01"),
          filter_double_greater = NULL,
          intensity_col_sep = NULL,
          intensity_col_id = NULL,
          na_string = "0",
          verbose = TRUE)
}
\arguments{
  \item{filename}{
A long-format tab-separated text file with a primary column of protein identification, secondary columns of fragment ions, a column of sample names, a column for quantitative intensities, and extra columns for annotation.
}
\item{primary_id}{Unique values in this column form the list of proteins to be quantified.}
\item{secondary_id}{One or more columns. A concatenation of these columns determines the fragment ions used for quantification.}
\item{sample_id}{Unique values in this column form the list of samples.}
\item{intensity_col}{The column for intensities, or several columns, one per fragment ion (e.g. \code{paste0("Fr.", 0:11, ".Quantity")} in DIA-NN 2).}
\item{annotation_col}{Annotation columns}
\item{filter_string_equal}{A named vector of strings. Only rows satisfying the condition are kept.}
\item{filter_string_not_equal}{A named vector of strings. Only rows satisfying the condition are kept.}
\item{filter_double_less}{A named vector of strings. Rows with a value less than the given value, or a missing value (empty, NaN or not a number), are kept. Default PG.Qvalue < 0.01 and EG.Qvalue < 0.01.}
\item{filter_double_greater}{A named vector of strings. Rows with a value greater than the given value, or a missing value (empty, NaN or not a number), are kept.}
\item{intensity_col_sep}{A separator character when entries in a single intensity column contain multiple values. It cannot be used with several intensity columns.}
\item{intensity_col_id}{Identities of multiple quantitative values: a single column split by \code{intensity_col_sep}, or one column per intensity column.}
\item{na_string}{The value considered as NA, in addition to empty values and NaN. It is used with \code{intensity_col_sep} or several intensity columns.}
\item{verbose}{A logical value. The default \code{TRUE} prints progress and summary messages; \code{FALSE} prints none of them. Warnings and errors are always shown.}
}
\details{
When there are multiple quantitative values per row (a single intensity column split by \code{intensity_col_sep}, or several intensity columns), this function will replicate entries in other columns and the \code{secondary_id} will be appended with corresponding entries in \code{intensity_col_id} when it is provided. Otherwise, integer values 1, 2, 3, etc... will be used.
}
\value{
A list is returned with following components
\item{protein}{A table of proteins in the first column followed by annotation columns.}
\item{sample}{A vector of samples.}
\item{ion}{A vector of fragment ions to be used for quantification: the \code{secondary_id} columns and the suffix of the intensity value (\code{1}, \code{2}, ... or its \code{intensity_col_id}), joined with \code{_}.}
\item{ion_prefix}{For each ion, its \code{secondary_id} part, without the suffix of the intensity value (used for \code{peptide_extractor} in \code{\link{process}}).}
\item{quant_table}{A list of four components: protein_list (index pointing to \code{protein}), sample_list (index pointing to \code{sample}), id (index pointing to \code{ion}), and quant (intensities).}

If the file cannot be read (e.g. a missing file or column) or no values are left after filtering, \code{fast_read} stops with an error.
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198-1203.

Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}

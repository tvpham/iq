\name{fast_read_parquet}
\alias{fast_read_parquet}
\title{
Reading data from a parquet file
}
\description{
Reads a long-format parquet file, for example a DIA-NN \code{report.parquet}, one row group at a time for iq processing. Only the needed columns are read, and only one row group is held in memory at a time. Requires the \pkg{arrow} package.
}
\usage{
fast_read_parquet(filename,
                  sample_id = "Run",
                  primary_id = "Protein.Group",
                  secondary_id = "Precursor.Id",
                  intensity_col = "Precursor.Normalised",
                  annotation_col = NULL,
                  filter_string_equal = NULL,
                  filter_string_not_equal = NULL,
                  filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.01"),
                  filter_double_greater = NULL,
                  intensity_col_sep = NULL,
                  intensity_col_id = NULL,
                  na_string = "0",
                  verbose = TRUE)
}
\arguments{
\item{filename}{A long-format parquet file.}
\item{sample_id}{Unique values in this column form the list of samples.}
\item{primary_id}{Unique values in this column form the list of proteins to be quantified.}
\item{secondary_id}{One or more columns. A concatenation of these columns determines the fragment ions used for quantification.}
\item{intensity_col}{The column for intensities, or several columns, one per fragment ion (e.g. \code{paste0("Fr.", 0:11, ".Quantity")} in DIA-NN 2).}
\item{annotation_col}{Annotation columns.}
\item{filter_string_equal}{A named vector of strings. Only rows satisfying the condition are kept.}
\item{filter_string_not_equal}{A named vector of strings. Only rows satisfying the condition are kept.}
\item{filter_double_less}{A named vector of numbers. Rows with a value less than the given value, or a missing value, are kept. Default Q.Value < 0.01 and PG.Q.Value < 0.01.}
\item{filter_double_greater}{A named vector of numbers. Rows with a value greater than the given value, or a missing value, are kept.}
\item{intensity_col_sep}{A separator character when a single text intensity column contains multiple values.}
\item{intensity_col_id}{Identities of multiple quantitative values: one column per intensity column, or a single text column split by \code{intensity_col_sep}.}
\item{na_string}{The value considered as NA, in addition to NA and NaN. For numeric intensity columns, a number (default 0).}
\item{verbose}{A logical value. The default \code{TRUE} prints progress and summary messages; \code{FALSE} prints none of them. Warnings and errors are always shown.}}
\details{
Filtering and missing values follow \code{\link{fast_read}}, and the result is the same as \code{\link{fast_read}} on the same data written as a tab-separated text file. A filter value is missing if it is NA, NaN, empty, or text that is not a number; such rows are kept.

When there are multiple quantitative values per row, the \code{secondary_id} is appended with the corresponding entries of \code{intensity_col_id} when it is provided. Otherwise, integer values 1, 2, 3, etc... are used.
}
\value{
A list as returned by \code{\link{fast_read}}.
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198-1203.
}
\author{
Thang V. Pham
}

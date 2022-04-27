\name{process_wide_format}
\alias{process_wide_format}
\title{
Merging rows with identical values in a particular column in a table
}
\description{
Collapses rows with identical values in a particular column in a table. When the values in each row are proportional such as intensities of multiple fragments of a protein, the MaxLFQ algorithm is recommended.
}
\usage{
process_wide_format(input_filename,
                    output_filename,
                    id_column,
                    quant_columns,
                    data_in_log_space = FALSE,
                    annotation_columns = NULL,
                    method = "maxLFQ")
}
\arguments{
  \item{input_filename}{Input filename of a tab-separated value text file.}
  \item{output_filename}{Output filename.}
  \item{id_column}{The column where unique values will be kept. Rows with identical values in this column are merged. Rows with empty values here are removed.}
  \item{quant_columns}{Columns containing numerical data to be merged.}
  \item{data_in_log_space}{A logical value. If \code{FALSE}, the numerical data will be log2-transformed.}
  \item{annotation_columns}{Columns in the input file apart from \code{id_column} and \code{quant_columns} that will be kept in the output.}
  \item{method}{Method for merging. Default value is \code{"maxLFQ"}. Possible values are \code{"maxLFQ"}, \code{"maxLFQ_R"}, \code{"median_polish"}, \code{"top3"}, \code{"top5"}, \code{"meanInt"}, \code{"maxInt"}, \code{"sum"}, \code{"least_na"} and any function for collapsing a numerical matrix to a row vector.}
}
\details{
Method \code{"maxLFQ_R"} implements the MaxLFQ algorithm pure R. It is slower than \code{"maxLFQ"}.

Method \code{"maxInt"} selects row with maximum intensity (top 1).

Method \code{"sum"} sum all intensities.

Method \code{"least_na"} selects row with the least number of missing values.

The value of \code{method} can be a function such as \code{function(x) log2(colSums(2^x, na.rm = TRUE))} for summing all intensities in the original space.
}
\value{
The result table is written to \code{output_filename}. A \code{NULL} value is returned.
}
\references{
Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}

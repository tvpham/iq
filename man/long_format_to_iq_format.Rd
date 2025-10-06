\name{long_format_to_iq_format}
\alias{long_format_to_iq_format}
\title{
Long format input data to iq format
}
\description{
This function reads long-format input data, filters it, and writes the results to the iq format.
}
\usage{
long_format_to_iq_format(input_data, 
                         output_data,
                         sample_id = "Run",
                         primary_id = "Protein.Group",
                         secondary_id = c("Precursor.Id"),
                         intensity_col = "Intensities",
                         annotation_col = NULL,
                         filter_string_equal = NULL,
                         filter_string_not_equal = NULL,
                         filter_double_less = NULL,
                         filter_double_greater = NULL,
                         intensity_col_sep = ";",
                         intensity_col_id = NULL, 
                         na_string = "0",
                         normalization = "median",
                         log2_intensity_cutoff = 0,
                         pdf_out = "qc-plots-iq.pdf",
                         pdf_width = 12,
                         pdf_height = 8,
                         show_boxplot = FALSE)
}

\arguments{
  \item{input_data}{See \code{filename} in \code{\link{fast_read}}.}
  \item{output_data}{Output data.}
  \item{sample_id}{See \code{sample_id} in \code{\link{fast_read}}.}
  \item{primary_id}{See \code{primary_id} in \code{\link{fast_read}}.}
  \item{secondary_id}{See \code{secondary_id} in \code{\link{fast_read}}.}
  \item{intensity_col}{See \code{intensity_col} in \code{\link{fast_read}}.}
  \item{annotation_col}{See \code{annotation_col} in \code{\link{fast_read}}.}
  \item{filter_string_equal}{See \code{filter_string_equal} in \code{\link{fast_read}}.}
  \item{filter_string_not_equal}{See \code{filter_string_not_equal} in \code{\link{fast_read}}.}
  \item{filter_double_less}{See \code{filter_double_less} in \code{\link{fast_read}}.}
  \item{filter_double_greater}{See \code{filter_double_greater} in \code{\link{fast_read}}.}
  \item{intensity_col_sep}{See \code{intensity_col_sep} in \code{\link{fast_read}}.}
  \item{intensity_col_id}{See \code{intensity_col_id} in \code{\link{fast_read}}.}
  \item{na_string}{See \code{intensity_col_id} in \code{\link{fast_read}}.}
  \item{normalization}{Normalization type. Possible values are \code{median} and \code{none}. The default value \code{median} is for median normalization in \code{\link{fast_preprocess}}.}
  \item{log2_intensity_cutoff}{See \code{log2_intensity_cutoff} in \code{\link{fast_preprocess}}.}
  \item{pdf_out}{See \code{pdf_out} in \code{\link{fast_preprocess}}.}
  \item{pdf_width}{See \code{pdf_width} in \code{\link{fast_preprocess}}.}
  \item{pdf_height}{See \code{pdf_height} in \code{\link{fast_preprocess}}.}
  \item{show_boxplot}{See \code{show_boxplot} in \code{\link{fast_preprocess}}.}
}
\value{
The output is written to a new directory named \code{output_data}.
}
\references{
Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}
\seealso{
\code{\link{fast_read}}, \code{\link{fast_preprocess}}
}


\name{process_long_format}
\alias{process_long_format}
\title{
Long format to a wide format table using the MaxLFQ algorithm
}
\description{
A convenient function combining multiple steps to process a long format table using the MaxLFQ algorithm.
}
\usage{
process_long_format(input_filename,
                    output_filename,
                    sample_id = "File.Name",
                    primary_id = "Protein.Group",
                    secondary_id = "Precursor.Id",
                    intensity_col = "Fragment.Quant.Corrected",
                    annotation_col = NULL,
                    filter_string_equal = NULL,
                    filter_string_not_equal = NULL,
                    filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.01"),
                    filter_double_greater = NULL,
                    intensity_col_sep = ";",
                    intensity_col_id = NULL,
                    na_string = "0",
                    normalization = "median",
                    log2_intensity_cutoff = 0,
                    pdf_out = "qc-plots.pdf",
                    pdf_width = 12,
                    pdf_height = 8,
                    show_boxplot = TRUE,
                    peptide_extractor = NULL)
}
\arguments{
  \item{input_filename}{See \code{filename} in \code{\link{fast_read}}.}
  \item{output_filename}{Output filename.}
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
  \item{peptide_extractor}{A function to parse peptides.}
}
\value{
After processing with \code{\link{fast_read}}, \code{\link{fast_preprocess}}, and \code{\link{fast_MaxLFQ}}, the result table is written to \code{output_filename}. A \code{NULL} value is returned. If \code{peptide_extractor} is not \code{NULL}, fragment statistics for each protein will be calculated based on the result of the extractor function. Counting the number of peptides contributing to a protein is possible using an appropriate extractor function. An example value for \code{peptide_extractor} is \code{function(x) gsub("[0-9].*$", "", x)}, which removes the charge state and fragment descriptors in an ion descriptor to obtain unique peptide sequences. One can examine the \code{ion} component returned by the \code{\link{fast_read}} function to derive a regular expression to be used in the \code{gsub} function above.
}
\references{
Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}
\seealso{
\code{\link{fast_read}}, \code{\link{fast_preprocess}}, \code{\link{fast_MaxLFQ}}
}


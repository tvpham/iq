\name{fast_preprocess}
\alias{fast_preprocess}
\title{
Data filtering and normalization
}
\description{
Filters out low intensities and performs median normalization.
}
\usage{
fast_preprocess(quant_table,
                median_normalization = TRUE,
                log2_intensity_cutoff = 0,
                pdf_out = "qc-plots-fast.pdf",
                pdf_width = 12,
                pdf_height = 8)
}
\arguments{
  \item{quant_table}{The \code{quant_table} component as returned by \code{fast_read}.
}
\item{median_normalization}{A logical value. The default \code{TRUE} value is to perform median normalization.}
\item{log2_intensity_cutoff}{Entries lower than this value in log2 space are ignored. Plot a histogram of all intensities to set this parameter.}
\item{pdf_out}{A character string specifying the name of the PDF output. A \code{NULL} value will suppress the PDF output.}
\item{pdf_width}{Width of the pdf output in inches.}
\item{pdf_height}{Height of the pdf output in inches.}
}
\value{
A list is returned with the same components as input data in which low intensities are filtered out and median normalization is performed if requested.
}
\references{
Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}
\seealso{
\code{\link{fast_read}}}

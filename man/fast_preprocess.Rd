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
                pdf_height = 8,
                show_boxplot = TRUE,
                boxplot_max_points = 10000,
                verbose = TRUE)
}
\arguments{
  \item{quant_table}{The \code{quant_table} component as returned by \code{fast_read}.
}
\item{median_normalization}{A logical value. The default \code{TRUE} value is to perform median normalization.}
\item{log2_intensity_cutoff}{Entries not above this value in log2 space are ignored; \code{NULL} for no cutoff. Plot a histogram of all intensities to set this parameter. Zero and negative intensities (not finite in log2 space) are always removed.}
\item{pdf_out}{A character string specifying the name of the PDF output. A \code{NULL} value will suppress the PDF output.}
\item{pdf_width}{Width of the pdf output in inches.}
\item{pdf_height}{Height of the pdf output in inches.}
\item{show_boxplot}{A logical value. The default \code{TRUE} value is to create boxplots of fragment intensities for each sample.}
\item{boxplot_max_points}{The maximum number of values per sample in the boxplots (evenly spaced values are used for a larger sample). The median normalization always uses all values. \code{NULL} for all values.}
\item{verbose}{A logical value. The default \code{TRUE} prints progress and summary messages; \code{FALSE} prints none of them. Warnings and errors are always shown.}
}
\value{
A list is returned with the same components as input data in which low intensities are filtered out and median normalization is performed if requested.
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198-1203.

Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}
\seealso{
\code{\link{fast_read}}}

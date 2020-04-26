\name{preprocess}
\alias{preprocess}
\title{
Data preprocessing for protein quantification
}
\description{
Prepares a long-format input including removing low-intensity ions and performing median normalization.
}
\usage{
preprocess(quant_table,
           primary_id = "PG.ProteinGroups",
           secondary_id = c("EG.ModifiedSequence", "FG.Charge", "F.FrgIon", "F.Charge"),
           sample_id = "R.Condition",
           intensity_col = "F.PeakArea",
           median_normalization = TRUE,
           log2_intensity_cutoff = 0,
           pdf_out = "qc-plots.pdf",
           pdf_width = 12,
           pdf_height = 8)
}
\arguments{
  \item{quant_table}{
A long-format table with a primary column of protein identification, secondary columns of fragment ions, a column of sample names, and a column for quantitative intensities.
}
\item{primary_id}{Unique values in this column form the list of proteins to be quantified.}
\item{secondary_id}{A concatenation of these columns determines the fragment ions used for quantification.}
\item{sample_id}{Unique values in this column form the list of samples.}
\item{intensity_col}{The column for intensities.}
\item{median_normalization}{A logical value. The default \code{TRUE} value is to perform median normalization.}
\item{log2_intensity_cutoff}{Entries lower than this value in log2 space are ignored. Plot a histogram of all intensities to set this parameter.}
\item{pdf_out}{A character string specifying the name of the PDF output. A NULL value will suppress the PDF output.}
\item{pdf_width}{Width of the pdf output in inches.}
\item{pdf_height}{Height of the pdf output in inches.}
}
\value{
A data frame is returned with following components
\item{protein_list}{A vector of proteins.}
\item{sample_list}{A vector of samples.}
\item{id}{A vector of fragment ions to be used for quantification.}
\item{quant}{A vector of log2 intensities.}
}
\examples{
\donttest{
data("spikeins")
head(spikeins)
# This example set of spike-in proteins has been 'median-normalized'.
norm_data <- iq::preprocess(spikeins, median_normalization = FALSE, pdf_out = NULL)
}
}
\references{
Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}

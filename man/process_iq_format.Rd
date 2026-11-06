\name{process_iq_format}
\alias{process_iq_format}
\title{
Process data in the iq format
}
\description{
Estimates protein abundance of all proteins in the iq format data.
}
\usage{
process_iq_format(iqformat_folder,
                  output_filename = NULL,
                  method = "maxlfq_bit", 
                  p1 = NULL, p2 = NULL, k = 1.345, min_M = 15,
                  rescale_method = "to_mean",
                  buffer_mb = 1024,
                  n_threads = -1,
                  verbose = TRUE)
}
\arguments{
  \item{iqformat_folder}{A folder with data in the iq format, as written by \code{\link{long_format_to_iq_format}}.}
  \item{output_filename}{Output file name. If \code{NULL}, a file name derived from \code{iqformat_folder} will be used.}
  \item{method}{See \code{method} in \code{\link{process_matrix}}, or a function. A function is called for each protein with the matrix of log2 intensities (ions in rows, with the ion indices of the iq format as row names; samples in columns, with the sample names as column names; \code{NA} for missing values) and must return a numeric vector with one value per sample.
Values that are not finite (\code{-Inf}, \code{Inf}, \code{NaN}) become \code{NA}.}
  \item{p1}{See \code{p1} in \code{\link{process_matrix}}.}
  \item{p2}{See \code{p2} in \code{\link{process_matrix}}.}
  \item{k}{See \code{k} in \code{\link{process_matrix}}.}
  \item{min_M}{See \code{min_M} in \code{\link{process_matrix}}.}
  \item{rescale_method}{Rescale the output intensities. See \code{method} in \code{\link{rescale}}. The default \code{"to_mean"} shifts each connected component of samples so that the mean of its estimates equals the mean of its data (in log2 space); \code{NULL} for no rescaling.}
  \item{buffer_mb}{Memory in megabytes for a batch of proteins quantified together (see Details), or \code{NULL} to quantify the proteins one by one in R.}
  \item{n_threads}{See \code{n_threads} in \code{\link{process_matrix}}.}
  \item{verbose}{A logical value. The default \code{TRUE} prints progress and summary messages; \code{FALSE} prints none of them. Warnings and errors are always shown.}
}
\details{
The protein files are read in batches of about \code{buffer_mb} megabytes of values, and the proteins of a batch are quantified in parallel (\code{n_threads}), as in \code{\link{process}}. With \code{buffer_mb = NULL}, or a function as \code{method}, each protein is read with \code{read.delim} and quantified with \code{\link{process_matrix}} in R. The log2 values are shifted by the normalization of each sample (the \code{normalization} column of \code{_sample.tsv}, if present. All columns of \code{_protein.tsv} are written before the samples: the protein annotation, and the fragment statistics if \code{\link{long_format_to_iq_format}} computed them.
}
\value{
\code{NULL}, invisibly. A tab-separated table is written to \code{output_filename} (by default \code{iqformat_folder} followed by \code{-quant.tsv}), with one row for each protein with values, columns for the protein annotation and statistics, and the samples.
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198-1203.
}
\author{
Thang V. Pham
}
\seealso{
\code{\link{long_format_to_iq_format}}, \code{\link{process}}, \code{\link{process_matrix}}, \code{\link{rescale}}
}

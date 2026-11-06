\name{iq_format_to_protein_list}
\alias{iq_format_to_protein_list}
\title{
Protein data from the iq format
}
\description{
Reads the data of one or more proteins from the iq format into a list of matrices, as returned by \code{\link{create_protein_list}}, for example to plot a protein with \code{\link{plot_protein}}.
}
\usage{
iq_format_to_protein_list(iqformat_folder, proteins, normalization = TRUE)
}
\arguments{
  \item{iqformat_folder}{A folder with data in the iq format, as written by \code{\link{long_format_to_iq_format}}.}
  \item{proteins}{A vector of protein identifiers, as in the first column of \code{_protein.tsv}.}
  \item{normalization}{A logical value. The default \code{TRUE} applies the normalization of each sample stored in \code{_sample.tsv}, as \code{\link{process_iq_format}} does; \code{FALSE} gives the intensities as read.}
}
\details{
Only the files of the requested proteins are read. The intensities of an ion in a sample are summed, as in \code{\link{process_iq_format}}, so that \code{\link{process_matrix}} on a matrix of the list gives the same result as the row of that protein in the output of \code{\link{process_iq_format}} with the same method. A protein without values (it has no protein file, e.g. all its values are below the intensity cut-off) gets a matrix without rows, with a warning. An identifier that is not in \code{_protein.tsv} stops with an error.
}
\value{
A list with one element for each protein in \code{proteins}, named by its identifier: a matrix of log2 intensities with the ions in rows (named as in \code{_ion.tsv}) and all samples in columns (named as in \code{_sample.tsv}); \code{NA} for a missing value.
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198-1203.
}
\author{
Thang V. Pham
}
\seealso{
\code{\link{long_format_to_iq_format}}, \code{\link{process_iq_format}}, \code{\link{plot_protein}}, \code{\link{process_matrix}}
}
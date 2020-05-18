\name{fast_MaxLFQ}
\alias{fast_MaxLFQ}
\title{
The MaxLFQ algorithm
}
\description{
A fast implementation of the MaxLFQ algorithm.
}
\usage{
fast_MaxLFQ(norm_data, row_names = NULL, col_names = NULL)
}
\arguments{
  \item{norm_data}{A list of four vectors with equal length \code{protein_list}, \code{sample_list}, \code{id} and \code{quant} as prepared by the \code{fast_preprocess} function or the \code{quant_table} component returned by the \code{fast_read} function. Note that \code{quant} should contain log2 intensities.}
  \item{row_names}{A vector of character strings for row names. If \code{NULL}, unique values in the \code{protein_list} component of \code{norm_data} will be used. Otherwise, it should be the \code{sample} component returned by the \code{fast_read}.}
  \item{col_names}{A vector of character strings for column names. If \code{NULL}, unique values in the \code{sample_list} component of \code{norm_data} will be used. Otherwise, it should be the \code{sample} component returned by the \code{fast_read}.}
}
\value{
A list is returned with two components
\item{estimate}{A quantification result table.}
\item{annotation}{A vector of strings indicating membership in case of multiple connected components for each row of \code{estimate}.}
}
\references{
Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}
\seealso{
\code{\link{fast_read}}, \code{\link{fast_preprocess}}}

\name{process_iq_format}
\alias{process_iq_format}
\title{
Process data in the iq format
}
\description{
Estimates protein abundance of all proteins in the iq format data.
}
\usage{
process_iq_format(input_data,
                  output_filename = NULL,
                  method = "maxlfq_bit", 
                  p1 = NULL, p2 = NULL, k = 1.345, min_M = 15, n_threads = -1,
                  rescale_method = NULL)
}
\arguments{
  \item{input_data}{Input data in the iq format.}
  \item{output_filename}{Output file name. If \code{NULL}, a file name derived from \code{input_data} will be used.}
  \item{method}{See \code{method} in \code{\link{process_matrix}}.}
  \item{p1}{See \code{p1} in \code{\link{process_matrix}}.}
  \item{p2}{See \code{p2} in \code{\link{process_matrix}}.}
  \item{k}{See \code{k} in \code{\link{process_matrix}}.}
  \item{min_M}{See \code{min_M} in \code{\link{process_matrix}}.}
  \item{n_threads}{See \code{n_threads} in \code{\link{process_matrix}}.}
  \item{rescale_method}{Rescale the output intensities. See \code{method} in \code{\link{rescale}}.}
}
\value{
A tab-separated table with rows representing proteins and columns representing protein annotations and samples.
}
\references{
Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}

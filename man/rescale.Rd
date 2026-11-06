\name{rescale}
\alias{rescale}
\title{
Rescale relative quantification
}
\description{
Each connected component of the input values can be rescaled to match the level of the overall input intensities.
}
\usage{
rescale(x, X, method)
}
\arguments{
  \item{x}{A row vector containing the quantification result.}
  \item{X}{An input matrix, or a data frame with numeric columns, of log2-transformed intensities. Columns represent samples, and rows represent fragment ions.}
  \item{method}{One of the methods for rescaling: \code{"median_mean"}, \code{"to_median"}, \code{"to_mean"} or \code{"to_sum"}.}
}
\value{
A row vector of rescaled values of the same size as \code{x}, with each connected component rescaled independently. 
For the \code{"median_mean"} method, single-sample
components are rescaled to the median of input data, while all other components are rescaled to the mean. For the \code{"to_median"} and \code{"to_mean"} methods, each component is rescaled to the median or
mean of the input data. For the \code{"to_sum"} method, each component is shifted so that its values, in linear scale (\code{2^x}), add up to the summed intensities of its input data (\code{2^X}) in the samples of the component; a single-sample component gets the log2 of the summed intensities of that sample.
For all methods, a sample without input data (an all-\code{NA} column of \code{X}) gets \code{NA}, and a component whose values in \code{x} are all \code{NA} is left unchanged.
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198-1203.
}
\author{
Thang V. Pham
}

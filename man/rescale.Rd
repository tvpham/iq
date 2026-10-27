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
  \item{X}{An input matrix of log2-transformed intensities. Columns represent samples, and rows represent fragment ions.}
  \item{method}{One of the methods for rescaling: \code{"median-mean"}, \code{"median-median"}, \code{"mean-mean"} or \code{"sum"}.}
}
\value{
A row vector of rescaled values of the same size as \code{x}, with each connected component rescaled independently. 
For the \code{"median-mean"} method, single-sample
components are rescaled to the median of input data, while all other components are rescaled to the mean. For other methods, each component is rescaled to the median, 
mean, or sum of the input data.
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198–1203.
}
\author{
Thang V. Pham
}

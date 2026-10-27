\name{huber_weight}
\alias{huber_weight}
\title{
Calculate Huber weights for an individual matrix 
}
\description{
Calculate Huber weights for an individual matrix 
}
\usage{
huber_weight(X, x, c = 1.4826, k = 1.345)

}
\arguments{
  \item{X}{A matrix of intensities in the log2 space. Columns are samples and rows are fragment ions.}
  \item{x}{An initial estimate}
  \item{c}{Weighting function parameter.}
  \item{k}{Weighting function parameter.}
}
\value{
A numeric weight vector of size equal to the number of rows of X.
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198–1203.
}
\author{
Thang V. Pham
}

\name{connected_component}
\alias{connected_component}
\title{
Group connected samples
}
\description{
Two samples are directly connected if they share at least one quantified fragment ion.
Two samples in a connected component are either directly connected or indirectly via other samples.
}
\usage{
connected_component(X)
}
\arguments{
  \item{X}{A matrix of intensities in the log2 space. Columns are samples and rows are fragment ions.}
}
\value{
A vector with length equal to the number of columns of the input containing membership of the connected components.
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198–1203.
}
\author{
Thang V. Pham
}

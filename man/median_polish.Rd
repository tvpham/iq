\name{median_polish}
\alias{median_polish}
\title{
A wrapper for the R implementation of the median polish algorithm
}
\description{
Estimates protein abundances using the Tukey median polish algorithm.
}
\usage{
median_polish(X)
}
\arguments{
  \item{X}{A matrix of ion intensities in log2 space. Columns are samples and rows are fragment ions.}
}
\value{
A list of two components is returned
\item{estimate}{A vector with length equal to the number of columns of the input containing the protein abundances in log2 space.}
\item{annotation}{Reserved, currently an empty string}
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198–1203.

Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.

Tukey JW. \emph{Exploratory Data Analysis}, Reading Massachusetts: Addison-Wesley,
1977.
}
\author{
Thang V. Pham
}

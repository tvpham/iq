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
\item{estimate}{A vector with length equal to the number of columns of the input containing the protein abundances.}
\item{annotation}{Reserved, currently an empty string}
}
\references{
Tukey JW. \emph{Exploratory Data Analysis}, Reading Massachusetts: Addison-Wesley,
1977.
}
\author{
Thang V. Pham
}

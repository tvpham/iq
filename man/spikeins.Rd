\name{spikeins}
\alias{spikeins}
\docType{data}
\title{
An example dataset of 12 spike-in proteins
}
\description{
A subset of the Bruderer 2015 dataset containing 12 spike-in proteins. The full dataset was exported from the Spectronaut software. The complete dataset has been median-normalized.
}
\usage{data("spikeins")}
\format{
  A data frame with 18189 observations on the following 9 variables.
  \describe{
    \item{\code{R.Condition}}{Sample names.}
    \item{\code{PG.ProteinGroups}}{Protein identifiers.}
    \item{\code{EG.ModifiedSequence}}{Sequence of the fragment ions.}
    \item{\code{FG.Charge}}{Fragment group charge.}
    \item{\code{F.FrgIon}}{Fragment ions.}
    \item{\code{F.Charge}}{Fragment charges.}
    \item{\code{F.PeakArea}}{Quantitative values.}
    \item{\code{PG.Genes}}{Gene names.}
    \item{\code{PG.ProteinNames}}{Protein names.}
  }
}
\examples{
data("spikeins")
head(spikeins)
}
\keyword{datasets}

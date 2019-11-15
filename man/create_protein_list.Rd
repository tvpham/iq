\name{create_protein_list}
\alias{create_protein_list}
\title{
Creating a list of matrices of fragment ion intensities for all proteins
}
\description{
For each protein, a numerical matrix is formed where the columns are samples and rows are fragment ions.
}
\usage{
create_protein_list(preprocessed_data)
}
\arguments{
  \item{preprocessed_data}{A data frame of four components as output of the \code{preprocess} function.}
}
\value{
A list where each element contains the quantitative data of a protein. The column names are sample names and the row names fragment ions.
}
\examples{
\donttest{
data("spikeins")
head(spikeins)
# This example set of spike-in proteins has been 'median-normalized'.
norm_data <- iq::preprocess(spikeins, median_normalization = FALSE)
protein_list <- iq::create_protein_list(norm_data)
}
}
\author{
Thang V. Pham
}
\seealso{
\code{\link{preprocess}}
}

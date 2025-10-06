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
norm_data <- iq::preprocess(spikeins, median_normalization = FALSE, pdf_out = NULL)
protein_list <- iq::create_protein_list(norm_data)
}
}
\references{
Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}
\seealso{
\code{\link{preprocess}}
}

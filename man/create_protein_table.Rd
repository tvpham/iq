\name{create_protein_table}
\alias{create_protein_table}
\title{
Protein quantification for a list of proteins
}
\description{
Travels through the input list and quantifies all proteins one by one.
}
\usage{
create_protein_table(protein_list, method = "maxLFQ", ...)
}
\arguments{
\item{protein_list}{The input protein list}
\item{method}{Possible values are "maxLFQ", "median_polish", "topN", and "meanInt".}
\item{...}{Additional parameters for individual quantitation methods.}
}
\value{
A list of two components is returned
\item{estimate}{A table of protein abundances for all samples.}
\item{annotation}{A vector of annotations, one for each protein.}
}
\examples{
data("spikeins")
# This example set of spike-in proteins has been 'median-normalized'.
norm_data <- iq::preprocess(spikeins, median_normalization = FALSE, pdf_out = NULL)
protein_list <- iq::create_protein_list(norm_data)
result <- iq::create_protein_table(protein_list)
head(result)
}
\references{
Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}
\seealso{
\code{\link{create_protein_list}}, \code{\link{maxLFQ}}, \code{\link{median_polish}}, \code{\link{topN}}, \code{\link{meanInt}}

}

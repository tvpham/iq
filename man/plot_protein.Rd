\name{plot_protein}
\alias{plot_protein}
\title{
Plotting the underlying quantitative data for a protein
}
\description{
Displays the underlying data for a protein.
}
\usage{
plot_protein(X, main = "", col = NULL, split = 0.6, ...)
}
\arguments{
\item{X}{Protein data matrix.}
\item{main}{Title of the plot.}
\item{col}{Colors of the rows of the data matrix.}
\item{split}{Fraction of the plotting area for the main figure. The remaining one is for legend. Set this parameter to \code{NULL} to ignore the legend area.}
\item{...}{Additional parameters for plotting.}
}
\value{
A \code{NULL} value is returned.
}
\examples{
data("spikeins")
head(spikeins)
# This example set of spike-in proteins has been 'median-normalized'.
norm_data <- iq::preprocess(spikeins, median_normalization = FALSE, pdf_out = NULL)
protein_list <- iq::create_protein_list(norm_data)
iq::plot_protein(protein_list$P00366, main = "Protein P00366", split = NULL)
}
\references{
Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}

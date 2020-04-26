\name{maxLFQ}
\alias{maxLFQ}
\title{
The MaxLFQ algorithm for protein quantification
}
\description{
Estimates protein abundances by aiming to maintain the fragment intensity ratios between samples.
}
\usage{
maxLFQ(X)
}
\arguments{
  \item{X}{A matrix of ion intensities in log2 space. Columns are samples and rows are fragment ions.}
}
\value{
A list of two components is returned
\item{estimate}{A vector with length equal to the number of columns of the input containing the protein abundances.}
\item{annotation}{An empty string if all quantified samples are connected. Otherwise, a string of membership of the connected components is returned.}
}
\references{
Cox J, Hein MY, Luber CA, et al. Accurate proteome-wide label-free quantification
by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ.
\emph{Mol Cell Proteomics}. 2014;13(9):2513–2526.

Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}

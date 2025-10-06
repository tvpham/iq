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
Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}

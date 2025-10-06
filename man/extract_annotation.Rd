\name{extract_annotation}
\alias{extract_annotation}
\title{
Protein annotation extraction
}
\description{
Extracts annotation columns from a long-format input
}
\usage{
extract_annotation(protein_ids, quant_table, primary_id = "PG.ProteinGroups",
                   annotation_columns = NULL)
}
\arguments{
  \item{protein_ids}{A vector of protein ids.}
  \item{quant_table}{A long-format input table. The input is typically the same as input to the \code{preprocess} function.}
  \item{primary_id}{The column containing protein ids.}
  \item{annotation_columns}{A vector of columns for annotation.}
}
\value{
A table of proteins and associated annotation extracted from the input.
}
\examples{
\donttest{
data("spikeins")
extra_names <- iq::extract_annotation(levels(spikeins$PG.ProteinGroups),
                                      spikeins,
                                      annotation_columns = c("PG.Genes", "PG.ProteinNames"))
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

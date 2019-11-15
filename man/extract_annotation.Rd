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
extra_names <- iq::extract_annotation(rownames(result$estimate),
                                      spikeins,
                                      annotation_columns = c("PG.Genes", "PG.ProteinNames"))
}
}
\author{
Thang V. Pham
}
\seealso{
\code{\link{preprocess}}
}

\name{read_fasta_file}
\alias{read_fasta_file}
\title{
Reading a FASTA file
}
\description{
Reads protein sequences from a FASTA file into a named list, for example for the \code{rfasta} argument of \code{\link{process}} and \code{\link{process_long_format}}.
}
\usage{
read_fasta_file(fasta_file_name, regex_id = "^>.*\\\\|(.*)\\\\|.*")
}
\arguments{
  \item{fasta_file_name}{A FASTA file.}
  \item{regex_id}{A regular expression applied to each header line (the line starting with \code{>}); the first group in parentheses is the name of the sequence. The default takes the accession between the two \code{|} of a UniProt header such as \code{>sp|P04406|G3P_HUMAN ...}.}
}
\details{
A sequence is the concatenation of the lines below its header, with all white space removed. A header that does not match \code{regex_id} gives the whole header line, without \code{>}, as the name. Lines before the first header are ignored.
}
\value{
A named list of sequences (character strings), in the order of the file.
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198-1203.
}
\author{
Thang V. Pham
}
\seealso{
\code{\link{process}}
}

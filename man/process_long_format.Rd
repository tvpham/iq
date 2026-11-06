\name{process_long_format}
\alias{process_long_format}
\title{
Long format to a wide format table using the MaxLFQ algorithm
}
\description{
A convenient function combining multiple steps to process a long format table using the MaxLFQ algorithm.
}
\usage{
process_long_format(input_data,
                    output_filename,
                    sample_id = "File.Name",
                    primary_id = "Protein.Group",
                    secondary_id = "Precursor.Id",
                    intensity_col = "Fragment.Quant.Corrected",
                    annotation_col = NULL,
                    filter_string_equal = NULL,
                    filter_string_not_equal = NULL,
                    filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.01"),
                    filter_double_greater = NULL,
                    intensity_col_sep = ";",
                    intensity_col_id = NULL,
                    na_string = "0",
                    normalization = "median",
                    log2_intensity_cutoff = 0,
                    pdf_out = "qc-plots.pdf",
                    pdf_width = 12,
                    pdf_height = 8,
                    show_boxplot = TRUE,
                    peptide_extractor = NULL,
                    rfasta = NULL,
                    verbose = TRUE)
}
\arguments{
  \item{input_data}{A data frame or a filename. See \code{filename} in \code{\link{fast_read}}.}
  \item{output_filename}{Output filename.}
  \item{sample_id}{See \code{sample_id} in \code{\link{fast_read}}.}
  \item{primary_id}{See \code{primary_id} in \code{\link{fast_read}}.}
  \item{secondary_id}{See \code{secondary_id} in \code{\link{fast_read}}.}
  \item{intensity_col}{See \code{intensity_col} in \code{\link{fast_read}}.}
  \item{annotation_col}{See \code{annotation_col} in \code{\link{fast_read}}.}
  \item{filter_string_equal}{See \code{filter_string_equal} in \code{\link{fast_read}}.}
  \item{filter_string_not_equal}{See \code{filter_string_not_equal} in \code{\link{fast_read}}.}
  \item{filter_double_less}{See \code{filter_double_less} in \code{\link{fast_read}}.}
  \item{filter_double_greater}{See \code{filter_double_greater} in \code{\link{fast_read}}.}
  \item{intensity_col_sep}{See \code{intensity_col_sep} in \code{\link{fast_read}}.}
  \item{intensity_col_id}{See \code{intensity_col_id} in \code{\link{fast_read}}.}
  \item{na_string}{See \code{na_string} in \code{\link{fast_read}}.}
  \item{normalization}{Normalization type. Possible values are \code{median} and \code{none}. The default value \code{median} is for median normalization in \code{\link{fast_preprocess}}.}
  \item{log2_intensity_cutoff}{See \code{log2_intensity_cutoff} in \code{\link{fast_preprocess}}.}
  \item{pdf_out}{See \code{pdf_out} in \code{\link{fast_preprocess}}.}
  \item{pdf_width}{See \code{pdf_width} in \code{\link{fast_preprocess}}.}
  \item{pdf_height}{See \code{pdf_height} in \code{\link{fast_preprocess}}.}
  \item{show_boxplot}{See \code{show_boxplot} in \code{\link{fast_preprocess}}.}
  \item{peptide_extractor}{A function that maps ion identifiers to peptides, to compute fragment statistics per protein, or \code{NULL} (the default) for none. It is called with the ions of one protein and must return one character value per ion. See Value.}
  \item{rfasta}{A named list of protein sequences, named by protein accession as in \code{primary_id} (a protein group is split at \code{;}), e.g. from \code{\link{read_fasta_file}}. With \code{peptide_extractor} and an input file, the sequence coverage of each protein is added.}
  \item{verbose}{A logical value. The default \code{TRUE} prints progress and summary messages; \code{FALSE} prints none of them. Warnings and errors are always shown.}
}
\value{
An input data frame is filtered and processed with \code{\link{preprocess}} and \code{\link{fast_MaxLFQ}}; an input file is processed with \code{\link{fast_read}}, \code{\link{fast_preprocess}}, and \code{\link{fast_MaxLFQ}}.
Subsequently, the result is written to \code{output_filename}. The quantification values are in log2 space.
A \code{NULL} value is returned. If \code{peptide_extractor} is not \code{NULL}, fragment statistics
for each protein will be calculated based on the result of the extractor function: the columns \code{n_fragments} and \code{n_peptides}, and for an input file also \code{peptides}, \code{XIC} (the summed intensities, before preprocessing) and, with \code{rfasta}, \code{coverage}. Counting the number of peptides contributing to a protein is possible using an appropriate extractor function. An example value for \code{peptide_extractor} is \code{function(x) gsub("\\\\(UniMod:\\\\d+\\\\)|\\\\d|_", "", x)}, which removes modifications, digit numbers, and underscores in an ion descriptor to obtain peptide sequences. One can examine the \code{ion} component returned by the \code{\link{fast_read}} function to derive a regular expression to be used in the \code{gsub} function above.
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198-1203.

Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}
\seealso{
\code{\link{fast_read}}, \code{\link{fast_preprocess}}, \code{\link{fast_MaxLFQ}}
}


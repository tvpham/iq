\name{long_format_to_iq_format}
\alias{long_format_to_iq_format}
\title{
Long format input data to iq format
}
\description{
This function reads long-format input data (a text file or a parquet file), filters it, and writes the results to the iq format. The input is read and written incrementally, so it is never held in memory as a whole.
}
\usage{
long_format_to_iq_format(input_filename, 
                         output_iqformat_folder,
                         sample_id = "Run",
                         primary_id = "Protein.Group",
                         secondary_id = c("Precursor.Id"),
                         intensity_col = "Precursor.Normalised",
                         annotation_col = NULL,
                         filter_string_equal = NULL,
                         filter_string_not_equal = NULL,
                         filter_double_less = c("Q.Value" = "0.01",
                                                "PG.Q.Value" = "0.01"),
                         filter_double_greater = NULL,
                         intensity_col_sep = NULL,
                         intensity_col_id = NULL, 
                         na_string = "0",
                         normalization = "median",
                         log2_intensity_cutoff = 0,
                         pdf_out = "qc-plots-iq-format.pdf",
                         pdf_width = 12,
                         pdf_height = 8,
                         show_boxplot = TRUE,
                         peptide_extractor = NULL,
                         rfasta = NULL,
                         buffer_mb = 1024,
                         verbose = TRUE)
}

\arguments{
  \item{input_filename}{The name of a parquet file (ending in \code{.parquet}; requires the \pkg{arrow} package) or of a tab-separated text file. The filtering is the same as in \code{\link{fast_read_parquet}} and \code{\link{fast_read}}. For a data frame, use \code{\link{process}}.}
  \item{output_iqformat_folder}{A folder for the iq format, new or empty.}
  \item{sample_id}{See \code{sample_id} in \code{\link{fast_read}}.}
  \item{primary_id}{See \code{primary_id} in \code{\link{fast_read}}.}
  \item{secondary_id}{See \code{secondary_id} in \code{\link{fast_read}}.}
  \item{intensity_col}{See \code{intensity_col} in \code{\link{fast_read}}. Several columns, one per fragment ion, can be given (e.g. \code{paste0("Fr.", 0:11, ".Quantity")} for a DIA-NN 2 parquet file, with \code{intensity_col_sep = NULL}).}
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
  \item{peptide_extractor}{See \code{peptide_extractor} in \code{\link{process}}: if not \code{NULL}, the fragment statistics of each protein are computed and stored in \code{_protein.tsv}.}
  \item{rfasta}{See \code{rfasta} in \code{\link{process}}: the sequence coverage, with \code{peptide_extractor}.}
  \item{buffer_mb}{Memory in megabytes for buffering data. It changes only the memory use and speed, not the result.}
  \item{verbose}{A logical value. The default \code{TRUE} prints progress and summary messages; \code{FALSE} prints none of them. Warnings and errors are always shown.}
}
\details{
The input is read in chunks: blocks of lines of a text file, or row groups of a parquet file (converted to R in slices of rows). The values of each chunk are buffered (\code{buffer_mb}) and appended to the files. Only the lists of proteins, samples and ions are kept in memory; for a parquet file, also one row group at a time, as decoded by \pkg{arrow}.

The defaults are those of \code{\link{process}} (DIA-NN precursor quantities \code{Precursor.Normalised} per \code{Run}, filtered on \code{Q.Value} and \code{PG.Q.Value} below 0.01, log2 intensity cut-off 0), so that the same arguments give the same proteins as \code{\link{process}}.

The preprocessing is that of \code{\link{fast_preprocess}}: values that are zero or negative are left out, as are values not above \code{log2_intensity_cutoff} (by default 0). The protein files contain the intensities as read; the normalization, a shift of each sample in log2 space, is stored in the \code{normalization} column of \code{_sample.tsv} and applied by \code{\link{process_iq_format}}. If \code{peptide_extractor} is given, the fragment statistics of \code{\link{process}} (\code{n_fragments}, \code{n_peptides}, \code{peptides}, \code{XIC} and, with \code{rfasta}, \code{coverage}) are computed while the input is read, from all values read as in \code{\link{process}} (also those left out by the preprocessing), and stored as further columns of \code{_protein.tsv}; \code{\link{process_iq_format}} writes them with the protein annotation. Results are the same as with \code{\link{process}} on the same input.
}
\value{
\code{NULL}, invisibly. The output is written to the folder \code{output_iqformat_folder}: \code{_protein.tsv} (proteins and annotation, and the fragment statistics if \code{peptide_extractor} is given), \code{_sample.tsv} (samples and their normalization), \code{_ion.tsv} (ions), and one file \code{<i>.tsv} for each protein with values (\code{i} its row in \code{_protein.tsv}), with the sample index, ion index and intensity of each value.
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198-1203.}
\author{
Thang V. Pham
}
\seealso{
\code{\link{fast_read}}, \code{\link{fast_read_parquet}}, \code{\link{fast_preprocess}}, \code{\link{process_iq_format}}, \code{\link{process}}
}


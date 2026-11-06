\name{process}
\alias{process}
\title{
Protein quantification from a long-format table
}
\description{
Reads a long-format table of ion intensities (a tab-separated text file, a parquet file or a data frame), filters and normalizes the data, quantifies all proteins with one of the methods of \code{\link{process_matrix}}, and writes a table with one row per protein and one column per sample.
}
\usage{
process(input_data,
        output_filename,
        sample_id = "Run",
        primary_id = "Protein.Group",
        secondary_id = "Precursor.Id",
        intensity_col = "Precursor.Normalised",
        annotation_col = NULL,
        filter_string_equal = NULL,
        filter_string_not_equal = NULL,
        filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.01"),
        filter_double_greater = NULL,
        intensity_col_sep = NULL,
        intensity_col_id = NULL,
        na_string = "0",
        normalization = "median",
        log2_intensity_cutoff = 0,
        pdf_out = "qc-plots.pdf", pdf_width = 12, pdf_height = 8,
        show_boxplot = TRUE, boxplot_max_points = 10000,
        peptide_extractor = NULL,
        rfasta = NULL,
        method = "maxlfq_bit",
        p1 = NULL, p2 = NULL, k = 1.345, min_M = 15,
        rescale_method = "to_mean",
        n_threads = -1,
        verbose = TRUE)
}
\arguments{
  \item{input_data}{A data frame, the name of a parquet file (ending in \code{.parquet}), or the name of a tab-separated text file. See \code{\link{fast_read_data_frame}}, \code{\link{fast_read_parquet}} and \code{\link{fast_read}}.}
  \item{output_filename}{Output filename.}
  \item{sample_id}{The column of sample names.}
  \item{primary_id}{The column of protein identifiers. Unique values form the list of proteins.}
  \item{secondary_id}{The column(s) identifying an ion within a protein, concatenated with \code{"_"}.}
  \item{intensity_col}{The column of intensities. Several columns can be given, one per fragment ion, e.g. \code{paste0("Fr.", 0:11, ".Quantity")} in DIA-NN 2.}
  \item{annotation_col}{Annotation columns, written after the protein identifier.}
  \item{filter_string_equal}{A named vector of strings. Only rows where the column equals the value are kept.}
  \item{filter_string_not_equal}{A named vector of strings. Only rows where the column differs from the value are kept.}
  \item{filter_double_less}{A named vector of numbers. Rows with a value less than the given value, or a missing value, are kept.}
  \item{filter_double_greater}{A named vector of numbers. Rows with a value greater than the given value, or a missing value, are kept.}
  \item{intensity_col_sep}{A separator character when a single text intensity column contains multiple values (e.g. \code{";"} for the fragment intensities in \code{Fragment.Quant.Corrected} of DIA-NN 1.x), or \code{NULL} (the default) for one value per column.}
  \item{intensity_col_id}{Identities of multiple intensity values: one column per intensity column, or a single text column split by \code{intensity_col_sep}. If \code{NULL}, the values are numbered 1, 2, 3, ...}
  \item{na_string}{An intensity value considered missing, in addition to NA, NaN and empty values. For a text file, it is used with \code{intensity_col_sep} or several intensity columns.}
  \item{normalization}{\code{"median"} for median normalization (see \code{\link{fast_preprocess}}) or \code{"none"}.}
  \item{log2_intensity_cutoff}{Values with a log2 intensity not above this cutoff are removed; the default 0 removes intensities of 1 or less. \code{NULL} for no cutoff. See \code{\link{fast_preprocess}}.}
  \item{pdf_out}{A PDF file for quality control plots, or \code{NULL} for none. See \code{\link{fast_preprocess}}.}
  \item{pdf_width}{Width of the PDF pages.}
  \item{pdf_height}{Height of the PDF pages.}
  \item{show_boxplot}{Whether to draw boxplots of the intensities per sample.}
  \item{boxplot_max_points}{The maximum number of values per sample in the boxplots (evenly spaced values are used for a larger sample), which keeps the plots fast and small for large data; the median normalization always uses all values. \code{NULL} for all values.}
  \item{peptide_extractor}{A function that maps ion identifiers to peptides, to compute fragment statistics per protein. It is called with the ions of one protein and must return one character value per ion. See Details.}
  \item{rfasta}{A named list or a named character vector of protein sequences, named by protein accession as in \code{primary_id} (a protein group is split at \code{;}), to compute the sequence coverage. Used with \code{peptide_extractor}. See \code{\link{read_fasta_file}}.}
  \item{method}{The quantification method: one of the methods of \code{\link{process_matrix}} (\code{"maxlfq"}, \code{"maxlfq_bit"}, \code{"weighted_maxlfq"}, \code{"median_polish"}, \code{"weighted_median_polish"}, \code{"rlm"}, \code{"weighted_rlm"}, \code{"sum"}, \code{"mean"} or \code{"topN"}), or a function. See Details.}
  \item{p1}{See \code{p1} in \code{\link{process_matrix}}.}
  \item{p2}{See \code{p2} in \code{\link{process_matrix}}.}
  \item{k}{See \code{k} in \code{\link{process_matrix}}.}
  \item{min_M}{See \code{min_M} in \code{\link{process_matrix}}.}
  \item{rescale_method}{See \code{method} in \code{\link{rescale}}. The default \code{"to_mean"} shifts each connected component of samples so that the mean of its estimates equals the mean of its data (in log2 space). \code{NULL} for no rescaling.}
  \item{n_threads}{The number of threads, used in parallel over proteins. See \code{n_threads} in \code{\link{process_matrix}}.}
  \item{verbose}{A logical value. The default \code{TRUE} prints progress and summary messages; \code{FALSE} prints none of them. Warnings and errors are always shown.}
}
\details{
\strong{Input.} The defaults are the same for all inputs: the precursor quantities \code{Precursor.Normalised} per \code{Run} of a DIA-NN report. A data frame is read in place, without copying (\code{\link{fast_read_data_frame}}). A parquet file is read one row group at a time, only the needed columns (\code{\link{fast_read_parquet}}; requires the \pkg{arrow} package). Any other input is read as a tab-separated text file (\code{\link{fast_read}}). Filtering is the same for all inputs; a missing filter value (NA, NaN, empty, or text that is not a number) keeps the row.

\strong{Preprocessing.} Intensities are log2-transformed; zero and negative intensities are removed, as are values not above \code{log2_intensity_cutoff} (by default 0, i.e. intensities of 1 or less); and samples are median normalized (\code{\link{fast_preprocess}}).

\strong{Quantification.} For each protein, a matrix of log2 intensities is formed with ions as rows (in input order) and all samples as columns; an ion with more than one value in a sample gets the log2 of the summed intensities. The matrix is quantified with \code{method} and the result rescaled with \code{rescale_method}, as in \code{\link{process_iq_format}}, without writing iq-format files. For the built-in methods, all proteins are processed in C++ in parallel. A function given as \code{method} is called in R for each protein with this matrix, with the ion identifiers as row names and the sample names as column names, and must return a numeric vector with one value per sample.
Values that are not finite (\code{-Inf}, \code{Inf}, \code{NaN}) become \code{NA}. For example, \code{function(X) log2(colSums(2^X, na.rm = TRUE))} with \code{rescale_method = NULL} gives the log2 summed intensities per sample, with \code{NA} for a sample without data, as the method \code{"sum"}.

\strong{Output.} A tab-separated file with one row for each protein that has data after preprocessing. The columns are \code{primary_id}, the \code{annotation_col} columns, the fragment statistics (if \code{peptide_extractor} is given), and one column per sample with the log2 estimates; a sample without data for a protein is \code{NA}.

\strong{Fragment statistics.} If \code{peptide_extractor} is not \code{NULL}, the columns \code{n_fragments} (the number of MS2 fragment ions or precursors depending on \code{intensity_col}), \code{n_peptides} (the number of unique values of \code{peptide_extractor} on the ions), \code{peptides} and \code{XIC} (the summed intensities, before preprocessing) are added, and \code{coverage} (sequence coverage in percent for each protein in the group) if \code{rfasta} is given. \code{peptide_extractor} is called with the \code{secondary_id} part of the ion identifiers, without the suffix of the intensity values (\code{1}, \code{2}, ... or the \code{intensity_col_id} values), e.g. a precursor id such as \code{(UniMod:1)AAAAAAGAGPEMVR2}, also when \code{intensity_col} has one column per fragment ion; the \code{ion_prefix} component returned by \code{\link{fast_read}}, \code{\link{fast_read_parquet}} or \code{\link{fast_read_data_frame}} shows these values. An example extractor for DIA-NN precursor ids is \code{function(x) gsub("\\\\(UniMod:\\\\d+\\\\)|\\\\d", "", x)}, which removes the modifications and the charge state: \code{(UniMod:1)AAAAAAGAGPEMVR2} becomes \code{AAAAAAGAGPEMVR}.
}
\value{
\code{NULL}, invisibly. The result is written to \code{output_filename}.
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198-1203.
}
\author{
Thang V. Pham
}
\seealso{
\code{\link{fast_read}}, \code{\link{fast_read_parquet}}, \code{\link{fast_read_data_frame}}, \code{\link{fast_preprocess}}, \code{\link{process_matrix}}, \code{\link{rescale}}, \code{\link{process_iq_format}}, \code{\link{read_fasta_file}}
}

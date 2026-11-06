#########################################################################
#
# Author: Thang V. Pham {t.pham@amsterdamumc.nl}
#
# All rights reserved.
#
# Citation:
#
# Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, 
# Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein 
# quantification algorithms in mass spectrometry-based proteomics.
# Journal of Proteome Research 2026 Feb 6; 25(2):1198-1203.
#
# Software version: 2.2.0
#
#########################################################################

# Read a FASTA file into a named list of sequences, e.g. for 'rfasta' in process().
# The name of each sequence is taken from its header line with regex_id (the first
# group); a header that does not match keeps the whole line without '>'.
read_fasta_file <- function(fasta_file_name, regex_id = "^>.*\\|(.*)\\|.*") {

    lines <- readLines(fasta_file_name, warn = FALSE)   # also without a final newline

    header <- startsWith(lines, ">")
    if (!any(header)) {
        stop("No sequence header (a line starting with '>') in ", fasta_file_name)
    }

    # sequence lines belong to the header above them; lines before the first header are ignored
    entry <- cumsum(header)
    sequence_line <- !header & entry > 0
    sequences <- vapply(split(gsub("[[:space:]]", "", lines[sequence_line]),
                              factor(entry[sequence_line], levels = seq_len(sum(header)))),
                        paste, character(1), collapse = "")

    ids <- sub(regex_id, "\\1", lines[header])
    ids <- sub("^>", "", ids)   # not matched (e.g. not UniProt): the whole header line

    ret <- as.list(unname(sequences))
    names(ret) <- ids
    ret
}

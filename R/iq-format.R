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

# quantification methods of process_matrix (method code in iq_quant)
.quant_methods <- c(
    "maxlfq" = 0,
    "maxlfq_bit" = 1,
    "weighted_maxlfq" = 2,
    "median_polish" = 3,
    "weighted_median_polish" = 4,
    "rlm" = 5,
    "weighted_rlm" = 6,
    "sum" = 7,
    "mean" = 8,
    "topN" = 9)

.rescale_methods <- c("median_mean", "to_median", "to_mean", "to_sum")

process <- function(input_data,
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
                    pdf_out = "qc-plots.pdf",
                    pdf_width = 12,
                    pdf_height = 8,
                    show_boxplot = TRUE,
                    boxplot_max_points = 10000,
                    peptide_extractor = NULL,
                    rfasta = NULL,
                    method = "maxlfq_bit", 
                    p1 = NULL, p2 = NULL, k = 1.345, min_M = 15,
                    rescale_method = "to_mean",
                    n_threads = -1,
                    verbose = TRUE) {


    .check_input_data(input_data)

    # output files
    if (missing(output_filename)) {
        stop("output_filename is missing: the name of the output file.")
    }
    .check_output_file(output_filename, "output_filename")
    if (!is.null(pdf_out)) {
        .check_output_file(pdf_out, "pdf_out")
    }

    # check the input file, output_filename and pdf_out must be different files
    if (is.character(input_data)) {
        if (.same_file(output_filename, input_data)) {
            stop("output_filename is the input file (", input_data, "); choose another name.")
        }
        if (!is.null(pdf_out) && .same_file(pdf_out, input_data)) {
            stop("pdf_out is the input file (", input_data, "); choose another name.")
        }
    }
    if (!is.null(pdf_out) && .same_file(pdf_out, output_filename)) {
        stop("pdf_out and output_filename are the same file (", output_filename, "); choose another name.")
    }

    median_normalization <- .median_normalization(normalization)

    # check the methods and their parameters before reading the data
    .check_methods(method, rescale_method)
    .check_parameters(p1, p2, k, min_M, n_threads)
    .check_id_columns(sample_id, primary_id, secondary_id)
    .check_filters(filter_string_equal, filter_string_not_equal, filter_double_less, filter_double_greater)
    .check_intensity_columns(intensity_col, intensity_col_sep, intensity_col_id)
    .check_cutoff(log2_intensity_cutoff)

    .check_fragment_options(peptide_extractor, rfasta)

    iq_dat <- .read_long_format(input_data,
                                primary_id = primary_id,
                                sample_id  = sample_id,
                                secondary_id = secondary_id,
                                intensity_col = intensity_col,
                                intensity_col_sep = intensity_col_sep,
                                intensity_col_id = intensity_col_id,
                                annotation_col = annotation_col,
                                filter_string_equal = filter_string_equal,
                                filter_string_not_equal = filter_string_not_equal,
                                filter_double_less = filter_double_less,
                                filter_double_greater = filter_double_greater,
                                na_string = na_string,
                                verbose = verbose)

    iq_dat <- .drop_unidentified(iq_dat, primary_id, sample_id)

    # samples in the input without any value (all missing or filtered out): columns of NA
    no_values <- setdiff(seq_along(iq_dat$sample), as.integer(iq_dat$quant_table$sample_list))
    if (length(no_values) > 0) {
        warning("No values for sample(s) ", .name_list(iq_dat$sample[no_values]),
                " (all missing, filtered out or without a protein identifier); their estimates are NA.", call. = FALSE)
    }

    # a sample left without values is named here (fast_preprocess only knows its index)
    iq_norm_data <- withCallingHandlers(
        fast_preprocess(iq_dat$quant_table,
                        median_normalization = median_normalization,
                        log2_intensity_cutoff = log2_intensity_cutoff,
                        pdf_out = pdf_out,
                        pdf_width = pdf_width,
                        pdf_height = pdf_height,
                        show_boxplot = show_boxplot,
                        boxplot_max_points = boxplot_max_points,
                        verbose = verbose),
        iq_empty_samples = function(w) {
            warning("No values left after removing zero, negative and low intensities in sample(s) ",
                    .name_list(iq_dat$sample[w$samples]),
                    "; their estimates are NA.", call. = FALSE)
            invokeRestart("muffleWarning")
        })

    n_proteins <- nrow(iq_dat$protein)
    samples <- as.character(iq_dat$sample)

    # the rows of each protein, in input order
    protein_rows <- split(seq_along(iq_norm_data$protein_list),
                          factor(iq_norm_data$protein_list, levels = seq_len(n_proteins)))

    # proteins with data after preprocessing
    quantified <- which(lengths(protein_rows) > 0)
    if (length(quantified) == 0) {
        stop("No quantitative values left after preprocessing.")
    }

    if (verbose) {
        message("Quantifying ", length(quantified), " proteins: ",
                if (is.function(method)) "user-defined method" else method,
                if (!is.null(rescale_method)) paste0(", rescaling: ", rescale_method))
    }

    if (is.function(method)) {
        estimate <- matrix(NA_real_, nrow = length(quantified), ncol = length(samples))
    } else {
        # all proteins in C++, in parallel, including the rescaling; matrices built as in process_iq_format
        parameters <- .method_parameters(method, p1, p2)
        estimate <- .Call("iq_process",
                          list(protein_index = as.integer(iq_norm_data$protein_list),
                               ion_index = as.integer(iq_norm_data$id),
                               sample_index = as.integer(iq_norm_data$sample_list),
                               quant = as.double(iq_norm_data$quant),
                               n_proteins = as.integer(n_proteins),
                               n_samples = length(samples),
                               method = as.integer(.quant_methods[method]),
                               p1 = as.double(parameters$p1),
                               p2 = as.integer(parameters$p2),
                               k = as.double(k),
                               min_M = as.integer(min_M),
                               n_threads = as.integer(n_threads),
                               rescale_method = if (is.null(rescale_method)) 0L else
                                   match(rescale_method, .rescale_methods),
                               verbose = verbose))
        if (is.null(estimate)) {
            stop("Cancelled.")
        }
        estimate <- estimate[quantified, , drop = FALSE]
    }
    colnames(estimate) <- samples

    # user-defined method: protein by protein in R
    if (is.function(method)) {

        n_display <- length(quantified) %/% 20

        for (r in seq_along(quantified)) {

            if (verbose && n_display > 0 && r %% n_display == 0) {
                message(format(r * 100 / length(quantified), digits = 2), "%")
            }

            ind <- protein_rows[[quantified[r]]]
            X <- .protein_matrix(iq_norm_data$id[ind], iq_norm_data$sample_list[ind],
                                 iq_norm_data$quant[ind], length(samples))
            # ion identifiers (rows, in the order of .protein_matrix) and sample names (columns)
            dimnames(X) <- list(as.character(iq_dat$ion[unique(iq_norm_data$id[ind])]), samples)

            estimate[r, ] <- .for_protein(process_matrix(X, method = method, rescale_method = rescale_method),
                                          iq_dat$protein[quantified[r], 1])
        }
        if (verbose) {
            message("Completed.")  # the built-in methods: printed by iq_process
        }
    }

    if (!is.null(peptide_extractor)) {
        # all values read, before preprocessing
        raw_rows <- split(seq_along(iq_dat$quant_table$protein_list),
                          factor(iq_dat$quant_table$protein_list, levels = seq_len(n_proteins)))
        s <- .fragment_statistics(iq_dat$protein[quantified, 1],
                                  lapply(quantified, function(i) unique(iq_dat$quant_table$id[raw_rows[[i]]])),
                                  vapply(quantified, function(i) sum(iq_dat$quant_table$quant[raw_rows[[i]]]), numeric(1)),
                                  iq_dat$ion_prefix, peptide_extractor, rfasta, verbose)
    }

    if (verbose) {
        message("Writing to: ", output_filename)
    }

    # protein id and annotation columns, then statistics, then estimates
    tab <- iq_dat$protein[quantified, , drop = FALSE]
    if (!is.null(peptide_extractor)) {
        tab <- cbind(tab, s)
    }
    tab <- cbind(tab, estimate)

    write.table(tab, output_filename, sep = "\t", row.names = FALSE, quote = FALSE)

    invisible(NULL)
}

# Read a long-format data frame (also a tibble or data.table) with the same rules as
# fast_read() and fast_read_parquet(). Returns the same structure as fast_read().
# The columns are read in place, without copying the data frame.
fast_read_data_frame <- function(input_data,
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
                                 verbose = TRUE) {

    if (!is.data.frame(input_data)) {
        stop("input_data must be a data frame")
    }

    options <- .stream_options(sample_id, primary_id, secondary_id, intensity_col, annotation_col,
                               filter_string_equal, filter_string_not_equal,
                               filter_double_less, filter_double_greater,
                               intensity_col_sep, intensity_col_id, na_string)

    columns <- .stream_columns(options)
    missing_columns <- setdiff(columns, names(input_data))
    if (length(missing_columns) > 0) {
        stop("The input data frame has no column(s): ", paste(missing_columns, collapse = ", "))
    }

    state <- .Call("iq_stream_begin", options)
    .Call("iq_stream_add", state, input_data)
    .stream_end(state, verbose)
}

# Read a long-format parquet file (e.g. DIA-NN report.parquet) one row group at a time.
# Returns the same structure as fast_read(). Needs the 'arrow' package.
# intensity_col: one column (split by intensity_col_sep if it is text), or several
#                columns, one per fragment (e.g. paste0("Fr.", 0:11, ".Quantity")).
# intensity_col_id: optional ion ids, one column per intensity column, or a single
#                text column split by intensity_col_sep.
fast_read_parquet <- function(filename,
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
                              verbose = TRUE) {

    if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("Reading parquet files requires the 'arrow' package: install.packages(\"arrow\")")
    }

    if (!file.exists(filename)) {
        stop("Cannot open: ", filename)
    }

    options <- .stream_options(sample_id, primary_id, secondary_id, intensity_col, annotation_col,
                               filter_string_equal, filter_string_not_equal,
                               filter_double_less, filter_double_greater,
                               intensity_col_sep, intensity_col_id, na_string)

    columns <- .stream_columns(options)

    reader <- arrow::ParquetFileReader$create(filename)
    column_indices <- match(columns, names(reader$GetSchema()))
    if (any(is.na(column_indices))) {
        stop("Cannot find column(s) in ", filename, ": ",
             paste(columns[is.na(column_indices)], collapse = ", "))
    }

    if (verbose) {
        message("Reading: ", filename)
    }

    # check the options before reading
    state <- .Call("iq_stream_begin", options)

    # one row group at a time, only the needed columns (0-based indices)
    for (g in seq_len(reader$num_row_groups) - 1L) {
        tab <- reader$ReadRowGroup(g, column_indices = column_indices - 1L)
        if (tab$num_rows > 0) {
            .Call("iq_stream_add", state, as.data.frame(tab))
        }
        rm(tab)
    }

    .stream_end(state, verbose)
}

long_format_to_iq_format <- function(input_filename,
                                     output_iqformat_folder,
                                     sample_id = "Run",
                                     primary_id = "Protein.Group",
                                     secondary_id = c("Precursor.Id"),
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
                                     pdf_out = "qc-plots-iq-format.pdf",
                                     pdf_width = 12,
                                     pdf_height = 8,
                                     show_boxplot = TRUE,
                                     peptide_extractor = NULL,
                                     rfasta = NULL,
                                     buffer_mb = 1024,
                                     verbose = TRUE) {

    # The input is read chunk by chunk and written to the protein files as it comes, so it is
    # never held in memory as a whole: the files hold the intensities as read (zero, negative
    # and low intensities left out), and the median normalization is stored per sample in
    # _sample.tsv, applied by process_iq_format. The sample medians are computed afterwards
    # from the protein files, reading the values of as many samples at a time as fit in
    # buffer_mb. buffer_mb is used in both steps (one after the other) and does not change
    # the result.

    if (is.data.frame(input_filename)) {
        stop("long_format_to_iq_format reads a text file or a parquet file; for a data frame, use process().")
    }
    if (!is.character(input_filename) || length(input_filename) != 1 || is.na(input_filename)) {
        stop("input_filename must be the name of one text or parquet file.")
    }
    .check_id_columns(sample_id, primary_id, secondary_id)
    .check_filters(filter_string_equal, filter_string_not_equal, filter_double_less, filter_double_greater)
    .check_intensity_columns(intensity_col, intensity_col_sep, intensity_col_id)

    median_normalization <- .median_normalization(normalization)

    .check_cutoff(log2_intensity_cutoff)
    if (!is.null(pdf_out)) {
        .check_output_file(pdf_out, "pdf_out")
    }
    if (!is.numeric(buffer_mb) || length(buffer_mb) != 1 || is.na(buffer_mb) || buffer_mb <= 0) {
        stop("buffer_mb must be a positive number (megabytes).")
    }
    .check_fragment_options(peptide_extractor, rfasta)

    # the protein files are appended to: a new or empty directory
    if (!is.character(output_iqformat_folder) || length(output_iqformat_folder) != 1 || is.na(output_iqformat_folder) || !nzchar(output_iqformat_folder)) {
        stop("output_iqformat_folder must be the name of a folder.")
    }
    if (file.exists(output_iqformat_folder)) {
        if (!dir.exists(output_iqformat_folder)) {
            stop("output_iqformat_folder exists and is not a folder: ", output_iqformat_folder)
        }
        if (length(list.files(output_iqformat_folder, all.files = TRUE, no.. = TRUE)) > 0) {
            stop("output_iqformat_folder is not empty: ", output_iqformat_folder, "; the iq format is written to a new or empty folder.")
        }
    } else if (!dir.create(output_iqformat_folder, recursive = TRUE)) {
        stop("Cannot create the folder ", output_iqformat_folder)
    }
    out_dir <- paste0(output_iqformat_folder, "/")

    writer <- .Call("iq_writer_begin",
                    list(dir = out_dir,
                         cutoff = if (is.null(log2_intensity_cutoff)) NULL else as.double(log2_intensity_cutoff),
                         buffer = as.double(buffer_mb),
                         statistics = !is.null(peptide_extractor)))

    iq_data <- .read_long_format(input_filename,
                                 primary_id = primary_id,
                                 sample_id  = sample_id,
                                 secondary_id = secondary_id,
                                 intensity_col = intensity_col,
                                 intensity_col_sep = intensity_col_sep,
                                 intensity_col_id = intensity_col_id,
                                 annotation_col = annotation_col,
                                 filter_string_equal = filter_string_equal,
                                 filter_string_not_equal = filter_string_not_equal,
                                 filter_double_less = filter_double_less,
                                 filter_double_greater = filter_double_greater,
                                 na_string = na_string,
                                 writer = writer,
                                 verbose = verbose)

    info <- .Call("iq_writer_end", writer)

    if (info$n_unidentified > 0) {
        warning("Left out ", info$n_unidentified, " value(s) without a protein identifier (", primary_id,
                ") or a sample identifier (", sample_id, ").", call. = FALSE)
    }
    if (info$n_written == 0) {
        stop("No quantitative values left: all intensities are ",
             if (is.null(log2_intensity_cutoff)) "missing, zero or negative."
             else "missing, zero, negative or not above the intensity cut-off.")
    }

    samples <- as.character(iq_data$sample)
    n_samples <- length(samples)
    count <- c(info$sample_count, rep(0, n_samples - length(info$sample_count)))

    if (verbose) {
        message("Values written: ", info$n_written,
                "\n  (left out: ", info$n_not_positive, " zero, negative or missing",
                if (!is.null(log2_intensity_cutoff)) paste0(",\n  ", info$n_below_cutoff, " not above the cut-off"),
                ");\n proteins with values: ", sum(info$protein_has_file), " of ", nrow(iq_data$protein), ".")
    }

    #--- sample medians, from the protein files: groups of samples whose values (8 bytes each)
    #    fit in buffer_mb (at least one sample), one pass over the files per group

    boxplots <- !is.null(pdf_out) && show_boxplot

    medians <- rep(NA_real_, n_samples)
    box_stats <- matrix(NA_real_, nrow = 5, ncol = n_samples)  # boxplot statistics of all values

    with_files <- which(info$protein_has_file)
    with_count <- which(count > 0)
    values_per_pass <- buffer_mb * 2^20 / 8
    group <- cumsum(count[with_count]) %/% values_per_pass   # consecutive samples
    if (verbose) {
        message("Computing sample medians (", length(unique(group)), " pass(es) over the protein files)...")
    }
    for (g in unique(group)) {
        wanted <- with_count[group == g]
        values <- .Call("iq_read_sample_values", out_dir, as.integer(with_files),
                        as.integer(wanted), as.double(count[wanted]))
        for (k in seq_along(wanted)) {
            s <- wanted[k]
            medians[s] <- median(values[[k]])  # the values written: above the cut-off, if any
            if (boxplots) {
                box_stats[, s] <- boxplot.stats(values[[k]])$stats
            }
            values[k] <- list(NULL)
        }
        rm(values)
    }

    with_values <- !is.na(medians)
    if (any(!with_values)) {
        warning("No values for sample(s) ", .name_list(samples[!with_values]),
                " (all missing, filtered out, zero, negative or not above the cut-off); ",
                "they are left out of the normalization.", call. = FALSE)
    }

    # the median normalization of fast_preprocess, a shift of each sample in log2 space
    normalization_shift <- rep(0, n_samples)
    if (median_normalization) {
        if (verbose) {
            message("Median normalization ...")
        }
        normalization_shift[with_values] <- mean(medians[with_values]) - medians[with_values]
    }

    if (!is.null(pdf_out)) {
        pdf(pdf_out, pdf_width, pdf_height)
        pdf_device <- dev.cur()
        on.exit(dev.off(pdf_device), add = TRUE)  # also on an error

        # histogram of all finite log2 values (also those not above the cut-off), about 100
        # bars made of the writer's bins of width 1 / hist_bins_per_unit
        per_unit <- info$hist_bins_per_unit
        merge <- max(1, ceiling((max(info$hist_bin) - min(info$hist_bin) + 1) / 100))
        bar <- floor(info$hist_bin / merge)
        bars <- seq(min(bar), max(bar))
        hist_counts <- as.numeric(tapply(info$hist_count, factor(bar, levels = bars), sum))
        hist_counts[is.na(hist_counts)] <- 0
        breaks <- c(bars, max(bars) + 1) * merge / per_unit
        h <- structure(list(breaks = breaks, counts = hist_counts,
                            density = hist_counts / (sum(hist_counts) * diff(breaks)),
                            mids = (breaks[-1] + breaks[-length(breaks)]) / 2,
                            xname = "log2 intensity", equidist = TRUE),
                       class = "histogram")
        plot(h, col = "steelblue", border = "steelblue", freq = FALSE,
             main = "Histogram of log2 intensities", xlab = "log2 intensity")
        if (!is.null(log2_intensity_cutoff)) {
            arrows(log2_intensity_cutoff, 0, log2_intensity_cutoff, max(h$density) / 2.0, col = "red", code = 1, lwd = 2)
        }

        # the boxplots of fast_preprocess, from the statistics of all values of each sample
        boxplot_samples <- function(shift) {
            bxp(list(stats = box_stats[, with_values, drop = FALSE] + rep(shift[with_values], each = 5),
                     n = count[with_values],
                     names = samples[with_values]),
                main = "Boxplot of fragment intensities per sample",
                ylab = "log2 intensity",
                outline = FALSE,
                boxfill = "steelblue",
                whisklty = 1,
                staplelty = 0,
                las = 2)
        }
        if (boxplots) {
            boxplot_samples(rep(0, n_samples))
            if (median_normalization) {
                boxplot_samples(normalization_shift)
            }
        }
    }

    #--- the lists of proteins, samples (with the normalization) and ions
    # (quoted text; a quote in a name is doubled, as read.delim in process_iq_format expects)

    # fragment statistics as in process(), for the proteins with values: further columns of
    # _protein.tsv, which process_iq_format writes with the protein annotation
    protein_table <- iq_data$protein
    if (!is.null(peptide_extractor)) {
        proteins_with_values <- which(c(info$protein_has_file, logical(nrow(protein_table)))[seq_len(nrow(protein_table))])
        stats <- .fragment_statistics(protein_table[proteins_with_values, 1], info$protein_ions[proteins_with_values],
                                      info$protein_xic[proteins_with_values], iq_data$ion_prefix,
                                      peptide_extractor, rfasta, verbose)
        all_stats <- matrix("", nrow = nrow(protein_table), ncol = ncol(stats),
                            dimnames = list(NULL, colnames(stats)))
        all_stats[proteins_with_values, ] <- stats
        protein_table <- cbind(protein_table, all_stats)
    }

    write.table(protein_table,
                file = paste0(out_dir, "_protein.tsv"), sep = "\t", row.names = FALSE, qmethod = "double")

    write.table(data.frame(sample = samples, normalization = normalization_shift),
                file = paste0(out_dir, "_sample.tsv"), sep = "\t", row.names = FALSE, qmethod = "double")

    ion <- iq_data$ion
    colnames(ion) <- "ion"
    write.table(ion,
                file = paste0(out_dir, "_ion.tsv"), sep = "\t", row.names = FALSE, qmethod = "double")

    if (verbose) {
        message("Done.")
    }

    invisible(NULL)
}

connected_component <- function(X) {

    X <- .as_matrix(X)

    if (ncol(X) == 1) {
        return(0)
    }
    
    .Call("iq_connected_components",
          list(X = as.double(X),
               M = as.integer(nrow(X)),
               N = as.integer(ncol(X))))
}

# Rescale the quantification x of one protein with data X (log2, ions x samples):
# each connected component of samples separately. 
rescale <- function(x, X, method) {

    code <- match(method, .rescale_methods)
    if (length(method) != 1 || is.na(code)) {
        stop("Unknown rescaling method. Use one of: ", paste(.rescale_methods, collapse = ", "))
    }

    X <- as.matrix(.as_matrix(X))
    if (length(x) != ncol(X)) {
        stop("x must have one value per column of X")
    }

    x_new <- .Call("iq_rescale",
                   list(X = as.double(X),
                        M = as.integer(nrow(X)),
                        N = as.integer(ncol(X)),
                        x = as.double(x),
                        method = as.integer(code)))
    names(x_new) <- names(x)
    x_new
}

huber_weight <- function(X, x, c = 1.4826, k = 1.345) {

    X <- .as_matrix(X)

    if (nrow(X) == 1) {
        return(1)
    }
    
    .Call("iq_huber_weight",
          list(X = as.double(X),
               M = as.integer(nrow(X)),
               N = as.integer(ncol(X)),
               x = as.double(x),
               c = as.double(c),
               k = as.double(k)))
}

process_matrix <- function(X, method = "maxlfq_bit",
                           p1 = NULL, p2 = NULL, k = 1.345, min_M = 15,
                           rescale_method = "to_mean",
                           n_threads = -1) {

    .check_methods(method, rescale_method)
    .check_parameters(p1, p2, k, min_M, n_threads)

    X <- .as_matrix(X)

    if (is.function(method)) {
        x <- .user_method(method, X)
    } else if (nrow(X) == 1) {
        x <- as.numeric(X[1, ])
    } else {
        parameters <- .method_parameters(method, p1, p2)

        x <- .Call("iq_quant",
                   list(X = as.double(X),
                        M = as.integer(nrow(X)),
                        N = as.integer(ncol(X)),
                        method = as.integer(.quant_methods[method]),
                        p1 = as.double(parameters$p1),
                        p2 = as.integer(parameters$p2),
                        k = as.double(k),
                        min_M = as.integer(min_M),
                        n_threads = as.integer(n_threads)))$x
    }

    if (!is.null(rescale_method)) {
        x <- rescale(x, X, method = rescale_method)
    }

    x
}

process_iq_format <- function(iqformat_folder,
                              output_filename = NULL,
                              method = "maxlfq_bit", 
                              p1 = NULL, p2 = NULL, k = 1.345, min_M = 15,
                              rescale_method = "to_mean",
                              buffer_mb = 1024,
                              n_threads = -1,
                              verbose = TRUE) {

    .check_methods(method, rescale_method)
    .check_parameters(p1, p2, k, min_M, n_threads)
    if (!is.null(buffer_mb) &&
        (!is.numeric(buffer_mb) || length(buffer_mb) != 1 || is.na(buffer_mb) || buffer_mb <= 0)) {
        stop("buffer_mb must be a positive number (megabytes), or NULL to quantify protein by protein in R.")
    }

    # without a trailing / or \ ("folder/" -> "folder"), for the default output name
    iqformat_folder <- .no_trailing_slash(iqformat_folder)
    dir_in <- paste0(iqformat_folder, "/")

    p <- read.delim(paste0(dir_in, "_protein.tsv"), colClasses = "character")
    s <- read.delim(paste0(dir_in, "_sample.tsv"), colClasses = "character")

    # the median normalization of each sample (log2), written by long_format_to_iq_format
    if (!all(c("sample", "normalization") %in% names(s))) {
        stop(dir_in, "_sample.tsv must have the columns sample and normalization, ",
             "as written by long_format_to_iq_format.")
    }
    shift <- as.numeric(s$normalization)

    res <- matrix(NA, nrow = nrow(p), ncol = nrow(s))
    colnames(res) <- s$sample

    # a protein without values after preprocessing has no file
    has_values <- file.exists(paste0(dir_in, seq_len(nrow(p)), ".tsv"))

    if (!is.null(buffer_mb) && !is.function(method)) {

        # batches of proteins in C++: the protein files are read until a batch holds about
        # buffer_mb of values (20 bytes each), and the proteins of a batch are quantified in
        # parallel by iq_process, as in process(); the same results as protein by protein below
        parameters <- .method_parameters(method, p1, p2)
        proteins <- which(has_values)
        max_values <- buffer_mb * 2^20 / 20
        start <- 1L
        n_batch <- 0L
        while (start <= length(proteins)) {
            n_batch <- n_batch + 1L
            if (verbose) {
                message("Batch ", n_batch, ": reading from protein ", start, " of ", length(proteins),
                        " (", format((start - 1) * 100 / length(proteins), digits = 2), "% done)")
            }
            batch <- .Call("iq_read_protein_batch", dir_in, as.integer(proteins), as.integer(start),
                           as.double(max_values), as.double(shift))
            rows <- proteins[start:(batch$next_start - 1L)]
            if (verbose) {
                message("Batch ", n_batch, ": quantifying proteins ", start, "-", batch$next_start - 1L)
            }
            estimate <- .Call("iq_process",
                              list(protein_index = batch$protein_index,
                                   ion_index = batch$ion_index,
                                   sample_index = batch$sample_index,
                                   quant = batch$quant,
                                   n_proteins = length(rows),
                                   n_samples = nrow(s),
                                   method = as.integer(.quant_methods[method]),
                                   p1 = as.double(parameters$p1),
                                   p2 = as.integer(parameters$p2),
                                   k = as.double(k),
                                   min_M = as.integer(min_M),
                                   n_threads = as.integer(n_threads),
                                   rescale_method = if (is.null(rescale_method)) 0L else
                                       match(rescale_method, .rescale_methods),
                                   verbose = FALSE))  # progress per batch above
            if (is.null(estimate)) {
                stop("Cancelled.")
            }
            res[rows, ] <- estimate
            start <- batch$next_start
            rm(batch, estimate)
        }
        if (verbose) {
            message("Completed.")
        }

    } else {

        # protein by protein in R (buffer_mb = NULL, or a user-defined method)

        n_display <- nrow(p) %/% 20

        for (j in 1:nrow(p)) {

            if (verbose && n_display > 0 && j %% n_display == 0) {
                message(format(j * 100 / nrow(p), digits = 2), "%")
            }

            # a protein without values after preprocessing has no file
            f <- paste0(dir_in, j, ".tsv")
            if (!has_values[j]) {
                next
            }

            q <- read.delim(f, header = FALSE)  # sample, ion, intensity

            # ions x samples: log2 of the summed intensities of an ion in a sample, NA without a value
            rows <- unique(q$V2)
            cell <- match(q$V2, rows) + (q$V1 - 1) * length(rows)
            summed <- rowsum(q$V3, cell)
            m <- matrix(0, nrow = length(rows), ncol = nrow(s))
            m[as.numeric(rownames(summed))] <- summed[, 1]
            rownames(m) <- as.character(rows)
            colnames(m) <- s$sample

            m <- log2(m)
            m[m == -Inf] <- NA
            m <- m + rep(shift, each = nrow(m))

            res[j, ] <- .for_protein(process_matrix(m, method = method, p1 = p1, p2 = p2, k = k, min_M = min_M,
                                                    n_threads = n_threads, rescale_method = rescale_method),
                                     p[j, 1])
        }

    }

    # proteins with values, and the samples (a nameless sample has no values)
    keep_samples <- s$sample != "" | colSums(!is.na(res[has_values, , drop = FALSE])) > 0

    write.table(cbind(p[has_values, , drop = FALSE], res[has_values, keep_samples, drop = FALSE]),
                file = ifelse(is.null(output_filename), paste0(iqformat_folder, "-quant.tsv"), output_filename),
                sep = "\t", row.names = FALSE, quote = FALSE)

    invisible(NULL)
}

iq_format_to_protein_list <- function(iqformat_folder, proteins, normalization = TRUE) {

    # The data of some proteins in the iq format, as create_protein_list returns them: a named
    # list with a matrix for each protein, ions in rows, samples in columns, log2 intensities
    # (an ion measured more than once in a sample summed, as in process_iq_format), NA if missing.

    if (!is.character(proteins) || length(proteins) == 0 || anyNA(proteins)) {
        stop("proteins must be one or more protein identifiers (the first column of _protein.tsv).")
    }
    if (!is.logical(normalization) || length(normalization) != 1 || is.na(normalization)) {
        stop("normalization must be TRUE or FALSE.")
    }
    if (!dir.exists(iqformat_folder)) {
        stop("Cannot find the iq format folder: ", iqformat_folder)
    }

    dir_in <- paste0(.no_trailing_slash(iqformat_folder), "/")

    # names as text, also a name "NA"
    read_text <- function(f) read.delim(paste0(dir_in, f), colClasses = "character", na.strings = character(0))
    p <- read_text("_protein.tsv")
    s <- read_text("_sample.tsv")
    ions <- read_text("_ion.tsv")[, 1]

    if (!all(c("sample", "normalization") %in% names(s))) {
        stop(dir_in, "_sample.tsv must have the columns sample and normalization, ",
             "as written by long_format_to_iq_format.")
    }
    # the median normalization of each sample (log2), written by long_format_to_iq_format
    shift <- if (normalization) as.numeric(s$normalization) else rep(0, nrow(s))

    rows_p <- match(proteins, p[, 1])
    if (anyNA(rows_p)) {
        stop("Not in ", dir_in, "_protein.tsv: ", .name_list(proteins[is.na(rows_p)]))
    }

    res <- vector("list", length(proteins))
    names(res) <- proteins

    for (k in seq_along(proteins)) {

        # a protein without values after preprocessing has no file
        f <- paste0(dir_in, rows_p[k], ".tsv")
        if (!file.exists(f)) {
            warning("No values for protein ", proteins[k], "; its matrix has no rows.", call. = FALSE)
            res[[k]] <- matrix(NA_real_, nrow = 0, ncol = nrow(s), dimnames = list(NULL, s$sample))
            next
        }

        q <- read.delim(f, header = FALSE)  # sample, ion, intensity

        # ions x samples: log2 of the summed intensities of an ion in a sample, NA without a value
        rows <- unique(q$V2)
        cell <- match(q$V2, rows) + (q$V1 - 1) * length(rows)
        summed <- rowsum(q$V3, cell)
        m <- matrix(0, nrow = length(rows), ncol = nrow(s))
        m[as.numeric(rownames(summed))] <- summed[, 1]

        m <- log2(m)
        m[m == -Inf] <- NA
        m <- m + rep(shift, each = nrow(m))

        dimnames(m) <- list(ions[rows], s$sample)
        res[[k]] <- m
    }

    res
}

#------------------------------ helpers ------------------------------------

# A folder name without trailing / or \ ("folder/" -> "folder"); a root such as "/" is kept
.no_trailing_slash <- function(path) {
    trimmed <- sub("[/\\\\]+$", "", path)
    ifelse(nzchar(trimmed), trimmed, path)
}

# X as a matrix: a data frame must have numeric columns (or logical ones, e.g. a column that
# read.delim read as all NA); other input is returned as it is
.as_matrix <- function(X) {
    if (is.data.frame(X)) {
        if (!all(vapply(X, function(v) is.numeric(v) || is.logical(v), logical(1)))) {
            stop("X must have numeric columns only.")
        }
        X <- as.matrix(X)
    }
    X
}

# The matrix of one protein for process_matrix, built as in process_iq_format:
# ions (rows) in input order, all samples (columns); log2 values, NA if missing.
# An ion measured more than once in a sample gets the sum of its intensities.
.protein_matrix <- function(ion, sample, quant, n_samples) {
    ions <- unique(ion)
    cell <- match(ion, ions) + (as.integer(sample) - 1L) * length(ions)

    X <- matrix(NA_real_, nrow = length(ions), ncol = n_samples)

    if (anyDuplicated(cell)) {
        total <- rowsum(2^quant, cell)
        X[as.integer(rownames(total))] <- log2(total[, 1])
    } else {
        X[cell] <- quant
    }

    X[!is.finite(X)] <- NA   # zero intensity (-Inf), as in process_iq_format
    X
}

# Call a user-defined quantification method on the matrix X of one protein and check
# that it returns one number per sample (column of X). Values that are not finite
# (-Inf, e.g. log2(0) for a sample without data, Inf, NaN) become NA.
.user_method <- function(method, X) {
    x <- method(X)
    if (length(x) != ncol(X)) {
        stop("The quantification method must return one value per sample (", ncol(X),
             "), but returned ", length(x), ".")
    }
    if (!is.numeric(x) && !all(is.na(x))) {
        stop("The quantification method must return numbers, but returned ", class(x)[1], ".")
    }
    x <- as.numeric(x)
    x[!is.finite(x)] <- NA
    x
}

# Evaluate the quantification of one protein; an error names the protein
.for_protein <- function(expr, protein) {
    tryCatch(expr, error = function(e) stop("Protein ", protein, ": ", conditionMessage(e), call. = FALSE))
}

# Check the method parameters: each a single finite number; p1 / p2 may also be NULL
# (the default of the method). An empty value (e.g. k = NULL) would be read past its end
# in C++.
.check_parameters <- function(p1, p2, k, min_M, n_threads) {
    values <- list(p1 = p1, p2 = p2, k = k, min_M = min_M, n_threads = n_threads)
    for (name in names(values)) {
        v <- values[[name]]
        if (is.null(v) && name %in% c("p1", "p2")) {
            next
        }
        if (!is.numeric(v) || length(v) != 1 || !is.finite(v)) {
            stop(name, " must be a single number", if (name %in% c("p1", "p2")) " or NULL", ".")
        }
    }
}

# One warning for the protein accessions without a sequence in rfasta (coverage 0%)
.warn_missing_fasta <- function(all, missing) {
    missing <- unique(missing)
    if (length(missing) > 0) {
        warning("No sequence in rfasta for ", length(missing), " of ", length(unique(all)),
                " protein accessions, e.g. ", paste(head(missing, 3), collapse = ", "),
                "; their coverage is reported as 0%.", call. = FALSE)
    }
}

# peptide_extractor (a function) and rfasta (named sequences) for the fragment statistics
.check_fragment_options <- function(peptide_extractor, rfasta) {
    if (!is.null(peptide_extractor) && !is.function(peptide_extractor)) {
        stop("peptide_extractor must be a function, e.g. function(x) gsub(\"\\\\(UniMod:\\\\d+\\\\)|\\\\d\", \"\", x).")
    }
    if (!is.null(rfasta)) {
        if (!(is.list(rfasta) || is.character(rfasta)) || is.null(names(rfasta))) {
            stop("rfasta must be a named list or a named character vector of protein sequences.")
        }
        if (is.null(peptide_extractor)) {
            warning("rfasta is used only with peptide_extractor (sequence coverage); it is ignored.", call. = FALSE)
        }
    }
}

# The unique peptides of some ions: peptide_extractor must return one text value per ion
.peptides <- function(peptide_extractor, ions) {
    peptides <- peptide_extractor(ions)
    if (!(is.character(peptides) || is.factor(peptides)) || length(peptides) != length(ions)) {
        stop("peptide_extractor must return one character value per ion: it was given ", length(ions),
             " ion(s) and returned ", length(peptides), " value(s) of type ", class(peptides)[1], ".", call. = FALSE)
    }
    unique(as.character(peptides))
}

# Fragment statistics of proteins (process, long_format_to_iq_format): for each protein, its
# identifier (a protein group is split at ";" for the coverage), its ions (indices into
# ion_prefix, in input order) and its XIC (summed intensities of all values read, before
# preprocessing). Returns a character matrix with the columns n_fragments, n_peptides,
# peptides and XIC, and coverage if rfasta is given.
.fragment_statistics <- function(protein_ids, ions, xic, ion_prefix, peptide_extractor, rfasta,
                                 verbose = TRUE) {

    if (verbose) {
        message("Calculating fragment statistics... ")
    }

    s <- matrix("", nrow = length(protein_ids),
                ncol = if (is.null(rfasta)) 4 else 5)

    colnames(s) <- if (is.null(rfasta)) c("n_fragments", "n_peptides", "peptides", "XIC") else c("n_fragments", "n_peptides", "peptides", "XIC", "coverage")

    # accessions looked up in rfasta, and those without a sequence there
    fasta_all <- character(0)
    fasta_missing <- character(0)

    n_display <- nrow(s) %/% 20

    for (r in seq_along(protein_ids)) {

        if (verbose && n_display > 0 && r %% n_display == 0) {
            message(format(r * 100 / nrow(s), digits = 2), "%")
        }

        s[r, 4] <- xic[r]

        ind <- ions[[r]]
        s[r, 1] <- length(ind)
        # the extractor gets the secondary_id part of each ion, without the suffix of the
        # intensity value (1, 2, ... or its intensity_col_id), e.g. a precursor id
        frags <- .peptides(peptide_extractor, ion_prefix[ind])
        s[r, 2] <- length(frags)
        s[r, 3] <- paste(frags, collapse = ";")

        if (!is.null(rfasta)) {
            aa <- unlist(strsplit(protein_ids[r], ";"))
            fasta_all <- c(fasta_all, aa)
            if (length(aa) >0) {
                coverages <- rep(0, length(aa))
                for (i2 in seq_along(aa)) {
                    if (aa[i2] %in% names(rfasta)) {
                        # one string, also for a sequence of single letters (seqinr::read.fasta)
                        seq <- paste(rfasta[[aa[i2]]], collapse = "")
                        tmp <- rep(0, nchar(seq))
                        for (i3 in seq_along(frags)) {  # no iteration if there are no peptides
                            # an NA or empty peptide is not a sequence ("" would match everywhere)
                            if (is.na(frags[i3]) || !nzchar(frags[i3])) {
                                next
                            }
                            bb <- gregexpr(frags[i3], seq, fixed = TRUE)[[1]]
                            # gregexpr returns -1 when the peptide is not in this sequence
                            if (is.na(bb[1]) || bb[1] < 0) {
                                next
                            }
                            for (i4 in seq_along(bb)) {
                                tmp[bb[i4] : (bb[i4] + nchar(frags[i3])-1)] <- 1
                            }
                        }
                        coverages[i2] <- round(sum(tmp) * 100 / length(tmp))
                    }
                    else {
                        fasta_missing <- c(fasta_missing, aa[i2])  # coverage stays 0%
                    }
                }
                s[r, 5] <- paste(paste0(coverages, "%"), collapse = ";")
            }
        }
    }
    .warn_missing_fasta(fasta_all, fasta_missing)
    if (verbose) {
        message("Completed. ")
    }

    s
}

# input_data: a data frame, or the name of one text or parquet file
.check_input_data <- function(input_data) {
    if (!is.data.frame(input_data) &&
        !(is.character(input_data) && length(input_data) == 1 && !is.na(input_data))) {
        stop("input_data must be a data frame or the name of one text or parquet file",
             if (inherits(input_data, "ArrowTabular")) "; for an arrow Table or RecordBatch use as.data.frame(input_data)",
             ".")
    }
}

# Read long-format data with the reader for its type: a data frame (fast_read_data_frame), a
# parquet file, the name ending in .parquet (fast_read_parquet), or a tab-separated text file
# (fast_read). The filtering is the same for all. Stops if nothing could be read.
#
# writer: an iq format writer (long_format_to_iq_format) that receives the values chunk by
# chunk; the result then has no values, only the lists of proteins, samples and ions.
.read_long_format <- function(input_data, primary_id, sample_id, secondary_id,
                              intensity_col, intensity_col_sep, intensity_col_id, annotation_col,
                              filter_string_equal, filter_string_not_equal,
                              filter_double_less, filter_double_greater, na_string,
                              writer = NULL, verbose = TRUE) {

    if (!is.null(writer)) {
        if (is.data.frame(input_data)) {
            stop("Only a text or parquet file is read incrementally.")
        }
        if (grepl("\\.parquet$", input_data, ignore.case = TRUE)) {
            options <- .stream_options(sample_id, primary_id, secondary_id, intensity_col, annotation_col,
                                       filter_string_equal, filter_string_not_equal,
                                       filter_double_less, filter_double_greater,
                                       intensity_col_sep, intensity_col_id, na_string)
            return(.write_parquet_values(input_data, options, writer, verbose = verbose))
        }
        iq_dat <- .Call("iq_filter_write",
                        .fast_read_command(input_data, sample_id, primary_id, secondary_id, intensity_col,
                                           annotation_col, filter_string_equal, filter_string_not_equal,
                                           filter_double_less, filter_double_greater,
                                           intensity_col_sep, intensity_col_id, na_string, verbose),
                        writer)
        return(iq_dat)
    }

    if (is.data.frame(input_data)) {
        iq_dat <- fast_read_data_frame(input_data,
                                       primary_id = primary_id,
                                       sample_id  = sample_id,
                                       secondary_id = secondary_id,
                                       intensity_col = intensity_col,
                                       intensity_col_sep = intensity_col_sep,
                                       intensity_col_id = intensity_col_id,
                                       annotation_col = annotation_col,
                                       filter_string_equal = filter_string_equal,
                                       filter_string_not_equal = filter_string_not_equal,
                                       filter_double_less = filter_double_less,
                                       filter_double_greater = filter_double_greater,
                                       na_string = na_string,
                                       verbose = verbose)
    } else if (grepl("\\.parquet$", input_data, ignore.case = TRUE)) {
        iq_dat <- fast_read_parquet(input_data,
                                    primary_id = primary_id,
                                    sample_id  = sample_id,
                                    secondary_id = secondary_id,
                                    intensity_col = intensity_col,
                                    intensity_col_sep = intensity_col_sep,
                                    intensity_col_id = intensity_col_id,
                                    annotation_col = annotation_col,
                                    filter_string_equal = filter_string_equal,
                                    filter_string_not_equal = filter_string_not_equal,
                                    filter_double_less = filter_double_less,
                                    filter_double_greater = filter_double_greater,
                                    na_string = na_string,
                                    verbose = verbose)
    } else {
        iq_dat <- fast_read(input_data,
                            primary_id = primary_id,
                            sample_id  = sample_id,
                            secondary_id = secondary_id,
                            intensity_col = intensity_col,
                            intensity_col_sep = intensity_col_sep,
                            intensity_col_id = intensity_col_id,
                            annotation_col = annotation_col,
                            filter_string_equal = filter_string_equal,
                            filter_string_not_equal = filter_string_not_equal,
                            filter_double_less = filter_double_less,
                            filter_double_greater = filter_double_greater,
                            na_string = na_string,
                            verbose = verbose)
    }

    iq_dat
}

# Leave out the values without a protein or sample identifier (NA or empty in the input; the
# readers turn both into ""): they would be merged into one protein or one sample named "".
.drop_unidentified <- function(iq_dat, primary_id, sample_id) {
    no_protein <- which(iq_dat$protein[, 1] == "")
    no_sample <- which(as.character(iq_dat$sample) == "")
    if (length(no_protein) + length(no_sample) == 0) {
        return(iq_dat)
    }

    qt <- iq_dat$quant_table
    in_no_protein <- qt$protein_list %in% no_protein
    in_no_sample <- qt$sample_list %in% no_sample
    warning("Left out ",
            paste(c(if (any(in_no_protein)) paste0(sum(in_no_protein), " value(s) without a protein identifier (", primary_id, ")"),
                    if (any(in_no_sample)) paste0(sum(in_no_sample), " value(s) without a sample identifier (", sample_id, ")")),
                  collapse = " and "),
            ".", call. = FALSE)

    keep <- which(!(in_no_protein | in_no_sample))
    for (nm in names(qt)) {
        qt[[nm]] <- qt[[nm]][keep]
    }
    if (length(qt$quant) == 0) {
        stop("No values left with a protein and a sample identifier.")
    }

    # the nameless sample is removed from the list of samples (no "" column in the output);
    # the nameless protein simply has no values left
    if (length(no_sample) > 0) {
        kept_samples <- setdiff(seq_along(iq_dat$sample), no_sample)
        qt$sample_list <- match(qt$sample_list, kept_samples)
        iq_dat$sample <- iq_dat$sample[kept_samples]
    }

    iq_dat$quant_table <- qt
    iq_dat
}

# Names for a message: at most max of them, then "... and n more" 
.name_list <- function(x, max = 10) {
    x <- as.character(x)
    if (length(x) <= max) {
        return(paste(x, collapse = ", "))
    }
    paste0(paste(x[seq_len(max)], collapse = ", "), " ... and ", length(x) - max, " more")
}

# Do two paths name the same file? An existing file is compared by its resolved path (also
# through a symbolic link), another by its resolved directory and name. Without case on Windows
# and macOS, whose file systems usually ignore case.
.same_file <- function(a, b) {
    resolve <- function(path) {
        if (file.exists(path)) {
            normalizePath(path)
        } else {
            file.path(normalizePath(dirname(path), mustWork = FALSE), basename(path))
        }
    }
    a <- resolve(a)
    b <- resolve(b)
    if (.Platform$OS.type == "windows" || Sys.info()[["sysname"]] == "Darwin") {
        a <- tolower(a)
        b <- tolower(b)
    }
    a == b
}

# The identifier columns, the same rules for all readers (text, parquet, data frame):
# sample_id and primary_id one column name each, secondary_id one or more
.check_id_columns <- function(sample_id, primary_id, secondary_id) {
    one_name <- function(x) is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x)
    if (!one_name(sample_id)) {
        stop("sample_id must be one column name.")
    }
    if (!one_name(primary_id)) {
        stop("primary_id must be one column name.")
    }
    if (!is.character(secondary_id) || length(secondary_id) < 1 || anyNA(secondary_id) || !all(nzchar(secondary_id))) {
        stop("secondary_id must be one or more column names.")
    }
}

# The filters, the same rules for all readers (text, parquet, data frame): named vectors (an
# unnamed filter would be ignored), and numbers for the double filters (the text reader would
# read "0,01" or "abc" as 0). Returns them as named character / numeric vectors (NULL if not given).
.check_filters <- function(filter_string_equal, filter_string_not_equal,
                           filter_double_less, filter_double_greater) {

    named <- function(x, what, numeric = FALSE) {
        if (is.null(x)) {
            return(NULL)
        }
        if (is.null(names(x)) || any(is.na(names(x)) | names(x) == "")) {
            stop(what, " must be a named vector, e.g. c(\"Q.Value\" = \"0.01\")")
        }
        if (numeric) {
            v <- suppressWarnings(as.numeric(x))
            if (any(is.na(v))) {
                stop(what, " must have numeric values")
            }
            names(v) <- names(x)
            return(v)
        }
        v <- as.character(x)
        names(v) <- names(x)
        v
    }

    list(filter_string_equal = named(filter_string_equal, "filter_string_equal"),
         filter_string_not_equal = named(filter_string_not_equal, "filter_string_not_equal"),
         filter_double_less = named(filter_double_less, "filter_double_less", TRUE),
         filter_double_greater = named(filter_double_greater, "filter_double_greater", TRUE))
}

# log2_intensity_cutoff: one number, or NULL for no cut-off (a text value would be compared with
# the log2 intensities as text, e.g. 6 > "12")
.check_cutoff <- function(log2_intensity_cutoff) {
    if (!is.null(log2_intensity_cutoff) &&
        (!is.numeric(log2_intensity_cutoff) || length(log2_intensity_cutoff) != 1 || is.na(log2_intensity_cutoff))) {
        stop("log2_intensity_cutoff must be a number, or NULL for no cut-off.")
    }
}

# The intensity columns, the same rules for all readers (text, parquet, data frame)
.check_intensity_columns <- function(intensity_col, intensity_col_sep, intensity_col_id) {
    if (!is.character(intensity_col) || length(intensity_col) < 1 || anyNA(intensity_col)) {
        stop("intensity_col must be one or more column names")
    }
    if (length(intensity_col) > 1 && !is.null(intensity_col_sep)) {
        stop("intensity_col_sep can only be used with a single intensity column")
    }
    if (!is.null(intensity_col_id) && length(intensity_col_id) != length(intensity_col)) {
        stop("intensity_col_id must have the same length as intensity_col")
    }
    if (length(intensity_col) == 1 && !is.null(intensity_col_id) && is.null(intensity_col_sep)) {
        stop("intensity_col_id with a single intensity column needs intensity_col_sep")
    }
}

# A file to be written: one non-empty name, in a directory that exists, not a directory itself
.check_output_file <- function(path, argument) {
    if (!is.character(path) || length(path) != 1 || is.na(path) || !nzchar(path)) {
        stop(argument, " must be one file name.")
    }
    if (dir.exists(path)) {
        stop(argument, ": ", path, " is a directory, not a file name.")
    }
    if (!dir.exists(dirname(path))) {
        stop(argument, ": the directory ", dirname(path), " does not exist.")
    }
}

# normalization = "median" or "none" -> TRUE / FALSE
.median_normalization <- function(normalization) {
    if (!is.character(normalization) || length(normalization) != 1 || is.na(normalization) ||
        !(normalization %in% c("median", "none"))) {
        stop("normalization must be \"median\" or \"none\".")
    }
    normalization == "median"
}

# Check the quantification and rescaling methods before reading any data
.check_methods <- function(method, rescale_method) {
    if (!is.function(method) && !(length(method) == 1 && method %in% names(.quant_methods))) {
        stop("Unknown quantitation method. Use a function or one of: ",
             paste(names(.quant_methods), collapse = ", "))
    }
    if (!is.null(rescale_method) && !(length(rescale_method) == 1 && rescale_method %in% .rescale_methods)) {
        stop("Unknown rescaling method. Use one of: ", paste(.rescale_methods, collapse = ", "))
    }
}

# Options for the shared reader (src/stream.cpp), checked in R
.stream_options <- function(sample_id, primary_id, secondary_id, intensity_col, annotation_col,
                            filter_string_equal, filter_string_not_equal,
                            filter_double_less, filter_double_greater,
                            intensity_col_sep, intensity_col_id, na_string) {

    .check_id_columns(sample_id, primary_id, secondary_id)
    filters <- .check_filters(filter_string_equal, filter_string_not_equal,
                              filter_double_less, filter_double_greater)

    chr <- function(x) if (is.null(x)) NULL else as.character(x)

    list(sample_id = chr(sample_id),
         primary_id = chr(primary_id),
         secondary_id = chr(secondary_id),
         intensity_col = chr(intensity_col),
         annotation_col = chr(annotation_col),
         filter_string_equal = filters$filter_string_equal,
         filter_string_not_equal = filters$filter_string_not_equal,
         filter_double_less = filters$filter_double_less,
         filter_double_greater = filters$filter_double_greater,
         intensity_col_sep = chr(intensity_col_sep),
         intensity_col_id = chr(intensity_col_id),
         na_string = chr(na_string))
}

# All columns the shared reader needs
.stream_columns <- function(options) {
    unique(c(options$sample_id, options$primary_id, options$secondary_id,
             options$intensity_col, options$intensity_col_id, options$annotation_col,
             names(options$filter_string_equal), names(options$filter_string_not_equal),
             names(options$filter_double_less), names(options$filter_double_greater)))
}

.stream_end <- function(state, verbose = TRUE) {
    ret <- .Call("iq_stream_end", state, verbose)
    if (is.null(ret)) {
        stop("No quantitative values left after filtering.")
    }
    ret
}

# Read a parquet file one row group at a time, and pass the values of each slice of rows to an
# iq format writer (long_format_to_iq_format). Arrow decodes a whole row group (only the needed
# columns), so the memory used is that of the largest row group, not of the file; the slices
# keep the R copy of it small. Returns the lists of proteins, samples and ions.
.write_parquet_values <- function(filename, options, writer, slice_rows = 100000, verbose = TRUE) {

    if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("Reading parquet files requires the 'arrow' package: install.packages(\"arrow\")")
    }

    if (!file.exists(filename)) {
        stop("Cannot open: ", filename)
    }

    columns <- .stream_columns(options)

    reader <- arrow::ParquetFileReader$create(filename, mmap = FALSE)  # read, not mapped: no growing mapped pages
    column_indices <- match(columns, names(reader$GetSchema()))
    if (any(is.na(column_indices))) {
        stop("Cannot find column(s) in ", filename, ": ",
             paste(columns[is.na(column_indices)], collapse = ", "))
    }

    if (verbose) {
        message("Reading: ", filename)
    }

    state <- .Call("iq_stream_begin", options)

    for (g in seq_len(reader$num_row_groups) - 1L) {
        tab <- reader$ReadRowGroup(g, column_indices = column_indices - 1L)
        n <- tab$num_rows
        start <- 0
        while (start < n) {
            size <- min(slice_rows, n - start)
            .Call("iq_stream_add", state, as.data.frame(tab$Slice(start, size)))
            .Call("iq_stream_write", state, writer)
            start <- start + size
        }
        # free the Arrow memory of the row group now: R does not see its size and would
        # otherwise collect it late, so that the memory would grow with the file
        rm(tab)
        invisible(gc(verbose = FALSE))
    }

    .stream_end(state, verbose)
}

# TRUE if the package was built with the original MaxLFQ implementation (Eigen); otherwise
# "maxlfq", "weighted_maxlfq" and fast_MaxLFQ() use the solver of "maxlfq_bit". See src/Makevars.
.maxlfq_eigen <- function() {
    .Call("iq_eigen")
}

# default p1 / p2 of a method (NULL if the method does not use it)
.method_parameters <- function(method, p1, p2) {

    if (method == "rlm" || method == "weighted_rlm") {
        if (is.null(p1)) {
            p1 <- 1.345
        }
        if (is.null(p2)) {
            p2 <- 1000      # max iterations
        }
    }

    if (method == "median_polish" || method == "weighted_median_polish") {
        if (is.null(p1)) {
            p1 <- 1e-4      # epsilon
        }
        if (is.null(p2)) {
            p2 <- 1000      # max iterations
        }
    }

    if (method == "maxlfq_bit" || method == "maxlfq" || method == "weighted_maxlfq") {
        if (is.null(p2)) {
            p2 <- 0         # memory level of the maxlfq_bit solver: middle
        }
    }

    if (method == "topN") {
        if (is.null(p2)) {
            p2 <- 3         # the number of ions
        }
        if (p2 < 1) {
            stop("p2 must be 1 or more for the method \"topN\" (the number of ions).")
        }
    }

    list(p1 = p1, p2 = p2)
}

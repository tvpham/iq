#########################################################################
#
# Author: Thang V. Pham {t.pham@amsterdamumc.nl}
#
# All rights reserved.
#
# Citations:
#
# Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, 
# Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein 
# quantification algorithms in mass spectrometry-based proteomics.
# Journal of Proteome Research 2026 Feb 6; 25(2):1198-1203.
#
# Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative
# protein abundances from ion quantification in DIA-MS-based proteomics,
# Bioinformatics 2020 Apr 15;36(8):2611-2613.
#
# Software version: 2.2.0
#
#########################################################################

fast_read <- function(filename,
                      sample_id = "R.Condition",
                      primary_id = "PG.ProteinGroups",
                      secondary_id = c("EG.ModifiedSequence", "FG.Charge", "F.FrgIon", "F.Charge"),
                      intensity_col = "F.PeakArea",
                      annotation_col = c("PG.Genes", "PG.ProteinNames"),
                      filter_string_equal = c("F.ExcludedFromQuantification" = "False"),
                      filter_string_not_equal = NULL,
                      filter_double_less = c("PG.Qvalue" = "0.01", "EG.Qvalue" = "0.01"),
                      filter_double_greater = NULL,
                      intensity_col_sep = NULL,
                      intensity_col_id = NULL,
                      na_string = "0",
                      verbose = TRUE) {

    cmd <- .fast_read_command(filename, sample_id, primary_id, secondary_id, intensity_col, annotation_col,
                              filter_string_equal, filter_string_not_equal,
                              filter_double_less, filter_double_greater,
                              intensity_col_sep, intensity_col_id, na_string, verbose)

    return(.Call("iq_filter", cmd))
}

fast_MaxLFQ <- function(norm_data, row_names = NULL, col_names = NULL, verbose = TRUE) {

    n <- length(norm_data$quant)
    if (n == 0) {
        stop("No quantitative values.")
    }
    if (length(norm_data$protein_list) != n || length(norm_data$sample_list) != n || length(norm_data$id) != n) {
        stop("$protein_list, $sample_list, $id and $quant must have the same length.")
    }

    # check for NA values
    if (any(is.na(norm_data$protein_list))) {
        stop("NA value in $protein_list.\n")
    }

    if (any(is.na(norm_data$sample_list))) {
        stop("NA value in $sample_list.\n")
    }

    if (any(is.na(norm_data$id))) {
        stop("NA value in $id.\n")
    }

    if (any(is.na(norm_data$quant))) {
        stop("NA value in $quant.\n")
    }

    if (is.null(row_names)) {
        proteins <- unique(as.character(norm_data$protein_list))
        p_list <- as.integer(factor(as.character(norm_data$protein_list), levels = proteins))
    } else {
        proteins <- as.character(row_names)
        p_list <- as.integer(norm_data$protein_list)
        if (any(p_list < 1 | p_list > length(proteins))) {
            stop("With row_names, $protein_list must be indices from 1 to length(row_names).")
        }
    }

    if (is.null(col_names)) {
        s <- unique(as.character(norm_data$sample_list))
        s_list <- as.integer(factor(as.character(norm_data$sample_list), levels = s))
    } else {
        s <- as.character(col_names)
        s_list <- as.integer(norm_data$sample_list)
        if (any(s_list < 1 | s_list > length(s))) {
            stop("With col_names, $sample_list must be indices from 1 to length(col_names).")
        }
    }

    ret <- .Call("iq_MaxLFQ",
                 list(protein_index = p_list,
                      ion_index = as.integer(factor(norm_data$id)),
                      sample_index = s_list,
                      quant = as.double(norm_data$quant),
                      verbose = verbose))

    if (is.null(ret)) {  # interrupted by the user (iq_MaxLFQ prints "Cancelled.")
        stop("Cancelled.")
    }

    rownames(ret$estimate) <- proteins[as.integer(rownames(ret$estimate))]

    colnames(ret$estimate) <- s[as.integer(colnames(ret$estimate))]

    return(ret)
}

fast_preprocess <- function(quant_table,
                            median_normalization = TRUE,
                            log2_intensity_cutoff = 0,
                            pdf_out = "qc-plots-fast.pdf",
                            pdf_width = 12,
                            pdf_height = 8,
                            show_boxplot = TRUE,
                            boxplot_max_points = 10000,
                            verbose = TRUE) {

    .check_cutoff(log2_intensity_cutoff)

    if (!is.null(boxplot_max_points) &&
        (!is.numeric(boxplot_max_points) || length(boxplot_max_points) != 1 ||
         is.na(boxplot_max_points) || boxplot_max_points < 1)) {
        stop("boxplot_max_points must be a positive number, or NULL for all values.")
    }

    if (!is.null(pdf_out)) {
        pdf(pdf_out, pdf_width, pdf_height)
        pdf_device <- dev.cur()
        on.exit(dev.off(pdf_device), add = TRUE)  # also on an error
    }

    d <- quant_table

    d$quant <- log2(d$quant)

    samples_read <- unique(d$sample_list)

    # always remove values that are not finite after log2 (-Inf for a zero intensity, NaN
    # for a negative one): one -Inf would make the median of its sample -Inf and every
    # normalized value NA or -Inf. Then the intensity cut-off, if given.
    keep <- is.finite(d$quant)
    if (!is.null(log2_intensity_cutoff)) {
        if (verbose) {
            message("Removing low intensities...")
        }
        keep <- keep & d$quant > log2_intensity_cutoff
    }

    if (!any(keep)) {
        stop("No quantitative values left: all intensities are ",
             if (is.null(log2_intensity_cutoff)) "zero or negative."
             else "zero, negative or not above the intensity cut-off.")
    }

    if (!is.null(pdf_out)) {
        a <- hist(d$quant[is.finite(d$quant)], 100, col = "steelblue", border = "steelblue", freq = FALSE,
                  main = "Histogram of log2 intensities",
                  xlab = "log2 intensity")
        if (!is.null(log2_intensity_cutoff)) {
            arrows(log2_intensity_cutoff, 0, log2_intensity_cutoff, max(a$density) / 2.0 , col = "red", code = 1, lwd = 2)
        }
    }

    selected <- which(keep)
    d$quant <- d$quant[selected]
    d$protein_list <- d$protein_list[selected]
    d$sample_list <- d$sample_list[selected]
    d$id <- d$id[selected]

    # samples with values after the cut-off (a sample without values would have median NA,
    # which would make every normalized value NA)
    samples <- unique(d$sample_list)

    if (length(samples) < length(samples_read)) {
        # a classed warning that carries the samples, so that process() can name them
        empty <- setdiff(samples_read, samples)
        warning(structure(class = c("iq_empty_samples", "warning", "condition"),
                          list(message = paste0("No values left after removing zero, negative and low intensities in sample(s) ",
                                                .name_list(empty),
                                                " (as in sample_list); these samples are left out of the normalization."),
                               call = sys.call(), samples = empty)))
    }

    # the values of each sample, found once (not a scan of all values for every sample)
    g <- match(d$sample_list, samples)
    rows <- split(seq_along(d$quant), factor(g, levels = seq_along(samples)))

    m <- vapply(rows, function(i) median(d$quant[i], na.rm = TRUE), numeric(1), USE.NAMES = FALSE)

    # boxplots of at most boxplot_max_points evenly spaced values per sample
    boxplot_samples <- function() {
        boxplot(lapply(rows, function(i) .thin_for_plot(d$quant[i], boxplot_max_points)),
                names = as.character(samples),
                main = "Boxplot of fragment intensities per sample",
                ylab = "log2 intensity",
                outline = FALSE,
                col = "steelblue",
                whisklty = 1,
                staplelty = 0,
                las = 2)
    }

    if (!is.null(pdf_out) && show_boxplot) {
        if (verbose) {
            message("Barplotting raw data ...")
        }
        boxplot_samples()
    }

    if (median_normalization) {

        if (verbose) {
            message("Median normalization ...")
        }

        f <- mean(m) - m
        d$quant <- d$quant + f[g]

        if (!is.null(pdf_out) && show_boxplot) {
            if (verbose) {
                message("Barplotting after normalization ...")
            }
            boxplot_samples()
        }
    }

    return(d)
}

process_long_format <- function(input_data,
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
                                verbose = TRUE) {


    median_normalization <- .median_normalization(normalization)
    .check_filters(filter_string_equal, filter_string_not_equal, filter_double_less, filter_double_greater)
    .check_cutoff(log2_intensity_cutoff)


    if (is.data.frame(input_data)) {

        required_columns <- c(sample_id, primary_id, secondary_id, intensity_col, annotation_col,
                              names(filter_string_equal),
                              names(filter_string_not_equal),
                              names(filter_double_less),
                              names(filter_double_greater))

        index <- match(required_columns, names(input_data))
        # check for NA values
        if (any(is.na(index))) {
            stop(paste0("The input data frame has no column: ", required_columns[which(is.na(index))[1]]))
        }

        d <- as.data.frame(input_data)

        if (!is.null(filter_string_equal)) {
            cmd <- paste0("d <- d[",
                          paste0("d$", names(filter_string_equal), " == ", '"', filter_string_equal, '"', collapse = " & "), ",]")
            eval(parse(text = cmd))
        }

        if (!is.null(filter_string_not_equal)) {
            cmd <- paste0("d <- d[",
                          paste0("d$", names(filter_string_not_equal), " != ", '"', filter_string_not_equal, '"', collapse = " & "), ",]")
            eval(parse(text = cmd))
        }

        # as fast_read: a row is kept if its value is less (greater) than the given value, or
        # missing (NA, NaN, empty or text that is not a number)
        for (f in names(filter_double_less)) {
            v <- suppressWarnings(as.numeric(as.character(d[[f]])))
            d <- d[is.na(v) | v < as.numeric(filter_double_less[[f]]), , drop = FALSE]
        }

        for (f in names(filter_double_greater)) {
            v <- suppressWarnings(as.numeric(as.character(d[[f]])))
            d <- d[is.na(v) | v > as.numeric(filter_double_greater[[f]]), , drop = FALSE]
        }

        norm_data <- preprocess(d,
                                primary_id = primary_id,
                                sample_id  = sample_id,
                                secondary_id = secondary_id,
                                intensity_col = intensity_col,
                                intensity_col_sep = NULL,
                                intensity_col_id = NULL,
                                median_normalization = median_normalization,
                                log2_intensity_cutoff = log2_intensity_cutoff,
                                pdf_out = pdf_out,
                                pdf_width = pdf_width,
                                pdf_height = pdf_height,
                                na_string = na_string,
                                show_boxplot = show_boxplot,
                                verbose = verbose)

        res <- fast_MaxLFQ(norm_data, verbose = verbose)

        if (!is.null(peptide_extractor)) {
            if (verbose) {
                message("Calculating fragment statistics... ")
            }
            s <- matrix(0, nrow = nrow(res$estimate), ncol = 2)
            colnames(s) <- c("n_fragments", "n_peptides")

            n_display <- nrow(s) %/% 20

            for (i in 1:nrow(s)) {

                if (verbose && n_display > 0 && i %% n_display == 0) {
                    message(format(i * 100 / nrow(s), digits = 2), "%")
                }

                #ind <- unique(iq_dat$quant_table$id[iq_dat$quant_table$protein_list == i])
                frags <- unique(norm_data$id[norm_data$protein_list == rownames(res$estimate)[i]])
                s[i, 1] <- length(frags)
                frags <- .peptides(peptide_extractor, frags)
                s[i, 2] <- length(frags)
            }
            if (verbose) {
                message("Completed. ")
            }
        }

        if (verbose) {
            message("Writing to: ", output_filename)
        }

        if (is.null(annotation_col)) {
            if (!is.null(peptide_extractor)) {
                tab <- cbind(rownames(res$estimate), s, res$estimate)
                colnames(tab)[1] <- primary_id
            }
            else {
                tab <- cbind(rownames(res$estimate), res$estimate)
                colnames(tab)[1] <- primary_id
            }
        } else {
            extra_annotation <- extract_annotation(rownames(res$estimate),
                                                   d,
                                                   primary_id = primary_id,
                                                   annotation_columns = annotation_col)

            if (!is.null(peptide_extractor)) {
                tab <- cbind(rownames(res$estimate),
                             extra_annotation[, annotation_col], s,
                             res$estimate)
                colnames(tab)[1:(length(annotation_col)+1)] <- c(primary_id, annotation_col)
            }
            else {
                tab <- cbind(rownames(res$estimate),
                             extra_annotation[, annotation_col],
                             res$estimate)
                colnames(tab)[1:(length(annotation_col)+1)] <- c(primary_id, annotation_col)
            }
        }

        write.table(tab, output_filename, sep = "\t", row.names = FALSE, quote = FALSE)


        return(invisible(NULL))

    }


    iq_dat <- fast_read(input_data,
                        primary_id = primary_id,
                        sample_id  = sample_id,
                        secondary_id = secondary_id,
                        intensity_col = intensity_col,
                        intensity_col_sep = intensity_col_sep,
                        annotation_col = annotation_col,
                        filter_string_equal = filter_string_equal,
                        filter_string_not_equal = filter_string_not_equal,
                        filter_double_less = filter_double_less,
                        filter_double_greater = filter_double_greater,
                        na_string = na_string,
                        verbose = verbose)

    iq_norm_data <- fast_preprocess(iq_dat$quant_table,
                                    median_normalization = median_normalization,
                                    log2_intensity_cutoff = log2_intensity_cutoff,
                                    pdf_out = pdf_out,
                                    pdf_width = pdf_width,
                                    pdf_height = pdf_height,
                                    show_boxplot = show_boxplot,
                                    verbose = verbose)

    res <- fast_MaxLFQ(iq_norm_data,
                       row_names = iq_dat$protein[, 1],
                       col_names = iq_dat$sample,
                       verbose = verbose)


    if (!is.null(peptide_extractor)) {
        if (verbose) {
            message("Calculating fragment statistics... ")
        }
        
        s <- matrix("", nrow = nrow(iq_dat$protein), 
                    ncol = if (is.null(rfasta)) 4 else 5)
        
        colnames(s) <- if (is.null(rfasta)) c("n_fragments", "n_peptides", "peptides", "XIC") else c("n_fragments", "n_peptides", "peptides", "XIC", "coverage")

        # accessions looked up in rfasta, and those without a sequence there
        fasta_all <- character(0)
        fasta_missing <- character(0)

        n_display <- nrow(s) %/% 20

        for (i in 1:nrow(s)) {

            if (verbose && n_display > 0 && i %% n_display == 0) {
                message(format(i * 100 / nrow(s), digits = 2), "%")
            }
            
            index <- iq_dat$quant_table$protein_list == i
            XIC <- sum(iq_dat$quant_table$quant[index])
            s[i, 4] <- XIC

            ind <- unique(iq_dat$quant_table$id[index])
            frags <- iq_dat$ion[ind]
            s[i, 1] <- length(frags)
            frags <- .peptides(peptide_extractor, frags)
            s[i, 2] <- length(frags)
            s[i, 3] <- paste(frags, collapse = ";")

            if (!is.null(rfasta)) {
                pname <- iq_dat$protein[i]
                aa <- unlist(strsplit(pname, ";"))
                fasta_all <- c(fasta_all, aa)
                if (length(aa) >0) {
                    coverages <- rep(0, length(aa))
                    for (i2 in 1:length(aa)) {
                        seq <- rfasta[[aa[i2]]]
                        if (!is.null(seq)) {
                            tmp <- rep(0, nchar(seq))
                            for (i3 in 1:length(frags)) {
                                bb <- gregexpr(frags[i3], seq, fixed = TRUE)[[1]]
                                # gregexpr returns -1 when the peptide is not in this sequence
                                if (is.na(bb[1]) || bb[1] < 0) {
                                    next
                                }
                                for (i4 in 1:length(bb)) {
                                    tmp[bb[i4] : (bb[i4] + nchar(frags[i3])-1)] <- 1
                                }
                            }
                            coverages[i2] <- round(sum(tmp) * 100 / length(tmp))
                        }
                        else {
                            fasta_missing <- c(fasta_missing, aa[i2])  # coverage stays 0%
                        }
                    }
                    s[i, 5] <- paste(paste0(coverages, "%"), collapse = ";")
                }
            }
        }
        .warn_missing_fasta(fasta_all, fasta_missing)
        if (verbose) {
            message("Completed. ")
        }
    }

    if (verbose) {
        message("Writing to: ", output_filename)
    }

    if (is.null(annotation_col)) {
        if (!is.null(peptide_extractor)) {
            tab <- cbind(rownames(res$estimate), s, res$estimate)
            colnames(tab)[1] <- primary_id
        }
        else {
            tab <- cbind(rownames(res$estimate), res$estimate)
            colnames(tab)[1] <- primary_id
        }
    } else {
        extra_annotation <- extract_annotation(rownames(res$estimate),
                                               iq_dat$protein,
                                               primary_id = primary_id,
                                               annotation_columns = annotation_col)

        if (!is.null(peptide_extractor)) {
            tab <- cbind(rownames(res$estimate),
                         extra_annotation[, annotation_col], s,
                         res$estimate)
            colnames(tab)[1:(length(annotation_col)+1)] <- c(primary_id, annotation_col)
        }
        else {
            tab <- cbind(rownames(res$estimate),
                         extra_annotation[, annotation_col],
                         res$estimate)
            colnames(tab)[1:(length(annotation_col)+1)] <- c(primary_id, annotation_col)
        }
    }

    write.table(tab, output_filename, sep = "\t", row.names = FALSE, quote = FALSE)

    invisible(NULL)
}

#------------------------------ helpers ------------------------------------

# The command for the text-file reader (iq_filter, iq_filter_write in src/iq.cpp)
.fast_read_command <- function(filename, sample_id, primary_id, secondary_id, intensity_col, annotation_col,
                               filter_string_equal, filter_string_not_equal,
                               filter_double_less, filter_double_greater,
                               intensity_col_sep, intensity_col_id, na_string, verbose = TRUE) {

    # an empty secondary_id would crash the reader; an unnamed filter would be ignored, and a
    # double filter that is not a number read as 0
    .check_id_columns(sample_id, primary_id, secondary_id)
    .check_filters(filter_string_equal, filter_string_not_equal, filter_double_less, filter_double_greater)

    sep <- rawToChar(as.raw(1))

    cmd <- paste("--sample", sample_id,
                 "--primary", primary_id,
                 "--secondary", paste(secondary_id, collapse = sep),
                 "--quant", paste(intensity_col, collapse = sep), sep = sep)

    if (!is.null(annotation_col)) {
        cmd <- paste(cmd, "--annotation", paste(annotation_col, collapse = sep), sep = sep)
    }

    if (!is.null(filter_string_equal)) {
        for (f in names(filter_string_equal)) {
            cmd <- paste(cmd, "--filter-string-equal", f, filter_string_equal[f], sep = sep)
        }
    }

    if (!is.null(filter_string_not_equal)) {
        for (f in names(filter_string_not_equal)) {
            cmd <- paste(cmd, "--filter-string-not-equal", f, filter_string_not_equal[f], sep = sep)
        }
    }

    if (!is.null(filter_double_less)) {
        for (f in names(filter_double_less)) {
            cmd <- paste(cmd, "--filter-double-less", f, filter_double_less[f], sep = sep)
        }
    }

    if (!is.null(filter_double_greater)) {
        for (f in names(filter_double_greater)) {
            cmd <- paste(cmd, "--filter-double-greater", f, filter_double_greater[f], sep = sep)
        }
    }

    if (!is.null(intensity_col_sep)) {
        cmd <- paste(cmd, "--intensity_col_sep", intensity_col_sep, sep = sep)
    }

    # na_string: for values split from one column, or one value per column (several columns)
    if ((!is.null(intensity_col_sep) || length(intensity_col) > 1) && !is.null(na_string)) {
        cmd <- paste(cmd, "--na_string", na_string, sep = sep)
    }

    if (!is.null(intensity_col_id)) {
        cmd <- paste(cmd, "--intensity_col_id", paste(intensity_col_id, collapse = sep), sep = sep)
    }

    # --quiet: the reader prints nothing but errors
    if (!verbose) {
        cmd <- paste(cmd, "--quiet", sep = sep)
    }

    as.character(paste(cmd, filename, sep = sep))
}

# At most n evenly spaced values of v, for a boxplot (all values if n is NULL or v is shorter).
.thin_for_plot <- function(v, n) {
    if (is.null(n) || length(v) <= n) {
        return(v)
    }
    v[round(seq(1, length(v), length.out = floor(n)))]
}

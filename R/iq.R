#########################################################################
#
# Author: Thang V. Pham, t.pham@amsterdamumc.nl
#
# All rights reserved.
#
# Citation:
#
# Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein abundances from ion quantification in DIA-MS-based proteomics, Bioinformatics 2020 Apr 15;36(8):2611-2613.
#
# Software version: 1.9.7
#
#########################################################################

preprocess <- function(quant_table,
                       primary_id = "PG.ProteinGroups",
                       secondary_id = c("EG.ModifiedSequence", "FG.Charge", "F.FrgIon", "F.Charge"),
                       sample_id = "R.Condition",
                       intensity_col = "F.PeakArea",
                       median_normalization = TRUE,
                       log2_intensity_cutoff = 0,
                       pdf_out = "qc-plots.pdf",
                       pdf_width = 12,
                       pdf_height = 8,
                       intensity_col_sep = NULL,
                       intensity_col_id = NULL,
                       na_string = "0") {

    if (!is.null(pdf_out)) {
        pdf(pdf_out, pdf_width, pdf_height)
    }

    if (is.null(intensity_col_sep)) {
        if (!is.numeric(quant_table[, intensity_col])) {
            stop("Intensity column must be numeric when 'intensity_col_sep' is NULL.")
        }

        # make secondary ids
        message("Concatenating secondary ids...\n")
        second_id <- as.character(quant_table[, secondary_id[1]])
        if (length(secondary_id) > 1) {
            for (i in 2:length(secondary_id)) {
                second_id <- paste(second_id, as.character(quant_table[, secondary_id[i]]), sep = "_")
            }
        }

        d <- data.frame(protein_list = as.character(quant_table[, primary_id]),
                        sample_list = as.character(quant_table[, sample_id]),
                        quant = log2(quant_table[, intensity_col]),
                        id = second_id)
    } else {
        # make secondary ids
        message("Concatenating secondary ids...\n")
        second_id <- as.character(quant_table[, secondary_id[1]])
        if (length(secondary_id) > 1) {
            for (i in 2:length(secondary_id)) {
                second_id <- paste(second_id, as.character(quant_table[, secondary_id[i]]), sep = "_")
            }
        }

        v <- strsplit(as.character(quant_table[, intensity_col]), intensity_col_sep)
        if (!is.null(intensity_col_id)) {
            v_ids <- strsplit(as.character(quant_table[, intensity_col_id]), intensity_col_sep)
        }

        v_lengths <- lengths(v)
        n_row <- sum(v_lengths)

        d <- data.frame(protein_list = rep("", n_row),
                        sample_list = rep("", n_row),
                        quant = rep(0, n_row),
                        id = rep("", n_row))
        cc <- 1
        for (r in 1:nrow(quant_table)) {
            a <- v[[r]]
            a[a == na_string] <- NA

            cc_end <- cc + length(a) - 1

            d$protein_list[cc:cc_end] <- as.character(quant_table[r, primary_id])
            d$sample_list[cc:cc_end] <- as.character(quant_table[r, sample_id])
            d$quant[cc:cc_end] <- log2(as.numeric(a))

            if (is.null(intensity_col_id)) {
                d$id[cc:cc_end] <- paste0(second_id[r], "_", 1:length(a))
            } else {
                d$id[cc:cc_end] <- paste0(second_id[r], "_", v_ids[[r]])
            }
            cc <- cc_end + 1
        }
    }

    d <- d[complete.cases(d),] # remove NaN quant

    samples <- unique(d$sample_list)

    # intensity cut off
    if (!is.null(log2_intensity_cutoff)) {
        message("Removing low intensities...\n")
        if (!is.null(pdf_out)) {
            a <- hist(d$quant, 100, col = "steelblue", border = "steelblue", freq = FALSE,
                      main = "Histogram of log2 intensities",
                      xlab = "log2 intensity")
            arrows(log2_intensity_cutoff, 0, log2_intensity_cutoff, max(a$density) / 2.0 , col = "red", code = 1, lwd = 2)
        }
        d <- d[d$quant > log2_intensity_cutoff,]
    }

    if (!is.null(pdf_out)) {
        dl <- list()
    }

    m <- rep(NA, length(samples))

    for (i in 1:length(samples)) {
        v <- d$quant[d$sample_list == samples[i]]
        m[i] <- median(v, na.rm = TRUE)

        if (!is.null(pdf_out)) {
            dl[i] <- list(v)
        }
    }

    if (!is.null(pdf_out)) {
        message("Barplotting raw data ...\n")

        boxplot(dl,
                names = as.character(samples),
                main = "Boxplot of fragment intensities per sample",
                ylab = "log2 intensity",
                outline = FALSE,
                col = "steelblue",
                whisklty = 1,
                staplelty = 0,
                las = 2)
    }

    if (median_normalization) {

        message("Median normalization ...\n")

        f <- mean(m) - m

        # normalization
        for (i in 1:length(samples)) {
            idx <- d$sample_list == samples[i]
            d$quant[idx] <- d$quant[idx] + f[i]
        }

        if (!is.null(pdf_out)) {
            message("Barplotting after normalization ...\n")

            dl <- list()
            for (i in 1:length(samples)) {
                dl[i] <- list(d$quant[d$sample_list == samples[i]])
            }

            boxplot(dl,
                    names = as.character(samples),
                    main = "Boxplot of fragment intensities per sample",
                    ylab = "log2 intensity",
                    outline = FALSE,
                    col = "steelblue",
                    whisklty = 1,
                    staplelty = 0,
                    las = 2)
        }
    }

    if (!is.null(pdf_out)) {
        dev.off()
    }

    return(d)
}

create_protein_list <- function(preprocessed_data) {

    # check for NA values
    if (any(is.na(preprocessed_data$protein_list))) {
        stop("NA value in $protein_list.\n")
    }

    if (any(is.na(preprocessed_data$sample_list))) {
        stop("NA value in $sample_list.\n")
    }

    if (any(is.na(preprocessed_data$id))) {
        stop("NA value in $id.\n")
    }

    if (any(is.na(preprocessed_data$quant))) {
        stop("NA value in $quant.\n")
    }

    proteins <- unique(preprocessed_data$protein_list)
    samples <- unique(preprocessed_data$sample_list)

    p_list <- list()

    message("# proteins = ", length(proteins), ", # samples = ", length(samples))

    thres_display <- length(proteins) / 20

    for (i in 1:length(proteins)) {

        if (i > thres_display) {
            message(format(i * 100 / length(proteins), digits = 2), "%")
            thres_display <- i + length(proteins) / 20
        }

        tmp <- preprocessed_data[preprocessed_data$protein_list == proteins[i],]

        if (nrow(tmp) > 0) { # in the filtering step, we might have excluded some proteins, but the levels remain

            char_r <- as.character(tmp$id)
            char_c <- as.character(tmp$sample_list)

            id <- unique(as.character(tmp$id))

            m <- matrix(NA, nrow = length(id), ncol = length(samples), dimnames = list(id, samples))

            for (j in 1:nrow(tmp)) {
                if (is.na(m[char_r[j], char_c[j]])) {
                    m[char_r[j], char_c[j]] <- tmp[j, "quant"]
                } else {
                    message(j, " : sample ", char_c[j], "; id ", char_r[j], " not unique.\n")
                    stop("duplicate entry.\n")
                }
            }
            p_list[[as.character(proteins[i])]] <- m
        }
    }

    message("Completed.")

    return(p_list)
}


plot_protein <- function(X, main = "", col = NULL, split = 0.6, ...) {

    if (is.null(col)) {
        col <- 1:nrow(X)
    }

    if (!is.null(split)) {
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par))

        par(fig = c(0, split, 0, 1), mar = c(5, 4, 4, 0) + 0.1)
        on.exit(par(old_par))
    }

    matplot(t(X),
            type = "o", lty = 1, lwd = 2, pch = 19,
            col = col,
            ylab = "Intensity",
            axes = FALSE,
            main = main,
            xlab = "Sample")

    axis(2)
    axis(side = 1, at = 1:ncol(X), labels = colnames(X), las = 2)

    if (!is.null(split)) {
        par(fig = c(split, 1, 0, 1), xpd = TRUE, new = FALSE)
        on.exit(par(old_par))

        legend("topleft",
               legend = rownames(X),
               col = col,
               pch = 19,
               bty = "n", ...)

        par(old_par)
    }

    invisible(NULL)
}

maxLFQ <- function(X) {

    if (all(is.na(X))) {
        return(list(estimate = NA, annotation = "NA"))
    }

    if (nrow(X) == 1) {
        return(list(estimate = as.numeric(X[1,]), annotation = ""))
    }

    N <- ncol(X)

    cc <- 0
    g <- rep(NA, N) # group label

    spread <- function(i) {
        g[i] <<- cc
        for (r in 1:nrow(X)) {
            if (!is.na(X[r, i])) {
                for (k in 1:ncol(X)) {
                    if (!is.na(X[r, k]) && is.na(g[k])) {
                        spread(k)
                    }
                }
            }
        }
    }

    maxLFQ.do <- function(X) {

        N <- ncol(X)

        AtA <- matrix(0, nrow = N, ncol = N)
        Atb <- matrix(0, nrow = N, ncol = 1)

        for (i in 1:(N-1)) {
            for (j in (i+1):N) {

                r_i_j <- median(-X[,i] + X[,j], na.rm = TRUE)

                if (!is.na(r_i_j)) {

                    AtA[i, j] <- AtA[j, i] <- -1
                    AtA[i, i] <- AtA[i, i] + 1
                    AtA[j, j] <- AtA[j, j] + 1

                    Atb[i] <- Atb[i] - r_i_j
                    Atb[j] <- Atb[j] + r_i_j
                }
            }
        }

        res <- lsfit(rbind(cbind(2*AtA , rep(1, N)),
                           c(rep(1, N), 0)),
                     c(2*Atb, mean(data.matrix(X), na.rm = TRUE) * N),
                     intercept = FALSE)

        return(res$coefficients[1:N])
    }

    for (i in 1:N) {
        if (is.na(g[i])) {
            cc <- cc + 1
            spread(i)
        }
    }

    w <- rep(NA, N)
    for (i in 1:cc) {
        ind <- g == i
        if (sum(ind) == 1) {
            w[ind] <- median(X[, ind], na.rm = TRUE)
        } else {
            w[ind] <- maxLFQ.do(X[, ind])
        }
    }

    if (all(is.na(w))) {
        return(list(estimate = w, annotation = "NA"))
    } else {
        quantified_samples <- which(!is.na(w))
        if (all(g[quantified_samples] == g[quantified_samples[1]])) {
            return(list(estimate = w, annotation = ""))
        } else {
            g[is.na(w)] <- NA
            return(list(estimate = w, annotation = paste(g, collapse = ";")))
        }
    }
}

median_polish <- function(X) {
    out  <-  medpolish(X, na.rm = TRUE, trace.iter = FALSE)
    return(list(estimate = out$overall + out$col, annotation = ""))
}


topN <- function(X, N = 3, aggregation_in_log_space = TRUE) {

    if (nrow(X) == 1) {
        return(list(estimate = X[1,], annotation = ""))
    }

    if (aggregation_in_log_space) {
        v <- rowMeans(X, na.rm = TRUE)
        v_sorted <- sort(v, decreasing = TRUE, na.last = TRUE, index.return = TRUE)
        out <- colMeans(X[v_sorted$ix[1:min(N, length(v))],], na.rm = TRUE)
    } else {
        XX <- 2^X
        v <- rowMeans(XX, na.rm = TRUE)
        v_sorted <- sort(v, decreasing = TRUE, na.last = TRUE, index.return = TRUE)
        out <- log2(colMeans(XX[v_sorted$ix[1:min(N, length(v))],], na.rm = TRUE))
    }

    out[is.nan(out)] <- NA

    return(list(estimate = out, annotation = ""))
}

meanInt <- function(X, aggregation_in_log_space = TRUE) {

    if (nrow(X)==1) {
        return(list(estimate = X[1,], annotation = ""))
    }

    if (aggregation_in_log_space) {
        out <- colMeans(X, na.rm = TRUE)
    } else {
        out <- log2(colMeans(2.0^X, na.rm = TRUE))
    }

    out[is.nan(out)] <- NA
    return(list(estimate = out, annotation = ""))
}

create_protein_table <- function(protein_list, method = "maxLFQ", ...) {

    if (!is.list(protein_list)) {
        stop("Input is not a list.")
    }

    if (length(protein_list) == 0) {
        return(NA)
    }

    tab <- matrix(NA,
                  nrow = length(protein_list),
                  ncol = ncol(protein_list[[1]]),
                  dimnames = list(names(protein_list), colnames(protein_list[[1]])))

    annotation <- rep(NA, length(protein_list))

    thres_display <- nrow(tab) / 20

    for (i in 1:nrow(tab)) {

        if (i > thres_display) {
            message(format(i * 100 / nrow(tab), digits = 2), "%")
            thres_display <- i + nrow(tab) / 20
        }

        if (method == "maxLFQ") {
            out <- maxLFQ(protein_list[[i]], ...)
        } else if (method == "median_polish") {
            out <- median_polish(protein_list[[i]], ...)
        } else if (method == "topN") {
            out <- topN(protein_list[[i]], ...)
        } else if (method == "meanInt") {
            out <- meanInt(protein_list[[i]], ...)
        } else {
            stop(paste0("Unknown method: ", method))
        }
        tab[i,] <- out$estimate
        annotation[i] <- out$annotation
    }

    message("Completed.")

    return(list(estimate = tab, annotation = annotation))
}

extract_annotation <- function(protein_ids,
                               quant_table,
                               primary_id = "PG.ProteinGroups",
                               annotation_columns = NULL) {

    all_columns <- c(primary_id, annotation_columns)

    index <- match(all_columns, colnames(quant_table))
    # check for NA values
    if (any(is.na(index))) {
        stop(paste0("The input table has no column: ", all_columns[which(is.na(index))[1]]))
    }

    index <- match(protein_ids, quant_table[, primary_id])
    if (any(is.na(index))) {
        stop(paste0("Cannot find ", protein_ids[which(is.na(index))[1]]))
    }

    tab <- quant_table[index, c(primary_id, annotation_columns)]
    rownames(tab) <- protein_ids

    return(tab)
}

process_wide_format <- function(input_filename,
                                output_filename,
                                id_column,
                                quant_columns,
                                data_in_log_space = FALSE,
                                annotation_columns = NULL,
                                method = "maxLFQ") {

    message("Reading file: ", input_filename)
    d <- read.delim(input_filename)
    message("# rows = ", nrow(d), "; # cols = ", ncol(d))

    if (!(id_column %in% colnames(d))) {
        stop("'id_column' not found")
    }

    if (is.numeric(quant_columns)) {
        if (max(quant_columns) > ncol(d)) {
            stop("'quant_columns' out of range")
        }
        if (min(quant_columns) < 1) {
            stop("'quant_columns' out of range")
        }
    } else {
        if (!all(quant_columns %in% colnames(d))) {
            stop("One or more column name in 'quant_columns' not found")
        }
    }

    if (!all(annotation_columns %in% colnames(d))) {
        stop("One or more column name in 'annotation_columns' not found")
    }

    message("Preparing data...")

    # collapse to unique 'id_column'
    d <- d[!is.na(d[, id_column]),]

    d <- d[d[, id_column] != "",] # bad for rownames

    d_long <- reshape(d,
                      direction = "long",
                      varying = quant_columns,
                      timevar = "sample_list",
                      times = colnames(d)[quant_columns],
                      v.names = "quant")

    d_long <- d_long[!is.na(d_long$quant) & !is.na(d_long[, id_column]),]

    unique_values <- unique(d_long[, id_column])

    if (!is.function(method) && method == "maxLFQ") {

        df <- data.frame(protein_list = d_long[, id_column],
                         sample_list = d_long$sample_list,
                         id = d_long$id,
                         quant = if (!data_in_log_space) log2(d_long$quant) else d_long$quant)

        message("Running MaxLFQ...")

        res <- fast_MaxLFQ(df)

        res$estimate <- res$estimate[unique_values, ]

    } else {

        # set up result tables
        res <- list(estimate = matrix(NA, nrow = length(unique_values), ncol = length(quant_columns)))

        if (is.numeric(quant_columns)) {
            colnames(res$estimate) <- colnames(d)[quant_columns]
        }
        else {
            colnames(res$estimate) <- quant_columns
        }
        rownames(res$estimate) <- unique_values

        thres_display <- length(unique_values)/20

        message("Processing...")
        for (r in 1:length(unique_values)) {

            if (r > thres_display) {
                message(format(r * 100/length(unique_values), digits = 2), "%")
                thres_display <- r + length(unique_values)/20
            }

            a <- d[d[, id_column] == unique_values[r], quant_columns]

            if (!data_in_log_space) {
                a <- log2(a)
            }

            if (is.function(method)) {
                out <- list(estimate = method(a))
            } else if (method == "maxLFQ_R") {
                out <- maxLFQ(a)
            }
            else if (method == "median_polish") {
                out <- median_polish(a)
            }
            else if (method == "top3") {
                out <- topN(a, N = 3)
            }
            else if (method == "top5") {
                out <- topN(a, N = 5)
            }
            else if (method == "maxInt") {
                out <- topN(a, N = 1)
            }
            else if (method == "meanInt") {
                out <- meanInt(a)
            }
            else if (method == "sum") {
                out <- list(estimate = log2(colSums(2^a, na.rm = TRUE)))
            }
            else if (method == "least_na") {
                out <- list(estimate = a[which.min(rowSums(is.na(a))),])
            }
            else {
                stop(paste0("Unknown method: ", method))
            }

            res$estimate[r,] <- as.numeric(out$estimate)
        }
    }

    message("Writing to: ", output_filename)

    if (is.null(annotation_columns)) {

        tab <- cbind(rownames(res$estimate), res$estimate)
        colnames(tab)[1] <- id_column

    } else {

        extra_annotation <- extract_annotation(rownames(res$estimate),
                                               d,
                                               primary_id = id_column,
                                               annotation_columns = annotation_columns)

        tab <- cbind(rownames(res$estimate),
                     extra_annotation[, annotation_columns],
                     res$estimate)
        colnames(tab)[1:(length(annotation_columns)+1)] <- c(id_column, annotation_columns)

    }

    write.table(tab, output_filename, sep = "\t", row.names = FALSE, quote = FALSE)

    invisible(NULL)
}

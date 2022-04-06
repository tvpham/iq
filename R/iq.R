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
# Software version: 1.9.3
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
                       pdf_height = 8) {

    if (!is.null(pdf_out)) {
        pdf(pdf_out, pdf_width, pdf_height)
    }

    if (class(quant_table[, intensity_col]) != "numeric") {
        stop("Intensity column must be numeric.")
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

    dl <- list()
    m <- rep(NA, length(samples))
    for (i in 1:length(samples)) {
        dl[i] <- list(d[d$sample_list == samples[i], "quant"])
        m[i] <- median(dl[[i]], na.rm = TRUE)
    }

    if (!is.null(pdf_out)) {
        message("Barplotting raw data ...\n")

        boxplot(dl,
                names = as.character(samples),
                col = "steelblue",
                whisklty = 1,
                staplelty = 0,
                outpch = ".",
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
            m <- rep(NA, length(samples))
            for (i in 1:length(samples)) {
                dl[i] <- list(d[d$sample_list == samples[i], "quant"])
                m[i] <- median(dl[[i]], na.rm = TRUE)
            }

            boxplot(dl,
                    names = as.character(samples),
                    col = "steelblue",
                    whisklty = 1,
                    staplelty = 0,
                    outpch = ".",
                    las = 2)
        }
    }

    if (!is.null(pdf_out)) {
        dev.off()
    }

    return(d)
}


create_protein_list <- function(preprocessed_data) {

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
                    stop("Error: duplicate entry.\n")
                }
            }
            p_list[[as.character(proteins[i])]] <- m
        }
    }

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
                     c(2*Atb, mean(X, na.rm = TRUE) * N),
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
            out <- meanInt(protein_list[[i]])
        } else {
            stop(paste0("Unknown method: ", method))
        }
        tab[i,] <- out$estimate
        annotation[i] <- out$annotation
    }
    return(list(estimate = tab, annotation = annotation))
}

extract_annotation <- function(protein_ids,
                               quant_table,
                               primary_id = "PG.ProteinGroups",
                               annotation_columns = NULL) {

    all_columns <- c(primary_id, annotation_columns)

    index <- match(all_columns, colnames(quant_table))

    if (any(index == 0)) {
        stop(paste0("The input table has no column: ", all_columns[which(index == 0)[1]]))
    }

    index <- match(protein_ids, quant_table[, primary_id])
    if (any(index == 0)) {
        stop(paste0("Cannot find ", protein_ids[which(index == 0)[1]]))
    }

    tab <- quant_table[index, c(primary_id, annotation_columns)]
    rownames(tab) <- protein_ids

    return(tab)
}

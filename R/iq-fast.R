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

fast_read <- function(filename,
                      sample_id = "R.Condition",
                      primary_id = "PG.ProteinGroups",
                      secondary_id = c("EG.ModifiedSequence", "FG.Charge", "F.FrgIon", "F.Charge"),
                      intensity_col = "F.PeakArea",
                      annotation_col = c("PG.Genes", "PG.ProteinNames"),
                      filter_string_equal = c("F.ExcludedFromQuantification" = "False"),
                      filter_double_less = c("PG.Qvalue" = "0.01", "EG.Qvalue" = "0.01")) {

    cmd <- paste0("--sample ", sample_id,
                  " --primary ", primary_id,
                  " --secondary ", paste(secondary_id, collapse = " "),
                  " --quant ", intensity_col)

    if (!is.null(annotation_col)) {
        cmd <- paste0(cmd, " --annotation ", paste(annotation_col, collapse = " "))
    }

    if (!is.null(filter_string_equal)) {
        for (f in names(filter_string_equal)) {
            cmd <- paste0(cmd, " --filter-string-equal ", f, " ", filter_string_equal[f])
        }
    }

    if (!is.null(filter_double_less)) {
        for (f in names(filter_double_less)) {
            cmd <- paste0(cmd, " --filter-double-less ", f, " ", filter_double_less[f])
        }
    }

    return(.Call("iq_filter", as.character(paste0(cmd, " ", filename))))
}

fast_MaxLFQ <- function(norm_data, row_names = NULL, col_names = NULL) {

    if (is.null(row_names)) {
        proteins <- unique(as.character(norm_data$protein_list))
        p_list <- as.integer(factor(as.character(norm_data$protein_list), levels = proteins))
    } else {
        proteins <- as.character(row_names)
        p_list <- as.integer(norm_data$protein_list)
    }

    if (is.null(col_names)) {
        s <- unique(as.character(norm_data$sample_list))
        s_list <- as.integer(factor(as.character(norm_data$sample_list), levels = s))
    } else {
        s <- as.character(col_names)
        s_list <- as.integer(norm_data$sample_list)
    }

    ret <- .Call("iq_MaxLFQ",
                 list(protein_index = p_list,
                      ion_index = as.integer(factor(norm_data$id)),
                      sample_index = s_list,
                      quant = as.double(norm_data$quant)))

    rownames(ret$estimate) <- proteins[as.integer(rownames(ret$estimate))]

    colnames(ret$estimate) <- s[as.integer(colnames(ret$estimate))]

    return(ret)

}

fast_preprocess <- function(quant_table,
                            median_normalization = TRUE,
                            log2_intensity_cutoff = 0,
                            pdf_out = "qc-plots-fast.pdf",
                            pdf_width = 12,
                            pdf_height = 8) {

    if (!is.null(pdf_out)) {
        pdf(pdf_out, pdf_width, pdf_height)
    }

    d <- quant_table

    d$quant <- log2(d$quant)

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
        selected <- d$quant > log2_intensity_cutoff
        d$quant <- d$quant[selected]
        d$protein_list <- d$protein_list[selected]
        d$sample_list <- d$sample_list[selected]
        d$id <- d$id[selected]
    }


    dl <- list()
    m <- rep(NA, length(samples))
    for (i in 1:length(samples)) {
      dl[i] <- list(d$quant[d$sample_list == samples[i]])
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

        for (i in 1:length(samples)) {
            idx <- d$sample_list == samples[i]
            d$quant[idx] <- d$quant[idx] + f[i]
        }

        if (!is.null(pdf_out)) {
            message("Barplotting after normalization ...\n")

            dl <- list()
            m <- rep(NA, length(samples))
            for (i in 1:length(samples)) {
                dl[i] <- list(d$quant[d$sample_list == samples[i]])
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

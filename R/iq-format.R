long_format_to_iq_format <- function(input_data, 
                                     output_data,
                                     sample_id = "Run",
                                     primary_id = "Protein.Group",
                                     secondary_id = c("Precursor.Id"),
                                     intensity_col = "Intensities",
                                     annotation_col = NULL,
                                     filter_string_equal = NULL,
                                     filter_string_not_equal = NULL,
                                     filter_double_less = NULL,
                                     filter_double_greater = NULL,
                                     intensity_col_sep = ";",
                                     intensity_col_id = NULL, 
                                     na_string = "0",
                                     normalization = "median",
                                     log2_intensity_cutoff = 0,
                                     pdf_out = "qc-plots-iq.pdf",
                                     pdf_width = 12,
                                     pdf_height = 8,
                                     show_boxplot = FALSE) {
    
    iq_data <- fast_read(input_data,
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
                         na_string = na_string)
    
    if (normalization == "median") {
        median_normalization <- TRUE
    } else if (normalization == "none") {
        median_normalization <- FALSE
    } else {
        stop("Unknown normalization method.")
    }
    
    if (median_normalization) {
        norm_data <- fast_preprocess(iq_data$quant_table,
                                     median_normalization = median_normalization,
                                     log2_intensity_cutoff = log2_intensity_cutoff,
                                     pdf_out = pdf_out,
                                     pdf_width = pdf_width,
                                     pdf_height = pdf_height,
                                     show_boxplot = show_boxplot)
    }
    
    out_dir <- paste0(output_data, "/")
    dir.create(out_dir)
    
    write.table(iq_data$protein, 
                file = paste0(out_dir, "_protein.tsv"), sep = "\t", row.names = FALSE)
    
    colnames(iq_data$sample)[1] <- "sample"
    write.table(iq_data$sample, 
                file = paste0(out_dir, "_sample.tsv"), sep = "\t", row.names = FALSE)
    
    colnames(iq_data$ion)[1] <- "ion"
    write.table(iq_data$ion, 
                file = paste0(out_dir, "_ion.tsv"), sep = "\t", row.names = FALSE)
    
    if (median_normalization) {
        n_display <- length(norm_data$protein_list) %/% 20
        
        for (i in 1:length(norm_data$protein_list)) {
            
            if (n_display > 0 && i %% n_display == 0) {
                message(format(i * 100 / length(norm_data$protein_list), digits = 2), "%")
            }
            
            j <- norm_data$protein_list[i]
            
            cat(paste0(c(norm_data$sample_list[i],
                         norm_data$id[i],
                         2^norm_data$quant[i]), collapse = "\t"),
                "\n",
                file = paste0(out_dir, j, ".tsv"),
                sep = "", append = TRUE)
        }
    }
    else {
        n_display <- length(iq_data$quant_table$protein_list) %/% 20
        
        for (i in 1:length(iq_data$quant_table$protein_list)) {
            
            if (n_display > 0 && i %% n_display == 0) {
                message(format(i * 100 / length(iq_data$quant_table$protein_list), digits = 2), "%")
            }
            
            j <- iq_data$quant_table$protein_list[i]
            
            cat(paste0(c(iq_data$quant_table$sample_list[i],
                         iq_data$quant_table$id[i],
                         iq_data$quant_table$quant[i]), collapse = "\t"),
                "\n",
                file = paste0(out_dir, j, ".tsv"),
                sep = "", append = TRUE)
        }
    }
    
    message("Done.")
    
    invisible(NULL)
}


connected_component <- function(X) {
    
    if (ncol(X) == 1) {
        return(0)
    }
    
    .Call("iq_connected_components",
          list(X = as.double(X),
               M = as.integer(nrow(X)),
               N = as.integer(ncol(X))))
}


rescale <- function(x, X, method) {
    
    g <- connected_component(X)
    
    x_new <- x
    for (i in 0:max(g)) {
        ind <- g == i
        
        if (all(is.na(x_new[ind]))) next
        
        if (method == "median-mean") {
            if (sum(ind) == 1) {
                x_new[ind] <- median(data.matrix(X[, ind]), na.rm = TRUE) 
            }
            else {
                m <- mean(data.matrix(X[, ind]), na.rm = TRUE) 
                x_new[ind] <- x_new[ind] - mean(x_new[ind], na.rm = TRUE) + m
            }
        }
        else if (method == "median-median") {
            if (sum(ind) == 1) {
                x_new[ind] <- median(data.matrix(X[, ind]), na.rm = TRUE) 
            }
            else {
                m <- median(data.matrix(X[, ind]), na.rm = TRUE) 
                x_new[ind] <- x_new[ind] - median(x_new[ind], na.rm = TRUE) + m
            }
        }
        else if (method == "mean-mean") {
            if (sum(ind) == 1) {
                x_new[ind] <- mean(data.matrix(X[, ind]), na.rm = TRUE) 
            }
            else {
                m <- mean(data.matrix(X[, ind]), na.rm = TRUE) 
                x_new[ind] <- x_new[ind] - mean(x_new[ind], na.rm = TRUE) + m
            }
        }
        else if (method == "sum") {
            if (sum(ind) == 1) {
                x_new[ind] <- sum(data.matrix(X[, ind]), na.rm = TRUE) 
            }
            else {
                m <- sum(data.matrix(X[, ind]), na.rm = TRUE) 
                x_new[ind] <- x_new[ind] + (m - sum(mean(x_new[ind], na.rm = TRUE))) / sum(!is.na(x_new[ind]))
            }
        }
        else {
            message("unknown rescaling method")
        }
    }
    return(x_new)
}


process_matrix <- function(X, method, 
                           p1 = NULL, p2 = NULL, k = 1.345, min_M = 15, n_threads = -1) {
    
    method_list <- c(
        "maxlfq" = 0,
        "maxlfq_bit" = 1,
        "weighted_maxlfq" = 2,
        "median_polish" = 3,
        "weighted_median_polish" = 4,
        "rlm" = 5,
        "weighted_rlm" = 6)
    
    if (!(method %in% names(method_list))) {
        stop("Unknown quantitation method.\n")
    }
    
    if (nrow(X) == 1) {
        #return(list(x = X[1,], s = 0))
        return(as.numeric(X[1,]))
    }
    
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
    
    if (method == "maxlfq_bit") {
        if (is.null(p2)) {
            p2 <- 0         # memory level: middle
        }
    }
    
    .Call("iq_quant",
          list(X = as.double(X),
               M = as.integer(nrow(X)),
               N = as.integer(ncol(X)),
               method = as.integer(method_list[method]),
               p1 = as.double(p1),
               p2 = as.integer(p2),
               k = as.double(k),
               min_M = as.integer(min_M),
               n_threads = as.integer(n_threads)))$x
}


process_iq_format <- function(input_data,
                              output_filename = NULL,
                              method = "maxlfq_bit", 
                              p1 = NULL, p2 = NULL, k = 1.345, min_M = 15,
                              n_threads = -1,
                              rescale_method = NULL) {
    
    long_to_wide <- function(q, rows, ncol) { # sample, ion, quant
        m <- matrix(0, nrow = length(rows), ncol = ncol)
        rownames(m) <- as.character(rows)
        for (k in 1:nrow(q)) {
            m[as.character(q[k, 2]), q[k, 1]] <- m[as.character(q[k, 2]), q[k, 1]] + q[k, 3]
        }
        return(m)
    }
    
    dir_in <- paste0(input_data, "/")
    
    p <- read.delim(paste0(dir_in, "_protein.tsv"))
    s <- read.delim(paste0(dir_in, "_sample.tsv"))
    
    res <- matrix(NA, nrow = nrow(p), ncol = nrow(s))
    
    colnames(res) <- s$sample
    
    n_display <- nrow(p) %/% 20
    
    for (j in 1:nrow(p)) {
        
        if (n_display > 0 && j %% n_display == 0) {
            message(format(j * 100 / nrow(p), digits = 2), "%")
        }
        
        q <- read.delim(paste0(dir_in, j, ".tsv"), header = FALSE)
        
        rows <- unique(q$V2)
        
        m <- long_to_wide(q, rows, nrow(s))
        
        colnames(m) <- s$sample
        
        m <- log2(m)
        m[m == -Inf] <- NA
        
        if (is.function(method)) {
            res[j,] <- method(m)
        }
        else {
            res[j,] <- process_matrix(m, method = method, p1 = p1, p2 = p2, k = k, n_threads = n_threads)
        }
        
        if (!is.null(rescale_method)) {
            res[j, ] <- rescale(res[j, ], m, method = rescale_method)
        }
    }
    
    
    write.table(cbind(p, res), 
                file = ifelse(is.null(output_filename), paste0(input_data, "-quant.tsv"), output_filename),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    invisible(NULL)
}

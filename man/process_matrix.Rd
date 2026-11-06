\name{process_matrix}
\alias{process_matrix}
\title{
Process individual matrix 
}
\description{
Estimates a row vector best representing all rows in the input numeric matrix.
}
\usage{
process_matrix(X, 
               method = "maxlfq_bit", 
               p1 = NULL, p2 = NULL, k = 1.345, min_M = 15,
               rescale_method = "to_mean",
               n_threads = -1)

}
\arguments{
  \item{X}{A matrix, or a data frame with numeric columns, of intensities in the log2 space. Columns are samples and rows are fragment ions.}
  \item{method}{One of the quantification methods: \code{"maxlfq"}, \code{"maxlfq_bit"}, 
        \code{"weighted_maxlfq"}, \code{"median_polish"}, \code{"weighted_median_polish"}, \code{"rlm"}, 
        \code{"weighted_rlm"}, \code{"sum"}, \code{"mean"} or \code{"topN"}; or a function. The default is \code{"maxlfq_bit"}, as in \code{\link{process}} and \code{\link{process_iq_format}}. See Details for \code{"sum"}, \code{"mean"} and \code{"topN"}.
        A function is called with \code{X} and must return a numeric vector with one value per sample (column of \code{X}).
        Values that are not finite (\code{-Inf}, \code{Inf}, \code{NaN}) become \code{NA}.
        For example, \code{function(X) log2(colSums(2^X, na.rm = TRUE))} with \code{rescale_method = NULL} gives the log2 summed intensities per sample, with \code{NA} for a sample without data, as the method \code{"sum"}.}
  \item{p1}{Method parameter. For \code{"rlm"} and \code{"weighted_rlm"}, it is the robust parameter. 
        For \code{"median_polish"} and \code{"weighted_median_polish"}, it is the convergence parameter. Unused in other cases.}
  \item{p2}{Method parameter. For \code{"maxlfq_bit"}, \code{"maxlfq"} and \code{"weighted_maxlfq"}, it specifies the memory level of the solver (default 0). For \code{"median_polish"}, \code{"weighted_median_polish"}, \code{"rlm"} and \code{"weighted_rlm"}, it is the maximum number of optimization rounds. For \code{"topN"}, it is the number of ions (default 3).}
  \item{k}{Weighting function parameter.}
  \item{min_M}{Minimum number of rows for the weighted methods.}
  \item{rescale_method}{The result is rescaled with \code{\link{rescale}} using this method. The default \code{"to_mean"} shifts each connected component of samples so that the mean of its estimates equals the mean of its data (in log2 space), as in \code{\link{process}} and \code{\link{process_iq_format}}. \code{NULL} for no rescaling.}
  \item{n_threads}{The number of threads to be used. When \code{n_threads} is 0, the maximal number of CPU cores is used. 
  When \code{n_threads} is -1 (default), one CPU core less than the maximum is used, and so on.}
}
\details{
The methods \code{"sum"}, \code{"mean"} and \code{"topN"} summarize the ions of each sample directly: the log2 of the summed intensities, the mean of the log2 intensities, and the mean of the log2 intensities of the \code{p2} ions with the highest mean log2 intensity over the samples (as \code{\link{topN}}). A sample without data gets \code{NA}. To keep these values as they are, use \code{rescale_method = NULL}; the default \code{rescale_method} shifts them like the estimates of the other methods.

The methods \code{"maxlfq"} and \code{"weighted_maxlfq"} solve the MaxLFQ least-squares problem with the iterative solver of \code{"maxlfq_bit"}. A package built with the original MaxLFQ implementation (using the Eigen library, see the file \code{src/Makevars} of the source package) solves it exactly instead.
}
\value{
A numeric row vector, one value per sample (column of \code{X}). For the built-in methods, a matrix with a single row gives that row (before rescaling; the \code{"to_mean"} and \code{"to_sum"} rescaling leave it unchanged).
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198-1203.
}
\author{
Thang V. Pham
}
\examples{
# a data matrix for a single protein: log2 intensities, ions in rows, samples in columns
X <- matrix(c(10, 12, 12, NA, NA, 10), nrow = 2,
            dimnames = list(c("ion1", "ion2"),
                            c("sample1", "sample2", "sample3")))
X

# relative quantification: sample2 is 4-fold higher than sample1 (a difference of 2 in
# log2 space), the difference of ion1, the only ion they share; similarly, sample1 is
# 4-fold higher than sample3 (ion2)
process_matrix(X, method = "maxlfq_bit", rescale_method = NULL)  # 11 13 9

# rescaled to the total intensity: the same differences between samples, and the
# intensities of the estimates add up to those of the data
x <- process_matrix(X, method = "maxlfq_bit", rescale_method = "to_sum")
x                       # 10.93 12.93 8.93
sum(2^X, na.rm = TRUE)  # 10240
sum(2^x)                # 10240
}

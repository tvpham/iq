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
               method, 
               p1 = NULL, p2 = NULL, k = 1.345, min_M = 15, n_threads = -1)

}
\arguments{
  \item{X}{A matrix of intensities in the log2 space. Columns are samples and rows are fragment ions.}
  \item{method}{One of the quantification methods: \code{"maxlfq"}, \code{"maxlfq_bit"}, 
        \code{"weighted_maxlfq"}, \code{"median_polish"}, \code{"weighted_median_polish"}, \code{"rlm"}, 
        or \code{"weighted_rlm"}.}
  \item{p1}{Method parameter. For \code{"rlm"} and \code{"weighted_rlm"}, it is the robust parameter. 
        For \code{"median_polish"} and \code{"weighted_median_polish"}, it is the convergence parameter. Unused in other cases.}
  \item{p2}{Method parameter. For \code{"maxlfq_bit", it specifies the memory level. For others, it is the maximum number of optimization rounds.}}
  \item{k}{Weighting function parameter.}
  \item{min_M}{Minimum number of rows for the weighted methods.}
  \item{n_threads}{The number of threads to be used. When \code{n_threads} is 0, the maximal number of CPU cores is used. 
  When \code{n_threads} is -1 (default), one CPU core less than the maximum is used, and so on.}
}
\value{
A numeric row vector.
}
\references{
Pham TV, Tran CTM, Henneman AA, Pham LHC, Le DG, Can AH, Bui PHL, Piersma SR, Jimenez CR. Boosting the speed and accuracy of protein quantification algorithms in mass spectrometry-based proteomics. \emph{Journal of Proteome Research} 2026 Feb 6; 25(2):1198–1203.
}
\author{
Thang V. Pham
}

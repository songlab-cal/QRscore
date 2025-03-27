#' QRscore: Nonparametric Quantile Rank Score Tests for Distributional Shifts 
#' in Gene Expression
#'
#' The QRscore package implements a family of nonparametric two-sample, and 
#' multi-sample tests for detecting shifts in central tendency (mean) and 
#' spread (variance/dispersion) in gene expression data.
#' The package includes functionality for negative binomial and zero-inflated 
#' negative binomial models, resampling-based and asymptotic approximations, 
#' and is designed to be robust, flexible and efficient.
#'
#' @section Main Functions:
#' - [`QRscoreTest()`]: Wrapper for one-, two-, and multi-sample QRscore tests.
#' - [`QRscoreGenetest()`]: DE pipeline over many genes.
#'
#' @docType package
#' @name QRscore
#' @aliases QRscore-package
#' @importFrom stats pnorm qnbinom dnbinom rnbinom
#' @importFrom pscl zeroinfl
#' @importFrom MASS glm.nb
#' @importFrom assertthat assert_that
#' @importFrom arrangements compositions
#' @importFrom hitandrun simplex.sample
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam
#' @importFrom dplyr %>%
#' @keywords package
"_PACKAGE"

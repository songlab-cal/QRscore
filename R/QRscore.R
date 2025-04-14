#' @keywords internal
.qrscore_flex_one_sample <- function(x, k, p, wList, alternative,
                                type, resamp_number) {
      # internal function for QRscore_Flex
    assertthat::assert_that(!any(x < 0))
    Sk <- x / sum(x)
    t <- sum(Sk^p * wList)
    message(date(), ": The test statistic for the data is ", t)
    message(date(), ": Resampling on simplex with ", resamp_number)
    res <- simplex.sample(n = k, N = resamp_number)$samples
    resampled_ts <- as.vector((res)^p %*% wList)
    cdf_t <- mean(resampled_ts < t)
    cdf_upp <- 1 - mean(c(resampled_ts, t) >= t)
    cdf_low <- mean(c(resampled_ts, t) <= t)
      
    if (alternative == "two.sided") {
        message(date(), ": Computing two-sided p-value")
        if (type == "unbiased") return(2 * min(cdf_t, 1 - cdf_t))
        if (type == "valid") return(2 * min(cdf_low, 1 - cdf_upp))
        return(setNames(c(2 * min(cdf_t, 1 - cdf_t),
                          2 * min(cdf_low, 1 - cdf_upp)),
                        c("unbiased", "valid")))
    }
  
    message(date(), ": Computing one-sided p-value with alternative = ",
            alternative)
    if (alternative == "greater") {
        if (type == "unbiased") return(1 - cdf_t)
        if (type == "valid") return(1 - cdf_upp)
        return(setNames(c(1 - cdf_t, 1 - cdf_upp),
                        c("unbiased", "valid")))
    }
    if (type == "unbiased") return(cdf_t)
    if (type == "valid") return(cdf_low)
    return(setNames(c(cdf_t, cdf_low), c("unbiased", "valid")))
}

#' @keywords internal
.qrscore_flex_two_sample <- function(x, y, k, p, wList, alternative,
                                approx, type, resamp_number) {
    # internal function for QRscore_Flex
    x_ordered <- sort(x)
    x_ordered <- c(-Inf, x_ordered, Inf)
    n <- length(y)
    
    Snk <- vapply(seq_len(k), function(i) {
        sum(y >= x_ordered[i] & y < x_ordered[i + 1])
    }, numeric(1))
    
    t <- if (approx == "resample" & !(n >= 100 & k >= 50)) {
        message(date(), ": Adjusting n and k for resampling")
        sum(((Snk + 1) / (n + k))^p * wList)
    } else {
        sum((Snk / n)^p * wList)
    }
  
    message(date(), ": The test statistic for the data is ", t)
    assertthat::assert_that(k / n >= 1e-3)
    
    if (approx == "asymptotic" | (n >= 100 & k >= 50)) {
        if (!(n >= 100 & k >= 50)) {
            warning("Sample sizes may be too small...")
        }
        message(date(), ": Applying Gaussian asymptotics...")
        moments <- .get_gaussion_asymptotics_moments(Snk, n, k, wList, p)
        z <- (t - moments$mean) / sqrt(moments$var)
        if (alternative == "two.sided") return(2 * min(pnorm(z), pnorm(-z)))
        if (alternative == "greater") return(pnorm(-z))
        return(pnorm(z))
    } else {
        message(date(), ": Using resampling approach with ", resamp_number)
        return(getCompositionPValue(
            t = t, n = n + k, k = k, p = p, wList = wList,
            alternative = alternative, type = type,
            resamp_number = resamp_number
        ))
    }
}

#' @keywords internal
.run_qrscore_one_sample <- function(samples, p, wList, alternative,
                                    approx, type, n_mom, resamp_number) {
  # internal function for QRscoreTest
  message(date(), ": Performing a one-sample test.")
  assertthat::assert_that(!is.null(p),
                          msg = "Exponent p must be provided."
  )
  assertthat::assert_that(!is.null(wList),
                          msg = "wList must be provided."
  )
  assertthat::assert_that(length(samples) == length(wList),
                          msg = "wList and sample lengths must match."
  )
  assertthat::assert_that(
    p > 0 &&
      abs(p - round(p)) < .Machine$double.eps^0.5,
    msg = "Only integer p supported."
  )
  assertthat::assert_that(approx == "resample",
                          msg = "Only 'resample' supported."
  )
  
  QRscore_Flex(
    samples, NULL, p, wList, alternative,
    approx, type, n_mom, resamp_number
  )
}


#' @keywords internal
.run_qrscore_one_sample <- function(samples, p, wList, alternative,
                               approx, type, n_mom, resamp_number) {
    # internal function for QRscoreTest
    message(date(), ": Performing a one-sample test.")
    assertthat::assert_that(!is.null(p),
        msg = "Exponent p must be provided."
    )
    assertthat::assert_that(!is.null(wList),
        msg = "wList must be provided."
    )
    assertthat::assert_that(length(samples) == length(wList),
        msg = "wList and sample lengths must match."
    )
    assertthat::assert_that(
        p > 0 &&
            abs(p - round(p)) < .Machine$double.eps^0.5,
        msg = "Only integer p supported."
    )
    assertthat::assert_that(approx == "resample",
        msg = "Only 'resample' supported."
    )

    QRscore_Flex(
        samples, NULL, p, wList, alternative,
        approx, type, n_mom, resamp_number
    )
}

#' @keywords internal
.run_qrscore_two_sample <- function(sample_list, unique_labels,
                               alternative, approx, resamp_number,
                               zero_inflation, LR.test, pi_threshold,
                               gene.name, measure, p, wList,
                               type, n_mom, use_base_r) {
    # internal function for QRscoreTest
    x <- sample_list[[1]]
    y <- sample_list[[2]]
    x_label <- unique_labels[1]
    y_label <- unique_labels[2]
    if (length(x) > length(y)) {
        tmp <- x; x <- y; y <- tmp
        tmp <- x_label; x_label <- y_label; y_label <- tmp
    }
    message(date(), ": Performing a two-sample test. '", x_label,
            "' is x and '", y_label, "' is y.")
    if (is.null(wList) || is.null(p)) {
        message(date(), ": Use ZINB or NB model to estimate weights.")
        assertthat::assert_that(all(c(x, y) >= 0),msg = 
                                  "Samples must not contain negative values.")
        return(QRscore_ZINB(x, y, zero_inflation, LR.test,
                            approx, alternative, resamp_number, pi_threshold,
                            gene.name, measure))
    }
    assertthat::assert_that(length(x) + 1 == length(wList),
                            msg = "wList must have length = length(x) + 1.")
    assertthat::assert_that(p > 0 &&
                              abs(p - round(p)) < .Machine$double.eps^0.5,
                            msg = "Only integer p supported.")
    assertthat::assert_that(approx == "resample",
                            msg = "Only 'resample' supported.")
    is_mann_whitney <- identical(wList, c(length(x):0)) && p == 1
    if (is_mann_whitney && use_base_r) {
        message(date(), ": Using base R Mann-Whitney test...")
        return(wilcox.test(x, y, alternative = alternative,
                           correct = FALSE)$p.value)
    }
    message(date(), ": Using QRscore_Flex with given weights.")
    return(QRscore_Flex(x, y, p, wList, alternative, approx,
                        type, n_mom, resamp_number))
  }




#' Flexible Non-Parametric One- and Two-Sample Tests (Native R version)
#'
#' Given data consisting of either a single sample
#' \eqn{\boldsymbol{x}=(x_1,\ldots,x_k)},
#' or two samples \eqn{\boldsymbol{x}=(x_1,\ldots,x_k)} and
#' \eqn{\boldsymbol{y}=(y_1,\ldots,y_n)},
#' this function uses summary statistics computed on weighted linear
#' combinations of powers of the spacing statistics
#' \eqn{S_k} (former) or \eqn{S_{n,k}} (latter).
#'
#' More precisely, this function does the following:
#'
#' For a single sample \eqn{x}, the function tests for uniformity of its
#' entries. When \eqn{p=2} and a particular choice of \eqn{\boldsymbol{w}} is
#' specified, we recover Greenwood's test.
#'
#' For two samples, the function tests the null of \eqn{\boldsymbol{x}} and
#' \eqn{\boldsymbol{y}} being drawn from the same distribution
#' (i.e., stochastic equality), against flexible alternatives that correspond
#' to specific choices of the test statistic parameters, \eqn{\boldsymbol{w}}
#' (weight vector) and \eqn{p} (power). These parameters not only determine the
#' test statistic \eqn{||S_k||_{p,\boldsymbol{w}}^p=\sum_{j=1}^k w_iS_{k}[j]^p}
#' (analogously defined for \eqn{||S_{n,k}||_{p,\boldsymbol{w}}^p}), but also
#' encode alternative hypotheses ranging from different populational means
#' (i.e., \eqn{\mu_x \neq \mu_y}), different populational spreads
#' (i.e., \eqn{\sigma^2_x \neq \sigma^2_y}), etc.
#'
#' Additional tuning parameters include:
#' (1) choice of p-value computation (one- or two-sided);
#' (2) approximation method
#' (3) number of moments accompanying the approximation chosen if using
#' moment-based approximation (recommended 200, typically at least 100); and
#' (4) in case of two samples, whether the user prefers to use exact discrete
#' moments (more accurate but slower) or to use continuous approximations of the
#' discrete moments (less accurate but faster).
#'
#' Currently, only resampling and Gaussian asymptotics are supported.
#' Both are efficient and well-calibrated. For \eqn{n\geqslant 100} and
#' \eqn{k\geqslant 50} such that \eqn{\frac{k}{n}\geqslant 0.001},
#' function automatically uses Gaussian approximation to the null.
#'
#' Dependencies: functions in `utils.R`
#' @param x First sample
#' @param y Second sample
#' @param p Exponent value in defining test statistic (must be integer)
#' @param wList Vector of weights. It should have length equal to \eqn{x} when
#' \eqn{y} is `NULL`, and one more than the length of \eqn{x} when \eqn{y} is
#' not `NULL`
#' @param alternative How p-value should be computed; i.e., a character
#' specifying the alternative hypothesis, must be one of "`two.sided`",
#' "`greater`" or "`less`"
#' @param approx Which approximation method to use (choose "`resample`",
#' "`asymptotic`")
#' @param type If using resampling approximation, either an unbiased estimate of
#' ("`unbiased`", default), or valid, but biased estimate of, ("`valid`")
#' p-value (see Hemerik and Goeman, 2018), or both ("`both`").
#' Default is "`unbiased`".
#' @param n_mom The number of moments to accompany the approximation
#' (recommended 200, if not at least 100)
#' @param resamp_number Number of \eqn{k}-compositions of \eqn{n} or simplex
#' vectors in \eqn{[0,1]^k}  to draw
#' @return Returns the p-value.
#' @export
#' @examples
#'
#' set.seed(1)
#' # One-sample examples
#' QRscore_Flex(
#'     x = abs(rnorm(10)), p = 2, wList = rep(1, 10),
#'     alternative = "two.sided", approx = "resample"
#' )
#'
#' # Two-sample examples
#' QRscore_Flex(
#'     x = abs(rnorm(30)), y = abs(rnorm(100)), p = 2,
#'     wList = rep(1, 31), alternative = "two.sided",
#'     approx = "resample", resamp_number = 5000
#' )
#'
#' QRscore_Flex(
#'     x = abs(rnorm(100)), y = abs(rnorm(100)), p = 1,
#'     wList = 0:100, alternative = "two.sided",
#'     approx = "asymptotic"
#' )
#' @export
QRscore_Flex <- function(x, y = NULL, p = 1, wList,
                         alternative = "two.sided",
                         approx = "resample",
                         type = "unbiased",
                         n_mom = NULL,
                         resamp_number = 5000) {
    k <- if (!is.null(y)) length(x) + 1 else length(x)
    if (!is.null(y)) assertthat::assert_that(length(wList) == k)
    assertthat::assert_that(p %in% c(1, 2))
    assertthat::assert_that(approx %in% c("resample", "asymptotic"))
    wList <- wList / max(wList)
    message(date(), ": Normalizing weight vector...")

    if (!is.null(y)) {
        return(.qrscore_flex_two_sample(
            x, y, k, p, wList, alternative,
            approx, type, resamp_number
        ))
    } else {
        return(.qrscore_flex_one_sample(
            x, k, p, wList, alternative,
            type, resamp_number
        ))
    }
}



#' Non-Parametric Two-Sample Tests Designed for Testing Differences in Mean or
#' Dispersion Parameters in (Zero-Inflated) Negative Binomial Distributions.
#'
#' This function evaluates the null hypothesis that two samples,
#' \eqn{\boldsymbol{x}} and \eqn{\boldsymbol{y}}, are drawn from the same
#' distribution, specifically designed for NB or ZINB models. It is particularly
#' effective in detecting shifts in either the mean or the dispersion
#' parameters.
#'
#' The function automatically computes optimal weights for the chosen model and
#' derives a p-value based on the selected test statistic and approximation
#' method.
#'
#' Additional tuning parameters include: (1) whether to use a likelihood ratio
#' test to determine which model (NB or ZINB) to fit, (2) the approximation
#' method (default is resampling, with asymptotic estimation for large 
#' samples),(3) choice of p-value computation (one- or two-sided), 
#' (4) threshold for estimated proportion of zeros in ZINB model 
#' (returns NA if exceeded).
#'
#' Dependencies: pscl::zeroinfl, MASS::glm.nb, and auxiliary functions from
#' `auxiliary.R`
#' @param x First sample
#' @param y Second sample
#' @param zero_inflation If TRUE, automatically chooses between ZINB and NB
#' models based on the data; if FALSE, applies NB model estimation.
#' @param LR.test Whether to use a likelihood ratio test to determine which
#' model (NB or ZINB) to fit
#' @param approx Which approximation method to use (default `resample`)
#' @param alternative How p-value should be computed; must be one of
#' "`two.sided`", "`greater`" or "`less`".
#' @param resamp_num Number of \eqn{k}-compositions of \eqn{n} or simplex
#' vectors in \eqn{[0,1]^k} to draw
#' @param pi_threshold Threshold for estimated proportion of zeros in ZINB model
#' @param gene.name Optional, name of the gene if applicable, used for
#' customized messages.
#' @param measure Specifies whether to test for shifts in "`mean`" or
#' "`dispersion`".
#' @param p_value If TRUE, returns a p-value, else returns test statistics and
#' weights.
#' @return p-value or test statistics depending on `p_value` parameter.
#' @examples
#'
#' set.seed(1)
#' # Two-sample example comparing mean shifts
#' QRscore_ZINB(
#'     x = rzinbinom(100, size = 2, mu = 20, pi = 0),
#'     y = rzinbinom(100, size = 2, mu = 30, pi = 0),
#'     zero_inflation = FALSE, LR.test = FALSE, alternative = "greater",
#'     approx = "asymptotic", measure = "mean"
#' )
#'
#' # Two-sample example comparing dispersion shifts
#' QRscore_ZINB(
#'     x = rzinbinom(100, size = 2, mu = 20, pi = 0.1),
#'     y = rzinbinom(100, size = 1, mu = 20, pi = 0.1),
#'     zero_inflation = TRUE, LR.test = TRUE, alternative = "two.sided",
#'     approx = "asymptotic", measure = "dispersion"
#' )
#'
#' # Two-sample example with significant zero inflation and variance shift
#' QRscore_ZINB(
#'     x = rzinbinom(30, size = 4, mu = 20, pi = 0.1),
#'     y = rzinbinom(30, size = 1, mu = 20, pi = 0.3),
#'     zero_inflation = TRUE, LR.test = FALSE, alternative = "two.sided",
#'     approx = "resample", resamp_num = 50000, measure = "dispersion"
#' )
#'
#' @export
#'
#'
QRscore_ZINB <- function(x, y, zero_inflation = TRUE, LR.test = FALSE,
                         approx = "resample", alternative = "two.sided",
                         resamp_num = 20000, pi_threshold = 0.95,
                         gene.name = NULL, measure = "mean", p_value = TRUE) {
    if (sum(round(x) != x) > 0 | sum(round(y) != y) > 0) {
        warning("Input values are not integers. Round input.")
        x <- round(x); y <- round(y)
    }
    assertthat::assert_that(!is.null(y), msg = "y should not be NULL.")
    assertthat::assert_that(all(x >= 0) & all(y >= 0),
        msg = "Input contains negative values."
    )
    combined <- c(x, y); n <- length(y); k <- length(x); alpha <- k / (n + k)
    computeweight <- ifelse(
        measure == "dispersion", computeweight_disp, computeweight_mean
        )
    if (sum(combined == 0) == 0 || !zero_inflation) {
        est <- .estimate_nb_parameters(combined)
    } else {
        est <- .fit_and_validate_zinb(
            combined, gene.name, pi_threshold, LR.test)
        if (is.null(est)) return(NA)
    }
    results <- computeweight(est$beta, est$mu, est$pi, n, k)
    if (!p_value) return(results)
    weights <- results$weight
    est_var <- alpha * (1 - alpha) * results$var
    zscore <- mean(weights[rank_x(x, y)])
    if (approx == "asymptotic" | (n >= 100 & k >= 50)) {
        zscore <- k / sqrt(n + k) * zscore
        p.value <- switch(alternative,
            "two.sided" = 2 * pnorm(
                abs(zscore), sd = sqrt(est_var), lower.tail = FALSE),
            "greater" = pnorm(
                zscore, sd = sqrt(est_var), lower.tail = FALSE),
            pnorm(zscore, sd = sqrt(est_var))
        )
    } else {
        message(date(), ": Using resampling with ", resamp_num)
        null_zscore <- replicate(
            resamp_num, mean(weights[sample(seq_len(n + k), k)])
        )
        p.permute <- mean(null_zscore < zscore)
        p.value <- switch(alternative,
            "two.sided" = 2 * min(p.permute, 1 - p.permute),
            "greater" = 1 - p.permute,
            p.permute)
    }
    return(p.value)
}

#' Multi-Sample Nonparametric Test for Mean or Dispersion Differences in
#' (Zero-Inflated) Negative Binomial Distributions.
#'
#' This function conducts statistical tests across multiple samples to evaluate
#' the null hypothesis that all groups are drawn from the same distribution. It
#' is optimized for data modeled by Negative Binomial (NB) or Zero-Inflated
#' Negative Binomial (ZINB) distributions and is capable of detecting shifts in
#' mean or dispersion parameters. The function can handle any number of groups
#' and automatically computes optimal weights for the specified measure (mean or
#' dispersion).
#'
#' The computation involves constructing a B matrix that transforms
#' group-specific scores into a space where independence among groups is
#' maximized. It then uses these transformed scores to calculate a test
#' statistic, which follows a chi-square distribution under the null 
#' hypothesis.
#'
#' Additional tuning parameters allow customization of the model fitting and
#' statistical testing, including:
#' - Selection between NB and ZINB models based on presence of zero inflation.
#' - Choice of approximation method for computing p-values - 'asymptotic' is
#'   recommended.
#' - Decision criteria for statistical tests (one-sided or two-sided).
#' - Threshold for the estimated proportion of zeros beyond which results are
#'   considered unreliable.
#'
#' Dependencies: Requires `pscl::zeroinfl` for zero-inflated models,
#' `MASS::glm.nb` for NB models, and other auxiliary functions as needed.
#'
#' @param samples Vector of all sample measurements.
#' @param labels Group labels for each sample.
#' @param zero_inflation Boolean, if TRUE, the function chooses between ZINB and
#'   NB models based on data; if FALSE, only NB model is applied.
#' @param LR.test Boolean, if TRUE, performs a likelihood ratio test to select
#'   between NB and ZINB models.
#' @param approx The method used for p-value approximation; "resample"
#'   (default) or "asymptotic".
#' @param resamp_num The number of resampling iterations used if `approx` is
#'   "resample".
#' @param pi_threshold Threshold for proportion of zeros at which to return NA,
#'   indicating unreliable results due to excessive zero inflation.
#' @param gene.name Optional, name of the gene if applicable, enhancing the
#'   relevance of output in genetic studies.
#' @param measure Specifies whether the test focuses on "mean" or "dispersion"
#'   differences.
#' @param perturb Boolean, if TRUE, adds small noise to data to avoid ties and
#'   improve model stability.
#' @return Returns the p-value of the test if `p_value` is TRUE, otherwise
#'   returns test statistics and weights.
#' @examples
#'
#' set.seed(1)
#' data <- c(
#'     rnbinom(100, size = 2, mu = 20), rnbinom(100, size = 2, mu = 25),
#'     rnbinom(100, size = 2, mu = 30)
#' )
#' labels <- factor(c(rep("Group1", 100), rep("Group2", 100), 
#'           rep("Group3", 100)))
#' QRscore_ZINB_nSamples(
#'     samples = data, labels = labels,
#'     zero_inflation = FALSE, LR.test = FALSE, approx = "resample",
#'     resamp_num = 5000, pi_threshold = 0.95, measure = "mean"
#' )
#'
#' # Example with zero inflation and dispersion shift detection
#' data_zi <- c(
#'     rzinbinom(100, size = 2, mu = 20, pi = 0.1),
#'     rzinbinom(100, size = 3, mu = 20, pi = 0.1),
#'     rzinbinom(100, size = 4, mu = 20, pi = 0.1)
#' )
#' labels_zi <- factor(c(rep("Group1", 100), rep("Group2", 100), 
#'               rep("Group3", 100)))
#' QRscore_ZINB_nSamples(
#'     samples = data_zi, labels = labels_zi,
#'     zero_inflation = TRUE, LR.test = TRUE, approx = "asymptotic",
#'     resamp_num = 2000, pi_threshold = 0.95, measure = "dispersion"
#' )
#' @export
QRscore_ZINB_nSamples <- function(samples, labels, zero_inflation = TRUE,
                                  LR.test = FALSE, approx = "resample",
                                  resamp_num = 20000, pi_threshold = 0.95,
                                  gene.name = NULL, measure = "mean",
                                  perturb = TRUE) {
    unique_labels <- unique(labels)
    n_groups <- length(unique_labels)
    sample_list <- lapply(unique_labels, function(l) {
        samples[labels == l]
    })
    group_sizes <- vapply(sample_list, length, integer(1))
    N_all <- sum(group_sizes)
    B_matrix <- .generateBMatrix(n_groups)
    group_proportion <- group_sizes / N_all
    AA_matrix <- .generateAAmatrix(group_proportion)

    x <- sample_list[[1]]
    y <- unlist(sample_list[-1], use.names = FALSE)
    results <- QRscore_ZINB(x, y,
        zero_inflation = zero_inflation,
        LR.test = LR.test, approx = approx,
        resamp_num = resamp_num,
        pi_threshold = pi_threshold, gene.name = gene.name,
        measure = measure, p_value = FALSE
    )

    if (!is.list(results) || is.null(results$weight)) {
        return(NA)
    }

    weights <- results$weight
    S_scores <- .computeSscores(sample_list, weights, N_all)
    Lambda_matrix <- diag(1 / sqrt(group_proportion))
    Test_vec <- t(B_matrix) %*% Lambda_matrix %*% S_scores
    H_new <- .get_H_new(results$var, B_matrix, AA_matrix, gene.name)
    H_inv <- .get_H_inv(H_new, gene.name)
    if (is.na(H_inv)[1]) {
        return(NA)
    }

    Q_all <- t(Test_vec) %*% H_inv %*% Test_vec
    p.value <- pchisq(Q_all, df = n_groups - 1, lower.tail = FALSE)
    return(p.value[1, 1])
}


#' QRscore Test
#'
#' This function performs statistical tests on data from one or more groups 
#' using summary statistics computed on weighted linear combinations of powers 
#' of spacing statistics. It is capable of conducting one-sample tests, 
#' two-sample tests, and multi-sample tests, utilizing either user-defined 
#' weights or automatically generated weights based on Negative Binomial (NB) 
#' or Zero-Inflated Negative Binomial (ZINB) models.
#'
#' For one-sample tests, the function assesses the uniformity of data entries.
#' For two-sample and multi-sample tests, it evaluates whether groups are drawn
#' from the same distribution, with alternative hypotheses considering
#' differences in means or dispersions.
#'
#' If the weights and \eqn{p} are given, the function calculates the test
#' statistic as: \deqn{||S_k||_{p,\boldsymbol{w}}^p=\sum_{j=1}^k
#' w_iS_{k}[j]^p} where \eqn{w_i} are weights, \eqn{x_i} are data points, and
#' \eqn{p} is the power specified.
#'
#' In two-sample and multi-sample settings without specified weights, the
#' function can automatically estimate weights using score function for a
#' Negative Binomial or a Zero-Inflated Negative Binomial model, optimizing for
#' dispersion or mean shifts.
#'
#' @param samples A numeric vector containing all sample measurements.
#' @param labels An optional vector of group labels corresponding to each entry
#'   in `samples`.
#' @param p The exponent used in the power sum test statistic, required if
#'   `wList` is not `NULL`.
#' @param wList An optional vector of weights; if `NULL`, weights are estimated
#'   using an NB or ZINB model for multiple groups.
#' @param alternative Specifies the alternative hypothesis; must be one of
#'   "two.sided", "greater", or "less".
#' @param approx The method used for p-value approximation, either "resample" or
#'   "asymptotic".
#' @param type Specifies if the estimation of the p-value should be "unbiased",
#'   "valid", or "both".
#' @param n_mom The number of moments to accompany the approximation, relevant
#'   for moment-based methods.
#' @param resamp_number The number of resampling iterations used if `approx` is
#'   "resample".
#' @param zero_inflation Indicates whether to account for zero inflation in
#'   model-based weight estimation.
#' @param LR.test Whether a likelihood ratio test is used to decide between NB
#'   and ZINB models.
#' @param pi_threshold Threshold for the proportion of zeros in ZINB models;
#'   results in NA if exceeded.
#' @param gene.name Optional identifier for a gene, used in output messages.
#' @param measure Specifies the statistical measure to be analyzed ("mean" or
#'   "dispersion") when weights are auto-generated.
#' @param perturb Boolean to indicate if data should be perturbed slightly to
#'   prevent ties.
#' @param use_base_r Boolean to decide whether to use base R functions for
#'   certain edge cases like Mann-Whitney tests.
#' @return Returns the p-value of the test.
#' @examples
#'
#' set.seed(1)
#' # One-sample test example with normally distributed data
#' data <- abs(rnorm(10))
#' QRscoreTest(data,
#'     p = 2, wList = rep(1, 10), alternative = "two.sided",
#'     approx = "resample"
#' )
#'
#' # Two-sample test with specified weights using normally distributed data
#' group1 <- rnorm(120, sd = 1)
#' group2 <- rnorm(120, sd = 2) # Different mean
#' data <- c(group1, group2)
#' labels <- c(rep("Group1", 120), rep("Group2", 120))
#' QRscoreTest(
#'     samples = data, labels = labels, p = 1, wList = c(60:0, seq_len(60)),
#'     alternative = "two.sided", approx = "resample"
#' )
#'
#' # Two-sample test with automatically estimated weights from NB model
#' group1 <- rzinbinom(120, size = 2, mu = 20, pi = 0)
#' group2 <- rzinbinom(100, size = 2, mu = 30, pi = 0) # Different mean
#' data <- c(group1, group2)
#' labels <- c(rep("Group1", 120), rep("Group2", 100))
#' QRscoreTest(
#'     samples = data, labels = labels,
#'     approx = "asymptotic", measure = "mean", zero_inflation = FALSE
#' )
#'
#' # Two-sample test with automatically estimated weights from ZINB model
#' group1 <- rzinbinom(100, size = 2, mu = 40, pi = 0.1)
#' group2 <- rzinbinom(200, size = 1, mu = 40, pi = 0.1)
#' data <- c(group1, group2)
#' labels <- c(rep("Group1", 100), rep("Group2", 200))
#' QRscoreTest(
#'     samples = data, labels = labels, alternative = "two.sided",
#'     approx = "asymptotic", measure = "dispersion"
#' )
#'
#' # Three-sample test with automatically estimated weights from NB model
#' group1 <- rzinbinom(150, size = 1, mu = 30, pi = 0.1)
#' group2 <- rzinbinom(100, size = 2, mu = 30, pi = 0.1)
#' group3 <- rzinbinom(30, size = 3, mu = 30, pi = 0.1)
#' data <- c(group1, group2, group3)
#' labels <- c(rep("Group1", 150), rep("Group2", 100), rep("Group3", 30))
#' QRscoreTest(
#'     samples = data, labels = labels, alternative = "two.sided",
#'     approx = "asymptotic", measure = "dispersion"
#' )
#' @export
#' 
QRscoreTest <- function(samples, labels = NULL, p = NULL,
                        wList = NULL, alternative = "two.sided",
                        approx = "resample", type = "unbiased",
                        n_mom = 100, resamp_number = 5000,
                        zero_inflation = TRUE, LR.test = FALSE,
                        pi_threshold = 0.95, gene.name = NULL,
                        measure = "mean", perturb = TRUE,
                        use_base_r = TRUE) {
    assertthat::assert_that(!is.null(samples),
        msg = "Samples cannot be NULL."
    )
    assertthat::assert_that(
        alternative %in% c("two.sided", "greater", "less"),
        msg = "Invalid alternative value."
    )
    assertthat::assert_that(
        approx %in% c("resample", "asymptotic"),
        msg = "Invalid approximation method."
    )
    assertthat::assert_that(
        type %in% c("unbiased", "valid", "both"), msg = "Invalid type value.")
    if (is.null(labels) || length(unique(labels)) == 1) {
        return(.run_qrscore_one_sample(
            samples, p, wList, alternative, approx, type, n_mom, resamp_number
        ))
    }
    assertthat::assert_that(length(labels) == length(samples),
        msg = "Sample and label lengths must match."
    )
    unique_labels <- sort(unique(labels))
    n_groups <- length(unique_labels)
    sample_list <- lapply(unique_labels, function(l) {
        samples[labels == l]
    })
    if (n_groups == 2) {
        return(.run_qrscore_two_sample(
            sample_list, unique_labels, alternative, approx, resamp_number,
            zero_inflation, LR.test, pi_threshold, gene.name, measure, p,
            wList, type, n_mom, use_base_r
        ))
    }
    message(
        date(), ": Performing multi-sample test. Use ZINB or NB model to",
        " estimate the weights."
    )
    return(QRscore_ZINB_nSamples(
        samples, labels, zero_inflation, LR.test, approx, resamp_number,
        pi_threshold, gene.name, measure, perturb
    ))
}
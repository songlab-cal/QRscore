#' @keywords internal
.dzinbinom <- function(x, mu, theta, size, pi, log = FALSE) {
    if (any(pi < 0) | any(pi > 1)) {
        warning("'pi' must be in [0, 1]")
    }
    if (!missing(theta) & !missing(size)) {
        stop("only 'theta' or 'size' may be specified")
    }
    if (!missing(size)) {
        theta <- size
    }
    rval <- log(1 - pi) + dnbinom(x, mu = mu, size = theta, log = TRUE)
    if (any(x0 <- (x == 0L))) {
        rval[x0] <- log(exp(rval) + pi)[x0]
    }
    if (log) {
        return(rval)
    } else {
        return(exp(rval))
    }
}

#' @keywords internal
.qzinbinom <- function(p, mu, theta, size, pi,
                       lower.tail = TRUE, log.p = FALSE) {
    if (any(pi < 0) | any(pi > 1)) {
        warning("'pi' must be in [0, 1]")
    }
    if (!missing(theta) & !missing(size)) {
        stop("only 'theta' or 'size' may be specified")
    }
    if (!missing(size)) {
        theta <- size
    }
    if (log.p) {
        p <- exp(p)
    }
    if (!lower.tail) {
        p <- 1 - p
    }
    p <- pmax(0, (p - pi) / (1 - pi))
    rval <- qnbinom(p,
        mu = mu, size = theta,
        lower.tail = TRUE, log.p = FALSE
    )
    return(rval)
}

#' @keywords internal
.get_gaussion_asymptotics_moments <- function(Snk, n, k, wList, p) {
    # internal function for .qrscore_flex_two_sample (QRscore_Flex)
    if (p == 1) {
        m1 <- sum(wList) / k
        m2 <- ((k / n + 1) / (k^2 * (k + 1))) *
            sum(wList * ((k * diag(k) - outer(
                rep(1, k),
                rep(1, k)
            )) %*% wList))
    } else {
        m1 <- ((2 + k / n - 1 / n) / ((k + 1) * k)) * sum(wList)
        w2 <- sum(wList^2)
        coef1 <- (k - 1) * (k / n + 1) * (2 + k / n - 1 / n) *
            (12 - 6 / n + k * (k / n + 10 - 5 / n)) /
            (k^2 * (1 + k)^2 * (2 + k) * (3 + k))
        offdiag <- sum(outer(wList, wList)) - w2
        coef2 <- (k / n + 1) * (6 / n^2 + k * (3 + k * (k - 2)) /
            n^2 - 24 / n + 8 * (k - 1) * k / n +
            8 * (3 + 2 * k)) /
            (k^2 * (1 + k)^2 * (2 + k) * (3 + k))
        m2 <- w2 * coef1 - offdiag * coef2
    }
    return(list(mean = m1, var = m2))
}

#' @keywords internal
.compute_exact_pval <- function(t, exact_ts, alternative) {
    # internal function for getCompositionPValue()
    if (alternative == "two.sided") {
        message(date(), ": Computing exact two-sided p-value")
        upper_tail <- mean(exact_ts >= t)
        cdf_at_t <- mean(exact_ts <= t)
        return(2 * min(cdf_at_t, upper_tail))
    } else if (alternative == "greater") {
        message(
            date(),
            ": Computing exact one-sided p-value with alternative `greater`"
        )
        return(mean(exact_ts >= t))
    } else {
        message(
            date(),
            ": Computing exact one-sided p-value with alternative `less`"
        )
        return(mean(exact_ts <= t))
    }
}

#' @keywords internal
.compute_resample_pval <- function(t, resampled_ts, alternative, type) {
    # internal function for getCompositionPValue()
    cdf_at_t <- mean(resampled_ts < t)
    cdf_at_t_upp_tail <- 1 - mean(c(resampled_ts, t) >= t)
    cdf_at_t_low_tail <- mean(c(resampled_ts, t) <= t)

    if (alternative == "two.sided") {
        message(date(), ": Computing two-sided p-value")
        if (type == "unbiased") {
            return(2 * min(cdf_at_t, 1 - cdf_at_t))
        } else if (type == "valid") {
            return(2 * min(cdf_at_t_low_tail, 1 - cdf_at_t_upp_tail))
        } else {
            return(c(
                unbiased = 2 * min(cdf_at_t, 1 - cdf_at_t),
                valid = 2 * min(cdf_at_t_low_tail, 1 - cdf_at_t_upp_tail)
            ))
        }
    } else if (alternative == "greater") {
        message(
            date(),
            ": Computing one-sided p-value with alternative `greater`"
        )
        if (type == "unbiased") {
            return(1 - cdf_at_t)
        } else if (type == "valid") {
            return(1 - cdf_at_t_upp_tail)
        } else {
            return(c(
                unbiased = 1 - cdf_at_t,
                valid = 1 - cdf_at_t_upp_tail
            ))
        }
    } else {
        message(date(), ": Computing one-sided
                p-value with alternative `less`")
        if (type == "unbiased") {
            return(cdf_at_t)
        } else if (type == "valid") {
            return(cdf_at_t_low_tail)
        } else {
            return(c(
                unbiased = cdf_at_t,
                valid = cdf_at_t_low_tail
            ))
        }
    }
}


#' @keywords internal
.Gamma_derivative_disp <- function(theta, x, beta) {
    # internal function for computeweight_disp
    x4 <- gamma(theta * beta)
    x1 <- beta * digamma(theta * beta + x)
    x2 <- (beta * (x4 * digamma(theta * beta))) / x4
    return(x1 - x2)
}

#' @keywords internal
.z_function_disp <- function(x, beta, mu, pi) {
    # internal function for computeweight_disp
    if (x == 0) {
        num <- (1 - pi) * (beta / (beta + mu))^(beta + 1) *
            ((beta + mu) * log(beta / (beta + mu)) + mu)
        denom <- (1 - pi) * (beta / (beta + mu))^beta + pi
        return(num / denom)
    } else {
        return((.Gamma_derivative_disp(1, x, beta) +
            beta * (mu - x) / (beta + mu) +
            beta * log(beta / (beta + mu))) * (1 - pi))
    }
}



#' @keywords internal
.estimate_nb_parameters <- function(combined_data) {
    # internal function for QRscore_ZINB
    message(date(), ": Estimating NB parameters.")
    fitNB <- MASS::glm.nb(combined_data ~ 1)
    return(list(
        beta = min(fitNB$theta, 100),
        mu = exp(fitNB$coefficients["(Intercept)"]),
        pi = 0
    ))
}

#' @keywords internal
.estimate_zinb_parameters <- function(combined_data, gene.name) {
    # internal function for .fit_and_validate_zinb (QRscore_ZINB)
    message(date(), ": Estimating ZINB parameters.")
    fitZINB <- pscl::zeroinfl(combined_data ~ 1 | 1, dist = "negbin")
    beta <- min(fitZINB$theta, 100)
    mu <- exp(fitZINB$coef$count[1])
    zcoef <- fitZINB$coef$zero[1]
    pi <- exp(zcoef) / (1 + exp(zcoef))
    return(list(fitZINB = fitZINB, beta = beta, mu = mu, pi = pi))
}

#' @keywords internal
.fit_and_validate_zinb <- function(combined, gene.name,
                                   pi_threshold, LR.test) {
    # internal function for QRscore_ZINB
    zinb <- .estimate_zinb_parameters(combined, gene.name)

    if (LR.test) {
        message(date(), ": Applying likelihood ratio test.")
        fitNB <- MASS::glm.nb(combined ~ 1)
        lrt_p <- as.numeric(pchisq(
            2 * (logLik(zinb$fitZINB) - logLik(fitNB)),
            df = 1, lower.tail = FALSE
        ))
        if (lrt_p >= 0.05) {
            return(.estimate_nb_parameters(combined))
        }
    }

    if (zinb$pi > pi_threshold) {
        msg <- if (is.null(gene.name)) {
            sprintf(
                "The estimated probability of zeros is > %.2f.",
                pi_threshold
            )
        } else {
            sprintf(
                "The estimated probability of zeros for gene %s is > %.2f.",
                gene.name, pi_threshold
            )
        }
        warning(sprintf("%s Returning NA.", msg))
        return(NULL)
    }

    return(zinb)
}

#' @keywords internal
.generateBMatrix <- function(d) {
    # internal function for QRscore_ZINB_nSmaples
    vectors <- matrix(0, nrow = d - 1, ncol = d)
    for (k in seq_len(d - 1)) {
        vector <- rep(0, d)
        vector[seq_len(k)] <- 1
        vector[k + 1] <- -k
        norm_factor <- sqrt(sum(vector^2))
        vectors[k, ] <- vector / norm_factor
    }
    return(t(vectors))
}

#' @keywords internal
.generateAAmatrix <- function(props) {
    # internal function for QRscore_ZINB_nSmaples
    m <- length(props)
    mat <- matrix(NA, m, m)
    for (i in seq_len(m)) {
        for (j in seq_len(m)) {
            mat[i, j] <- if (i == j) {
                1 - props[i]
            } else {
                -sqrt(props[i] * props[j])
            }
        }
    }
    return(mat)
}

#' @keywords internal
.computeSscores <- function(sample_list, weights, N_all) {
    # internal function for QRscore_ZINB_nSmaples
    scores <- numeric(length(sample_list))
    for (i in seq_along(sample_list)) {
        x <- sample_list[[i]]
        y <- unlist(sample_list[-i], use.names = FALSE)
        idx <- rank_x(x, y)
        scores[i] <- sum(weights[idx]) / sqrt(N_all)
    }
    return(scores)
}

#' @keywords internal
.get_H_new <- function(H_w, B_matrix, AA_matrix, gene.name) {
    # internal function for QRscore_ZINB_nSmaples
    H_modifier <- t(B_matrix) %*% AA_matrix %*% B_matrix
    H_modifier <- (H_modifier + t(H_modifier)) / 2
    H_new <- H_w * H_modifier
    if (any(is.na(H_new)) || any(is.infinite(H_new))) {
        msg <- "Covariance matrix contains NA or Inf values"
        if (!is.null(gene.name)) {
            msg <- sprintf("%s for gene %s", msg, gene.name)
        }
        warning(sprintf("%s, returning NA.", msg))
        return(NULL)
    }
    return(H_new)
}

#' @keywords internal
.get_H_inv <- function(H_new, gene.name) {
    # internal function for QRscore_ZINB_nSmaples
    if (isSymmetric(H_new) && all(eigen(H_new)$values > 0)) {
        return(chol2inv(chol(H_new)))
    }
    tryCatch(
        {
            solve(H_new)
        },
        error = function(e) {
            msg <- "Covariance matrix is not invertible"
            if (!is.null(gene.name)) {
                msg <- sprintf("%s for gene %s", msg, gene.name)
            }
            warning(sprintf("%s, returning NA.", msg))
            return(NA)
        }
    )
}

#' Generate random samples from a zero-inflated negative binomial distribution
#'
#' @param n Integer. Number of samples to generate.
#' @param mu Numeric. Mean of the negative binomial distribution.
#' @param theta Numeric. Dispersion parameter (size) of the negative binomial
#' distribution.
#' @param size Numeric. Alternative name for the dispersion parameter
#' (used interchangeably with `theta`).
#' @param pi Numeric. Zero-inflation probability; must be in the range
#' \eqn{[0,1]}.
#' @return A vector of random samples from a Zero-Inflated Negative Binomial
#' (ZINB) distribution.
#' @examples
#' set.seed(1)
#' rzinbinom(n = 5, mu = 10, theta = 5, pi = 0.2)
#' @export
rzinbinom <- function(n, mu, theta, size, pi) {
    if (any(pi < 0) | any(pi > 1)) {
        warning("'pi' must be in [0, 1]")
    }
    if (!missing(theta) &
        !missing(size)) {
        stop("only 'theta' or 'size' may be specified")
    }
    if (!missing(size)) {
        theta <- size
    }
    rval <- stats::rnbinom(n, mu = mu, size = theta)
    rval[stats::runif(n) < pi] <- 0
    return(rval)
}







#' Approximate p-value by Resampling Integer Compositions
#'
#' Given the value of the test statistic \eqn{t}, the sample sizes \eqn{n} and
#' \eqn{k},power exponent \eqn{p} and vector of weights that together determine
#' the test statistic (by default \eqn{n\geqslant k}), as well as the
#' user-specified resampling number (by default this is \eqn{5000}),
#' performs resampling from the collection of integer compositions
#' to approximate the p-value of the observed test statistic.
#'
#' For \eqn{n} and \eqn{k} small enough (\eqn{n\leqslant 40, k\leqslant 10}),
#' computes exact p-value by computing the test statistic across all
#' \eqn{k}-compositions of \eqn{n} under the uniform distribution.
#'
#' The function returns a two-sided p-value by default, which is more
#' conservative. Users can choose other p-values corresponding to different
#' alternatives; see documentation on `alternative`. Note that the
#' interpretation of the choice of `alternative` depends on the choice of
#' weight vector. For example, a weight vector that is a quadratic kernel will
#' upweight the extreme components of the weight vector. For this choice,
#' setting `alternative` to be `greater` translates into an alternative
#' hypothesis of a bigger spread in the larger sample
#' (the one with sample size \eqn{n}).
#'
#' Dependencies: arrangements::compositions
#' @param t Value of test statistic \eqn{||S_{n,k}(D)/n||_{p,\boldsymbol{w}}^p}
#' computed from data \eqn{D}
#' @param n Sample size of \eqn{y}
#' @param k Sample size of \eqn{x}
#' @param p Power exponent of test statistic
#' @param wList Weight vector
#' @param alternative Character string that should be one of "`two.sided`"
#' (default), "`greater`" or "`less`"
#' @param type If using resampling approximation, either an unbiased estimate
#' of ("`unbiased`", default), or valid, but biased estimate of, ("`valid`")
#' p-value (see Hemerik and Goeman, 2018), or both ("`both`").
#' Default is "`unbiased`".
#' @param resamp_number Number of compositions of \eqn{n} to draw
#' (default is 5000)
#' @return p-value (scalar)
#' @examples
#'
#' getCompositionPValue(
#'     t = 0.5, n = 50, k = 11, p = 1, wList = (10:0) / 10,
#'     alternative = "two.sided", type = "unbiased", resamp_number = 5000
#' )
#'
#' @export
getCompositionPValue <- function(t, n, k, p, wList,
                                 alternative, type, resamp_number) {
    if (n <= 40 & k <= 10) {
        message(
            date(),
            ": n and k are small enough, computing exact p-value..."
        )
        exact_ts <- ((arrangements::compositions(n = n, k = k) / n)^p %*%
            wList) %>% as.vector()
        return(.compute_exact_pval(t, exact_ts, alternative))
    } else {
        resampled_ts <- ((
            arrangements::compositions(
                n = n, k = k,
                nsample = resamp_number
            ) / n
        )^p %*% wList) %>% as.vector()
        return(
            .compute_resample_pval(t, resampled_ts, alternative, type)
        )
    }
}

#' Retrieve indices of \eqn{x_i}'s after merging \eqn{\boldsymbol{x}} and
#' \eqn{\boldsymbol{y}} in ascending order.
#'
#' Given data consisting of either a single sample
#' \eqn{\boldsymbol{x}=(x_1,\ldots,x_k)},
#' or two samples \eqn{\boldsymbol{x}=(x_1,\ldots,x_k)} and
#' \eqn{\boldsymbol{y}=(y_1,\ldots,y_n)},
#' this function obtains the indices of \eqn{x_i}'s after merging
#' \eqn{\boldsymbol{x}} and \eqn{\boldsymbol{y}} in ascending order.
#'
#' Dependencies: None
#' @param x First sample
#' @param y Second sample
#' @param ties.break Whether to break the ties when ordering `x` and `y`.
#' Default is `'TRUE'`.
#' @return Ranks of \eqn{x_i}'s after merging `x` and `y` in ascending order
#' @examples
#'
#' set.seed(1)
#' rank_x(x = abs(rnorm(10)))
#' rank_x(x = abs(rnorm(10)), y = abs(rnorm(100)))
#' rank_x(
#'     x = rnbinom(10, size = 5, prob = 0.3),
#'     y = rnbinom(20, size = 2, prob = 0.3), ties.break = TRUE
#' )
#' @export
#'
rank_x <- function(x, y = NULL, ties.break = TRUE) {
    k <- length(x)
    n <- length(y)
    if (ties.break) {
        x_y_sorted <- sort(c(x, y))
        differences <- diff(x_y_sorted)
        # Find the minimum difference that greater than 0
        minimal_gap <- min(differences[differences > 0])
        # Perturb x and y to avoid ties
        x <- runif(k, 0, minimal_gap / 2) + x
        y <- runif(n, 0, minimal_gap / 2) + y
    }
    x_ordered <- sort(x)
    x <- sort(x)
    x_rank <- c()
    # Obtain the rank of x_i after merging x and y in ascending order
    for (i in seq_len(k)) {
        x_rank[i] <- sum(c(x, y) < x_ordered[i]) + 1
    }
    return(x_rank)
}

#' Compute weights for mean shift analysis
#'
#' Given the value of the estimated parameters for zero-inflated negative
#' binomial distribution (\eqn{\beta}, \eqn{\mu}, \eqn{\pi}),the sample sizes
#' \eqn{n} and \eqn{k} (by default \eqn{n\geqslant k}), this function computes
#' weights that target changes in the mean parameter of the zero-inflated
#' negative binomial model.
#'
#' The function generates a vector of weights and an estimated variance that
#' are used to construct test statistics. The weight vector is computed by
#' considering how each potential count value from the combined distribution
#' contributes to the expected shift in mean, considering both inflated and
#' uninflated components.
#'
#' @param beta Dispersion parameter of the negative binomial component. Must be
#' strictly positive.
#' @param mu Mean of the uninflated negative binomial distribution. Should be
#' non-negative.
#' @param pi Probability of zeros due to zero inflation.
#' @param n Sample size of \eqn{y}
#' @param k Sample size of \eqn{x}
#' @return A list containing a weight vector (`$weight`) and
#' estimated variance (`$var`) useful for constructing test statistics.
#' @examples
#'
#' computeweight_mean(beta = 5, mu = 40, pi = 0.1, n = 100, k = 60)
#'
#' @export
computeweight_mean <- function(beta, mu, pi, n, k) {
    p <- beta / (beta + mu)
    Gdict <- .qzinbinom(seq_len(n + k) / (n + k + 1),
        size = beta,
        mu = mu,
        pi = pi
    )
    z_function <- function(x) {
        if (x == 0) {
            mu * (pi - 1) * (beta / (beta + mu))^(beta + 1) /
                (pi + (1 - pi) * (beta / (beta + mu))^beta)
        } else {
            (beta * x / (beta + mu) - beta * mu / (beta + mu)) * (1 - pi)
        }
    }
    term1 <- beta * mu / (beta + mu)
    term2 <- (beta / (beta + mu))^beta
    term3 <- (1 - pi)^2
    term4 <- beta^2 * term3 * (beta / (beta + mu))^(2 * beta + 2)
    term5 <- pi + (1 - pi) * term2

    H_value <- (term1 - term1^2 * term2) * term3 + term4 / term5

    test_var <- H_value # Put alpha*(1-alpha)* in QRscore_ZINB
    results <- list()
    zvals <- vapply(Gdict, z_function, numeric(1))
    results$weight <- zvals - mean(zvals)
    if (any(is.na(test_var)) ||
        any(is.infinite(test_var)) || all(test_var == 0)) {
        results$var <- var(zvals)
        message("Using weight dictionary to calculate variance.")
    } else {
        results$var <- test_var
    }
    results$var <- test_var
    return(results)
}




#' Compute weights for dispersion shift analysis
#'
#' Given the value of the estimated parameters for zero-inflated negative
#' binomial distribution (\eqn{\beta}, \eqn{\mu}, \eqn{\pi}), the sample sizes
#' \eqn{n} and \eqn{k} (by default \eqn{n\geqslant k}), this function computes
#' weights that target on changes of the dispersion parameter in zero-inflated
#' negative binomial model.
#'
#' The function returns a list containing a weight vector and estimated
#' variance for constructing test statistics.
#'
#' @param beta Target for number of successful trials, or dispersion parameter.
#' Must be strictly positive.
#' @param mu Non-negative mean of the uninflated negative binomial
#' distribution.
#' @param pi Zero inflation probability for structural zeros.
#' @param n Sample size of \eqn{y}
#' @param k Sample size of \eqn{x}
#' @param tail Distribution tail trimmed for numerically calculate expectation.
#' Default is \eqn{10^{-4}}
#' @param bigN Sampling numbers for numerically calculate expectation, this
#' will only be applied when mean parameter is extremely large.
#' Default is \eqn{10^6}
#' @return A list containing a weight vector (`$weight`) and estimated variance
#' (`$var`) useful for constructing test statistics.
#' @examples
#'
#' set.seed(1)
#' computeweight_disp(
#'     beta = 5, mu = 40, pi = 0.1, n = 100, k = 60,
#'     tail = 10^(-4), bigN = 10^6
#' )
#'
#' @export
#'
computeweight_disp <- function(beta, mu, pi, n, k,
                               tail = 10^(-4), bigN = 10^6) {
    p <- beta / (beta + mu)
    Gdict <- .qzinbinom(seq_len(n + k) / (n + k + 1),
        size = beta, mu = mu, pi = pi
    )

    if (mu > 10^4) {
        ## This is for extremely large mus to avoid large number issue
        data <- rzinbinom(bigN, size = beta, mu = mu, pi = pi)
        x_prob <- table(data) / bigN
        x_values <- as.numeric(names(x_prob))
    } else {
        ## For smaller mu's, calculate the expectation without sampling
        tail.q <- .qzinbinom(1 - tail, size = beta, mu = mu, pi = pi)
        x_prob <- .dzinbinom(0:tail.q, size = beta, mu = mu, pi = pi)
        x_values <- 0:tail.q
    }

    z_x_values <- vapply(x_values, function(x) {
        .z_function_disp(x, beta, mu, pi)
    }, numeric(1))

    E_z <- sum(z_x_values * x_prob)
    H_value <- sum((z_x_values - E_z)^2 * x_prob)
    test_var <- H_value

    zvals <- vapply(Gdict, function(x) {
        .z_function_disp(x, beta, mu, pi)
    }, numeric(1))

    results <- list()
    # Adjust the weights to mean 0.
    results$weight <- zvals - mean(zvals)
    assertthat::assert_that(sum(is.na(results$weight)) == 0,
        msg = "Weight cannot be calculated.
                          Check if estimated dispersion parameter beta is
                          less than 100."
    )

    if (any(is.na(test_var)) ||
        any(is.infinite(test_var)) || all(test_var == 0)) {
        results$var <- var(zvals)
        message("Using weight dictionary to calculate variance.")
    } else {
        results$var <- test_var
    }
    return(results)
}

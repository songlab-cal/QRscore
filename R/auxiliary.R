#' @keywords internal
dzinbinom <- function(x, mu, theta, size, pi, log = FALSE) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- log(1 - pi) + dnbinom(x, mu = mu, size = theta, log = TRUE)
  if(any(x0 <- (x == 0L))) rval[x0] <- log(exp(rval) + pi)[x0]
  if(log) rval else exp(rval)
}

#' @keywords internal
qzinbinom <- function(p, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  p <- pmax(0, (p - pi)/(1 - pi))
  rval <- qnbinom(p, mu = mu, size = theta, lower.tail = TRUE, log.p = FALSE)
  rval
}

#' Generate random samples from a zero-inflated negative binomial distribution
#'
#' @param n Integer. Number of samples to generate.
#' @param mu Numeric. Mean of the negative binomial distribution.
#' @param theta Numeric. Dispersion parameter (size) of the negative binomial distribution.
#' @param size Numeric. Alternative name for the dispersion parameter (used interchangeably with 'theta').
#' @param pi Numeric. Zero-inflation probability; must be in the range [0, 1].
#' @return A vector of random samples from a zero-inflated negative binomial distribution.
#' @export
rzinbinom <- function(n, mu, theta, size, pi) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- rnbinom(n, mu = mu, size = theta)
  rval[runif(n) < pi] <- 0
  rval
}



#' Approximate p-value by Resampling Integer Compositions
#'
#' Given the value of the test statistic \eqn{t}, the sample sizes \eqn{n} and \eqn{k},
#' power exponent \eqn{p} and vector of weights that together determine the test statistic
#' (by default \eqn{n\geqslant k}), as well as the user-specified resampling number
#' (by default this is \eqn{5000}), performs resampling from the collection of integer compositions
#' to approximate the p-value of the observed test statistic.
#'
#' For \eqn{n} and \eqn{k} small enough (\eqn{n\leqslant 40, k\leqslant 10}), computes exact p-value by
#' computing the test statistic across all \eqn{k}-compositions of \eqn{n} under the uniform distribution.
#'
#' The function returns a two-sided p-value by default, which is more conservative. Users can
#' choose other p-values corresponding to different alternatives; see documentation on `alternative`.
#' Note that the interpretation of the choice of `alternative` depends on the choice of weight vector.
#' For example, a weight vector that is a quadratic kernel will upweight the extreme components of
#' the weight vector. For this choice, setting `alternative` to be `greater` translates into an alternative
#' hypothesis of a bigger spread in the larger sample (the one with sample size \eqn{n}).
#'
#' Dependencies: arrangements::compositions
#' @param t Value of test statistic \eqn{||S_{n,k}(D)/n||_{p,\boldsymbol{w}}^p} computed from data \eqn{D}
#' @param n Sample size of \eqn{y}
#' @param k Sample size of \eqn{x}
#' @param p Power exponent of test statistic
#' @param wList Weight vector
#' @param alternative Character string that should be one of "`two.sided`" (default), "`greater`" or "`less`"
#' @param type If using resampling approximation, either an unbiased estimate of (`'unbiased'`, default),
#' or valid, but biased estimate of, (`'valid'`) p-value (see Hemerik and Goeman, 2018), or both (`'both'`). Default is `'unbiased'`.
#' @param resamp_number Number of compositions of \eqn{n} to draw (default is 5000)
#' @return p-value (scalar)
#' @examples
#'
#' getCompositionPValue(t = 0.5, n = 50, k = 11, p = 1, wList = (10:0)/10, alternative = "two.sided", type = "unbiased", resamp_number = 5000)
#'
#' @export
getCompositionPValue <- function(t, n, k, p, wList, alternative, type, resamp_number) {
  # if n and k are small, compute exact probability by enumerating all k-compositions of n
  if (n <= 40 & k <= 10) {
    message(date(), ": n and k are small enough, computing exact p-value...")
    exact_ts <- ((arrangements::compositions(n = n, k = k)/n)^p %*% wList) %>% as.vector()

    if (alternative == "two.sided") {
      message(date(), ": Computing exact two-sided p-value")
      upper_tail <- mean(exact_ts >= t)
      cdf_at_t <- mean(exact_ts <= t)
      return(2*(min(cdf_at_t, upper_tail)))
    } else if (alternative == "greater") {
      message(date(), ": Computing exact one-sided p-value with alternative set to `greater`")
      return(mean(exact_ts >= t))
    } else {
      message(date(), ": Computing exact one-sided p-value with alternative set to `less`")
      return(mean(exact_ts <= t))
    }
  } else {
    # otherwise, sample test statistic and compute empirical CDF at t
    resampled_ts <- ((arrangements::compositions(n = n, k = k, nsample = resamp_number)/n)^p %*% wList) %>% as.vector()
    cdf_at_t <- mean(resampled_ts < t)
    cdf_at_t_upp_tail <- 1 - mean(c(resampled_ts,t) >= t)
    cdf_at_t_low_tail <- mean(c(resampled_ts,t) <= t)

    if (alternative == "two.sided") {
      message(date(), ": Computing two-sided p-value")
      if (type == "unbiased") {
        return(2*min(cdf_at_t, 1 - cdf_at_t))
      } else if (type == "valid") {
        return(2*min(cdf_at_t_low_tail, 1 - cdf_at_t_upp_tail))
      } else {
        unbiased <- 2*min(cdf_at_t, 1 - cdf_at_t)
        valid <- 2*min(cdf_at_t_low_tail, 1 - cdf_at_t_upp_tail)
        to_return <- c(unbiased,valid)
        names(to_return) <- c("unbiased","valid")
        return(to_return)
      }
    } else if (alternative == "greater") {
      message(date(), ": Computing one-sided p-value with alternative set to `greater`")
      if (type == "unbiased") {
        return(1 - cdf_at_t)
      } else if (type == "valid") {
        return(1 - cdf_at_t_upp_tail)
      } else {
        unbiased <- 1 - cdf_at_t
        valid <- 1 - cdf_at_t_upp_tail
        to_return <- c(unbiased,valid)
        names(to_return) <- c("unbiased","valid")
        return(to_return)
      }
    } else {
      message(date(), ": Computing one-sided p-value with alternative set to `less`")
      if (type == "unbiased") {
        return(cdf_at_t)
      } else if (type == "valid") {
        return(cdf_at_t_low_tail)
      } else {
        unbiased <- cdf_at_t
        valid <- cdf_at_t_low_tail
        to_return <- c(unbiased,valid)
        names(to_return) <- c("unbiased","valid")
        return(to_return)
      }
    }
  }
}


#' Retrieve the indices of \eqn{x_i}s after merging x and y in ascending order.
#'
#' Given data consisting of either a single sample \eqn{\boldsymbol{x}=(x_1,\ldots,x_k)},
#' or two samples \eqn{\boldsymbol{x}=(x_1,\ldots,x_k)} and \eqn{\boldsymbol{y}=(y_1,\ldots,y_n)},
#' this function obtains the indices of \eqn{x_i}s after merging x and y in ascending order.
#'
#' Dependencies: None
#' @param x First sample
#' @param y Second sample
#' @param ties.break Whether to break the ties when ordering x and y. Default is `'TRUE'`.
#' @param seed Random seed for tie breaking.
#' @return Ranks of \eqn{x_i}s after merging x and y in ascending order
#' @examples
#'
#' rank_x (x = abs(rnorm(10)))
#' rank_x (x = abs(rnorm(10)), y = abs(rnorm(100))) 
#' rank_x (x = rnbinom(10, size = 5, prob = 0.3), y = rnbinom(20, size = 2, prob = 0.3), ties.break = TRUE)
#' @export
#' 
rank_x = function(x,
                  y = NULL,
                  ties.break = TRUE,
                  seed = 1) {
  k = length(x)
  n = length(y)
  if (ties.break) {
    set.seed(seed)
    # Calculate the differences between adjacent elements after merging the vector
    x_y_sorted = sort(c(x, y))
    differences <- diff(x_y_sorted)
    # Find the minimum difference that greater than 0
    minimal_gap <- min(differences[differences > 0])
    # Perturb x and y to avoid ties
    x = runif(k, 0, minimal_gap / 2) + x
    y = runif(n, 0, minimal_gap / 2) + y
  }
  x_ordered <- sort(x)
  x = sort(x)
  x_rank <- c()
  # Obtain the rank of x_i after merging x and y in ascending order
  for (i in 1:k) {
    x_rank[i] =  sum(c(x, y) < x_ordered[i]) + 1
  }
  return(x_rank)
}

#' Compute weights for mean shift analysis
#'
#' Given the value of the estimated parameters for zero-inflated negative binomial distribution (\eqn{beta}, \eqn{mu}, \eqn{pi}),
#' the sample sizes \eqn{n} and \eqn{k} (by default \eqn{n\geqslant k}), this function computes weights that target changes in the 
#' mean parameter of the zero-inflated negative binomial model.
#'
#' The function generates a vector of weights and an estimated variance that are used to construct test statistics. 
#' The weight vector is computed by considering how each potential count value from the combined distribution contributes 
#' to the expected shift in mean, considering both inflated and uninflated components.
#'
#' @param beta Dispersion parameter of the negative binomial component. Must be strictly positive.
#' @param mu Mean of the uninflated negative binomial distribution. Should be non-negative.
#' @param pi Probability of zeros due to zero inflation.
#' @param n Sample size of \eqn{y}
#' @param k Sample size of \eqn{x}
#' @return A list containing a weight vector 'weight' and estimated variance 'var' useful for constructing test statistics.
#' @examples
#'
#' computeweight_mean(beta = 5, mu = 40, pi = 0.1, n = 100, k = 60)
#'
#' @export
computeweight_mean = function(beta,mu,pi, n,k){
  p = beta/(beta + mu)
  Gdict = qzinbinom((1:(n+k))/(n+k+1), size = beta, mu = mu,pi = pi) 
  z_function = function(x){
    if(x ==0){
      mu*(pi-1)*(beta/(beta+mu))^(beta+1)/(pi+(1-pi)*(beta/(beta+mu))^beta)
    }else{
      (beta*x/(beta+mu) - beta*mu/(beta+mu))*(1-pi)
    }
    
  }
  H_value = (beta*mu/(beta+mu)-(beta*mu/(beta+mu))^2*(beta/(beta+mu))^beta)*(1-pi)^2 + beta^2*(1-pi)^2*(beta/(beta+mu))^(2*beta+2)/((pi+(1-pi)*(beta/(beta+mu))^beta))
  test_var = H_value # Put alpha*(1-alpha)* in QRscore_ZINB
  results = list()
  results$weight = sapply(Gdict,z_function) - mean(sapply(Gdict, z_function)) # weight vector is a dictionary for all possible ranks
  if (any(is.na(test_var)) || any(is.infinite(test_var)) || all(test_var == 0)) {
    results$var <- var(sapply(Gdict, z_function))
    message(paste0("Using weight dictionary to calcualte variance."))
  } else{
    results$var = test_var
  }
  results$var = test_var
  return(results)
}




#' Compute weights for dispersion shift analysis
#'
#' Given the value of the estimated parameters for zero-inflated negative binomial distribution (\eqn{beta}, \eqn{mu}, \eqn{pi}),
#' the sample sizes \eqn{n} and \eqn{k} (by default \eqn{n\geqslant k}), this function computes weights that target on changes of the 
#' dispersion parameter in zero-inflated negative binomial model.
#' 
#' The function returns a list containing a weight vector 'weight' and estimated variance 'var' for constructing test statistics.
#'
#' @param beta Target for number of successful trials, or dispersion parameter. Must be strictly positive.
#' @param mu Non-negative mean of the uninflated negative binomial distribution.
#' @param pi Zero inflation probability for structural zeros.
#' @param n Sample size of \eqn{y}
#' @param k Sample size of \eqn{x}
#' @param tail Distribution tail trimmed for numerically calculate expectation. Default is \eqn{10^(-4)}
#' @param bigN Sampling numbers for numerically calculate expectation, this will only be applied when mean parameter is extremely large. Default is \eqn{10^6}
#' @param seed Random seed for sampling.
#' @return A list containing a weight vector 'weight' and estimated variance 'var' useful for constructing test statistics.
#' @examples
#'
#' computeweight_disp(beta = 5, mu = 40, pi = 0.1, n = 100, k = 60, tail = 10 ^ (-4), bigN = 10 ^ 6)
#'
#' @export
#' 
computeweight_disp = function(beta,
                             mu,
                             pi,
                             n,
                             k,
                             tail = 10 ^ (-4),
                             bigN = 10 ^ 6,
                             seed = 1) {
  p = beta / (beta + mu)
  Gdict = qzinbinom((1:(n + k)) / (n + k+1),
                              size = beta,
                              mu = mu,
                              pi = pi)
  if (mu > 10 ^ 4) {
    ## This is for extremely large mus to avoid large number issue
    set.seed(seed)
    data = rzinbinom(bigN,
                               size = beta,
                               mu = mu,
                               pi = pi)
    x_prob = table(data) / bigN
    x_values = as.numeric(names(x_prob))
  } else{
    ## For smaller mus, we can calculate numerically calculate the expectation without sampling
    tail.quantile = qzinbinom(1 - tail,
                                        size = beta,
                                        mu = mu,
                                        pi = pi)
    x_prob = dzinbinom(0:tail.quantile,
                                 size = beta,
                                 mu = mu,
                                 pi = pi)
    x_values = 0:tail.quantile
  }
  
  # The original form of this function:
  # Gamma_derivative = function(theta,x){
  #   Gamma_expr = expression(log(gamma(theta*beta+x)/gamma(theta*beta)))
  #   eval(D(Gamma_expr, 'theta'))
  # }
  # Simlified version of the derivative is shown below, which help to avoid large number issue.

  Gamma_derivative = function(theta, x) {
    x4 = gamma(theta * beta)
    x1 = beta  * digamma(theta * beta + x)
    x2 =  (beta * (x4 * digamma(theta * beta))) / x4
    return(x1 - x2)
  }
  
  z_function = function(x) {
    if (x == 0) {
      (1 - pi) * (beta / (beta + mu)) ^ (beta + 1) * ((beta + mu) * log(beta /
                                                                          (beta + mu)) + mu) / ((1 - pi) * (beta / (beta + mu)) ^ beta + pi)
    } else{
      (Gamma_derivative(1, x) + beta * (mu - x) / (beta + mu) + beta * log(beta /
                                                                             (beta + mu))) * (1 - pi)
    }
  }
  z_x_values = sapply(x_values, z_function)
  E_z = sum(z_x_values * x_prob)
  H_value = sum((z_x_values - E_z) ^ 2 * x_prob)
  test_var = H_value # Put alpha*(1-alpha)* in QRscore_ZINB; another: var(sapply(Gdict, z_function))
  results = list()
  # weight vector is a dictionary for all possible ranks, adjust the weights to mean 0.
  results$weight = sapply(Gdict, z_function) - mean(sapply(Gdict, z_function))
  assertthat::assert_that(sum(is.na(results$weight)) == 0,
                          msg = "Weight cannot be calculated. Check if estimated dispersion parameter beta is less than 100.")
  if (any(is.na(test_var)) || any(is.infinite(test_var)) || all(test_var == 0)) {
    results$var <- var(sapply(Gdict, z_function))
    message(paste0("Using weight dictionary to calcualte variance."))
  } else{
    results$var = test_var
  }
  
  return(results)
}



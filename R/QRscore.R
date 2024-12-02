#' Flexible Non-Parametric One- and Two-Sample Tests (Native R version)
#'
#' Given data consisting of either a single sample \eqn{\boldsymbol{x}=(x_1,\ldots,x_k)},
#' or two samples \eqn{\boldsymbol{x}=(x_1,\ldots,x_k)} and \eqn{\boldsymbol{y}=(y_1,\ldots,y_n)},
#' this function uses summary statistics computed on weighted linear combinations of powers of
#' the spacing statistics \eqn{S_k} (former) or \eqn{S_{n,k}} (latter).
#'
#' More precisely, this function does the following:
#'
#' For a single sample \eqn{x}, the function tests for uniformity of its entries. When \eqn{p=2}
#' and a particular choice of \eqn{\boldsymbol{w}} is specified, we recover Greenwood's test.
#'
#' For two samples, the function tests the null of \eqn{\boldsymbol{x}} and \eqn{\boldsymbol{y}}
#' being drawn from the same distribution (i.e., stochastic equality), against flexible alternatives
#' that correspond to specific choices of the test statistic parameters, \eqn{\boldsymbol{w}} (weight vector)
#' and \eqn{p} (power). These parameters not only determine the test statistic
#' \eqn{||S_k||_{p,\boldsymbol{w}}^p=\sum_{j=1}^k w_iS_{k}[j]^p} (analogously defined for
#' \eqn{||S_{n,k}||_{p,\boldsymbol{w}}^p}), but also encode alternative hypotheses
#' ranging from different populational means (i.e., \eqn{\mu_x \neq \mu_y}), different
#' populational spreads (i.e., \eqn{\sigma^2_x \neq \sigma^2_y}), etc.
#'
#' Additional tuning parameters include (1) choice of p-value computation (one- or two-sided);
#' (2) approximation method (moment-based such as Bernstein, Chebyshev or Jacobi, or resampling-based);
#' (3) number of moments accompanying the approximation chosen if using moment-based approximation
#' (recommended 200, typically at least 100); and (4) in case of two samples,
#' whether the user prefers to use exact discrete moments (more accurate but slower) or to use
#' continuous approximations of the discrete moments (less accurate but faster).
#'
#' (4/21/22) Currently, only resampling and Gaussian asymptotics are supported. Both are efficient and well-calibrated.
#'
#' (4/14/22) Currently, for \eqn{n\geqslant 100} and \eqn{k\geqslant 50} such that \eqn{\frac{k}{n}\geqslant 0.001}, function
#' automatically uses Gaussian approximation to the null.
#'
#' Dependencies: functions in `auxiliary.R`
#' @param x First sample
#' @param y Second sample
#' @param p Exponent value in defining test statistic (must be integer)
#' @param wList Vector of weights. It should have length equal to \eqn{x} when \eqn{y} is `NULL`,
#' and one more than the length of \eqn{x} when \eqn{y} is not `NULL`
#' @param alternative How p-value should be computed; i.e., a character specifying the alternative hypothesis,
#' must be one of "`two.sided`", "`greater`" or "`less`"
#' @param approx Which approximation method to use (choose "`resample`", "`asymptotic`")
#' @param type If using resampling approximation, either an unbiased estimate of ("`unbiased`", default),
#' or valid, but biased estimate of, ("`valid`") p-value (see Hemerik and Goeman, 2018), or both ("`both`"). Default is "`unbiased`".
#' @param n_mom The number of moments to accompany the approximation (recommended 200, if not at least 100)
#' @param resamp_number Number of \eqn{k}-compositions of \eqn{n} or simplex vectors in \eqn{[0,1]^k}  to draw
#' @export
#' @examples
#'
#' set.seed(1)
#' # One-sample examples
#' QRscore_Flex(x = abs(rnorm(10)), p = 2, wList = rep(1,10), 
#'               alternative = "two.sided", approx = "resample")
#'
#' # Two-sample examples
#' QRscore_Flex(x = abs(rnorm(30)), y = abs(rnorm(100)), p = 2, 
#'              wList = rep(1,31), alternative = "two.sided", 
#'               approx = "resample", resamp_number = 5000)
#' 
#' QRscore_Flex(x = abs(rnorm(100)), y = abs(rnorm(100)), p = 1, 
#'              wList = 0:100, alternative = "two.sided", 
#'               approx = "asymptotic")
#' @export
#' 
#' 
QRscore_Flex <- function(x, y = NULL,
                        p = 1,
                        wList,
                        alternative = "two.sided",
                        approx = "resample",
                        type = "unbiased",
                        n_mom = NULL,
                        resamp_number = 5000
                        ) {
  # 1. Get number of bins
  if (!is.null(y)) {
    k <- length(x) + 1
    assertthat::assert_that(length(wList)==k, msg = "Length of wList must be length(x) + 1.")
  } else {
    k <- length(x)
  }
  # 2. Assert conditions
  assertthat::assert_that(p == 1 | p == 2, 
                          msg = "Currently only p = 1 or p = 2 is supported.")
  assertthat::assert_that(approx == "resample" | approx == "asymptotic", 
                          msg= "Currently only approx = 'resample' or approx = 'asymptotic' is supported")
  # 3. Normalize weights
  message(date(), ": Normalizing weight vector...")
  wList <- wList / max(wList)
  
  # 4. Compute test statistic t and return its p-value
  # 4.1. Case 1: y is not NULL
  if (!is.null(y)) {

    # construct ordering of x_i's
    x_ordered <- x[order(x)]
    x_ordered <- c(-Inf, x_ordered, Inf)
    n <- length(y) # get sample size / number of balls

    # construct Snk
    Snk <- c()
    for (i in 1:k) {
      # count number of y_j's between x_i and x_{i+1}
      Snk <- c(Snk, sum(y >= x_ordered[i] & y < x_ordered[i+1]))
    }
    
    # construct t
    if (approx == "resample" & !(n >= 100 & k >= 50)) {
      t <- sum(((Snk + 1) / (n + k))^p * wList)
      message(date(), ": Adjusting n and k for resampling")
    } else {
      t <- sum((Snk / n)^p * wList)
    }
    message(date(), ": The test statistic for the data is ", t)
    assertthat::assert_that(k/n >= 1e-3, msg = "Condition k/n >= 1e-3 must be satisfied.")
    # decide on an approximation:
    if (approx == "asymptotic" | (n >= 100 & k >= 50)) {
      if (!(n >= 100 & k >= 50)) {
        warning("Sample sizes may be too small for Gaussian asymptotics to be accurate. Recommend to use approx = \"resample\" instead.")
      }
      message(
        date(),
        ": Applying Gaussian asymptotics...")
      # compute analytical mean and variance
      if (p == 1) {
        # p = 1
        first_moment <- sum(wList) / k
        second_moment <- ((k/n+1)/(k^2*(k+1))) *
          sum(wList *((k*diag(k) - outer(rep(1,k),rep(1,k))) %*% wList))
      } else {
        # p = 2
        first_moment <- ((2+k/n-1/n) / ((k+1)*k)) * sum(wList)
        sum_of_wj2s <- sum(wList^2)
        coeff_sum_of_wj2s <- (k-1) * (k/n+1) * (2+k/n-1/n) * (12-6/n+k*(k/n+10-5/n)) / (k^2*(1+k)^2*(2+k)*(3+k))
        offdiag_sum <- sum(outer(wList, wList)) - sum(wList^2)
        coeff_offdiag_sum <- (k/n+1) * (6/n^2+k*(3+k*(k-2))/n^2-24/n+8*(k-1)*k/n+8*(3+2*k)) / (k^2*(1+k)^2*(2+k)*(3+k))
        second_moment <- sum_of_wj2s * coeff_sum_of_wj2s - offdiag_sum * coeff_offdiag_sum
      }
      
      z_score <- (t - first_moment) / sqrt(second_moment)
      
      if (alternative == "two.sided") {
        return(2 * min(pnorm(z_score), pnorm(-z_score)))
      } else if (alternative == "greater") {
        return(pnorm(-z_score))
      } else {
        return(pnorm(z_score))
      }
      
    } else {
      message(date(),
              ": Using resampling approach, with resampling number ", resamp_number, ", to approximate p-value...")
      return(getCompositionPValue(t = t, n = n + k, k = k, p = p,
                                  wList = wList,
                                  alternative = alternative,
                                  type = type,
                                  resamp_number = resamp_number))
      
    } 
    
  } else {
    # 4.2. Case 2: y is NULL => use continuous moments
    assertthat::assert_that(approx == "resample", msg = "Currently only approx = \"resample\" is supported for one-sample test.")
    assertthat::assert_that(!any(x < 0),
                            msg = "Observations should be non-negative")
    Sk <- x / sum(x)
    
    # construct t
    t <- sum(Sk^p * wList)
    message(date(), ": The test statistic for the data is ", t)
    
    # decide on approximation
    
    message(date(),
            ": Using resampling approach on continuous simplex, with resampling number ",
            resamp_number, ", to approximate p-value...")
    resampled_ts <- as.vector((simplex.sample(n = k, N = resamp_number)$samples)^p %*% wList)
    
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
        unbiased <- 2 * min(cdf_at_t, 1 - cdf_at_t)
        valid <- 2 * min(cdf_at_t_low_tail, 1 - cdf_at_t_upp_tail)
        to_return <- c(unbiased, valid)
        names(to_return) <- c("unbiased", "valid")
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
        to_return <- c(unbiased, valid)
        names(to_return) <- c("unbiased", "valid")
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
        to_return <- c(unbiased, valid)
        names(to_return) <- c("unbiased", "valid")
        return(to_return)
      }
    }
    
  }
}


#' Non-Parametric Two-Sample Tests Designed for Testing Differences in Mean or Dispersion Parameters in 
#' (Zero-Inflated) Negative Binomial Distributions.
#'
#' This function evaluates the null hypothesis that two samples, \eqn{\boldsymbol{x}} and \eqn{\boldsymbol{y}},
#' are drawn from the same distribution, specifically designed for NB or ZINB models. It is particularly effective
#' in detecting shifts in either the mean or the dispersion parameters.
#'
#' The function automatically computes optimal weights for the chosen model and derives a p-value based on
#' the selected test statistic and approximation method.
#'
#' Additional tuning parameters include:
#' (1) whether to use a likelihood ratio test to determine which model (NB or ZINB) to fit,
#' (2) the approximation method (default is resampling, with asymptotic estimation for large samples),
#' (3) choice of p-value computation (one- or two-sided),
#' (4) threshold for estimated proportion of zeros in ZINB model (returns NA if exceeded).
#'
#' Dependencies: pscl::zeroinfl, MASS::glm.nb, and auxiliary functions from `auxiliary.R`
#' @param x First sample
#' @param y Second sample
#' @param zero_inflation If TRUE, automatically chooses between ZINB and NB models based on the data; if FALSE, applies NB model estimation.
#' @param LR.test Whether to use a likelihood ratio test to determine which model (NB or ZINB) to fit
#' @param approx Which approximation method to use (default `resample`)
#' @param alternative How p-value should be computed; must be one of "`two.sided`", "`greater`" or "`less`". 
#' @param resamp_num Number of \eqn{k}-compositions of \eqn{n} or simplex vectors in \eqn{[0,1]^k} to draw
#' @param pi_threshold Threshold for estimated proportion of zeros in ZINB model.
#' @param gene.name Optional, name of the gene if applicable, used for customized messages.
#' @param measure Specifies whether to test for shifts in "`mean`" or "`dispersion`".
#' @param p_value If TRUE, returns a p-value, else returns test statistics and weights.
#' @param seed Random seed for sampling, ensures reproducibility in resampling methods.
#' @return p-value or test statistics depending on `p_value` parameter
#' @export
#' @examples
#'
#' # Two-sample example comparing mean shifts
#' QRscore_ZINB(x = rzinbinom(100, size = 2, mu = 20, pi = 0),
#'            y = rzinbinom(100, size = 2, mu = 30, pi = 0), 
#'            zero_inflation = FALSE, LR.test = FALSE, alternative = "greater", 
#'            approx = "asymptotic", measure = "mean")
#'
#' # Two-sample example comparing dispersion shifts
#' QRscore_ZINB(x = rzinbinom(100, size = 2, mu = 20, pi = 0.1),
#'            y = rzinbinom(100, size = 1, mu = 20, pi = 0.1), 
#'            zero_inflation = TRUE, LR.test = TRUE, alternative = "two.sided", 
#'            approx = "asymptotic", measure = "dispersion")
#'
#' # Two-sample example with significant zero inflation and variance shift
#' QRscore_ZINB(x = rzinbinom(30, size = 4, mu = 20, pi = 0.1),
#'            y = rzinbinom(30, size = 1, mu = 20, pi = 0.3), 
#'            zero_inflation = TRUE, LR.test = FALSE, alternative = "two.sided", 
#'            approx = "resample", resamp_num = 50000, measure = "dispersion")
#' 
#' @export
#' 
QRscore_ZINB <- function(x, y, zero_inflation = TRUE, LR.test = FALSE, approx = "resample",
                         alternative = "two.sided", resamp_num = 20000, pi_threshold = 0.95, gene.name = NULL,
                         measure = "mean", p_value = TRUE, seed = 1) {
  
  if (sum(round(x) != x) > 0 | sum(round(y) != y) > 0) {
    warning("Input values are not integers. Automatically rounding the input.")
    x <- round(x)
    y <- round(y)
  }
  
  assertthat::assert_that(!is.null(y), msg = "Please provide a non-null value for y.")
  assertthat::assert_that(all(x >= 0) & all(y >= 0), msg = "Input contains negative values, which are not allowed.")
  
  combined <- c(x, y)
  n <- length(y)
  k <- length(x)
  alpha = k/(n + k)
  
  # Decide which weight function to use based on the 'measure' parameter
  computeweight <- ifelse(measure == "dispersion", computeweight_disp, computeweight_mean)
  
  if (sum(combined == 0) == 0 || !zero_inflation) {
    message(date(), ": Estimating negative binomial parameters.")
    fitNB = MASS::glm.nb(combined ~ 1)
    beta = min(fitNB$theta, 100)  # Set a lower bound for dispersion estimation.
    mu = exp(fitNB$coefficients["(Intercept)"])
    pi = 0 
  } else {
    message(date(), ": Estimating zero inflated negative binomial parameters.")
    fitZINB <- pscl::zeroinfl(combined ~ 1 | 1, dist = "negbin")
    beta = min(fitZINB$theta, 100)
    mu = exp(fitZINB$coef$count[1])
    pi = exp(fitZINB$coef$zero[1])/(1+exp(fitZINB$coef$zero[1]))
    
    if (LR.test) {
      message(date(), ": Applying likelihood ratio test to choose model.")
      fitNB <- MASS::glm.nb(combined ~ 1)
      lrt_pvalue =  as.numeric(pchisq(2 * (logLik(fitZINB) - logLik(fitNB)), df = 1, lower.tail = FALSE))
      
      # If the LRT supports the NB model, update the parameters
      if (lrt_pvalue >= 0.05) {
        beta = min(fitNB$theta, 100)
        mu = exp(fitNB$coefficients["(Intercept)"])
        pi = 0  
      }
    }
    
    if (pi > pi_threshold) {
      if(is.null(gene.name)){
        warning(paste0("The estimated probability of zeros is greater than ", pi_threshold, ", returning NA."))
        
      }else{
        warning(paste0("The estimated probability of zeros for gene ", gene.name," is greater than ", pi_threshold, ", returning NA."))
      }
      return(NA)
    }
  }
  
  results <- computeweight(beta, mu, pi, n, k)
  if (!p_value) {
    return(results)
  }
  
  weights = results$weight
  est_var = alpha*(1-alpha)*results$var
  if (approx == "asymptotic" | (n >= 100 & k >= 50)) {
    zscore = k / sqrt(n + k) * mean(weights[rank_x(x, y)])
    if (alternative == "two.sided") {
      
      #zscore = (abs(n / sqrt(n + k) * mean(weights[rank_x(y, x)]))+abs(zscore))/2
      p.value = 2 * pnorm(abs(zscore), mean = 0, sd = sqrt(est_var), lower.tail = FALSE) 
    } else if (alternative == "greater") { # "greater" in dispersion detection means greater size parameter!!
      p.value = pnorm(zscore, mean = 0, sd = sqrt(est_var), lower.tail = FALSE) 
    } else {
      p.value = pnorm(zscore, mean = 0, sd = sqrt(est_var))
    }
    
  } else{
    message(date(),
            ": Using resampling approach, with resampling number ",
            resamp_num,", to approximate p-value...")
    set.seed(seed)
    zscore = mean(weights[rank_x(x, y)])
    null_zscore = c()
    for (num in 1:resamp_num) {
      null_zscore[num] = mean(weights[sample(1:(n + k), k)])
    }
    p.permute = sum(null_zscore < zscore) / resamp_num
    if (alternative == "two.sided") {
      p.value = 2 * min(p.permute, 1 - p.permute)
    } else if (alternative == "greater") {
      p.value = 1 - p.permute
    } else {
      p.value = p.permute
    }
  }
  
  return(p.value)
}



#' Multi-Sample Nonparametric Test for Mean or Dispersion Differences in 
#' (Zero-Inflated) Negative Binomial Distributions.
#'
#' This function conducts statistical tests across multiple samples to evaluate the null hypothesis
#' that all groups are drawn from the same distribution. It is optimized for data modeled by
#' Negative Binomial (NB) or Zero-Inflated Negative Binomial (ZINB) distributions and is capable
#' of detecting shifts in mean or dispersion parameters. The function can handle any number of groups
#' and automatically computes optimal weights for the specified measure (mean or dispersion).
#'
#' The computation involves constructing a B matrix that transforms group-specific scores into
#' a space where independence among groups is maximized. It then uses these transformed scores
#' to calculate a test statistic, which follows a chi-square distribution under the null hypothesis.
#'
#' Additional tuning parameters allow customization of the model fitting and statistical testing,
#' including:
#' - Selection between NB and ZINB models based on presence of zero inflation.
#' - Choice of approximation method for computing p-values - 'asymptotic' is recommended. 
#' - Decision criteria for statistical tests (one-sided or two-sided).
#' - Threshold for the estimated proportion of zeros beyond which results are considered unreliable.
#'
#' Dependencies: Requires `pscl::zeroinfl` for zero-inflated models, `MASS::glm.nb` for NB models,
#' and other auxiliary functions as needed.
#'
#' @param samples Vector of all sample measurements.
#' @param labels Group labels for each sample.
#' @param zero_inflation Boolean, if TRUE, the function chooses between ZINB and NB models based on data; if FALSE, only NB model is applied.
#' @param LR.test Boolean, if TRUE, performs a likelihood ratio test to select between NB and ZINB models.
#' @param approx The method used for p-value approximation; "`resample`" (default) or "`asymptotic`".
#' @param resamp_num The number of resampling iterations used if `approx` is "`resample`".
#' @param pi_threshold Threshold for proportion of zeros at which to return NA, indicating unreliable results due to excessive zero inflation.
#' @param gene.name Optional, name of the gene if applicable, enhancing the relevance of output in genetic studies.
#' @param measure Specifies whether the test focuses on "`mean`" or "`dispersion`" differences.
#' @param perturb Boolean, if TRUE, adds small noise to data to avoid ties and improve model stability.
#' @param seed Seed for random number generation, ensuring reproducibility.
#' @return Returns the p-value of the test if `p_value` is TRUE, otherwise returns test statistics and weights.
#' @export
#' @examples
#'
#' set.seed(123)
#' data <- c(rnbinom(100, size = 2, mu = 20), rnbinom(100, size = 2, mu = 25), 
#'           rnbinom(100, size = 2, mu = 30))
#' labels <- factor(c(rep('Group1', 100), rep('Group2', 100), rep('Group3', 100)))
#' QRscore_ZINB_nSamples(samples = data, labels = labels,
#'                               zero_inflation = FALSE, LR.test = FALSE, approx = "resample",
#'                               resamp_num = 5000, pi_threshold = 0.95, measure = "mean")
#'
#' # Example with zero inflation and dispersion shift detection
#' data_zi <- c(rzinbinom(100, size = 2, mu = 20, pi = 0.1), 
#'              rzinbinom(100, size = 3, mu = 20, pi = 0.1),
#'              rzinbinom(100, size = 4, mu = 20, pi = 0.1))
#' labels_zi <- factor(c(rep('Group1', 100), rep('Group2', 100), rep('Group3', 100)))
#' QRscore_ZINB_nSamples(samples = data_zi, labels = labels_zi,
#'                                  zero_inflation = TRUE, LR.test = TRUE, approx = "asymptotic",
#'                                  resamp_num = 2000, pi_threshold = 0.95, measure = "dispersion")
#' @export
#' 
QRscore_ZINB_nSamples <- function(samples, labels, zero_inflation = T, LR.test = F, approx = "resample", resamp_num = 20000, pi_threshold = 0.95, gene.name = NULL, measure = "mean", perturb = TRUE, seed = 1) {
  unique_labels <- unique(labels)
  n_groups <- length(unique_labels)
  sample_list <- lapply(unique_labels, function(l) samples[labels == l])
  group_sizes <- sapply(sample_list, length)
  N_all <- sum(group_sizes)
  
  generateBMatrix <- function(d) {
    vectors <- matrix(0, nrow = d-1, ncol = d) # Initialize matrix with zeros
    for (k in 1:(d-1)) {
      vector <- rep(0, d) # Start with a vector of zeros
      vector[1:k] <- 1 # Set the first k elements to 1
      vector[k+1] <- -k # Set the (k+1)th element to -k
      norm_factor <- sqrt(sum(vector^2)) # Calculate the normalization factor
      vectors[k, ] <- vector / norm_factor # Assign the normalized vector to the matrix
    }
    return(t(vectors))
  }
  
  
  B_matrix <- generateBMatrix(n_groups)
  #print(B_matrix)
  
  # Prepare AA_matrix
  group_proportion <- group_sizes / N_all
  AA_matrix <- matrix(NA, n_groups, n_groups)
  for (i in 1:n_groups) {
    for (j in 1:n_groups) {
      if (i == j) {
        AA_matrix[i,j] <- 1 - group_proportion[i]
      } else {
        AA_matrix[i,j] <- -sqrt(group_proportion[i] * group_proportion[j])
      }
    }
  }
  #print(AA_matrix)
  
  # Weight calculation 
  x = sample_list[[1]]
  y = unlist(sample_list[-1], use.names = FALSE)
  results = QRscore_ZINB(x,y,zero_inflation = zero_inflation,LR.test = LR.test, 
                         approx = approx, resamp_num = resamp_num, pi_threshold = pi_threshold, 
                         gene.name = gene.name,measure = measure,p_value = FALSE,seed = seed)
  if(!is.list(results)){
    return(NA)
  }
  # S scores calculation
  S_scores <- numeric(n_groups)
  
  for (i in 1:n_groups) {
    # Extract the current group 'x' and other groups 'y'
    x <- sample_list[[i]]
    y <- unlist(sample_list[-i], use.names = FALSE)
    if(!is.list(results) || is.null(results$weight)) {
      S_scores[i] <- NA  
    } else {
      ranked_indices <- rank_x(x, y)  
      S_scores[i] <- sum(results$weight[ranked_indices]) / sqrt(N_all)
    }
  }
  
  Lambda_matrix <- diag(1/sqrt(group_proportion))
  Test_vec <- t(B_matrix) %*% Lambda_matrix %*% S_scores
  
  # Hessian matrix adjustments
  H_w <- results$var
  H_modifier <- t(B_matrix) %*% AA_matrix %*% B_matrix
  H_modifier <- (H_modifier + t(H_modifier)) / 2
  
  H_new <- H_w * H_modifier
  
  if (any(is.na(H_new)) || any(is.infinite(H_new))) {
    if (is.null(gene.name)) {
      warning("Covariance matrix contains NA or Inf values, returning NA.")
    } else {
      warning(paste0("Covariance matrix contains NA or Inf values for gene ", 
                     gene.name, ", returning NA."))
    }
    return(NA)
  }
  
  if (isSymmetric(H_new) && all(eigen(H_new)$values > 0)) {
    H_inv <- chol2inv(chol(H_new))
  } else {
    H_inv <- tryCatch({
      solve(H_new)
    }, error = function(e) {
      if (is.null(gene.name)) {
        warning("Covariance matrix is not invertible, returning NA.")
      } else {
        warning(paste0("Covariance matrix is not invertible for gene ", 
                       gene.name, ", returning NA."))
      }
      return(NA)
    })
  }
  
  # Chi-square test
  Q_all <- t(Test_vec) %*% H_inv %*% Test_vec
  p.value <- pchisq(Q_all, df = n_groups-1, lower.tail = FALSE)
  
  
  return(p.value[1,1])
}

#' QRscore Test
#'
#' This function performs statistical tests on data from one or more groups using summary statistics
#' computed on weighted linear combinations of powers of spacing statistics. It is capable of conducting
#' one-sample tests, two-sample tests, and multi-sample tests, utilizing either user-defined weights
#' or automatically generated weights based on Negative Binomial (NB) or Zero-Inflated Negative Binomial (ZINB) models.
#'
#' For one-sample tests, the function assesses the uniformity of data entries. For two-sample and multi-sample tests,
#' it evaluates whether groups are drawn from the same distribution, with alternative hypotheses considering
#' differences in means or dispersions.
#'
#' If the weights and \eqn{p} are given, the function calculates the test statistic as:
#' \deqn{||S_k||_{p,\boldsymbol{w}}^p=\sum_{j=1}^k w_iS_{k}[j]^p}
#' where \eqn{w_i} are weights, \eqn{x_i} are data points, and \eqn{p} is the power specified.
#'
#' In two-sample and multi-sample settings without specified weights, the function can automatically estimate weights
#' using score function for a Negative Binomial or a Zero-Inflated Negative Binomial model, optimizing for dispersion or mean shifts.
#'
#' @param samples A numeric vector containing all sample measurements.
#' @param labels An optional vector of group labels corresponding to each entry in `samples`.
#' @param p The exponent used in the power sum test statistic, required if `wList` is not `NULL`.
#' @param wList An optional vector of weights; if `NULL`, weights are estimated using an NB or ZINB model for multiple groups.
#' @param alternative Specifies the alternative hypothesis; must be one of "`two.sided`", "`greater`", or "`less`".
#' @param approx The method used for p-value approximation, either "`resample`" or "`asymptotic`".
#' @param type Specifies if the estimation of the p-value should be "`unbiased`", "`valid`", or both.
#' @param n_mom The number of moments to accompany the approximation, relevant for moment-based methods.
#' @param resamp_number The number of resampling iterations used if `approx` is "`resample`".
#' @param zero_inflation Indicates whether to account for zero inflation in model-based weight estimation.
#' @param LR.test Whether a likelihood ratio test is used to decide between NB and ZINB models.
#' @param pi_threshold Threshold for the proportion of zeros in ZINB models; results in NA if exceeded.
#' @param gene.name Optional identifier for a gene, used in output messages.
#' @param measure Specifies the statistical measure to be analyzed ("`mean`" or "`dispersion`") when weights are auto-generated.
#' @param perturb Boolean to indicate if data should be perturbed slightly to prevent ties.
#' @param use_base_r Boolean to decide whether to use base R functions for certain edge cases like Mann-Whitney tests.
#' @param seed Integer to set the seed for reproducibility in resampling methods.
#' @export
#' @examples
#'
#' set.seed(1)
#' # One-sample test example with normally distributed data
#' data <- abs(rnorm(10))
#' QRscore.test(data, p = 2, wList = rep(1,10), alternative = "two.sided", 
#'               approx = "resample")
#'
#' # Two-sample test with specified weights using normally distributed data
#' group1 <- rnorm(120, sd = 1)
#' group2 <- rnorm(120, sd = 2) # Different mean
#' data <- c(group1, group2)
#' labels <- c(rep("Group1", 120), rep("Group2", 120))
#' QRscore.test(samples = data, labels = labels, p = 1, wList = c(60:0,1:60), 
#'               alternative = "two.sided", approx = "resample")
#'
#' # Two-sample test with automatically estimated weights from NB model
#' group1 <- rzinbinom(120, size = 2, mu = 20, pi = 0)
#' group2 <- rzinbinom(100, size = 2, mu = 30, pi = 0) # Different mean
#' data <- c(group1, group2)
#' labels <- c(rep("Group1", 120), rep("Group2", 100))
#' QRscore.test(samples = data, labels = labels, 
#'             approx = "asymptotic", measure = "mean", zero_inflation = FALSE)
#'
#' # Two-sample test with automatically estimated weights from ZINB model
#' group1 <- rzinbinom(100, size = 2, mu = 40, pi = 0.1)
#' group2 <- rzinbinom(200, size = 1, mu = 40, pi = 0.1) # Different size and zero-inflation
#' data <- c(group1, group2)
#' labels <- c(rep("Group1", 100), rep("Group2", 200))
#' QRscore.test(samples = data, labels = labels, alternative = "two.sided", 
#'               approx = "asymptotic", measure = "dispersion")
#'
#' # Three-sample test with automatically estimated weights from NB model
#' group1 <- rzinbinom(150, size = 1, mu = 30, pi = 0.1)
#' group2 <- rzinbinom(100, size = 2, mu = 30, pi = 0.1)
#' group3 <- rzinbinom(30, size = 3, mu = 30, pi = 0.1)
#' data <- c(group1, group2, group3)
#' labels <- c(rep("Group1", 150), rep("Group2", 100), rep("Group3", 30))
#' QRscore.test(samples = data, labels = labels, alternative = "two.sided", 
#'               approx = "asymptotic", measure = "dispersion")
#' @export
#' 
QRscore.test <- function(samples, labels = NULL,
                         p = NULL, wList = NULL,
                         alternative = "two.sided", approx = "resample", type = "unbiased",
                         n_mom = 100, resamp_number = 5000,
                         zero_inflation = TRUE, LR.test = FALSE, pi_threshold = 0.95,
                         gene.name = NULL, measure = "mean", perturb = TRUE,
                         use_base_r = TRUE, seed = 1) {
  
  assertthat::assert_that(!is.null(samples), msg = "Samples cannot be NULL.")
  assertthat::assert_that((alternative == "two.sided" | alternative == "greater" | alternative == "less"),
                          msg = "Please specify a valid alternative (`two.sided`, `greater`, or `less`)")
  assertthat::assert_that((approx == "resample" | approx == "asymptotic"),
                          msg = "Please specify a valid approximation (`resample`, `asymptotic`)")
  assertthat::assert_that((type == "unbiased" | type == "valid" | type == "both"),
                          msg = "Please specify a valid type (`unbiased`, `valid`, or `both`)")
  
  
  if (is.null(labels) || length(unique(labels)) == 1) {
    message(date(), ": Performing a one-sample test.")
    assertthat::assert_that(!is.null(p), msg = "Exponent p must be provided for the one-sample test.")
    assertthat::assert_that(!is.null(wList), msg = "wList must be provided for one-sample test.")
    assertthat::assert_that(length(samples) == length(wList), msg = "Length of wList must be same as length of the sample.")
    assertthat::assert_that(p > 0 & (abs(p - round(p)) < .Machine$double.eps^0.5),
                            msg = "Currently only supporting positive integer values of p. Check that p is a positive integer.")
    assertthat::assert_that((approx == "resample"),
                            msg = "Currently only supporting `resample`.")
    return(QRscore_Flex(samples,NULL, p, wList, alternative, approx, type, n_mom, resamp_number))
  }
  assertthat::assert_that(length(labels)==length(samples), msg = "Sample length should be the same as label length")
  unique_labels <- sort(unique(labels))
  n_groups <- length(unique_labels)
  sample_list <- lapply(unique_labels, function(l) samples[labels == l])
  
  if (n_groups == 2) {
    x <- sample_list[[1]]
    y <- sample_list[[2]]
    x_label <- unique_labels[1]
    y_label <- unique_labels[2]
    
    if (length(x) > length(y)) {
      x <- sample_list[[2]]
      y <- sample_list[[1]]
      x_label <- unique_labels[2]
      y_label <- unique_labels[1]
    }
    
    
    
    message(date(), ": Performing a two-sample test. '", x_label, 
            "' is x and '", y_label, "' is y.")
    if(is.null(wList) || is.null(p)){
      message(date(), ": Use ZINB or NB model to estimate the weights since p and wList are not specified.")
      assertthat::assert_that(all(samples >= 0), msg = "Samples must not contain negative values.")
      return(QRscore_ZINB(x, y, zero_inflation, LR.test, approx, alternative, resamp_number,
                          pi_threshold, gene.name, measure, seed))
      
      
    } else{
      assertthat::assert_that(length(x) + 1 == length(wList), msg = "Length of wList must be one more than length of sample x.")
      assertthat::assert_that(p > 0 & (abs(p - round(p)) < .Machine$double.eps^0.5),
                              msg = "Currently only supporting positive integer values of p. Check that p is a positive integer.")
      assertthat::assert_that((approx == "resample"),
                              msg = "Currently only supporting `resample`.")
      # edge case
      if (identical(wList,(length(x):0)) & p == 1 & use_base_r) {
        message(date(), ": Settings correspond to Mann-Whitney, using base R...")
        return(wilcox.test(x, y, alternative = alternative, correct=FALSE)$p.value)
      } else {
        message(date(), ": Using native R to compute p-value...")
        return(QRscore_Flex(x, y, p, wList, alternative, approx, type, n_mom, resamp_number))
      }
    }
  } else if (n_groups > 2) {
    message(date(), ": Performing a multi-sample test.")
    message(date(), ": Use ZINB or NB model to estimate the weights.")
    return(QRscore_ZINB_nSamples(samples, labels, zero_inflation, LR.test, approx,
                                 resamp_number, pi_threshold, gene.name, measure, perturb, seed))
  }
}

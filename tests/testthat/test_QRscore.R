# Test QRscore_Flex function for one-sample case
test_that("QRscore_Flex handles one-sample case correctly", {
    set.seed(1)
    x <- abs(rnorm(10))
    wList <- rep(1, 10)
    p_value <- QRscore_Flex(
        x = x, p = 2, wList = wList,
        alternative = "two.sided", approx = "resample"
    )
    expect_true(is.numeric(p_value))
    expect_true(p_value <= 1 & p_value >= 0)
})

# Test QRscore_Flex function for two-sample case
test_that("QRscore_Flex handles two-sample case correctly", {
    set.seed(1)
    x <- abs(rnorm(30))
    y <- abs(rnorm(100))
    wList <- rep(1, 31)
    p_value <- QRscore_Flex(
        x = x, y = y, p = 2, wList = wList,
        alternative = "two.sided", approx = "resample",
        resamp_number = 5000
    )
    expect_true(is.numeric(p_value))
    expect_true(p_value <= 1 & p_value >= 0)
})

# Test QRscore_Flex function for two-sample case with asymptotic approximation
test_that("QRscore_Flex handles two-sample approximation correctly", {
    set.seed(1)
    x <- abs(rnorm(100))
    y <- abs(rnorm(100))
    wList <- 0:100
    p_value <- QRscore_Flex(
        x = x, y = y, p = 1, wList = wList,
        alternative = "two.sided",
        approx = "asymptotic"
    )
    expect_true(is.numeric(p_value))
    expect_true(p_value <= 1 & p_value >= 0)
})


# Test QRscore_ZINB
test_that("QRscore_ZINB runs without error and returns numeric", {
    set.seed(1)
    x <- rzinbinom(30, size = 2, mu = 20, pi = 0.1)
    y <- rzinbinom(30, size = 1, mu = 20, pi = 0.1)

    p_value <- QRscore_ZINB(x, y,
        zero_inflation = TRUE, LR.test = FALSE,
        alternative = "two.sided", approx = "asymptotic", 
        measure = "dispersion"
    )
    expect_true(p_value <= 1 & p_value >= 0)
    expect_true(is.numeric(p_value))
})


test_that("QRscore_ZINB identifies distribution differences correctly", {
    n_sim <- 500 # Number of simulations
    pvalues_null <- numeric(n_sim)
    pvalues_alt <- numeric(n_sim)

    for (i in seq_len(n_sim)) {
        # Null hypothesis: both samples from the same distribution
        x_null <- rnbinom(100, size = 2, mu = 20)
        y_null <- rnbinom(100, size = 2, mu = 20)
        pvalues_null[i] <- QRscore_ZINB(x_null, y_null,
            zero_inflation = FALSE, LR.test = FALSE,
            alternative = "two.sided", approx = "asymptotic", measure = "mean"
        )

        # Alternative hypothesis: samples from different distributions
        x_alt <- rnbinom(100, size = 2, mu = 20)
        y_alt <- rnbinom(100, size = 10, mu = 20)
        pvalues_alt[i] <- QRscore_ZINB(x_alt, y_alt,
            zero_inflation = FALSE, LR.test = FALSE,
            alternative = "two.sided", approx = "asymptotic", 
            measure = "dispersion"
        )
    }

    prop_null <- mean(pvalues_null < 0.05, na.rm = TRUE)
    prop_alt <- mean(pvalues_alt < 0.05, na.rm = TRUE)

    expect_true(prop_null <= 0.07)
    expect_true(prop_alt > 0.05)
})


test_that("size factors stay aligned when the first group is larger", {
    set.seed(9)
    sfA <- runif(120, 0.3, 3)
    sfB <- runif(60, 0.3, 3)
    a <- rnbinom(120, size = 3, mu = 20 * sfA)
    b <- rnbinom(60, size = 3, mu = 20 * sfB)
    set.seed(1)
    p_test <- QRscoreTest(
        samples = c(a, b), labels = rep(c("A", "B"), c(120, 60)),
        size_factors = c(sfA, sfB), measure = "mean",
        zero_inflation = FALSE, approx = "asymptotic"
    )
    # Group B is smaller, so it becomes x and must keep its own size factors.
    set.seed(1)
    p_direct <- QRscore_ZINB(
        b, a, size_factors = c(sfB, sfA), measure = "mean",
        zero_inflation = FALSE, approx = "asymptotic"
    )
    expect_equal(p_test, p_direct)
})


test_that("NB dispersion estimate is the MLE when glm.nb diverges", {
    set.seed(3)
    # Zero-heavy counts with a large mean: theta.ml's Newton iterations
    # start far from the optimum and run off to theta ~ 3e5 on this data.
    y <- rnbinom(80, mu = 60, size = 1)
    y[runif(80) < 0.4] <- 0
    est <- QRscore:::.estimate_nb_parameters(y)
    fit <- QRscore:::.fit_nb_model(y)
    glm_ll <- suppressWarnings(as.numeric(logLik(MASS::glm.nb(y ~ 1))))
    expect_true(fit$loglik >= glm_ll - 1e-6)
    expect_true(est$beta < 5)
})

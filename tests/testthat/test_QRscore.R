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

# Test dzinbinom function
test_that("dzinbinom handles invalid 'pi' values", {
    expect_warning(.dzinbinom(x = 1, mu = 10, theta = 5, pi = -0.1))
    expect_warning(.dzinbinom(x = 1, mu = 10, theta = 5, pi = 1.1))
})

test_that("dzinbinom returns expected probabilities", {
    x_values <- 0:5
    expected_prob <- log(0.8) + dnbinom(x_values, mu = 3, size = 2, log = TRUE)
    expected_prob[1] <- log(exp(expected_prob[1]) + 0.2)
    calculated_prob <- .dzinbinom(
        x = x_values, mu = 3, theta = 2, pi = 0.2,
        log = TRUE
    )
    expect_true(all.equal(calculated_prob, expected_prob))
})


# Test qzinbinom function
test_that("qzinbinom handles invalid 'pi' values", {
    expect_warning(.qzinbinom(p = 0.5, mu = 10, theta = 5, pi = -0.1))
    expect_warning(.qzinbinom(p = 0.5, mu = 10, theta = 5, pi = 1.1))
})

test_that("qzinbinom returns correct quantiles", {
    expect_equal(
        .qzinbinom(p = 0.5, mu = 10, theta = 5, pi = 0),
        qnbinom(0.5, mu = 10, size = 5)
    )
})

# Test rzinbinom function
test_that("rzinbinom handles invalid 'pi' values", {
    expect_warning(rzinbinom(n = 10, mu = 10, theta = 5, pi = -0.1))
    expect_warning(rzinbinom(n = 10, mu = 10, theta = 5, pi = 1.1))
})

test_that("rzinbinom generates correct number of samples", {
    set.seed(1)
    samples <- rzinbinom(n = 100, mu = 10, theta = 5, pi = 0.2)
    expect_equal(length(samples), 100)
    expect_true(sum(samples == 0) > 0) # Check for zero inflation
})

# Test getCompositionPValue function
test_that("getCompositionPValue computes correct exact p-values for small n and
          k", {
    expect_true(getCompositionPValue(
        t = 0.5, n = 30, k = 5, p = 1,
        wList = rep(1, 5),
        alternative = "two.sided"
    ) >= 0)
})

# Test rank_x function
test_that("rank_x returns correct indices with tie-breaking", {
    set.seed(1)
    x <- rnorm(5)
    y <- rnorm(5)

    computed_ranks <- rank_x(x, y)
    expect_true(all(computed_ranks >= 1 & computed_ranks <= length(c(x, y))))
    expect_equal(length(unique(computed_ranks)), length(computed_ranks))
})


# Test computeweight_mean function
test_that("computeweight_mean computes weights correctly and contains no NAs", {
    results <- computeweight_mean(beta = 2, mu = 10, pi = 0.1, n = 50, k = 25)
    expect_true(is.list(results))
    expect_true("weight" %in% names(results))
    expect_true("var" %in% names(results))
    expect_true(all(!is.na(results$weight)), "No NAs in weight vector")
    expect_true(!is.na(results$var), "Variance is not NA")
})

# Test computeweight_disp function
test_that("computeweight_disp computes weights correctly and contains no NAs", {
    results <- computeweight_disp(beta = 2, mu = 10, pi = 0.1, n = 50, k = 25)
    expect_true(is.list(results))
    expect_true("weight" %in% names(results))
    expect_true("var" %in% names(results))
    expect_true(all(!is.na(results$weight)), "No NAs in weight vector")
    expect_true(!is.na(results$var), "Variance is not NA")
})


test_that("computeweight_mean variance matches the variance of its weights", {
    for (p in list(c(2, 10, 0), c(5, 3, 0.2), c(100, 5, 0), c(1e4, 3, 0.1))) {
        res <- computeweight_mean(p[1], p[2], p[3], n = 20000, k = 20000)
        expect_equal(res$var, var(res$weight), tolerance = 0.01)
    }
})

test_that("weights are finite and approach the Poisson limit for large beta", {
    x <- 0:10
    for (beta in c(1e3, 1e6)) {
        expect_false(anyNA(computeweight_disp(beta, 4, 0.1, 50, 25)$weight))
        expect_false(anyNA(computeweight_mean(beta, 4, 0.1, 50, 25)$weight))
    }
    z <- QRscore:::.z_function_disp(x, 1e8, 4, 0)
    expect_equal(1e8 * z, -((x - 4)^2 - x) / 2, tolerance = 1e-5)
})


test_that("weights stay finite when P(X = 0) underflows (large mu, large beta)", {
    for (p in list(c(1e6, 1000, 0), c(1e6, 5000, 0.1), c(100, 3e5, 0))) {
        expect_silent(d <- computeweight_disp(p[1], p[2], p[3], 2000, 2000))
        m <- computeweight_mean(p[1], p[2], p[3], 2000, 2000)
        expect_false(anyNA(c(d$weight, m$weight, d$var, m$var)))
        expect_equal(d$var, var(d$weight), tolerance = 0.05)
    }
})

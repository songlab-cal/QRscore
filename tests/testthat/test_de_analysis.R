### Test QRscoreGenetest

test_that("QRscoreGenetest runs without error on simple input", {
    set.seed(1)
    data <- matrix(rnbinom(2000, mu = 20, size = 2), nrow = 10)
    rownames(data) <- paste0("gene", seq_len(10))
    labels <- rep(c("A", "B"), each = 100)
    result <- QRscoreGenetest(data, labels, test_mean = TRUE, test_dispersion = TRUE, pairwise_test = TRUE, pairwise_logFC = TRUE)
    expect_type(result, "list")
    expect_true(all(c("mean_test", "var_test") %in% names(result)))
    expect_s3_class(result$mean_test, "data.frame")
    expect_s3_class(result$var_test, "data.frame")
})


test_that("QRscoreGenetest output rownames set is the same as input gene names
          set", {
    set.seed(1)
    data <- matrix(rnbinom(2000, mu = 20, size = 2), nrow = 10)
    rownames(data) <- paste0("gene", seq_len(10))
    labels <- rep(c("A", "B"), each = 100)
    result <- QRscoreGenetest(data, labels, pairwise_test = TRUE)
    expect_true(setequal(rownames(result$mean_test), rownames(data)))
})

test_that("QRscoreGenetest get pairwise columns appear when requested", {
    set.seed(2)
    data <- matrix(rnbinom(1500, mu = 20, size = 2), nrow = 5)
    rownames(data) <- paste0("g", seq_len(5))
    labels <- rep(c("X", "Y", "Z"), each = 100)

    result <- QRscoreGenetest(data, labels,
        pairwise_test = TRUE,
        pairwise_logFC = TRUE
    )

    expect_true(any(grepl("Pairwise_Test", colnames(result$mean_test))))
    expect_true(any(grepl("Log_FC_Mean", colnames(result$mean_test))))
})


test_that("QRscoreTest performs one-sample test correctly", {
    set.seed(1)
    x <- abs(rnorm(10))
    w <- rep(1, 10)
    pval <- QRscoreTest(
        samples = x, p = 2, wList = w,
        alternative = "two.sided", approx = "resample"
    )
    expect_true(is.numeric(pval))
    expect_true(pval >= 0 && pval <= 1)
})

test_that("QRscoreTest performs two-sample test with weights", {
    set.seed(1)
    group1 <- rnorm(50)
    group2 <- rnorm(50, mean = 1)
    data <- c(group1, group2)
    labels <- rep(c("A", "B"), each = 50)
    wList <- c(50:0)
    pval <- QRscoreTest(
        samples = data, labels = labels, p = 1, wList = wList,
        alternative = "two.sided", approx = "resample"
    )
    expect_true(is.numeric(pval))
    expect_true(pval >= 0 && pval <= 1)
})

test_that("QRscoreTest performs two-sample test using NB weights", {
    set.seed(1)
    group1 <- rzinbinom(100, size = 2, mu = 20, pi = 0)
    group2 <- rzinbinom(100, size = 2, mu = 30, pi = 0)
    data <- c(group1, group2)
    labels <- rep(c("A", "B"), each = 100)
    pval <- QRscoreTest(
        samples = data, labels = labels, approx = "asymptotic",
        zero_inflation = FALSE, measure = "mean"
    )
    expect_true(is.numeric(pval))
    expect_true(pval >= 0 && pval <= 1)
})

test_that("QRscoreTest performs multi-sample test with NB model", {
    set.seed(1)
    g1 <- rnbinom(100, size = 1, mu = 20)
    g2 <- rnbinom(100, size = 2, mu = 20)
    g3 <- rnbinom(100, size = 3, mu = 20)
    data <- c(g1, g2, g3)
    labels <- rep(c("G1", "G2", "G3"), each = 100)
    pval <- QRscoreTest(
        samples = data, labels = labels, approx = "asymptotic",
        zero_inflation = FALSE, measure = "dispersion"
    )
    expect_true(is.numeric(pval))
    expect_true(pval >= 0 && pval <= 1)
})

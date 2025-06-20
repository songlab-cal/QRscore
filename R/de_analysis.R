#' @keywords internal
.prepare_gene_input <- function(expr_matrix) {
    expr_list <- split(
        t(expr_matrix),
        rep(seq_len(nrow(expr_matrix)),
            each = ncol(expr_matrix)
        )
    )
    gene_names <- rownames(expr_matrix)
    gene_names_list <- split(gene_names,
        f = seq_len(length(gene_names))
    )
    Map(function(x, y) list(x, y), expr_list, gene_names_list)
}

#' @keywords internal
.extract_and_format_results <- function(test_type, test_results,
                                        pval_name, rank = TRUE) {
    if (!is.null(test_results)) {
        test_results[[paste0(test_type, "_adj_p_value")]] <-
            p.adjust(test_results[[pval_name]], method = "BH")
        if (rank) {
            test_results <-
                test_results[order(test_results[[pval_name]]), ]
        }
        pvals <- c(pval_name, paste0(test_type, "_adj_p_value"))
        test_results <- test_results[, c(
            pvals,
            setdiff(names(test_results), pvals)
        )]
    }
    return(test_results)
}

# Updated .pairwise_tests to support size_factors
.pairwise_tests <- function(
    gene_vec, labels, unique_labels, test_mean, test_dispersion, pairwise_test, 
    pairwise_logFC, gene.name, mean_test_list, var_test_list, 
    size_factors = NULL, ...) {
    for (i in seq_len(length(unique_labels) - 1)) {
        for (j in (i + 1):length(unique_labels)) {
            idx_i <- labels == unique_labels[i]
            idx_j <- labels == unique_labels[j]
            x <- gene_vec[idx_i]
            y <- gene_vec[idx_j]
            sf_pair <- if (!is.null(size_factors)) 
              c(size_factors[idx_i], size_factors[idx_j]) else NULL
            lbl <- paste0(unique_labels[j], "_vs_", unique_labels[i])

            if (pairwise_test && test_mean) {
                mean_p <- QRscoreTest(
                    samples = c(x, y),
                    labels = rep(unique_labels[c(i, j)], c(length(x), length(y))),
                    measure = "mean",
                    gene.name = gene.name,
                    size_factors = sf_pair,
                    ...
                )
                mean_test_list[[paste0("Pairwise_Test_Mean_", lbl)]] <- mean_p
            }
            if (pairwise_test && test_dispersion) {
                var_p <- QRscoreTest(
                    samples = c(x, y),
                    labels = rep(unique_labels[c(i, j)], c(length(x), length(y))),
                    measure = "dispersion",
                    gene.name = gene.name,
                    size_factors = sf_pair,
                    ...
                )
                var_test_list[[paste0("Pairwise_Test_Var_", lbl)]] <- var_p
            }
            if (pairwise_logFC && test_mean) {
                logfc_mean <- if (!is.null(sf_pair)) {
                    log2(mean(y / size_factors[idx_j]) / 
                        mean(x / size_factors[idx_i]))
                } else {
                    log2(mean(y) / mean(x))
                }
                mean_test_list[[paste0("Log_FC_Mean_", lbl)]] <- logfc_mean
            }
            if (pairwise_logFC && test_dispersion) {
                logfc_var <- if (!is.null(sf_pair)) {
                    log2(var(y / size_factors[idx_j]) / 
                        var(x / size_factors[idx_i]))
                } else {
                    log2(var(y) / var(x))
                }
                var_test_list[[paste0("Log_FC_Var_", lbl)]] <- logfc_var
            }
        }
    }
    return(list(mean = mean_test_list, var = var_test_list))
}



#' @keywords internal
.add_test_result <- function(x, y, label1, label2,
                             gene_vec, gene.name, measure,
                             do_test, do_logfc,
                             test_list, suffix, size_factors = NULL,...) {
    if (do_test) {
        pval <- QRscoreTest(
            samples = c(x, y),
            labels = rep(
                c(label1, label2),
                c(length(x), length(y))
            ),
            measure = measure,
            gene.name = gene.name, size_factors = size_factors,...
        )
        test_list[[paste0("QRscore_", suffix, "_p_value")]] <- pval
    }

    if (do_logfc) {
        stat_fn <- if (measure == "mean") mean else var
        if (!is.null(size_factors)) {
            idx_x <- seq_len(length(x))
            idx_y <- (length(x)+1):(length(x)+length(y))
            x <- x / size_factors[idx_x]
            y <- y / size_factors[idx_y]
        }
        logfc <- log2(stat_fn(y) / stat_fn(x))
        test_list[[paste0(
            "Log_FC_", suffix, "_",
            label2, "_vs_", label1
        )]] <- logfc
    }
    return(test_list)
}

#' @keywords internal
.QRscore_one_gene_test <- function(gene_list, labels, size_factors = NULL,
                                   pairwise_test = TRUE, pairwise_logFC = TRUE, 
                                   test_mean = TRUE,
                                   test_dispersion = FALSE, ...) {
    gene_vec <- gene_list[[1]]; gene.name <- gene_list[[2]]
    unique_labels <- sort(unique(labels))
    mean_test_list <- list(); var_test_list <- list()
    if (length(unique_labels) == 2) {
        x <- gene_vec[labels == unique_labels[1]]
        y <- gene_vec[labels == unique_labels[2]]
        if (!is.null(size_factors)) {
            idx_x <- labels == unique_labels[1]
            idx_y <- labels == unique_labels[2]
            sf_combined <- c(size_factors[idx_x], size_factors[idx_y])
        } else {
            sf_combined <- NULL
        }
        
        if (test_mean) {
            mean_test_list <- .add_test_result(
                x, y, unique_labels[1], unique_labels[2], gene_vec, gene.name, 
                "mean", pairwise_test, pairwise_logFC, mean_test_list,"Mean",
                size_factors = sf_combined, ...
                )
        }
        if (test_dispersion) {
            var_test_list <- .add_test_result(
              x, y, unique_labels[1], unique_labels[2], gene_vec, gene.name, 
              "dispersion", pairwise_test, pairwise_logFC, var_test_list, 
              "Var", size_factors = sf_combined, ...)
        }
    } else if (length(unique_labels) > 2) {
        if (test_mean) {
            mean_test_list[["QRscore_Mean_p_value"]] <- QRscoreTest(
                samples = gene_vec, labels = labels, measure = "mean", 
                gene.name = gene.name, size_factors = size_factors,...)
        } 
        if (test_dispersion) {
            var_test_list[["QRscore_Var_p_value"]] <- QRscoreTest(
                samples = gene_vec, labels = labels, measure = "dispersion", 
                gene.name = gene.name,  size_factors = size_factors,...
            )
        }
        if (pairwise_test || pairwise_logFC) {
            out <- .pairwise_tests(
                gene_vec, labels, unique_labels, test_mean, 
                test_dispersion, pairwise_test, pairwise_logFC,
                gene.name, mean_test_list, var_test_list,  
                size_factors = size_factors,...
            )
            mean_test_list <- out$mean; var_test_list <- out$var
        }
    }
    mean_test_df <- if (length(mean_test_list) > 0) {
        as.data.frame(mean_test_list, check.names = FALSE) } else NULL
    var_test_df <- if (length(var_test_list) > 0) {
        as.data.frame(var_test_list, check.names = FALSE)} else NULL
    return(list(mean_test = mean_test_df, var_test = var_test_df))
}


#' Perform Differential Expression Analysis using QRscore
#'
#' This function performs differential expression analysis using the QRscore 
#' method.
#'
#' @param expr_matrix A prefiltered raw count matrix (gene \eqn{\times} 
#' sample).
#' @param labels A vector of labels corresponding to the samples.
#' @param size_factors A vector of size factors estimated from DESeq2
#' @param pairwise_test Logical, whether to perform pairwise test statistics 
#' for more than two groups.
#' @param pairwise_logFC Logical, whether to calculate pairwise log fold 
#' changes for mean and variance.
#' @param test_mean Logical, whether to test the mean.
#' @param test_dispersion Logical, whether to test the dispersion.
#' @param num_cores Integer, number of cores to use for parallel processing.
#' @param verbose Logical, whether to suppress all messages except warnings and 
#' errors.
#' @param rank Logical. If TRUE, returns a ranked list of DEGs or DDGs sorted 
#' by p-values.
#' @param seed Integer. Set seed in parallel computing.
#' @param ... Additional arguments passed to the `QRscoreTest` function.
#' @return A list with two data frames: `mean_test` and `var_test`.
#' @export
#'
#' @examples
#' data(example_dataset)
#' results <- QRscoreGenetest(
#'   expr_matrix = example_dataset$example_data1,
#'   labels = example_dataset$labels1,
#'   size_factors = example_dataset$size_factors1,
#'   pairwise_test = FALSE, pairwise_logFC = TRUE,
#'   test_mean = TRUE, test_dispersion = TRUE, num_cores = 2,
#'   approx = "asymptotic"
#' )
#' head(results$var_test)
#'
#' results2 <- QRscoreGenetest(
#'   expr_matrix = example_dataset$example_data2,
#'   labels = example_dataset$labels2,
#'   size_factors = example_dataset$size_factors2,
#'   pairwise_test = TRUE, pairwise_logFC = TRUE,
#'   test_mean = TRUE, test_dispersion = FALSE, num_cores = 2,
#'   approx = "asymptotic"
#' )
#' head(results2$mean_test)
QRscoreGenetest <- function(
    expr_matrix, labels, size_factors = NULL, pairwise_test = FALSE,
    pairwise_logFC = FALSE, test_mean = TRUE, test_dispersion = FALSE,
    num_cores = 1, verbose = FALSE, rank = TRUE, seed = 1, ...) {
    bpparam <- if (num_cores > 1 && 
                   BiocParallel::multicoreWorkers() > 1) {
        BiocParallel::MulticoreParam(workers = num_cores, 
                                     RNGseed = seed)
    } else {
        BiocParallel::SerialParam(RNGseed = seed)
    }
    if (!is.null(size_factors)) {
        if (any(size_factors <= 0)) stop("All size factors must be > 0.")
    }

    combined_list <- .prepare_gene_input(expr_matrix)

    run_test <- function() {
        BiocParallel::bplapply(
            combined_list, .QRscore_one_gene_test,
            labels = labels,
            size_factors = size_factors,
            pairwise_test = pairwise_test,
            pairwise_logFC = pairwise_logFC,
            test_mean = test_mean,
            test_dispersion = test_dispersion,
            BPPARAM = bpparam, ...
        )
    }

    results_list <- if (!verbose) suppressMessages(run_test()) else run_test()

    mean_test_results <- do.call(
      rbind, lapply(results_list, function(x) x$mean_test))
    var_test_results  <- do.call(
      rbind, lapply(results_list, function(x) x$var_test))

    gene_names <- rownames(expr_matrix)
    if (!is.null(mean_test_results)) {
        rownames(mean_test_results) <- gene_names[seq_len(nrow(mean_test_results))]
    }
    if (!is.null(var_test_results)) {
        rownames(var_test_results) <- gene_names[seq_len(nrow(var_test_results))]
    }

    mean_test_results <- .extract_and_format_results(
        "QRscore_Mean", mean_test_results, "QRscore_Mean_p_value", rank)
    var_test_results <- .extract_and_format_results(
        "QRscore_Var", var_test_results, "QRscore_Var_p_value", rank)

    return(list(mean_test = mean_test_results, 
                var_test = var_test_results))
}




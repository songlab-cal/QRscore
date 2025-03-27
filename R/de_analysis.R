#' @keywords internal
.prepare_gene_input <- function(normalized_matrix) {
    normalized_list <- split(
        t(normalized_matrix),
        rep(seq_len(nrow(normalized_matrix)), 
            each = ncol(normalized_matrix))
    )
    gene_names <- rownames(normalized_matrix)
    gene_names_list <- split(gene_names, 
                             f = seq_len(length(gene_names)))
    Map(function(x, y) list(x, y), normalized_list, gene_names_list)
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
      test_results <- test_results[, c(pvals, 
                                       setdiff(names(test_results), pvals))]
    }
    return(test_results)
}

#' @keywords internal
.pairwise_tests <- function(gene_vec, labels, unique_labels,
                            test_mean, test_dispersion,
                            pairwise_test, pairwise_logFC,
                            gene.name, mean_test_list,
                            var_test_list, ...) {
    for (i in seq_len(length(unique_labels) - 1)) {
        for (j in (i + 1):length(unique_labels)) {
            x <- gene_vec[labels == unique_labels[i]]
            y <- gene_vec[labels == unique_labels[j]]
            lbl <- paste0(unique_labels[j], "_vs_", unique_labels[i])
            if (pairwise_test && test_mean) {
                mean_p <- QRscoreTest(
                    samples = c(x, y),
                    labels = rep(unique_labels[c(i, j)],
                                  c(length(x), length(y))),
                    measure = "mean",
                    gene.name = gene.name, ...
                    )
                mean_test_list[[paste0("Pairwise_Test_Mean_", 
                                        lbl)]] <- mean_p
            }
            if (pairwise_test && test_dispersion) {
                var_p <- QRscoreTest(
                    samples = c(x, y),
                    labels = rep(unique_labels[c(i, j)],
                                  c(length(x), length(y))),
                    measure = "dispersion",
                    gene.name = gene.name, ...
                    )
                var_test_list[[paste0("Pairwise_Test_Var_", lbl)]] <- var_p
            }
            if (pairwise_logFC && test_mean) {
                logfc_mean <- log2(mean(y) / mean(x))
                mean_test_list[[paste0("Log_FC_Mean_", lbl)]] <- logfc_mean
            } 
            if (pairwise_logFC && test_dispersion) {
                logfc_var <- log2(var(y) / var(x))
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
                             test_list, suffix, ...) {
    if (do_test) {
        pval <- QRscoreTest(samples = c(x, y),
                            labels = rep(c(label1, label2),
                                         c(length(x), length(y))),
                            measure = measure,
                            gene.name = gene.name, ...)
        test_list[[paste0("QRscore_", suffix, "_p_value")]] <- pval
    }
    
    if (do_logfc) {
        stat_fn <- if (measure == "mean") mean else var
        logfc <- log2(stat_fn(y) / stat_fn(x))
        test_list[[paste0("Log_FC_", suffix, "_",
                          label2, "_vs_", label1)]] <- logfc
    }
    return(test_list)
}

#' @keywords internal
.QRscore_one_gene_test <- function(gene_list, labels, pairwise_test = TRUE,
                                   pairwise_logFC = TRUE, test_mean = TRUE,
                                   test_dispersion = FALSE, ...) {
    gene_vec <- gene_list[[1]]; gene.name <- gene_list[[2]]
    unique_labels <- sort(unique(labels))
    mean_test_list <- list(); var_test_list <- list()
    if (length(unique_labels) == 2) {
        x <- gene_vec[labels == unique_labels[1]]
        y <- gene_vec[labels == unique_labels[2]]
        if (test_mean) {
            mean_test_list <- .add_test_result(x, y, unique_labels[1],
                                               unique_labels[2], gene_vec,
                                               gene.name, "mean", pairwise_test,
                                               pairwise_logFC, mean_test_list,
                                               "Mean",...)
        }
        if (test_dispersion) {
            var_test_list <- .add_test_result(x, y, unique_labels[1],
                                              unique_labels[2], gene_vec,
                                              gene.name, "dispersion", 
                                              pairwise_test, pairwise_logFC, 
                                              var_test_list, "Var",...)
        }
    } else if (length(unique_labels) > 2) {
        if (test_mean) {
            mean_test_list[["QRscore_Mean_p_value"]] <- QRscoreTest(
                samples = gene_vec, labels = labels,
                measure = "mean", gene.name = gene.name, ...)
        } 
        if (test_dispersion) {
            var_test_list[["QRscore_Var_p_value"]] <- QRscoreTest(
                samples = gene_vec, labels = labels,
                measure = "dispersion", gene.name = gene.name, ...
            )
        }
        if (pairwise_test || pairwise_logFC) {
            out <- .pairwise_tests(
                gene_vec, labels, unique_labels, test_mean, 
                test_dispersion, pairwise_test, pairwise_logFC,
                gene.name, mean_test_list, var_test_list, ...
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
#' @param normalized_matrix A prefiltered normalized matrix (gene \eqn{\times} 
#' sample).
#' @param labels A vector of labels corresponding to the samples.
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
#' results <- QRscoreGenetest(example_dataset$example_data1, 
#'   example_dataset$labels1,
#'   pairwise_test = FALSE, pairwise_logFC = TRUE,
#'   test_mean = TRUE, test_dispersion = TRUE, num_cores = 2,
#'   approx = "asymptotic"
#' )
#'
#' head(results$var_test)
#'
#' results2 <- QRscoreGenetest(example_dataset$example_data2, 
#'   example_dataset$labels2,
#'   pairwise_test = TRUE, pairwise_logFC = TRUE,
#'   test_mean = TRUE, test_dispersion = FALSE, num_cores = 2,
#'   approx = "asymptotic"
#' )
#' head(results2$mean_test)
QRscoreGenetest <- function(normalized_matrix, labels,
                            pairwise_test = FALSE,
                            pairwise_logFC = FALSE,
                            test_mean = TRUE,
                            test_dispersion = FALSE,
                            num_cores = 1,
                            verbose = FALSE,
                            rank = TRUE,
                            seed = 1, ...) {
    bpparam <- if (num_cores > 1 && 
                   BiocParallel::multicoreWorkers() > 1) {
        BiocParallel::MulticoreParam(workers = num_cores, 
                                     RNGseed = seed)
    } else {
        BiocParallel::SerialParam(RNGseed = seed)
    }
    combined_list <- .prepare_gene_input(normalized_matrix)
    run_test <- function() {
        BiocParallel::bplapply(
            combined_list, .QRscore_one_gene_test,
            labels = labels,
            pairwise_test = pairwise_test,
            pairwise_logFC = pairwise_logFC,
            test_mean = test_mean,
            test_dispersion = test_dispersion,
            BPPARAM = bpparam, ...)
    }
    results_list <- if (!verbose) suppressMessages(run_test()) else run_test()
    mean_test_results <- do.call(rbind, lapply(results_list, function(x) {
        x$mean_test
    }))
    var_test_results <- do.call(rbind, lapply(results_list, function(x) {
        x$var_test
    }))
    gene_names <- rownames(normalized_matrix)
    if (!is.null(mean_test_results)) {
        rownames(mean_test_results) <- gene_names[
            seq_len(nrow(mean_test_results))]
    }
    if (!is.null(var_test_results)) {
        rownames(var_test_results) <- gene_names[
            seq_len(nrow(var_test_results))]
    }
    mean_test_results <- .extract_and_format_results(
        "QRscore_Mean", mean_test_results, "QRscore_Mean_p_value", rank)
    var_test_results <- .extract_and_format_results(
        "QRscore_Var", var_test_results, "QRscore_Var_p_value", rank)
    return(list(mean_test = mean_test_results, 
                var_test = var_test_results))
}


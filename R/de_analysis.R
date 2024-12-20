#' Perform Differential Expression Analysis using QRscore
#'
#' This function performs differential expression analysis using the QRscore method.
#'
#' @param normalized_matrix A prefiltered normalized matrix (gene \eqn{\times} sample).
#' @param labels A vector of labels corresponding to the samples.
#' @param pairwise_test Logical, whether to perform pairwise test statistics for more than two groups.
#' @param pairwise_logFC Logical, whether to calculate pairwise log fold changes for mean and variance.
#' @param test_mean Logical, whether to test the mean.
#' @param test_dispersion Logical, whether to test the dispersion.
#' @param num_cores Integer, number of cores to use for parallel processing.
#' @param silent Logical, whether to suppress all messages except warnings and errors.
#' @param rank Logical. If TRUE, returns a ranked list of DEGs or DDGs sorted by p-values.
#' @param seed Integer. Set seed in parallel computing.
#' @param ... Additional arguments passed to the `QRscore.test` function.
#' @return A list with two data frames: `mean_test` and `var_test`. 
#' @export
#'
#' @examples
#' data(example_dataset)
#' results <- QRscore.genetest(example_dataset$example_data1, example_dataset$labels1, 
#'                                pairwise_test = FALSE, pairwise_logFC = TRUE, 
#'                                test_mean = TRUE, test_dispersion = TRUE, num_cores = 2,
#'                                approx = "asymptotic")
#'                              
#' head(results$var_test)
#'
#' results2 <- QRscore.genetest(example_dataset$example_data2, example_dataset$labels2, 
#'                                pairwise_test = TRUE, pairwise_logFC = TRUE, 
#'                                test_mean = TRUE, test_dispersion = FALSE, num_cores = 2,
#'                                approx = "asymptotic")
#' head(results2$mean_test)
QRscore.genetest <- function(normalized_matrix, labels, pairwise_test = FALSE, 
                             pairwise_logFC = FALSE, test_mean = TRUE, 
                             test_dispersion = FALSE, num_cores = 1, 
                             silent = TRUE, rank = TRUE, seed = 1,...) {
  if (num_cores > 1 && BiocParallel::multicoreWorkers() > 1) {
    bpparam <- BiocParallel::MulticoreParam(workers = num_cores, RNGseed = seed) 
  } else {
    bpparam <- BiocParallel::SerialParam(RNGseed = seed) 
  }
  
  adjusted_counts <- normalized_matrix
  normalized_list <- split(t(adjusted_counts), rep(1:nrow(adjusted_counts), each = ncol(adjusted_counts)))
  
  gene_names <- rownames(normalized_matrix)
  gene_names_list <- split(gene_names, f = 1:length(gene_names))
  combined_list <- Map(function(x, y) list(x, y), normalized_list, gene_names_list)
  
  if (silent) {
    results_list <- suppressMessages(BiocParallel::bplapply(
      combined_list, QRscore_one_gene_test, labels = labels, pairwise_test = pairwise_test, 
      pairwise_logFC = pairwise_logFC, test_mean = test_mean, test_dispersion = test_dispersion, BPPARAM = bpparam, ...))
  } else {
    results_list <- BiocParallel::bplapply(
      combined_list, QRscore_one_gene_test, labels = labels, pairwise_test = pairwise_test, 
      pairwise_logFC = pairwise_logFC, test_mean = test_mean, test_dispersion = test_dispersion, BPPARAM = bpparam, ...)
  }
  # Extract and combine the results into two data frames
  mean_test_results <- do.call(rbind, lapply(results_list, function(x) x$mean_test))
  var_test_results <- do.call(rbind, lapply(results_list, function(x) x$var_test))
  
  # Set the row names as gene names only if results are not null
  if (!is.null(mean_test_results)) rownames(mean_test_results) <- gene_names[1:nrow(mean_test_results)]
  if (!is.null(var_test_results)) rownames(var_test_results) <- gene_names[1:nrow(var_test_results)]
  
  # Adjust p-values and place them in the second column
  if (!is.null(mean_test_results)) {
    mean_test_results$QRscore_Mean_adj_p_value <- p.adjust(mean_test_results$QRscore_Mean_p_value, method = "BH")
    if (rank) {
      mean_test_results <- mean_test_results[order(mean_test_results$QRscore_Mean_p_value), ]
    }
    mean_test_results <- mean_test_results[, c("QRscore_Mean_p_value", "QRscore_Mean_adj_p_value", setdiff(names(mean_test_results), c("QRscore_Mean_p_value", "QRscore_Mean_adj_p_value")))]
  }
  if (!is.null(var_test_results)) {
    var_test_results$QRscore_Var_adj_p_value <- p.adjust(var_test_results$QRscore_Var_p_value, method = "BH")
    if (rank) {
      var_test_results <- var_test_results[order(var_test_results$QRscore_Var_p_value), ]
    }
    var_test_results <- var_test_results[, c("QRscore_Var_p_value", "QRscore_Var_adj_p_value", setdiff(names(var_test_results), c("QRscore_Var_p_value", "QRscore_Var_adj_p_value")))]
  }
  
  result <- list(mean_test = mean_test_results, var_test = var_test_results)
  return(result)
}

#' @keywords internal
QRscore_one_gene_test <- function(gene_list, labels, pairwise_test = TRUE, pairwise_logFC = TRUE, test_mean = TRUE, test_dispersion = FALSE, ...) {
  gene_vec <- gene_list[[1]]
  gene.name <- gene_list[[2]]
  unique_labels <- sort(unique(labels))
  
  mean_test_list <- list()
  var_test_list <- list()
  
  if (length(unique_labels) == 2) {
    x <- gene_vec[labels == unique_labels[1]]
    y <- gene_vec[labels == unique_labels[2]]
    
    if (test_mean) {
      QRscore_test_mean <- QRscore.test(samples = gene_vec, labels = labels, 
                                        measure = "mean", gene.name = gene.name...)
      log2_mean_ratio <- log2(mean(y) / mean(x))
      mean_test_list[["QRscore_Mean_p_value"]] <- QRscore_test_mean
      mean_test_list[[paste0("Log_FC_Mean_", unique_labels[2], "_vs_", unique_labels[1])]] <- log2_mean_ratio
    }
    
    if (test_dispersion) {
      QRscore_test_var <- QRscore.test(samples = gene_vec, labels = labels, 
                                       measure = "dispersion", gene.name=gene.name,...)
      log2_var_ratio <- log2(var(y) / var(x))
      var_test_list[["QRscore_Var_p_value"]] <- QRscore_test_var
      var_test_list[[paste0("Log_FC_Var_", unique_labels[2], "_vs_", unique_labels[1])]] <- log2_var_ratio
    }
    
  } else if (length(unique_labels) > 2) {
    if (test_mean) {
      QRscore_test_mean <- QRscore.test(samples = gene_vec, labels = labels, 
                                        measure = "mean", gene.name=gene.name,...)
      mean_test_list[["QRscore_Mean_p_value"]] <- QRscore_test_mean
    }
    if (test_dispersion) {
      QRscore_test_var <- QRscore.test(samples = gene_vec, labels = labels, 
                                       measure = "dispersion", gene.name=gene.name,...)
      var_test_list[["QRscore_Var_p_value"]] <- QRscore_test_var
    }
    
    if (pairwise_test || pairwise_logFC) {
      for (i in 1:(length(unique_labels) - 1)) {
        for (j in (i + 1):length(unique_labels)) {
          x <- gene_vec[labels == unique_labels[i]]
          y <- gene_vec[labels == unique_labels[j]]
          
          label_combination <- paste0(unique_labels[j], "_vs_", unique_labels[i])
          
          if (pairwise_test) {
            if (test_mean) {
              QRscore_test_mean_pw <- QRscore.test(samples = c(x, y), 
                                                   labels = rep(unique_labels[c(i, j)], c(length(x), length(y))), 
                                                   measure = "mean", 
                                                   gene.name = gene.name,...)
              mean_test_list[[paste0("Pairwise_Test_Mean_", label_combination)]] <- QRscore_test_mean_pw
            }
            
            if (test_dispersion) {
              QRscore_test_var_pw <- QRscore.test(samples = c(x, y), 
                                                  labels = rep(unique_labels[c(i, j)], c(length(x), length(y))), 
                                                  measure = "dispersion", 
                                                  gene.name = gene.name, ...)
              var_test_list[[paste0("Pairwise_Test_Var_", label_combination)]] <- QRscore_test_var_pw
            }
          }
          
          if (pairwise_logFC) {
            if (test_mean) {
              log2_mean_ratio <- log2(mean(y) / mean(x))
              mean_test_list[[paste0("Log_FC_Mean_", label_combination)]] <- log2_mean_ratio
            }
            
            if (test_dispersion) {
              log2_var_ratio <- log2(var(y) / var(x))
              var_test_list[[paste0("Log_FC_Var_", label_combination)]] <- log2_var_ratio
            }
          }
        }
      }
    }
  }
  
  mean_test_df <- if (length(mean_test_list) > 0) as.data.frame(mean_test_list, check.names = FALSE) else NULL
  var_test_df <- if (length(var_test_list) > 0) as.data.frame(var_test_list, check.names = FALSE) else NULL
  
  result <- list(mean_test = mean_test_df, var_test = var_test_df)
  return(result)
}



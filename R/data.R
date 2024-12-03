#' Example Dataset for MOCHIS DEG Analysis
#'
#' This dataset contains a subset of rounded normalized gene expression counts from GTEx whole blood tissue and associated labels for demonstration purposes.
#'
#' @format This dataset contains the following objects:
#' \describe{
#'   \item{example_data1}{A numeric matrix of normalized counts for all age groups.}
#'   \item{labels1}{A vector of sample labels for all age groups.}
#'   \item{example_data2}{A numeric matrix of normalized counts for age groups 20-39 and 60-79.}
#'   \item{labels2}{A vector of sample labels for age groups 20-39 and 60-79.}
#' }
#' @source Normalized counts from GTEx for package demonstration purposes.
#' @examples
#' data(example_dataset)
"example_dataset"

#' Example Dataset for MOCHIS DEG Analysis with Raw Counts
#'
#' This dataset contains raw gene expression counts from 3000 genes in GTEx whole blood tissue for demonstration purposes.
#'
#' @format This dataset contains a list with the following components
#' \describe{
#'   \item{COUNTS}{A numeric matrix of raw counts for 3000 genes.}
#'   \item{METADATA}{The metadata for the corresponding raw counts.}
#' }
#' @source Raw counts from GTEx for package demonstration purposes.
#' @examples
#' data(example_dataset_raw_3000_genes)
"example_dataset_raw_3000_genes"

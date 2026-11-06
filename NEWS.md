# QRscore 1.5.2

* Weights for zero-inflated models now follow the exact score of the ZINB
  likelihood.

# QRscore 1.5.1

* Size factors now stay aligned with their samples in two-group tests.
* Improved numerical stability of NB/ZINB weight estimation, including
  `glm.nb` dispersion fits and large dispersion values.
* Corrected the variance used in the mean test, improving power for
  low-count genes.
* Updated the citation to the published article (Cell Reports Methods, 2025).
* Restored the package vignette, which was no longer being built.

# QRscore 0.99.10 (07/19/2025)

Revise on the algorithm - taking size factors into account when estimating ZINB/NB parameters.

# QRscore 0.99.9 (04/13/2025)

Our package is now on Bioconductor! See the link https://bioconductor.org/packages/devel/bioc/html/QRscore.html

# QRscore 0.99.8 (03/28/2025)

### General

* Added a `NEWS.md` file to document package changes, improving transparency 
and user engagement.

* Enhanced quality assurance by incorporating unit tests using the `testthat` 
package.

### DESCRIPTION

* Optimized dependency management by categorizing packages in Depends, Imports, 
and Suggests according to Bioconductor guidelines. 

### README

* Updated `README.md` to include installation instructions directly from 
Bioconductor, ensuring consistency with Bioconductor standards.

### Man pages

* Introduced a package-level man page to enhance documentation.

* Added a runnable example for the exported function `rzinbinom`.

### R code

* Updated function naming to adhere to Bioconductor guidelines, removing periods 
from function names and prefixing non-exported functions with a dot
(e.g., `.dzinbinom`).

* Replaced `sapply()` with `vapply()` for more robust and predictable outputs.

* Changed sequence generation from `1:...` to `seq_len()` or `seq_along()` and 
replaced assignment `=` with `<-`.

* Minimized the use of `paste` in condition signals and discouraged the default 
use of `suppressWarnings()`.

### Check

* Addressed R CMD check warnings and notes, including mismatches in 
documentation and code. Updated R version dependency from 3.5.0 to 4.4.0.

* Handled top-level file issues and improved several other aspects based on 
BiocCheck recommendations.

### Vignette

* Renamed the vignette to `QRScore.Rmd` for easier access and updated with 
BiocStyle for improved formatting.

* Consolidated settings in `knitr::opts_chunk$set` and adhered to R conventions 
for assignments.

### Citation

* A citation file is added
<!-- badges: start -->

[![R-CMD-check](https://github.com/imbs-hl/survivalsvm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/imbs-hl/survivalsvm/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
Stable](https://img.shields.io/badge/lifecycle-Stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html#Stable)
[![CRAN Status](https://img.shields.io/badge/CRAN-survivalsvm-blue)](https://cran.r-project.org/package=survivalsvm)

<!-- badges: end -->
## Survival Support Vector Analysis
Cesaire J. K. Fouodo

### Introduction
This package performs support vectors analysis for data sets with survival outcome. Three approaches are available in the package: The regression approach takes censoring into account when formulating the inequality constraints of the support vector problem. In the ranking approach, the inequality constraints set the objective to maximize the concordance index for comparable pairs of observations. The hybrid approach combines the regression and ranking constraints in the same model.

### Installation
Installation from Github:
```R
devtools::install_github("imbs-hl/survivalsvm")
```

CRAN release coming soon.

### Usage
For usage in R, see ?survivalsvm in R. Most importantly, see the Examples section. As a first example you could try 

```R  
survivalsvm(Surv(time, status) ~ ., veteran, gamma.mu = 0.1)
```

### References
* Van Belle, V., Pelcmans, K., Van Huffel S. and Suykens J. A.K. (2011a). Improved performance on high-dimensional survival data by application of Survival-SVM. Bioinformatics (Oxford, England) 27, 87-94.
* Van Belle, V., Pelcmans, K., Van Huffel S. and Suykens J. A.K. (2011b). Support vector methods for survival analysis: a comparaison between ranking and regression approaches. Artificial Intelligence in medecine 53, 107-118.
* Césaire J. K. Fouodo and Inke R. König and Claus Weihs and Andreas Ziegler and Marvin N. Wright (2018) Support Vector Machines for Survival Analysis with R. The R Journal 10, 412-423.

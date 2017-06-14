[![Travis Build Status](https://travis-ci.org/imbs-hl/survivalsvm.svg?branch=master)](https://travis-ci.org/imbs-hl/survivalsvm)
[![Coverage Status](https://coveralls.io/repos/github/imbs-hl/survivalsvm/badge.svg?branch=master)](https://coveralls.io/github/imbs-hl/survivalsvm?branch=master)
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

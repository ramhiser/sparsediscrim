# sparsediscrim

[![Build Status](https://travis-ci.org/ramhiser/sparsediscrim.svg?branch=master)](https://travis-ci.org/ramhiser/sparsediscrim)

The R package `sparsediscrim` provides a collection of sparse and regularized discriminant
analysis classifiers that are especially useful for when applied to
small-sample, high-dimensional data sets.

## Installation

You can install the stable version on [CRAN](https://cran.r-project.org/package=sparsediscrim):

```r
install.packages('sparsediscrim', dependencies = TRUE)
```

If you prefer to download the latest version, instead type:

```r
library(devtools)
install_github('ramhiser/sparsediscrim')
```

## Classifiers

The `sparsediscrim` package features the following classifier (the R function
is included within parentheses):

* [High-Dimensional Regularized Discriminant Analysis](https://arxiv.org/abs/1602.01182) (`rda_high_dim`) from Ramey et al. (2015)

The `sparsediscrim` package also includes a variety of additional classifiers
intended for small-sample, high-dimensional data sets. These include:

| Classifier                                                    | Author                                                                                             | R Function |
|---------------------------------------------------------------|----------------------------------------------------------------------------------------------------|------------|
| Diagonal Linear Discriminant Analysis                         | [Dudoit et al. (2002)](http://www.tandfonline.com/doi/abs/10.1198/016214502753479248)              | `lda_diag`     |
| Diagonal Quadratic Discriminant Analysis                      | [Dudoit et al. (2002)](http://www.tandfonline.com/doi/abs/10.1198/016214502753479248)              | `qda_diag`     |
| Shrinkage-based Diagonal Linear Discriminant Analysis         | [Pang et al. (2009)](http://onlinelibrary.wiley.com/doi/10.1111/j.1541-0420.2009.01200.x/abstract) | `sdlda`    |
| Shrinkage-based Diagonal Quadratic Discriminant Analysis      | [Pang et al. (2009)](http://onlinelibrary.wiley.com/doi/10.1111/j.1541-0420.2009.01200.x/abstract) | `qda_shrink_cov`    |
| Shrinkage-mean-based Diagonal Linear Discriminant Analysis    | [Tong et al. (2012)](http://bioinformatics.oxfordjournals.org/content/28/4/531.long)               | `lda_shrink_mean`   |
| Shrinkage-mean-based Diagonal Quadratic Discriminant Analysis | [Tong et al. (2012)](http://bioinformatics.oxfordjournals.org/content/28/4/531.long)               | `qda_shrink_mean`   |
| Minimum Distance Empirical Bayesian Estimator (MDEB)          | [Srivistava and Kubokawa (2007)](http://www.utstat.utoronto.ca/~srivasta/exp1.pdf)                 | `lda_emp_bayes`     |
| Minimum Distance Rule using Modified Empirical Bayes (MDMEB)  | [Srivistava and Kubokawa (2007)](http://www.utstat.utoronto.ca/~srivasta/exp1.pdf)                 | `lda_emp_bayes_eigen`    |
| Minimum Distance Rule using Moore-Penrose Inverse (MDMP)      | [Srivistava and Kubokawa (2007)](http://www.utstat.utoronto.ca/~srivasta/exp1.pdf)                 | `lda_eigen`     |

We also include modifications to Linear Discriminant Analysis (LDA) with
regularized covariance-matrix estimators:

* Moore-Penrose Pseudo-Inverse (`lda_pseudo`)
* Schafer-Strimmer estimator (`lda_schafer`)
* Thomaz-Kitani-Gillies estimator (`lda_thomaz`)

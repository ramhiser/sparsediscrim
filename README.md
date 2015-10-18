# sparsediscrim

The R package `sparsediscrim` provides a collection of sparse and regularized discriminant
analysis classifiers that are especially useful for when applied to
small-sample, high-dimensional data sets.

## Installation

You can install the stable version on [CRAN](http://cran.r-project.org/package=sparsediscrim):

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

* High-Dimensional Regularized Discriminant Analysis (`hdrda`) from Ramey et al. (2015)

The `sparsediscrim` package also includes a variety of additional classifiers
intended for small-sample, high-dimensional data sets. These include:

| Classifier                                                    | Author                         | R Function |
|---------------------------------------------------------------|--------------------------------|------------|
| Diagonal Linear Discriminant Analysis                         | Dudoit et al. (2002)           | `dlda`     |
| Diagonal Quadratic Discriminant Analysis                      | Dudoit et al. (2002)           | `dqda`     |
| Shrinkage-based Diagonal Linear Discriminant Analysis         | Pang et al. (2009)             | `sdlda`    |
| Shrinkage-based Diagonal Quadratic Discriminant Analysis      | Pang et al. (2009)             | `sdqda`    |
| Shrinkage-mean-based Diagonal Linear Discriminant Analysis    | Tong et al. (2012)             | `smdlda`   |
| Shrinkage-mean-based Diagonal Quadratic Discriminant Analysis | Tong et al. (2012)             | `smdqda`   |
| Minimum Distance Empirical Bayesian Estimator (MDEB)          | Srivistava and Kubokawa (2007) | `mdeb`     |
| Minimum Distance Rule using Modified Empirical Bayes (MDMEB)  | Srivistava and Kubokawa (2007) | `mdmeb`    |
| Minimum Distance Rule using Moore-Penrose Inverse (MDMP)      | Srivistava and Kubokawa (2007) | `mdmp`     |

We also include modifications to Linear Discriminant Analysis (LDA) with
regularized covariance-matrix estimators:

* Moore-Penrose Pseudo-Inverse (`lda_pseudo`)
* Schafer-Strimmer estimator (`lda_schafer`)
* Thomaz-Kitani-Gillies estimator (`lda_thomaz`)

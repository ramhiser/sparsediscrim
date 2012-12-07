# sparsediscrim

The R package `sparsediscrim` provides a collection of sparse discriminant analysis methods that are especially useful for applying supervised learning to small-sample, high-dimensional data sets. We include implementations of the following classifiers (the R functions are included in parentheses):

* SimDiag (`simdiag`)
* Diagonal Linear Discriminant Analysis (`dlda`)
* Diagonal Quadratic Discriminant Analysis (`dqda`)
* Shrinkage-based Diagonal Linear Discriminant Analysis (`sdlda`)
* Shrinkage-based Diagonal Quadratic Discriminant Analysis (`sdqda`)
* Regularized Shrinkage-based Diagonal Discriminant Analysis (`rsdda`)


## Installation

You can install the stable version on [CRAN](http://cran.r-project.org/package=sparsediscrim):

```r
install.packages('sparsediscrim', dependencies = TRUE)
```

If you prefer to download the latest version, instead type:

```r
library(devtools)
install_github('sparsediscrim', 'ramey')
```

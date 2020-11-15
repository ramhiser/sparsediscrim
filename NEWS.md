# sparsediscrim 0.2.5

* Fixed issue where classifiers had errors when only 1 feature used. #41
* Added arXiv reference for HDRDA to DESCRIPTION #40
* Added doi references for other classifiers to DESCRIPTION #40

# sparsediscrim 0.2.4

* Added "Suggests: testthat, caret" to DESCRIPTION to pass `R CMD CHECK`. #38
* Added [arXiv reference](https://arxiv.org/abs/1602.01182) to HDRDA paper

# sparsediscrim 0.2.3

* Fixed an issue with the HDRDA classifier's `predict` function. The posterior
  probabilities did not sum to 1 because they were unnormalized. #34

* Fixed another issue with the HDRDA classifier's `predict` function, where the
  class names were incorrect when predicting a single observation. #34

* Improved docs throughout the package to pass `R CMD CHECK`. #35

# sparsediscrim 0.2.2

* The `predict` function now returns posterior-probability estimates for each
  classifier.

* The object returned by `cv_hdrda()` can be plotted. A heatmap is produced
  using `ggplot2` to illustrate the cross-validation error rate for each
  tuning-parameter pair considered.

* The `predict` function for the HDRDA classifier is now substantially faster
  when classifying a large number of observations. #33

# sparsediscrim 0.2.1


* The cross-validation helper function `cv_hdrda()` for the HDRDA classifier now
  returns a trained classifier rather than the optimal model.

* `cv_hdrda()` also has an optional `verbose` argument to dump summary
  information while the cross-validation is running.

* Fixed issue with classifiers' documentation not appearing in help index. #26

* Better handling of HDRDA when its tuning parameters are both 0.

* Corrected calculation of W_k and Q_k in HDRDA classifier.

* Added unit tests for HDRDA.

* Can now specify population means in `generate_blockdiag()`.

* Added unit tests for `generate_blockdiag()`.

* Updated **man** docs with roxygen2 4.0.

* Added `log_determinant()` helper function to calculate the log-determinant of
  a matrix.

# sparsediscrim 0.2


* The High-Dimensional Regularized Discriminant Analysis (HDRDA) classifier from
  Ramey, Stein, and Young (2014) implemented in `hdrda()` has been revamped to
  improve its computational performance.

* `lda_pseudo()` is an implementation of Linear Discriminant Analysis (LDA) with
  the Moore-Penrose Pseudo-Inverse

* `lda_schafer()` is an implementation of Linear Discriminant Analysis (LDA)
  using the covariance matrix estimator from Schafer and Strimmer (2005)

* `lda_thomaz()` is an implementation of Linear Discriminant Analysis (LDA)
  using the covariance matrix estimator from Thomaz, Kitani, and Gillies (2006)

* `mdeb()` is an implementation of the Minimum Distance Empirical Bayesian
  Estimator (MDEB) classifier from Srivistava and Kubokawa (2007)

* `mdmeb()` is an implementation of the Minimum Distance Rule using Modified
  Empirical Bayes (MDMEB) classifier from Srivistava and Kubokawa (2007)

* `mdmp()` is an implementation of the Minimum Distance Rule using Moore-Penrose
  Inverse (MDMP) classifier from Srivistava and Kubokawa (2007)

* `smdlda()` is an implementation of the Shrinkage-mean-based Diagonal Linear
  Discriminant Analysis (SmDLDA) from Tong, Chen, and Zhao (2012)

* `smdqda()` is an implementation of the Shrinkage-mean-based Diagonal Quadratic
  Discriminant Analysis (SmDQDA) from Tong, Chen, and Zhao (2012)

* Added a summary function for `hdrda` classifiers

# sparsediscrim 0.1

* First version of the `sparsediscrim` package. With this package, we aim to
  provide a large collection of regularized and sparse discriminant analysis
  classifiers intended for high-dimensional classification.

* `hdrda()` is an implementation of the High-Dimensional Regularized
  Discriminant Analysis classifier from Ramey, Stein, and Young (2014).

* `dlda()` is an implementation of the Diagonal Linear Discriminant Analysis
  classifier from Dudoit, Fridlyand, and Speed (2002).

* `dqda()` is an implementation of the Diagonal Quadratic Discriminant Analysis
  classifier from Dudoit, Fridlyand, and Speed (2002).

* `sdlda()` is an implementation of the Shrinkage-based Diagonal Linear
  Discriminant Analysis classifier from Pang, Tong, and Zhao (2009).

* `sdqda()` is an implementation of the Shrinkage-based Diagonal Quadratic
  Discriminant Analysis classifier from Pang, Tong, and Zhao (2009).

* `generate_blockdiag()` generates random variates from K multivariate normal
  populations, where each class is generated with a constant mean vector and a
  covariance matrix consisting of block-diagonal autocorrelation matrices.

* `generate_intraclass()` generates random variates from K multivariate normal
  populations, where class is generated with a constant mean vector and an
  intraclass covariance matrix.

* `cv_partition()` randomly partitions data for cross-validation.

* `no_intercept()` removes the intercept term from a formula if it is included.

* `cov_mle()` computes the maximum likelihood estimator for the sample
  covariance matrix under the assumption of multivariate normality.

* `cov_pool()` computes the pooled maximum likelihood estimator for the common
  covariance matrix under the assumption of multivariate normality.

* `cov_eigen()` computes the eigenvalue decomposition of the maximum likelihood
  estimators of the covariance matrices for the given data matrix. We provide an
  option to calculate the eigenvalue decomposition using the Fast Singular Value
  Decomposition, which can greatly expedite the eigenvalue decomposition for
  very tall data (large n, small p) or very wide data (small n, large p).
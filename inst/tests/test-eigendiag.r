library('mvtnorm')

p <- 500
x1 <- rmvnorm(8, rep(0,p))
x2 <- rmvnorm(8, rep(1,p))

sig1 <- cov(x1)
sig2 <- cov(x2)

sig.list <- list(sig1 = sig1, sig2 = sig2)

Q <- eigendiag(sig.list)$Q

norm(Q %*% sig1 %*% t(Q) - diag(Q %*% sig1 %*% t(Q)), type = "F")
norm(Q %*% sig2 %*% t(Q) - diag(Q %*% sig2 %*% t(Q)), type = "F")


# Code is downloaded from: http://bioinformatics.med.yale.edu/rsdda/rsdda.htm
#.packageName <- "RSDDA"
RSDDA <- function(infile, genenum, sampsize, topgenes, numofruns) {

numofg <- genenum; myfilename <- infile;
library(MASS);

if (sampsize <= 5) {
stop("Warning you are trying to perform RSDDA with five or fewer subjects per group")
}

err.rate6 <- NULL

# We assume that the first column contains the class label.
# Sorts the training data frame by class.
#training.df <- training.df[order(training.df[,1]),]

#y1 <- training

  # class 1
  n1 <- round(sampsize*0.6)
  y1 <- rep(0,n1)
  testn1 <- sampsize - n1
  testy1 <- rep(0,testn1)

  # class 2
  n2 <- round(sampsize*0.6)
  y2 <- rep(1,n2)
  testn2 <- sampsize - n2
  testy2 <- rep(1,testn2)

  n <- n1 + n2
  testn <- testn1 + testn2
  y <- c(y1,y2)
  testy1 <- c(testy1,testy2)

  x <- as.matrix(read.table(myfilename, header=F))
 # x <- t(as.matrix(x))
  dimnames(x) <- NULL

  lambda <- seq(0,1, length=101)
  rsvar1mi.est <- rsvar2mi.est <- matrix(data=NA, nrow=101, ncol=topgenes, dimnames=NULL)
  err.rate.lamb <- matrix(data=NA, nrow=101, ncol=n1*n2, dimnames=NULL)
  lambvec <- NULL

for (sim in 1:numofruns) {
  print(sim)
  mysamindex1 <- sample(1:(n1+testn1),(n1+testn1))
  mysamindex2 <- sample(1:(n2+testn2),(n2+testn2))
  mysamindex2 <- mysamindex2+n1+testn1
  x1 <- x[mysamindex1[1:n1],]
  x2 <- x[mysamindex2[1:n2],]

  trainx <- t(rbind(x1,x2))

  # numofg controls the number of genes kept after examining the ratio of between sum of squares (BSS)
  #		and the within sum of squares (WSS)
  # We keep the "numofg" genes with the largest ratio of BSS to WSS.
  myscore <- 1:numofg
  for (i in 1:numofg) {
   rmean <- mean(trainx[i,])
   rc1mean <- mean(trainx[i,1:n1])
   rc2mean <- mean(trainx[i,((n1+1):n)])
   BSS <- (n1*(rc1mean - rmean)^2 + n2*(rc2mean - rmean)^2)
   WSS <- 0
   WSS <- WSS + sum((trainx[i,(1:n1)] - rc1mean)^2) + sum((trainx[i,((n1+1):n)] - rc2mean)^2)
   myscore[i] <- BSS/WSS
   }
  sorttrainx <- sort.list(myscore, decreasing=TRUE)
  sorttrainx <- sorttrainx[1:topgenes]

  x1 <- x[mysamindex1[1:n1],sorttrainx]
  x2 <- x[mysamindex2[1:n2],sorttrainx]
  testx1 <- x[mysamindex1[(n1+1):(n1+testn1)],sorttrainx]
  testx2 <- x[mysamindex2[(n2+1):(n2+testn2)],sorttrainx]

  truth <- truth2 <- NULL
  predict6 <- NULL
  predict.lamb <- matrix(data=NA, nrow=101, ncol=2, dimnames=NULL)

      mu1.est  <- apply(x1,2,mean)
      var1.est <- apply(x1,2,var)
      mu2.est  <- apply(x2,2,mean)
      var2.est <- apply(x2,2,var)

      mylen1 <- n1-1
      mylen2 <- n2-1
      var12.est <- (mylen1*var1.est + mylen2*var2.est)/(n-2)

      var1.shrink3 <- ishrink(var1.est,n1-1)
      var2.shrink3 <- ishrink(var2.est,n2-1)
      var12.shrink3 <- ishrink(var12.est,n-2)

    for (grid in 1:101) {    
       rsvar1mi.est[grid,] <- var1.shrink3^(1-lambda[grid])*var12.shrink3^lambda[grid]
       rsvar2mi.est[grid,] <- var2.shrink3^(1-lambda[grid])*var12.shrink3^lambda[grid]
    }

    numofcv <- n1*n1
    cvy <- cvy.pred <- NULL
    cvx.test <- matrix(data=NA, nrow=2, ncol=topgenes, dimnames=NULL)
    cvrsvar1mi.est <- cvrsvar2mi.est <- matrix(data=NA, nrow=101, ncol=topgenes, dimnames=NULL)

    tmpcv1 <- tmpcv2 <- NULL
 
    for (i in 1:n1) {
      tmpcv1 <- c(tmpcv1,rep(i,n1))
    }
    tmpcv2 <- rep(1:n1,n1)

    for (cvsim in 1:numofcv) {
        cvx1 <- x1[-tmpcv1[cvsim],]
        cvx2 <- x2[-tmpcv2[cvsim],]
    

        cvmu1.est  <- apply(cvx1,2,mean)
        cvvar1.est <- apply(cvx1,2,var)
        cvmu2.est  <- apply(cvx2,2,mean)
        cvvar2.est <- apply(cvx2,2,var) 
        mylen1 <- n1-2
        mylen2 <- n2-2
        cvvar12.est <- (mylen1*cvvar1.est + mylen2*cvvar2.est)/(n-4)
 
         cvvar1.shrink3 <- ishrink(cvvar1.est,n1-2)
         cvvar2.shrink3 <- ishrink(cvvar2.est,n2-2)
 
        cvvar12.shrink3 <- ishrink(cvvar12.est,n-4)
 
      for (grid in 1:101) {
         cvrsvar1mi.est[grid,] <- cvvar1.shrink3^(1-lambda[grid])*cvvar12.shrink3^lambda[grid]
         cvrsvar2mi.est[grid,] <- cvvar2.shrink3^(1-lambda[grid])*cvvar12.shrink3^lambda[grid]
      }

      cvy[1] <- 0
      cvy[2] <- 1
      truth2 <- NULL
      cvx.test[1,] <- x1[tmpcv1[cvsim],]
      cvx.test[2,] <- x2[tmpcv2[cvsim],]
      for (u in 1:2) {
        cvx.pred <- cvx.test[u]
        cvy.pred <- cvy[u]
        truth2 <- c(truth2,cvy.pred)
       for (grid in 1:101) {    
          D1 <- diss(cvx.pred,cvmu1.est,cvrsvar1mi.est[grid,], n1-1, n-2)
          D2 <- diss(cvx.pred,cvmu2.est,cvrsvar2mi.est[grid,], n2-1, n-2)
          predict.lamb[grid,u] <- 1-sum(D1<=D2)
       }
      }

      for (grid in 1:101) {
       err.rate.lamb[grid,cvsim] <- sum(predict.lamb[grid,]!=truth2)/2
      }

    }

      cv.err.rate <- NULL
      for (grid in 1:101) {
        cv.err.rate[grid] <- mean(err.rate.lamb[grid,])
      }

  err.lamb.ind <- which.min(cv.err.rate)
  err.lamb.ind2 <- 102-which.min(rev(cv.err.rate))
  lambnum <- floor((err.lamb.ind+err.lamb.ind2)/2)

  if (err.lamb.ind2 == 101) {
    lambnum <- err.lamb.ind2
    if (err.lamb.ind < 50) {
      lambnum <- floor((err.lamb.ind+err.lamb.ind2)/2)
    }
  }

  lambvec[sim] <-  lambda[lambnum]

  for (v in 1:testn1) {
      x.predict <- testx1[v,]
      y.predict <- testy1[v]
      truth <- c(truth,y.predict)
       D1 <- diss(x.predict,mu1.est,rsvar1mi.est[lambnum,], n1, n)
       D2 <- diss(x.predict,mu2.est,rsvar2mi.est[lambnum,], n2, n)
      predict6 <- c(predict6,1-sum(D1<=D2))
  }

  for (v in 1:testn2) {
      x.predict <- testx2[v,]
      y.predict <- testy2[v]
      truth <- c(truth,y.predict)
       D1 <- diss(x.predict,mu1.est,rsvar1mi.est[lambnum,], n1, n)
       D2 <- diss(x.predict,mu2.est,rsvar2mi.est[lambnum,], n2, n)
      predict6 <- c(predict6,1-sum(D1<=D2))
  }

    err.rate6 <- c(err.rate6,sum(predict6!=truth)/testn)

}

print("Mean Misclassification of RSDDA:");
print(mean(err.rate6));
print("SE of Misclassification of RSDDA:");
print(sd(err.rate6));
}
diss <- function(x,mu,sigma2,n,n.total){
  D <- 0
  G <- length(x)
  for (i in 1:G) {
   D <- D + (x[i]-mu[i])^2/sigma2[i] + log(sigma2[i])
   D <- D - 2*log(n/n.total)  
  }
  return(D)
}
ishrink <- function(s2, df,trun=0){
    len <- length(s2)

    C <- 2/df * (gamma(df/2)/gamma(df/2-1/len))^len

    spool.g <- 1
    power <- 1/len
    for (i in 1:len) 
    spool.g <- spool.g*(s2[i])^power
    G <- length(s2)
    G1 <- trunc(G*trun)
    G2 <- G-G1
    s2.sort1 <- sort(s2)
    s2.trunc1 <- s2.sort1[(G1+1):G2]

    variance <- s2.trunc1
    n <- df+1
    counts <- 1000

    alpha <- (0:counts)/counts
    len <- length(variance) #len=G
     vpool <- 1
     power <- 1/len
     for (i in 1:len) { vpool <- vpool*(variance[i])^power }
     k <- (n-1)/2
     C <- (gamma(k)/gamma(k-1/len))^len /k 
     tmp1 <- (gamma(k-alpha/len)/gamma(k))^(len-1)*gamma(k-(1-alpha+alpha/len))/gamma(k)
     tmp2 <- 0
     for (i in 1:len) { tmp2 <- tmp2 + (variance[i])^alpha }
     f <- C^alpha *k *((n-3)/(n-1))^(1-alpha) * tmp1 * vpool^(-alpha) * tmp2 - alpha*len*log(C) - (1-alpha)*len*log((n-3)/(n-1))
     A <- f
    min.order <- order(A)[1]
    alpha1 <- alpha[min.order]
    inv.s2.shrink1 <- (C/spool.g)^alpha1 * ((df-2)/(df*s2))^(1-alpha1)
    s2.shrink1 <- 1/inv.s2.shrink1

    s2.sort2 <- sort(s2.shrink1)
    s2.trunc2 <- s2.sort2[(G1+1):G2]

    variance <- s2.trunc2
    n <- df+1
    counts <- 1000

    alpha <- (0:counts)/counts
    len <- length(variance) #len=G
     vpool <- 1
     power <- 1/len
     for (i in 1:len) { vpool <- vpool*(variance[i])^power }
     k <- (n-1)/2
     C <- (gamma(k)/gamma(k-1/len))^len /k 
     tmp1 <- (gamma(k-alpha/len)/gamma(k))^(len-1)*gamma(k-(1-alpha+alpha/len))/gamma(k)
     tmp2 <- 0
     for (i in 1:len) { tmp2 <- tmp2 + (variance[i])^alpha }
     f <- C^alpha *k *((n-3)/(n-1))^(1-alpha) * tmp1 * vpool^(-alpha) * tmp2 - alpha*len*log(C) - (1-alpha)*len*log((n-3)/(n-1))
     A <- f
    min.order <- order(A)[1]
    alpha2 <- alpha[min.order]

    inv.s2.shrink2 <- (C/spool.g)^alpha2 * ((df-2)/(df*s2))^(1-alpha2)
    s2.shrink2 <- 1/inv.s2.shrink2

    return(s2.shrink2)
}

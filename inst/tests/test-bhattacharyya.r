library('sparsediscrim')
library('datamicroarray')
library('multtest')

data('chiaretti')
# table(chiaretti$y)

#ALL1/AF4  BCR/ABL E2A/PBX1      NEG   NUP-98  p15/p16
#      10       37        5       74        1        1

# We restrict the Chiaretti data set to the two classes with the largest sample
# sizes: 'NEG' and 'BCR/ABL'
two_classes <- c("NEG", "BCR/ABL")
idx <- which(chiaretti$y %in% two_classes)
chiaretti$x <- chiaretti$x[idx, ]
chiaretti$y <- factor(chiaretti$y[idx])


# Alon Data Set
# data('alon', package = 'datamicroarray')
# x <- alon$x
# y <- alon$y
# rm(alon)

# Singh
# data('singh', package = 'datamicroarray')
# x <- singh$x
# y <- singh$y
# rm(singh)


q <- 500

# Applies Dudoit et al. variable selection to the data set to reduce to 'q' genes
F_stat <- with(chiaretti, mt.teststat(t(x), y, test = "f"))
top_q_genes <- rev(order(F_stat))[seq_len(q)]
x <- chiaretti$x[, top_q_genes]
y <- chiaretti$y

simdiag_out <- simdiag_cov(x = x, y = y)
bhatta_out <-  with(simdiag_out, bhatta_simdiag(x = x, y = y))
bhatta_pool_out <-  with(simdiag_out, bhatta_simdiag(x = x, y = y, pool_cov = TRUE))

# Plot cumulative proportions
with(bhatta_out, plot(seq_along(cumprop), cumprop, type = "l"))
quartz()
with(bhatta_pool_out, plot(seq_along(cumprop), cumprop, type = "l"))

# Plot distances -- could be used for elbow criterion
with(bhatta_out, plot(seq_along(sort(dist)), sort(dist, dec = TRUE), type = "l"))
quartz()
with(bhatta_pool_out, plot(seq_along(sort(dist)), sort(dist, dec = TRUE), type = "l"))

train <- sample(seq_along(y), 75)
simdiag_out <- simdiag(x = x[train, ], y = y[train])
mean(predict(simdiag_out, x[-train, ])$class != y[-train])

simdiag_out <- simdiag(x = x[train, ], y = y[train], pool_cov = TRUE)
mean(predict(simdiag_out, x[-train, ])$class != y[-train])

simdiag_out <- simdiag(x = x[train, ], y = y[train], classifier = "quadratic")
mean(predict(simdiag_out, x[-train, ])$class != y[-train])

simdiag_out <- simdiag(x = x[train, ], y = y[train], classifier = "quadratic", pool_cov = TRUE)
mean(predict(simdiag_out, x[-train, ])$class != y[-train])


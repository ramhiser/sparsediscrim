library('diagdiscrim')
library('mvtnorm')
library('testthat')

# These unit tests are an attempt to replicate the results found in Simulation A of Pang et al. (2009).
# Results used in these unit tests come from the supplementary web materials found here:
# http://biometrics.tibs.org/datasets/071021M_FINAL_supp_materials.pdf

context("Pang et al.'s (2009) Simulation A")

generate_data_A <- function(n_k, p) {
	mu1 <- 0
	mu2 <- 0.5
	
	test_k <- 2 * n_k
	
	x1 <- rmvnorm(n_k + test_k, rep(mu1, p))
	x2 <- rmvnorm(n_k + test_k, rep(mu2, p))
	df <- rbind.data.frame(cbind(1, x1), cbind(2, x2))
	names(df) <- c("labels", paste("X", seq_len(p), sep = ""))
	df$labels <- factor(df$labels)
	
	which_are_training <- sapply(levels(df$labels), function(label) sample(which(df$labels == label), n_k))
	which_are_training <- as.vector(which_are_training)

	list(training = df[which_are_training,], test = df[-which_are_training,])
}

tol <- 0.04
B <- 500

n_k <- c(4, 5, 8, 10, 15)
p <- c(30, 50, 100)

sim_configs <- expand.grid(n_k = n_k, p = p)

test_that("DLDA's matches performance of Simulation A given in Pang (2009)", {
	set.seed(42)

	message("Pang et al. (2009) Simulation A - DLDA")
	dlda_error_A <<- ddply(sim_configs, .(n_k, p), function(sim_config) {
		print(sim_config)
		error_rates <- replicate(B, {
			data <- generate_data_A(n_k = sim_config$n_k, p = sim_config$p)
			train_df <- data$training
			test_df <- data$test

			dlda_out <- dlda(train_df)
			predictions <- predict(dlda_out, data.matrix(test_df[, -1]))
			test_error <- mean(predictions != test_df$labels)
			test_error
		})
		mean(error_rates)
	})
	names(dlda_error_A) <<- c("n_k", "p", "error")
	
	# n_k = 4. p = 30, 50, 100
	expect_that(dlda_error_A[1,3], equals(0.2806, tolerance = tol, scale = 1))
	expect_that(dlda_error_A[2,3], equals(0.2300, tolerance = tol, scale = 1))
	expect_that(dlda_error_A[3,3], equals(0.1475, tolerance = tol, scale = 1))
	
	# n_k = 5. p = 30, 50, 100
	expect_that(dlda_error_A[4,3], equals(0.2415, tolerance = tol, scale = 1))
	expect_that(dlda_error_A[5,3], equals(0.2015, tolerance = tol, scale = 1))
	expect_that(dlda_error_A[6,3], equals(0.1160, tolerance = tol, scale = 1))
	
	# n_k = 8. p = 30, 50, 100
	expect_that(dlda_error_A[7,3], equals(0.1928, tolerance = tol, scale = 1))
	expect_that(dlda_error_A[8,3], equals(0.1294, tolerance = tol, scale = 1))
	expect_that(dlda_error_A[9,3], equals(0.0619, tolerance = tol, scale = 1))
	
	# n_k = 10. p = 30, 50, 100
	expect_that(dlda_error_A[10,3], equals(0.1803, tolerance = tol, scale = 1))
	expect_that(dlda_error_A[11,3], equals(0.1030, tolerance = tol, scale = 1))
	expect_that(dlda_error_A[12,3], equals(0.0470, tolerance = tol, scale = 1))
	
	# n_k = 15. p = 30, 50, 100
	expect_that(dlda_error_A[13,3], equals(0.1505, tolerance = tol, scale = 1))
	expect_that(dlda_error_A[14,3], equals(0.0858, tolerance = tol, scale = 1))
	expect_that(dlda_error_A[15,3], equals(0.0272, tolerance = tol, scale = 1))
})

test_that("DQDA's matches performance of Simulation A given in Pang (2009)", {
	set.seed(42)

	message("Pang et al. (2009) Simulation A - DQDA")
	dqda_error_A <<- ddply(sim_configs, .(n_k, p), function(sim_config) {
		print(sim_config)
		error_rates <- replicate(B, {
			data <- generate_data_A(n_k = sim_config$n_k, p = sim_config$p)
			train_df <- data$training
			test_df <- data$test

			dqda_out <- dqda(train_df)
			predictions <- predict(dqda_out, data.matrix(test_df[, -1]))
			test_error <- mean(predictions != test_df$labels)
			test_error
		})
		mean(error_rates)
	})
	names(dqda_error_A) <<- c("n_k", "p", "error")
	
	# n_k = 4. p = 30, 50, 100
	expect_that(dqda_error_A[1,3], equals(0.3844, tolerance = tol, scale = 1))
	expect_that(dqda_error_A[2,3], equals(0.3631, tolerance = tol, scale = 1))
	expect_that(dqda_error_A[3,3], equals(0.3144, tolerance = tol, scale = 1))
	
	# n_k = 5. p = 30, 50, 100
	expect_that(dqda_error_A[4,3], equals(0.3410, tolerance = tol, scale = 1))
	expect_that(dqda_error_A[5,3], equals(0.3240, tolerance = tol, scale = 1))
	expect_that(dqda_error_A[6,3], equals(0.2710, tolerance = tol, scale = 1))
	
	# n_k = 8. p = 30, 50, 100
	expect_that(dqda_error_A[7,3], equals(0.2703, tolerance = tol, scale = 1))
	expect_that(dqda_error_A[8,3], equals(0.2194, tolerance = tol, scale = 1))
	expect_that(dqda_error_A[9,3], equals(0.1522, tolerance = tol, scale = 1))
	
	# n_k = 10. p = 30, 50, 100
	expect_that(dqda_error_A[10,3], equals(0.2443, tolerance = tol, scale = 1))
	expect_that(dqda_error_A[11,3], equals(0.1735, tolerance = tol, scale = 1))
	expect_that(dqda_error_A[12,3], equals(0.1078, tolerance = tol, scale = 1))
	
	# n_k = 15. p = 30, 50, 100
	expect_that(dqda_error_A[13,3], equals(0.2013, tolerance = tol, scale = 1))
	expect_that(dqda_error_A[14,3], equals(0.1418, tolerance = tol, scale = 1))
	expect_that(dqda_error_A[15,3], equals(0.0595, tolerance = tol, scale = 1))
})

test_that("SDLDA's matches performance of Simulation A given in Pang (2009)", {
	set.seed(42)

	message("Pang et al. (2009) Simulation A - SDLDA")
	sdlda_error_A <<- ddply(sim_configs, .(n_k, p), function(sim_config) {
		print(sim_config)
		error_rates <- replicate(B, {
			data <- generate_data_A(n_k = sim_config$n_k, p = sim_config$p)
			train_df <- data$training
			test_df <- data$test

			sdlda_out <- sdlda(train_df, num_alphas = 101)
			predictions <- predict(sdlda_out, data.matrix(test_df[, -1]))
			test_error <- mean(predictions != test_df$labels)
			test_error
		})
		mean(error_rates)
	})
	names(sdlda_error_A) <<- c("n_k", "p", "error")
	
	# n_k = 4. p = 30, 50, 100
	expect_that(sdlda_error_A[1,3], equals(0.2388, tolerance = tol, scale = 1))
	expect_that(sdlda_error_A[2,3], equals(0.1775, tolerance = tol, scale = 1))
	expect_that(sdlda_error_A[3,3], equals(0.0881, tolerance = tol, scale = 1))
	
	# n_k = 5. p = 30, 50, 100
	expect_that(sdlda_error_A[4,3], equals(0.2260, tolerance = tol, scale = 1))
	expect_that(sdlda_error_A[5,3], equals(0.1705, tolerance = tol, scale = 1))
	expect_that(sdlda_error_A[6,3], equals(0.0695, tolerance = tol, scale = 1))
	
	# n_k = 8. p = 30, 50, 100
	expect_that(sdlda_error_A[7,3], equals(0.1791, tolerance = tol, scale = 1))
	expect_that(sdlda_error_A[8,3], equals(0.1069, tolerance = tol, scale = 1))
	expect_that(sdlda_error_A[9,3], equals(0.0481, tolerance = tol, scale = 1))
	
	# n_k = 10. p = 30, 50, 100
	expect_that(sdlda_error_A[10,3], equals(0.1633, tolerance = tol, scale = 1))
	expect_that(sdlda_error_A[11,3], equals(0.0953, tolerance = tol, scale = 1))
	expect_that(sdlda_error_A[12,3], equals(0.0400, tolerance = tol, scale = 1))
	
	# n_k = 15. p = 30, 50, 100
	expect_that(sdlda_error_A[13,3], equals(0.1423, tolerance = tol, scale = 1))
	expect_that(sdlda_error_A[14,3], equals(0.0800, tolerance = tol, scale = 1))
	expect_that(sdlda_error_A[15,3], equals(0.0230, tolerance = tol, scale = 1))
})

test_that("SDQDA's matches performance of Simulation A given in Pang (2009)", {
	set.seed(42)

	message("Pang et al. (2009) Simulation A - SDQDA")
	sdqda_error_A <<- ddply(sim_configs, .(n_k, p), function(sim_config) {
		print(sim_config)
		error_rates <- replicate(B, {
			data <- generate_data_A(n_k = sim_config$n_k, p = sim_config$p)
			train_df <- data$training
			test_df <- data$test

			sdqda_out <- sdqda(train_df, num_alphas = 101)
			predictions <- predict(sdqda_out, data.matrix(test_df[, -1]))
			test_error <- mean(predictions != test_df$labels)
			test_error
		})
		mean(error_rates)
	})
	names(sdqda_error_A) <<- c("n_k", "p", "error")
	
	# n_k = 4. p = 30, 50, 100
	expect_that(sdqda_error_A[1,3], equals(0.2475, tolerance = tol, scale = 1))
	expect_that(sdqda_error_A[2,3], equals(0.1863, tolerance = tol, scale = 1))
	expect_that(sdqda_error_A[3,3], equals(0.0856, tolerance = tol, scale = 1))
	
	# n_k = 5. p = 30, 50, 100
	expect_that(sdqda_error_A[4,3], equals(0.2395, tolerance = tol, scale = 1))
	expect_that(sdqda_error_A[5,3], equals(0.1720, tolerance = tol, scale = 1))
	expect_that(sdqda_error_A[6,3], equals(0.0760, tolerance = tol, scale = 1))
	
	# n_k = 8. p = 30, 50, 100
	expect_that(sdqda_error_A[7,3], equals(0.1784, tolerance = tol, scale = 1))
	expect_that(sdqda_error_A[8,3], equals(0.1109, tolerance = tol, scale = 1))
	expect_that(sdqda_error_A[9,3], equals(0.0491, tolerance = tol, scale = 1))
	
	# n_k = 10. p = 30, 50, 100
	expect_that(sdqda_error_A[10,3], equals(0.1663, tolerance = tol, scale = 1))
	expect_that(sdqda_error_A[11,3], equals(0.0983, tolerance = tol, scale = 1))
	expect_that(sdqda_error_A[12,3], equals(0.0400, tolerance = tol, scale = 1))
	
	# n_k = 15. p = 30, 50, 100
	expect_that(sdqda_error_A[13,3], equals(0.1480, tolerance = tol, scale = 1))
	expect_that(sdqda_error_A[14,3], equals(0.0812, tolerance = tol, scale = 1))
	expect_that(sdqda_error_A[15,3], equals(0.0255, tolerance = tol, scale = 1))
})

test_that("RSDDA's matches performance of Simulation A given in Pang (2009)", {
	set.seed(42)

	message("Pang et al. (2009) Simulation A - RSDDA")
	rsdda_error_A <<- ddply(sim_configs, .(n_k, p), function(sim_config) {
		print(sim_config)
		error_rates <- replicate(B, {
			data <- generate_data_A(n_k = sim_config$n_k, p = sim_config$p)
			train_df <- data$training
			test_df <- data$test

			rsdda_out <- rsdda(train_df, num_alphas = 101)
			predictions <- predict(rsdda_out, data.matrix(test_df[, -1]), num_lambdas = 101)
			test_error <- mean(predictions != test_df$labels)
			test_error
		})
		mean(error_rates)
	})
	names(rsdda_error_A) <<- c("n_k", "p", "error")
	
	# n_k = 4. p = 30, 50, 100
	expect_that(rsdda_error_A[1,3], equals(0.2388, tolerance = tol, scale = 1))
	expect_that(rsdda_error_A[2,3], equals(0.1775, tolerance = tol, scale = 1))
	expect_that(rsdda_error_A[3,3], equals(0.0856, tolerance = tol, scale = 1))
	
	# n_k = 5. p = 30, 50, 100
	expect_that(rsdda_error_A[4,3], equals(0.2285, tolerance = tol, scale = 1))
	expect_that(rsdda_error_A[5,3], equals(0.1690, tolerance = tol, scale = 1))
	expect_that(rsdda_error_A[6,3], equals(0.0695, tolerance = tol, scale = 1))
	
	# n_k = 8. p = 30, 50, 100
	expect_that(rsdda_error_A[7,3], equals(0.1769, tolerance = tol, scale = 1))
	expect_that(rsdda_error_A[8,3], equals(0.1078, tolerance = tol, scale = 1))
	expect_that(rsdda_error_A[9,3], equals(0.0478, tolerance = tol, scale = 1))
	
	# n_k = 10. p = 30, 50, 100
	expect_that(rsdda_error_A[10,3], equals(0.1638, tolerance = tol, scale = 1))
	expect_that(rsdda_error_A[11,3], equals(0.0968, tolerance = tol, scale = 1))
	expect_that(rsdda_error_A[12,3], equals(0.0405, tolerance = tol, scale = 1))
	
	# n_k = 15. p = 30, 50, 100
	expect_that(rsdda_error_A[13,3], equals(0.1425, tolerance = tol, scale = 1))
	expect_that(rsdda_error_A[14,3], equals(0.0802, tolerance = tol, scale = 1))
	expect_that(rsdda_error_A[15,3], equals(0.0237, tolerance = tol, scale = 1))
})

library('diagdiscrim')
library('mvtnorm')
library('testthat')

B <- 3
rsdda_error_A <<- ddply(sim_configs, .(n_k, p), function(sim_config) {
	print(sim_config)
	error_rates <- replicate(B, {
		data <- generate_data_A(n_k = sim_config$n_k, p = sim_config$p)
		train_df <- data$training
		test_df <- data$test
		message("Training classifier")
		rsdda_out <- rsdda(train_df, num_alphas = 101)
		message("Classifier trained!")
		message("Making predictions")
		predictions <- predict(rsdda_out, data.matrix(test_df[, -1]), num_lambdas = 101)
		message("Made predictions")
		print(predictions)
		test_error <- mean(predictions != test_df$labels)
		print(test_error)
		test_error
	})
	mean(error_rates)
})
names(rsdda_error_A) <<- c("n_k", "p", "error")
library('diagdiscrim')
library('mvtnorm')
library('testthat')

# These unit tests are an attempt to replicate the results found in Simulation B of Pang et al. (2009).
# Results used in these unit tests come from the supplementary web materials found here:
# http://biometrics.tibs.org/datasets/071021M_FINBL_supp_materials.pdf

context("Pang et al.'s (2009) Simulation B")

generate_data_B <- function(n_k, p) {
	mu1 <- 0
	mu2 <- 1
	
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

tolerance <- 0.02
B <- 500

n_k <- c(4, 5, 8, 10, 15)
p <- c(30)

sim_configs <- expand.grid(n_k = n_k, p = p)

test_that("DLDA's matches performance of Simulation B given in Pang (2009)", {
	set.seed(42)

	message("Pang et al. (2009) Simulation B - DLDA")
	dlda_error_B <<- ddply(sim_configs, .(n_k, p), function(sim_config) {
		print(sim_config)
		error_rates <- replicate(B, {
			data <- generate_data_B(n_k = sim_config$n_k, p = sim_config$p)
			train_df <- data$training
			test_df <- data$test

			dlda_out <- dlda(train_df)
			predictions <- predict(dlda_out, data.matrix(test_df[, -1]))
			test_error <- mean(predictions != test_df$labels)
			test_error
		})
		mean(error_rates)
	})
	names(dlda_error_B) <- c("n_k", "p", "error")
	
	# n_k = 4. p = 30
	expect_that(dlda_error_B$error[1,], equals(0.0575, tolerance = tol))
	
	# n_k = 5. p = 30
	expect_that(dlda_error_B$error[2,], equals(0.0300, tolerance = tol))
	
	# n_k = 8. p = 30
	expect_that(dlda_error_B$error[3,], equals(0.0144, tolerance = tol))
	
	# n_k = 10. p = 30
	expect_that(dlda_error_B$error[4,], equals(0.0120, tolerance = tol))
	
	# n_k = 15. p = 30
	expect_that(dlda_error_B$error[5,], equals(0.0095, tolerance = tol))
})

test_that("DQDA's matches performance of Simulation B given in Pang (2009)", {
	set.seed(42)

	message("Pang et al. (2009) Simulation B - DQDA")
	dqda_error_B <<- ddply(sim_configs, .(n_k, p), function(sim_config) {
		print(sim_config)
		error_rates <- replicate(B, {
			data <- generate_data_B(n_k = sim_config$n_k, p = sim_config$p)
			train_df <- data$training
			test_df <- data$test

			dqda_out <- dqda(train_df)
			predictions <- predict(dqda_out, data.matrix(test_df[, -1]))
			test_error <- mean(predictions != test_df$labels)
			test_error
		})
		mean(error_rates)
	})
	names(dqda_error_B) <- c("n_k", "p", "error")
	
	# n_k = 4. p = 30
	expect_that(dqda_error_B$error[1,], equals(0.2025, tolerance = tol))
	
	# n_k = 5. p = 30
	expect_that(dqda_error_B$error[2,], equals(0.0995, tolerance = tol))
	
	# n_k = 8. p = 30
	expect_that(dqda_error_B$error[3,], equals(0.0416, tolerance = tol))
	
	# n_k = 10. p = 30
	expect_that(dqda_error_B$error[4,], equals(0.0285, tolerance = tol))
	
	# n_k = 15. p = 30
	expect_that(dqda_error_B$error[5,], equals(0.0158, tolerance = tol))
})

test_that("SDLDA's matches performance of Simulation B given in Pang (2009)", {
	set.seed(42)

	message("Pang et al. (2009) Simulation B - SDLDA")
	sdlda_error_B <<- ddply(sim_configs, .(n_k, p), function(sim_config) {
		print(sim_config)
		error_rates <- replicate(B, {
			data <- generate_data_B(n_k = sim_config$n_k, p = sim_config$p)
			train_df <- data$training
			test_df <- data$test

			sdlda_out <- sdlda(train_df, num.alphas = 101)
			predictions <- predict(sdlda_out, data.matrix(test_df[, -1]))
			test_error <- mean(predictions != test_df$labels)
			test_error
		})
		mean(error_rates)
	})
	names(sdlda_error_B) <- c("n_k", "p", "error")
	
	# n_k = 4. p = 30
	expect_that(sdlda_error_B$error[1,], equals(0.0244, tolerance = tol))
	
	# n_k = 5. p = 30
	expect_that(sdlda_error_B$error[2,], equals(0.0175, tolerance = tol))
	
	# n_k = 8. p = 30
	expect_that(sdlda_error_B$error[3,], equals(0.0081, tolerance = tol))
	
	# n_k = 10. p = 30
	expect_that(sdlda_error_B$error[4,], equals(0.0073, tolerance = tol))
	
	# n_k = 15. p = 30
	expect_that(sdlda_error_B$error[5,], equals(0.0072, tolerance = tol))
})

test_that("SDQDA's matches performance of Simulation B given in Pang (2009)", {
	set.seed(42)

	message("Pang et al. (2009) Simulation B - SDQDA")
	sdqda_error_B <<- ddply(sim_configs, .(n_k, p), function(sim_config) {
		print(sim_config)
		error_rates <- replicate(B, {
			data <- generate_data_B(n_k = sim_config$n_k, p = sim_config$p)
			train_df <- data$training
			test_df <- data$test

			sdqda_out <- sdqda(train_df, num.alphas = 101)
			predictions <- predict(sdqda_out, data.matrix(test_df[, -1]))
			test_error <- mean(predictions != test_df$labels)
			test_error
		})
		mean(error_rates)
	})
	names(sdqda_error_B) <- c("n_k", "p", "error")
	
	# n_k = 4. p = 30
	expect_that(sdqda_error_B$error[1,], equals(0.0294, tolerance = tol))
	
	# n_k = 5. p = 30
	expect_that(sdqda_error_B$error[2,], equals(0.0190, tolerance = tol))
	
	# n_k = 8. p = 30
	expect_that(sdqda_error_B$error[3,], equals(0.0094, tolerance = tol))
	
	# n_k = 10. p = 30
	expect_that(sdqda_error_B$error[4,], equals(0.0085, tolerance = tol))
	
	# n_k = 15. p = 30
	expect_that(sdqda_error_B$error[5,], equals(0.0083, tolerance = tol))
})

test_that("RSDDA's matches performance of Simulation B given in Pang (2009)", {
	set.seed(42)

	message("Pang et al. (2009) Simulation B - RSDDA")
	rsdda_error_B <<- ddply(sim_configs, .(n_k, p), function(sim_config) {
		print(sim_config)
		error_rates <- replicate(B, {
			data <- generate_data_B(n_k = sim_config$n_k, p = sim_config$p)
			train_df <- data$training
			test_df <- data$test

			rsdda_out <- rsdda(train_df, num.alphas = 101)
			predictions <- predict(rsdda_out, data.matrix(test_df[, -1]), num.lambdas = 101)
			test_error <- mean(predictions != test_df$labels)
			test_error
		})
		mean(error_rates)
	})
	names(rsdda_error_B) <- c("n_k", "p", "error")
	
	# n_k = 4. p = 30
	expect_that(rsdda_error_B$error[1,], equals(0.0231, tolerance = tol))
	
	# n_k = 5. p = 30
	expect_that(rsdda_error_B$error[2,], equals(0.0160, tolerance = tol))
	
	# n_k = 8. p = 30
	expect_that(rsdda_error_B$error[3,], equals(0.0081, tolerance = tol))
	
	# n_k = 10. p = 30
	expect_that(rsdda_error_B$error[4,], equals(0.0075, tolerance = tol))
	
	# n_k = 15. p = 30
	expect_that(rsdda_error_B$error[5,], equals(0.0075, tolerance = tol))
})
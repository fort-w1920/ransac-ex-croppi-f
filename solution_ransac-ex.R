ransaclm <- function(formula, data, error_threshold, inlier_threshold, iterations,
                     fraction, strata) {

  # removing NA
  data <- na.omit(data)
  # assertions & check
  checkmate::assert_data_frame(data)
  checkmate::assert_logical(data$inlier, any.missing = FALSE)
  # n<p
  if ((ncol(data) - 2 - nrow(data)) > 0) { # ncoef -y -inlier -nobs >0 bzw. n < p
    stop("The model ceofficients are not uniquely defined. Please check your data.")
  }
  # dealing with factors
  for (i in seq_along(data)) {
    if (class(data[, i]) == "character") {
      data[, i] <- as.factor(data[, i])
      droplevels(data[, i])
    }
  }

  # define a starting model with an error & an empty consensus_set
  best_model <- NULL
  best_error <- 100

  for (i in 1:iterations) {

    # Step 1: fit a model with a random subset of hypothetical inliers, the
    #         remaining data are considered as hypothetical outliers
    hypo_inliers <- resample(data, strata = strata, fraction = fraction)
    hypo_outliers <- dplyr::setdiff(data, hypo_inliers)
    rownames(hypo_outliers) <- dplyr::setdiff(rownames(data), rownames(hypo_inliers))
    hypo_model <- get_lm(formula, hypo_inliers) # fit the model with the formula

    # Step2: use get_consistent_data to find the inliers of the model fitted in
    #        Step 1 and exclude outliers
    inliers <- get_consistent_data(
      hypo_inliers, hypo_outliers, hypo_model,
      error_threshold
    )

    # Step3: fit the model again with get_lm using the inliers found in Step2
    better_model <- get_lm(formula, inliers)
    # evaluate model fit by the sd of the residuals
    this_error <- summary(better_model)$sigma

    # Step 4:update best_model and the model error, if inlier_threshold is fine &
    #        the fit of better_model has improved
    if (nrow(inliers) > inlier_threshold & this_error <= best_error) {
      best_model <- better_model
      best_error <- this_error
    } else {
      best_model <- best_model
      best_error <- best_error
    }
  }

  # create a new column in the data set named "consensus_set", that shows which
  # obs. were used to fit the best model
  data[, "consensus_set"] <- get_consensus_set(data, inliers)
  # Output
  best_fit <- list(model = best_model, data = data)
  best_fit
}


###################### get consistent data ########################################
get_consistent_data <- function(hypo_inliers, hypo_outliers,
                                hypo_model, error_threshold) {
  also_inliers <- data.frame()

  predicted_values <- predict(hypo_model, newdata = hypo_outliers)
  predicted_values <- unname(predicted_values)

  # hypo_outliers is a df, having in the 1.column the Response values
  for (i in 1:nrow(hypo_outliers)) {
    # use L1 loss to measure the goodness of the fit of each obs, if the loss
    # is less than the error threshold we will mark it as also_inlier
    if (abs(hypo_outliers[i, 1] - predicted_values[i]) <= error_threshold) {
      also_inliers <- rbind.data.frame(also_inliers, hypo_outliers[i, ])
    } else {
      also_inliers <- also_inliers
    }
  }
  # also_inlier & hypo_inliers are merged together
  inliers <- rbind.data.frame(hypo_inliers, also_inliers)
}

##################### get lm ##################################################
get_lm <- function(formula, data) {
  lm(formula = formula, data = data)
}

##################### get consensus set ########################################
get_consensus_set <- function(data, inliers) {
  consensus_set <- rep(NA, times = nrow(data))

  for (i in 1:nrow(data)) {
    if (any(rownames(data)[i] == rownames(inliers))) {
      consensus_set[i] <- TRUE
    } else {
      consensus_set[i] <- FALSE
    }
  }
  consensus_set
}

#################### sampling the hypo inliers ################################
# I took the solution of the stability selection ex. and adjusted it, in order
# to use in the ransac-ex

resample <- function(data, strata = strata, fraction = fraction) {
  nrows <- nrow(data)
  rows <- sample_without_replacement(nrows, strata, fraction)
  data[rows, ]
}

sample_without_replacement <- function(nrows, strata = strata, fraction = fraction) {
  if (is.null(strata)) {
    return(sample(1:nrows, size = ceiling(nrows * fraction), replace = FALSE))
  }
  # vectorize subsampling over:
  # <rows belonging to the different strata>,
  strata_index <- split(seq_len(nrows), strata)
  # <subsample sizes in different strata>
  # round up so that subsample has at least 1 member from each strata
  strata_sizes <- ceiling(table(strata) * fraction)
  rows <- mapply(sample,
    x = strata_index, size = strata_sizes,
    replace = FALSE
  )
  unlist(rows)
}

debugonce(ransaclm)
ransaclm(y ~ . - inlier,
  data = data_simple, error_threshold = 2,
  inlier_threshold = 50, iterations = 10,
  fraction = 0.4, strata = NULL
)

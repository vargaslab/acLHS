#' Set parameters for computing a Variogram.
#'
#' Sets specific parameters for computing Variograms within the acLHS 1D or 2D
#' function calls. Note that the lag value computed for Variograms will always
#' be the 'minimum' of the independent data (i.e., for 1D minimum time between
#' points and for 2D minimum distance between points).
#'
#' @param num_lags The number of lags
#' @param dir The direction
#' @param tol The tolerance
#' @param min_pairs The minimum number of pairs
#' @return A list of the set parameters
#' @examples
#' ## Store the parameters into a variable
#' v_params <- aclhs.vario_params(num_lags=10, dir=0, tol=90, min_pairs=1)
#'
#' ## Access one of the the set parameters
#' v_params$num_lags
#' @export
aclhs.vario_params <- function(num_lags = 8, dir = 0, tol = 90, min_pairs = 1) {
  # Return Variogram parameters as a list
  return (list(num_lags = num_lags,
               dir = dir,
               tol = tol,
               min_pairs = min_pairs))
}

#' Print out different correlation values.
#'
#' Prints out pearson, spearman, and kendall correlations
#' for the given dataset and it's acLHS subsampled version.
#'
#' @param df The original data in dataframe format
#' @param aclhs_samples The acLHS-derived sample indices
#' @examples
#' ## Get the data of interest and get the acLHS sample indices
#' data(ex_data_2D)
#' input2D <- ex_data_2D
#' aclhs_sam <- aclhs(df=input2D, num_samples=50, weights=c(1,1,1), iter=100)
#'
#' ## Print out correlations
#' aclhs.get_correlations(df=input2D, aclhs_samples=aclhs_sam)
#' @export
aclhs.get_correlations <- function(df, aclhs_samples) {
  # Subsample the original dataset
  df_sub <- df[aclhs_samples,]
  ncols <- ncol(df)

  cat("Original Data Correlations\n")
  cat("Pearson:", round(stats::cor(df[ncols-1], df[ncols], method="pearson"), 3), "\n")
  cat("Spearman:", round(stats::cor(df[ncols-1], df[ncols], method="spearman"), 3), "\n")
  cat("Kendall:", round(stats::cor(df[ncols-1], df[ncols], method="kendall"), 3), "\n")
  cat("\n")

  cat("Subsampled Data Correlations\n")
  cat("Pearson:", round(stats::cor(df_sub[ncols-1], df_sub[ncols], method="pearson"), 3), "\n")
  cat("Spearman:", round(stats::cor(df_sub[ncols-1], df_sub[ncols], method="spearman"), 3), "\n")
  cat("Kendall:", round(stats::cor(df_sub[ncols-1], df_sub[ncols], method="kendall"), 3), "\n")
}

#' Computes a score from three objective functions.
#'
#' Computes a score from the sum of three objective functions multiplied by
#' their respective weights. The score is used to determine the best set of
#' indices subsampled by the acLHS algorithm, where lower is better.
#'
#' @param var_samples Subsampled indices to test
#' @param df A dataframe with three columns of data
#' @param num_samples The number of subsamples
#' @param quantile_ind The quantile of the independent variable in `df`
#' @param corrs A vector of three correlations of the two variables in `df`
#' @param min_val The minimum time or distance between two points in `df`
#' @param vario_dep The computed Variogram of the data
#' @param vario_params The parameters to set for computing a Variogram
#' @param weights A vector of three weights for each objective function
#' @return Returns the summed score of the weighted objective functions
score_samples <- function(var_samples, df, num_samples, quantile_ind, corrs,
                          min_val, vario_dep, vario_params, weights) {
  # Extract valid samples
  pos_samples <- NULL
  for (i in 1:num_samples) {
    pos_samples <- rbind(pos_samples, which.min(abs(df$Dep-var_samples[i])))
  }

  # Objective function 1
  ind_samples <- df$Ind[pos_samples]
  ind_hist_count <- graphics::hist(ind_samples, breaks = quantile_ind,
                                   plot = FALSE)$counts
  frequency <- rep(1, num_samples)

  obj_result1 <- sum(abs(ind_hist_count-frequency))

  # Objective function 2
  pos_corr_p <- round(stats::cor(df$Ind[pos_samples], df$Dep[pos_samples],
                                 method="pearson"), 3)
  pos_corr_s <- round(stats::cor(df$Ind[pos_samples], df$Dep[pos_samples],
                                 method="spearman"), 3)
  pos_corr_k <- round(stats::cor(df$Ind[pos_samples], df$Dep[pos_samples],
                                 method="kendall"), 3)

  obj_result2 <- (abs(pos_corr_p-corrs[1]) + abs(pos_corr_s-corrs[2])
                  + abs(pos_corr_k-corrs[3]))

  # Objective function 3
  pre_geodata <- NULL
  if (ncol(df) == 3) {
    pre_geodata <- cbind(df$Time[pos_samples], numeric(length(df$Time[pos_samples])), var_samples)
  }
  if (ncol(df) == 4) {
    pre_geodata <- cbind(df$X[pos_samples], df$Y[pos_samples], var_samples)
  }

  sample_geodata <- geoR::as.geodata(
    pre_geodata,
    coords.col = 1:2,
    data.col = 3)

  sample_variogram <- geoR::variog(
    sample_geodata,
    breaks = c(seq(
      min_val/2,
      min_val*(vario_params$num_lags+1),
      min_val)),
    trend = "cte",
    lambda = 1,
    estimator.type = "classical",
    nugget.tolerance = 0,
    direction = vario_params$dir,
    tolerance = vario_params$tol,
    unit.angle = "degrees",
    pairs.min = vario_params$min_pairs,
    messages = FALSE,
    quiet = TRUE)

  obj_result3 <- sum(abs(sample_variogram$v-vario_dep$v))

  return ((weights[1]*obj_result1) + (weights[2]*obj_result2)
          + (weights[3]*obj_result3))
}

#' Get subsample indices using the acLHS algorithm.
#'
#' This function extracts a desired number of subsample indices from a dataframe
#' using the acLHS algorithm. The function works for either 1D or 2D data, where
#' it is assumed the last two columns of data are the independent and dependent
#' variables, respectively.
#'
#' @param df A dataframe with three columns of data
#' @param num_samples The number of desired subsamples
#' @param weights A vector of three weights for each objective function
#' @param iter The max number of iterations to perform to find optimized indices
#' @param vario_params A list of parameters to use when computing Variograms
#' @param export_file The name of a CSV to export subsampled rows to
#' @param seed A seed to use for randomization in the optimization process
#' @return A numeric vector of subsample indices of the original data
#' @examples
#' ## acLHS sampling example
#' data(ex_data_2D)
#' input2D <- ex_data_2D
#'
#' # Set Variogram parameters
#' v_params <- aclhs.vario_params(num_lags=10, dir=0, tol=90, min_pairs=1)
#'
#' ## Set weights for each objective function, respectively
#' w <- c(10, 1000, 0.001)
#'
#' ## Run the sampling algorithm
#' aclhs_samples <- aclhs(df=input2D, num_samples=50, weights=w, iter=100,
#'                        vario_params=v_params,
#'                        export_file=tempfile(fileext=".csv"))
#'
#' ## Subsample original data
#' df_sampled <- input2D[aclhs_samples,]
#' @export
aclhs <- function(df, num_samples, weights, iter = 1000,
                  vario_params = aclhs.vario_params(),
                  export_file = NULL, seed = NULL) {
  # Validate provided number of columns in dataframe is correct
  num_cols <- ncol(df)
  if (num_cols != 3 && num_cols != 4) {
    stop("Provided dataframe should only have three or four columns.\n")
  }

  # Set seed for randomization if one specified
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Initialize some variables
  min_val <- NULL
  geodata_cols <- NULL
  aclhs_samples <- NULL
  pos_opt_samples <- NULL
  prev_colnames <- colnames(df)

  # Do some computations based on dimensions of data
  if (num_cols == 3) {
    # Set column names
    colnames(df) <- c("Time", "Ind", "Dep")

    # Set columns of data to use for the Variogram
    geodata_cols <- cbind(df$Time, numeric(length(df$Time)), df$Dep)

    # Set minimum time
    min_val <- min(stats::dist(df$Time))
  }
  else {
    # Set column names
    colnames(df) <- c("X", "Y", "Ind", "Dep")

    # Set columns of data to use for the Variogram
    geodata_cols <- cbind(df$X, df$Y, df$Dep)

    # Set minimum distance
    min_val <- min(stats::dist(df[,c("X","Y")]))
  }

  # Compute the variogram of the depdendent variable
  geodata_dep <- geoR::as.geodata(
    geodata_cols,
    header = TRUE,
    coords.col = 1:2,
    data.col = 3)

  invisible(utils::capture.output(vario_dep <- geoR::variog(
    geodata_dep,
    breaks = c(seq(
      min_val/2,
      min_val*(vario_params$num_lags+1),
      min_val)),
    trend = "cte",
    lambda = 1,
    estimator.type = "classical",
    nugget.tolerance = 0,
    direction = vario_params$dir,
    tolerance = vario_params$tol,
    unit.angle = "degrees",
    pairs.min = vario_params$min_pairs,
    messages = FALSE,
    quiet = TRUE)))

  # Get the quantile for the independent variable
  fn_limit <- seq(0, 1, length.out = num_samples+1)
  quantile_ind <- stats::quantile(df$Ind, fn_limit)

  # Compute statistical correlations of the dependent and independent variable
  corr_pearson <- round(stats::cor(df$Ind, df$Dep, method="pearson"), 3)
  corr_spearman <- round(stats::cor(df$Ind, df$Dep, method="spearman"), 3)
  corr_kendall <- round(stats::cor(df$Ind, df$Dep, method="kendall"), 3)
  corrs <- c(corr_pearson, corr_spearman, corr_kendall)

  # Stage some data for use by DEoptim optimization
  fn_opt_limit <- seq(1/num_samples, 1, length.out = num_samples)
  quantile_opt_dep <- unname(stats::quantile(df$Dep, fn_opt_limit))
  ordered_dep <- df$Dep[order(df$Dep)]

  for (i in 1:num_samples) {
    pos_opt_samples <- rbind(pos_opt_samples,
                             which.min(abs(ordered_dep-quantile_opt_dep[i])))
  }

  total_samples <- round(length(df$Dep)/num_samples)
  initial_pop <- matrix(data = NA, nrow = total_samples, ncol = num_samples)

  for (i in 1:num_samples) {
    for (j in 1:total_samples) {
      position <- pos_opt_samples[i]
      initial_pop[j,i] <- ordered_dep[position-j]
    }
  }

  fn_sample_limit <- seq(0, 1, 1/num_samples)
  quantile_limit_dep <- unname(stats::quantile(df$Dep, probs = fn_sample_limit))
  subvector_min = quantile_limit_dep[-(num_samples+1)]
  subvector_max = quantile_limit_dep[-(1)]

  # Initialize control options for DEoptim function
  deoptim_controls <- DEoptim::DEoptim.control(
    VTR = 0.000001, strategy = 3,
    itermax = iter,
    reltol = 1e-8,
    CR = 0.5,
    F = 0.8,
    trace = FALSE,
    NP = nrow(initial_pop),
    initialpop = initial_pop)

  # Find optimal samples using DEoptim
  invisible(utils::capture.output(suppressWarnings(out_DEoptim <- DEoptim::DEoptim(
    fn = score_samples,
    lower = subvector_min,
    upper = subvector_max,
    df = df,
    num_samples = num_samples,
    quantile_ind = quantile_ind,
    corrs = corrs,
    min_val = min_val,
    vario_dep = vario_dep,
    vario_params = vario_params,
    weights = weights,
    control = deoptim_controls))))

  # Extract best samples indices according to DEoptim optimization
  best_samples <- as.numeric(out_DEoptim$optim$bestmem)
  for (i in 1:num_samples) {
    aclhs_samples <- rbind(aclhs_samples, which.min(abs(df$Dep-best_samples[i])))
  }

  # Export subsampled rows of original data to a CSV if requested
  if (!is.null(export_file)) {
    out_df <- df[aclhs_samples,]
    colnames(out_df) <- prev_colnames
    utils::write.csv(out_df, export_file, row.names = FALSE)
  }

  # Remove randomized seed
  if (!is.null(seed)) {
    rm(.Random.seed, envir=globalenv())
  }

  # Return the indices
  return(aclhs_samples)
}

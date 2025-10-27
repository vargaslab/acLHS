#' Set parameters for plotting.
#'
#' Sets various parameters for plotting including plot title, axis labels,
#' plot dimensions and resolution, and whether to add a legend to the plot.
#' By default, a plot will not be created, and the location of where the
#' legend should be placed on the plot should be passed (e.g., "topright").
#'
#' @param file_name The name of the file to store the plot in (should end with '.png')
#' @param plot_title The title of the plot (default is blank)
#' @param xlab The label for the x axis of the plot (default is blank)
#' @param ylab The label for the y axis of the plot (default is blank)
#' @param width The width of the plot (default is 1000)
#' @param height The height of the plot (default is 1000)
#' @param res The resolution of the plot (default is 150)
#' @param legend The location of the legend on the plot (default is NULL)
#' @return A list of the set plotting parameters
#' @examples
#' ## Set the parameters
#' p_params <- aclhs.plot_params(file_name=tempfile(fileext=".png"),
#'                               plot_title=expression(bold("Sample Distribution")),
#'                               xlab=expression(bold("X [km]")),
#'                               ylab=expression(bold("Y [km]")),
#'                               legend="topright")
#'
#' ## Access one of the the set parameters
#' p_params$plot_title
#' @export
aclhs.plot_params <- function(file_name, plot_title = "", xlab = "",
                              ylab = "", width = 1000, height = 1000,
                              res = 150, legend = NULL) {
  # Return plotting parameters as a list
  return (list(file_name = file_name,
               plot_title = plot_title,
               xlab = xlab,
               ylab = ylab,
               width = width,
               height = height,
               res = res,
               legend = legend))
}

#' Plots the acLHS samples distribution.
#'
#' Plots the acLHS sample distribution for either 1D or 2D data.
#' acLHS samples will be overlayed over the original data points in blue.
#'
#' @param df The original data in dataframe format
#' @param aclhs_samples The acLHS-derived sample indices
#' @param plot_params The plotting parameters to use
#' @examples
#' ## Get the data of interest and get the acLHS sample indices
#' data(ex_data_2D)
#' input2D <- ex_data_2D
#' aclhs_sam <- aclhs(df=input2D, num_samples=50, weights=c(1,1,1), iter=100)
#'
#' ## Set plotting parameters
#' p_params <- aclhs.plot_params(file_name=tempfile(fileext=".png"),
#'                               xlab=expression(bold("X [km]")),
#'                               ylab=expression(bold("Y [km]")))
#'
#' ## Create plot
#' aclhs.plot_sampling_distribution(df=input2D, aclhs_samples=aclhs_sam,
#'                                  plot_params=p_params)
#' @export
aclhs.plot_sampling_distribution <- function(df, aclhs_samples, plot_params) {
  # Create the PNG
  grDevices::png(plot_params$file_name, bg = "white", width = plot_params$width,
                 height = plot_params$height, res = plot_params$res)
  graphics::par(mfrow = c(1, 1), mar = c(4, 6, 2, 2), cex.lab = 1.0, cex.axis = 1.0)

  # Do correct plotting based on dimensions of data
  df_temp <- df
  if (ncol(df) == 3) {
    # Update column names and subsample original data
    colnames(df_temp) <- c("Time", "Ind", "Dep")
    df_sub <- df_temp[aclhs_samples,]

    # Plot all original points
    plot(df_temp$Time, df_temp$Dep, type = "p", lwd = 2, pch = 1, bty="o", col = "lightgray",
         xlim = c(min(df_temp$Time), max(df_temp$Time)), ylim = c(min(df_temp$Dep), max(df_temp$Dep)),
         xlab = "", ylab = "")
    graphics::grid(col = "lightgray", lty = "dashed", lwd = graphics::par("lwd"), equilogs = TRUE)

    # Plot subsampled points
    graphics::par(new = TRUE)
    plot(df_sub$Time, df_sub$Dep, pch = 18, cex = 1.5, bty = "o", col = "blue",
         xlim = c(min(df_temp$Time), max(df_temp$Time)), ylim = c(min(df_temp$Dep), max(df_temp$Dep)),
         xlab = plot_params$xlab, ylab = plot_params$ylab)
  }
  else {
    # Update column names and subsample original data
    colnames(df_temp) <- c("X", "Y", "Ind", "Dep")
    df_sub <- df_temp[aclhs_samples,]

    # Plot all original points
    plot(df_temp$X, df_temp$Y, pch = 0, xlim = c(min(df_temp$X), max(df_temp$X)),
         ylim = c(min(df_temp$Y), max(df_temp$Y)), xlab = "", ylab = "")

    # Plot subsampled points
    graphics::par(new = TRUE)
    plot(df_sub$X, df_sub$Y, pch = 18, col = "blue", main = plot_params$plot_title,
         xlim = c(min(df_temp$X), max(df_temp$X)), ylim = c(min(df_temp$Y), max(df_temp$Y)),
         xlab = plot_params$xlab, ylab = plot_params$ylab)
  }

  # Add legend if requested
  if (!is.null(plot_params$legend)) {
    graphics::legend(plot_params$legend, legend = c("Original", "acLHS"),
                     pch = c(19, 15, 17, 18), col = c("lightgray", "blue"), box.lty = 0)
  }

  graphics::box()
  grDevices::dev.off()
}

#' Plot the univariate PDF for a column of acLHS-derived samples.
#'
#' Plots the univariate PDF of acLHS-sampled points over the original
#' univariate PDF data. The PDF can be plotted for either the dependent
#' or independent variable of the original data.
#'
#' @param df The original data in dataframe format
#' @param aclhs_samples The acLHS-derived sample indices
#' @param col The column of data to plot
#' @param plot_params The plotting parameters to use
#' @examples
#' ## Get the data of interest and get the acLHS sample indices
#' data(ex_data_2D)
#' input2D <- ex_data_2D
#' aclhs_sam <- aclhs(df=input2D, num_samples=50, weights=c(1,1,1), iter=100)
#'
#' ## Set plotting parameters
#' p_params <- aclhs.plot_params(file_name=tempfile(fileext=".png"),
#'                               xlab=expression(bold("Temperature [Celsius]")),
#'                               ylab=expression(bold("Fn(Temperature)")))
#'
#' ## Create plot
#' aclhs.plot_univariate_pdf(df=input2D, aclhs_samples=aclhs_sam, col=3,
#'                           plot_params=p_params)
#' @export
aclhs.plot_univariate_pdf <- function(df, aclhs_samples, col, plot_params) {
  # Create the PNG
  grDevices::png(plot_params$file_name, bg = "white", width = plot_params$width,
                 height = plot_params$height, res = plot_params$res)
  graphics::par(mfrow = c(1, 1), mar = c(4, 6, 2, 2), cex.lab = 1.2, cex.axis = 1.2)

  # Subsample the data
  df_sub <- df[aclhs_samples,]

  # Plot the original points
  plot(stats::ecdf(df[,col]), pch = 1, yaxt = "n", col = "lightgray",
       xlim = c(min(df[,col]), max(df[,col])), ylim = c(0,1),
       main = "", xlab = "", ylab = "")
  graphics::axis(2, at = seq(0, 1, by = 0.2))
  graphics::grid(col = "lightgray", lty = "dashed", lwd = graphics::par("lwd"), equilogs = TRUE)

  # Plot the subsamples points
  graphics::par(new = TRUE)
  plot(stats::ecdf(df_sub[,col]), pch = 18, cex = 1.5, yaxt = "n", col = "blue",
       xlim = c(min(df[,col]), max(df[,col])), ylim = c(0,1), main = plot_params$plot_title,
       xlab = plot_params$xlab, ylab = plot_params$ylab)

  # Add legend if requested
  if (!is.null(plot_params$legend)) {
    graphics::legend(plot_params$legend, legend = c("Original", "acLHS"),
                     pch = c(19, 15, 17, 18), col = c("lightgray", "blue"), box.lty = 0)
  }

  graphics::box()
  grDevices::dev.off()
}

#' Plot the scatterplot of the acLHS subsamples.
#'
#' Plots the acLHS-sampled points of independent and dependent
#' variables of the data as a scatterplot over the original points.
#'
#' @param df The original data in dataframe format
#' @param aclhs_samples The acLHS-derived sample indices
#' @param plot_params The plotting parameters to use
#' @examples
#' #' ## Get the data of interest and get the acLHS sample indices
#' data(ex_data_2D)
#' input2D <- ex_data_2D
#' aclhs_sam <- aclhs(df=input2D, num_samples=50, weights=c(1,1,1), iter=100)
#'
#' ## Set plotting parameters
#' p_params <- aclhs.plot_params(file_name=tempfile(fileext=".png"),
#'                               xlab=expression(bold("Temperature")),
#'                               ylab=expression(bold("CO2 Efflux")))
#'
#' ## Create plot
#' aclhs.plot_scatterplot(df=input2D, aclhs_samples=aclhs_sam,
#'                        plot_params=p_params)
#' @export
aclhs.plot_scatterplot <- function(df, aclhs_samples, plot_params) {
  # Create the PNG
  grDevices::png(plot_params$file_name, bg = "white", width = plot_params$width,
                 height = plot_params$height, res = plot_params$res)
  graphics::par(mfrow = c(1, 1), mar = c(4, 6, 2, 2), cex.lab = 1.2, cex.axis = 1.2)

  # Update the column names based on the dimensions of the data
  df_temp <- df
  if (ncol(df) == 3) {
    colnames(df_temp) <- c("Time", "Ind", "Dep")
  }
  else {
    colnames(df_temp) <- c("X", "Y", "Ind", "Dep")
  }

  # Subsample the data
  df_sub <- df_temp[aclhs_samples,]

  # Plot the original points
  plot(df_temp$Ind, df_temp$Dep, pch = 1, col = "lightgray",
       xlim = c(min(df_temp$Ind), max(df_temp$Ind)),
       ylim = c(min(df_temp$Dep), max(df_temp$Dep)), xlab = "", ylab = "")
  graphics::grid(col = "lightgray", lty = "dashed", lwd = graphics::par("lwd"), equilogs = TRUE)

  # Plot the subsampled points
  graphics::par(new = TRUE)
  plot(df_sub$Ind, df_sub$Dep, pch = 18, cex = 1.5, col = "blue",
       xlim = c(min(df_temp$Ind), max(df_temp$Ind)), ylim = c(min(df_temp$Dep), max(df_temp$Dep)),
       main = plot_params$plot_title, xlab = plot_params$xlab,
       ylab = plot_params$ylab)

  # Add legend if requested
  if (!is.null(plot_params$legend)) {
    graphics::legend(plot_params$legend, legend = c("Original", "acLHS"),
                     pch = c(19, 15, 17, 18), col = c("lightgray", "blue"), box.lty = 0)
  }

  graphics::box()
  grDevices::dev.off()
}

#' Plot the Variogram comparison of the acLHS subsamples.
#'
#' Plots the acLHS-sampled Variogram against the Variogram of
#' the original data. A best-fit curve of the original Variogram
#' is added for clearer comparison.
#'
#' @param df The original dataframe
#' @param aclhs_samples The acLHS-derived sample indices
#' @param vario_params The parameters to set for computing a Variogram
#' @param plot_params The plotting parameters to use
#' @examples
#' #' ## Get the data of interest and get the acLHS sample indices
#' data(ex_data_2D)
#' input2D <- ex_data_2D
#' v_params <- aclhs.vario_params(num_lags=10, dir=0, tol=90, min_pairs=1)
#' aclhs_sam <- aclhs(df=input2D, num_samples=50, weights=c(1,1,1),
#'                    iter=100, vario_params=v_params)
#'
#' ## Set plotting parameters
#' p_params <- aclhs.plot_params(file_name=tempfile(fileext=".png"),
#'                               xlab=expression(bold("Distance [km]")),
#'                               ylab=expression(bold("Semivariance")))
#'
#' ## Create plot
#' aclhs.plot_variogram_comparison(df=input2D, aclhs_samples=aclhs_sam,
#'                                 vario_params=v_params, plot_params=p_params)
#' @export
aclhs.plot_variogram_comparison <- function(df, aclhs_samples, vario_params, plot_params) {
  # Initialize some variables
  min_val <- NULL
  cols_org <- NULL
  cols_sub <- NULL

  # Stage some data for the Variograms based on data dimensions
  df_temp <- df
  if (ncol(df) == 3) {
    # Update column names and subsample original data
    colnames(df_temp) <- c("Time", "Ind", "Dep")
    df_sub <- df_temp[aclhs_samples,]

    # Stage the geodata columns and calculate minimum time
    cols_org <- cbind(df_temp$Time, numeric(length(df_temp$Time)), df_temp$Dep)
    cols_sub <- cbind(df_sub$Time, numeric(length(df_sub$Time)), df_sub$Dep)
    min_val <- min(stats::dist(df_temp$Time))
  }
  else {
    # Update column names and subsample original data
    colnames(df_temp) <- c("X", "Y", "Ind", "Dep")
    df_sub <- df_temp[aclhs_samples,]

    # Stage the geodata columns and calculate the minimum distance
    cols_org <- cbind(df_temp$X, df_temp$Y, df_temp$Dep)
    cols_sub <- cbind(df_sub$X, df_sub$Y, df_sub$Dep)
    min_val <- min(stats::dist(df_temp[,c("X","Y")]))
  }

  # Compute the original Variogram
  geodata_org <- geoR::as.geodata(
    cols_org,
    header = TRUE,
    coords.col = 1:2,
    data.col = 3)

  invisible(utils::capture.output(vario_org <- geoR::variog(
    geodata_org,
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

  # Compute the subsampled Variogram
  geodata_sub <- geoR::as.geodata(
    cols_sub,
    header = TRUE,
    coords.col = 1:2,
    data.col = 3)

  invisible(utils::capture.output(vario_sub <- geoR::variog(
    geodata_sub,
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

  # Create the PNG
  grDevices::png(plot_params$file_name, bg = "white", width = plot_params$width,
                 height = plot_params$height, res = plot_params$res)
  graphics::par(mfrow = c(1, 1), mar = c(4, 6, 2, 2), cex.lab = 1.2, cex.axis = 1.2)

  # Plot the original points
  zoom_val <- 1.2
  plot(vario_org$u, vario_org$v, pch = 1, col = "black", xlim = c(0, max(vario_org$u)),
       ylim = c(0, max(vario_org$v)*zoom_val), xlab = "", ylab = "")
  graphics::grid(col = "lightgray", lty = "dashed", lwd = graphics::par("lwd"), equilogs = TRUE)

  # Plot the sampled points
  graphics::par(new = TRUE)
  plot(vario_sub$u, vario_sub$v, pch = 18, cex = 1.5, col = "blue", xlim = c(0, max(vario_org$u)),
       ylim = c(0, max(vario_org$v)*zoom_val), xlab = plot_params$xlab, ylab = plot_params$ylab)
  graphics::par(new = TRUE)

  # Plot the best fit curve for original Variogram
  suppressWarnings(vario_fit <- geoR::variofit(vario_org, messages = FALSE))
  graphics::lines(vario_fit, lwd = 3, col="black")

  # Add legend if requested
  if (!is.null(plot_params$legend)) {
    graphics::legend(plot_params$legend, legend = c("Original", "acLHS"),
                     pch = c(19, 15, 17, 18), col = c("gray", "blue"), box.lty = 0)
  }

  graphics::box()
  grDevices::dev.off()
}

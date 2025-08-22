.density_plot <- function(data, x, fill, fill_scale = NULL, color_scale = NULL,
                         title = NULL, subtitle = NULL,
                         xlab = x, fill_lab = fill) {
  p <- ggplot(data, aes(.data[[x]], fill = .data[[fill]], 
                        color = .data[[fill]])) +
    geom_density(alpha = 0.2) +
    fill_scale +
    labs(title = title, subtitle = subtitle, x = xlab, 
         fill = fill_lab, color = NULL) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    color_scale

  p
}

#' Plot distance density
#'
#' Plot density of distances between samples in QC samples and actual samples.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param all_features logical, should all features be used? 
#' If FALSE (the default), flagged features are removed before visualization.
#' @param dist_method method for calculating the distances, passed to 
#' \code{\link[stats]{dist}}
#' @param center logical, should the data be centered?
#' @param scale scaling used, as in \code{\link[pcaMethods]{prep}} 
#' Default is "uv" for unit variance
#' @param color_scale a scale for the color of the edge of density curves, as 
#' returned by a ggplot function
#' @param fill_scale a scale for the fill of the density curves, as returned by 
#' a ggplot function
#' @param title the plot title
#' @param subtitle the plot subtitle
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @return A ggplot object.
#'
#' @examples
#' data(example_set, package = "notame")
#' plot_dist_density(example_set)
#' # Drift correction tightens QCs together
#' plot_dist_density(notame::correct_drift(example_set))
#'
#' @seealso \code{\link[stats]{dist}}
#'
#' @export
plot_dist_density <- function(object, all_features = FALSE, 
                              dist_method = "euclidean", center = TRUE, 
                              scale = "uv", 
                              color_scale = getOption("notame.color_scale_dis"),
                              fill_scale = getOption("notame.fill_scale_dis"),
                              title = paste("Density plot of", dist_method,
                                            "distances between samples"),
                              subtitle = NULL, assay.type = NULL) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
    stop("Package \'pcaMethods\' needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_QC = TRUE, assay.type = from)

  assay <- pcaMethods::prep(t(assay(object, from)), 
                            center = center, scale = scale)

  qc_data <- assay[object$QC == "QC", ]
  sample_data <- assay[!object$QC == "QC", ]

  qc_dist <- stats::dist(qc_data, method = dist_method) %>% as.numeric()
  sample_dist <- stats::dist(sample_data, method = dist_method) %>% as.numeric()
  qc <- rep(c("QC", "Sample"), times = c(length(qc_dist), length(sample_dist)))
  qc <- rep(c("QC", "Sample"), times = c(length(qc_dist), length(sample_dist)))
  distances <- data.frame(dist = c(qc_dist, sample_dist), qc = qc)

  .density_plot(distances, x = "dist", fill = "qc", fill_scale = fill_scale,
               color_scale = color_scale, xlab = "Distance", fill_lab = NULL,
               title = title, subtitle = subtitle)
}

#' Estimate the magnitude of drift
#'
#' Plots histograms of p-values from linear regression model, where each 
#' feature is predicted
#' by injection order alone. The expected uniform distribution is represented 
#' by a dashed red line.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param all_features logical, should all features be used? 
#' If FALSE (the default), flagged features are removed before visualization.
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @return A ggplot object.
#'
#' @seealso \code{\link{plot_p_histogram}}
#'
#' @examples
#' data(example_set, package = "notame")
#' plot_injection_lm(example_set)
#'
#' @export
plot_injection_lm <- function(object, all_features = FALSE, assay.type = NULL) {
  if (!requireNamespace("notameStats", quietly = TRUE)) {
    stop("Package \'notameStats\' needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_injection = TRUE, pheno_QC = TRUE,
                         assay.type = from)

  # Apply linear model to QC samples and biological samples separately
  lm_all <- notameStats::perform_lm(object, "Feature ~ Injection_order", 
                                    assay.type = from)
  lm_sample <- notameStats::perform_lm(object[, object$QC != "QC"], 
                                       "Feature ~ Injection_order",
                                       assay.type = from)
  lm_qc <- notameStats::perform_lm(object[, object$QC == "QC"],
                                   "Feature ~ Injection_order",
                                   assay.type = from)

  # Only interested in the p_values
  p_values <- list("All samples" = lm_all$Injection_order_P,
                   "Biological samples" = lm_sample$Injection_order_P,
                   "QC samples" = lm_qc$Injection_order_P)
  # Plotting
  plot_p_histogram(p_values)
}

#' Histogram of p-values
#'
#' Draws histograms of p-values with expected uniform distribution represented 
#' by a dashed red line.
#'
#' @param p_values list or data frame, each element/column is a vector of p-
#' values. The list names are used as plot titles
#' @param hline logical, whether a horizontal line representing uniform 
#' distribution should be plotted
#' @param combine logical, whether plots of individual p-value vectors should 
#' be combined into a single object.
#' Set to FALSE if you want to add other plots to the list before plotting
#' @param x_label the x-axis label
#'
#' @examples 
#' data(example_set, package = "notame")
#' lm_sample <- notameStats::perform_lm(notame::drop_qcs(example_set),
#'   "Feature ~ Injection_order")
#' p_values <- list("Biological samples" = lm_sample$Injection_order_P)
#' plot_p_histogram(p_values)
#'
#' @return If combine = TRUE, a ggplot object. Otherwise a list of ggplot 
#' objects.
#'
#' @export
plot_p_histogram <- function(p_values, hline = TRUE, combine = TRUE, 
                             x_label = "p-value") {
  # Custom breaks for the x-axis
  breaks <- seq(0, 1, by = 0.05)

  # THree separate histograms
  plots <- list()
  for (i in seq_along(p_values)) {
    p <- ggplot(data.frame(P = p_values[[i]]), aes(.data$P)) +
      geom_histogram(breaks = breaks, col = "grey50", 
                     fill = "grey80", size = 1) +
      labs(x = x_label, y = "Frequency") +
      ggtitle(names(p_values)[i]) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

    if (hline) {
      # Compute the position of the expected line
      finite_count <- sum(is.finite(p_values[[i]]))
      h_line <- finite_count / (length(breaks) - 1)
      p <- p + geom_hline(yintercept = h_line, color = "red", 
                          linetype = "dashed", size = 1)
    }

    plots <- c(plots, list(p))
  }

  if (combine) {
    return(cowplot::plot_grid(plotlist = plots, ncol = 1))
  } else {
    return(plots)
  }
}


#' Plot quality metrics
#'
#' Plots distribution of each quality metric, and a distribution of the flags.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param plot_flags logical, should the distribution of flags be added as a 
#' barplot?
#' @param assay.type character, assay to be used in case of multiple assays and 
#' no quality metrics are present in feature data
#'
#' @return A ggplot object.
#'
#' @examples
#' data(example_set, package = "notame")
#' plot_quality(example_set)
#'
#' @export
plot_quality <- function(object, all_features = FALSE, plot_flags = TRUE,
                         assay.type = NULL) {
  # Drop flagged features
  object <- drop_flagged(object, all_features = all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, feature_flag = TRUE)
  if (plot_flags) {
    # Plot bar plot of flags
    flags <- flag(object)
    flags[is.na(flags)] <- "Good"
    flags <- factor(flags) %>% stats::relevel(ref = "Good")

    fp <- ggplot(data.frame(flags), aes(x = flags)) +
      geom_bar(col = "grey50", fill = "grey80", size = 1) +
      scale_y_continuous(sec.axis = sec_axis(~ . * 100 / length(flags), 
                         name = "Percentage")) +
      theme_minimal() +
      labs(x = "Flag")
  }

  if (is.null(quality(object))) {
    message("\n", "Quality metrics not found, computing them now")
    object <- assess_quality(object, assay.type = from)
  }

  # Distribution of quality metrics
  qps <- plot_p_histogram(quality(object)[, -1], hline = FALSE, 
                          combine = FALSE, x_label = "")

  if (plot_flags) {
    p <- cowplot::plot_grid(plotlist = c(qps, list(fp)), ncol = 1)
  } else {
    p <- cowplot::plot_grid(plotlist = qps, ncol = 1)
  }

  p
}


#' Plot a boxplot for each sample
#'
#' Plots a boxplot of the distribution of the metabolite values for each 
#' sample. The boxplots can be ordered and filled by any combination of columns 
#' in the pheno data. By default, order and fill are both determined by the 
#' combination of group and time columns.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param all_features logical, should all features be used? If FALSE (the 
#' default), flagged features are removed before visualization.
#' @param order_by character vector, names of columns used to order the samples
#' @param fill_by character vector, names of columns used to fill the boxplots
#' @param title,subtitle character, title and subtitle of the plot
#' @param fill_scale a scale for the fill of the boxplots, as returned by a 
#' ggplot function
#' @param zoom_boxplot logical, whether outliers should be left outside the 
#' plot and only the boxplots shown. Defaults to TRUE.
#' @param assay.type character, assay to be used in case of multiple assays
#'
#' @return A ggplot object.
#'
#' @examples
#' data(example_set, package = "notame")
#' plot_sample_boxplots(example_set, order_by = "Group", fill_by = "Group")
#'
#' @export
plot_sample_boxplots <- function(object, all_features = FALSE, order_by, 
                                 fill_by, title = "Boxplot of samples", 
                                 subtitle = NULL,
                                 fill_scale = 
                                 getOption("notame.fill_scale_dis"), 
                                 zoom_boxplot = TRUE,  assay.type = NULL) {
  # Drop flagged compounds if not told otherwise
  object <- drop_flagged(object, all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = c(order_by, fill_by),
                         assay.type = from)

  data <- combined_data(object, from)

  if (length(order_by) == 1) {
    data$order_by <- data[, order_by]
  } else {
    data <- tidyr::unite(data, "order_by", order_by, remove = FALSE)
  }

  if (length(fill_by) == 1) {
    data$fill_by <- data[, fill_by]
  } else {
    data <- tidyr::unite(data, "fill_by", fill_by, remove = FALSE)
  }

  data <- data %>% dplyr::arrange(order_by)

  data$Sample_ID <- factor(data$Sample_ID, levels = data$Sample_ID)

  data <- tidyr::gather(data, "Variable", "Value", rownames(object))

  p <- ggplot(data, aes(x = .data$Sample_ID, y = .data$Value, fill = fill_by))

  ## Zooming outliers out of view
  if (zoom_boxplot) {
    # compute lower and upper whiskers
    ylimits <- data %>%
      dplyr::group_by(.data$Sample_ID) %>%
      dplyr::summarise(low = grDevices::boxplot.stats(.data$Value)$stats[1],
                       high = grDevices::boxplot.stats(.data$Value)$stats[5])


    ylimits <- c(0, max(ylimits$high))
    # scale y limits based on ylim1
    p <- p +
      geom_boxplot(outlier.shape = NA) +
      coord_cartesian(ylim = ylimits)
    # add text to main title
    subtitle <- paste(subtitle, 
                      "(zoomed in boxplot: outliers out of view)", sep = " ")
  } else {
    p <- p + geom_boxplot()
  }

  p <- p +
    labs(x = paste(order_by, collapse = "_"),
         y = "Abundance of metabolites",
         fill = paste(fill_by, collapse = "_"),
         title = title, subtitle = subtitle) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 90, vjust = 0.3)) +
    fill_scale

  p
}

#' Save batch correction plots
#'
#' Saves plots of each feature showing the effect of batch correction.
#' Plots show QC samples and regular samples inside each batch, plus the
#' batch mean for biological samples and QC samples as a horizontal line.
#' The dashed line represents QC mean, the filled line represents biological
#' sample mean.
#' NOTE: if you change the shape variable, be sure to set a shape scale as well,
#' the default scale only has 2 values, so it can only accomodate 2 shapes.
#'
#' @param orig,corrected \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' objects before and after batch effect correction
#' @param file path to the PDF file where the plots will be saved
#' @param width,height width and height of the plots in inches
#' @param batch,color,shape column names of pheno data for batch labels,
#' and column used for coloring and shaping points (by default batch and QC)
#' @param color_scale,shape_scale scales for color and scale as returned by 
#' ggplot functions.
#' @param assay.type1 character, assay of orig to be used in case of 
#' multiple assays.
#' @param assay.type2 character, assay of corrected to be used in case of
#' multiple assays. If corrected is not supplied, this argument selects
#' another assay from orig.
#' @return None, the function is invoked for its plot-saving side effect.
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' data(example_set, package = "notame")
#' # Batch correction
#' batch_corrected <- batchCorr::normalizeBatches(example_set, assay.type = 1,
#'   batches = "Batch", sampleGroup = "Group", refGroup = "QC", 
#'   population = "all", name = "normalized")
#' # Plots of each feature
#' save_batch_plots(
#'   orig = example_set[1:10], corrected = batch_corrected[1:10],
#'   file = "batch_plots.pdf", assay.type2 = "normalized"
#' )
#' \dontshow{setwd(.old_wd)}
#' @export
save_batch_plots <- function(orig, corrected, file, width = 14, 
                             height = 10, batch = "Batch", color = "Batch", 
                             shape = "QC",
                             color_scale = getOption("notame.color_scale_dis"),
                             shape_scale = 
                             scale_shape_manual(values = c(15, 21)),
                             assay.type1 = NULL, assay.type2 = NULL) {
                               
  if (missing(corrected) && is(orig, "SummarizedExperiment")) {
    from1 <- .get_from_name(orig, assay.type1)
    from2 <- .get_from_name(orig, assay.type2)
    orig <- .check_object(orig, pheno_factors = c(batch, shape),
                          pheno_cols = color)
    data_orig <- combined_data(orig, from1)
    data_corr <-  combined_data(orig, from2)
  } else {
    from1 <- .get_from_name(orig, assay.type1)
    from2 <- .get_from_name(corrected, assay.type2)
    orig <- .check_object(orig, pheno_factors = c(batch, shape),
                          pheno_cols = color)
    corrected <- .check_object(corrected, pheno_factors = c(batch, shape),
                               pheno_cols = color)
    data_orig <- combined_data(orig, assay.type = from1)
    data_corr <- combined_data(corrected, from2)
  }
  
  # Prepare data.frame for batch means with batch and injection order range
  batch_injections <- data_orig %>%
    dplyr::group_by(!!dplyr::sym(batch)) %>%
    dplyr::summarise(start = min(.data$Injection_order), 
                     end = max(.data$Injection_order))

  batch_mean_helper <- function(data) {
    data %>%
      dplyr::group_by(!!dplyr::sym(batch)) %>%
      dplyr::summarise_at(rownames(orig), finite_mean) %>%
      dplyr::left_join(batch_injections, ., by = batch)
  }
  
  get_batch_means <- function(data) {
    batch_means <- batch_mean_helper(data) %>%
      dplyr::mutate(QC = "Sample")
    batch_means_qc <- data %>%
      dplyr::filter(.data$QC == "QC") %>%
      batch_mean_helper() %>%
      dplyr::mutate(QC = "QC")
    rbind(batch_means, batch_means_qc)
  }
  # Get batch means for QC and biological samples of the original data
  batch_means_orig <- get_batch_means(data_orig)
  # Get batch means for QC and biological samples of the corrected data
  batch_means_corr <- get_batch_means(data_corr)
  
  batch_plot_helper <- function(data, fname, batch_means) {
    p <- ggplot() +
      geom_point(data = data, 
                 mapping = aes(x = .data[["Injection_order"]], 
                               y = .data[[fname]],
                               color = .data[[color]], 
                               shape = .data[[shape]])) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      color_scale +
      shape_scale
    p <- p +
      geom_segment(data = batch_means, 
                   mapping = aes(x = .data[["start"]], xend = .data[["end"]],
                                 y = .data[[fname]], yend = .data[[fname]],
                                 color = .data[[color]], 
                                 linetype = .data[["QC"]]),
                             size = 1) +
      scale_linetype(guide = "none")
      
    p
  }
  # Save plots with original and corrected data to pdf
  grDevices::pdf(file, width = width, height = height)
  for (feature in rownames(orig)) {
    p1 <- batch_plot_helper(data_orig, feature, batch_means_orig)
    p2 <- batch_plot_helper(data_corr, feature, batch_means_corr)
    p <- cowplot::plot_grid(p1, p2, nrow = 2)
    plot(p)
  }
  grDevices::dev.off()
}

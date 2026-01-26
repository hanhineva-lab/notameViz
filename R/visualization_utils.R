#' Save plot to chosen format
#'
#' Saves the given plot to a file. Supports pdf, svg, emf, png and tiff formats.
#' If an error occurs with the plot, an empty file is created.
#'
#' @param p a ggplot object
#' @param file the file path
#' @param ... other arguments to plot function, like width and height
#'
#' @return None, the function is invoked for its plot-saving side effect.
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' data(toy_notame_set, package = "notame")
#'
#' p <- plot_sample_heatmap(toy_notame_set, group = "Group")
#'
#' save_plot(p, file = "test.pdf")
#' \dontshow{setwd(.old_wd)}
#'
#' @seealso \code{\link[grDevices]{pdf}},
#' \code{\link[devEMF]{emf}},
#' \code{\link[grDevices]{svg}},
#' \code{\link[grDevices]{png}},
#' \code{\link[grDevices]{tiff}}
#'
#' @export
save_plot <- function(p, file, ...) {
  # Create folder automatically
  folder <- dirname(file)
  if (!file.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }

  format <- utils::tail(unlist(strsplit(basename(file), split = "\\.")), 1)
  switch(format,
    "emf" = devEMF::emf(file, ...),
    "pdf" = grDevices::pdf(file, ...),
    "svg" = grDevices::svg(file, ...),
    "png" = grDevices::png(file, ...),
    "tiff" = grDevices::tiff(file, ...),
    stop("File format '", format, "' is not valid, saving failed"))
  tryCatch(
    {
      plot(p)
      grDevices::dev.off()
      log_text(paste("Saved to:", file))
    },
    error = function(e) {
      grDevices::dev.off()
      stop(e$message, call. = FALSE)
    }
  )
}

# Helper function for handling errors and keeping track of file names
.save_name <- function(object, prefix, format, fun, name, file_names,
                       width = 7, height = 7, ...) {
  save_seed <- .Random.seed
  p <- NULL
  tryCatch(
    {
      p <- fun(object, ...)
    },
    error = function(e) {
      message("Problem with plot named ", name, ":\n", e$message)
    }
  )
  
  if (!is.null(p)) {
    file_name <- paste0(prefix, "_", name, ".", format) 
    file_names <- c(file_names, file_name)
    save_plot(p, file = file_name, width = width, height = height)
  }
  assign(".Random.seed", save_seed, envir = .GlobalEnv)

  file_names
}

.merge_to_pdf <- function(prefix, file_names, remove_singles) {
  log_text("Merging plots to a single pdf")
  output_file <- paste0(prefix, ".pdf")
  qpdf::pdf_combine(file_names, output = output_file)
  log_text(paste("Saved to:", output_file))
  # Remove single files if wanted
  if (remove_singles) {
    file.remove(unlist(file_names))
    log_text("Removed single plot files")
  }
}

#' Write all relevant pretreatment visualizations to pdf
#'
#' A wrapper around all the major visualization functions, used for visualizing 
#' data between major steps of data preprocessing. Saves all visualizations as 
#' PDFs with a set prefix on filenames.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param prefix character, a file path prefix added to the file paths
#' @param format character, format in which the plots should be saved, DOES NOT 
#' support raster formats
#' @param perplexity perplexity for t-SNE plots
#' @param merge logical, whether the files should be merged to a single PDF, 
#' see Details
#' @param remove_singles logical, whether to remove single plot files after 
#' merging. Only used if \code{merge = TRUE}
#' @param group character, name of pheno data column containing the group labels
#' @param time character, name of pheno data column containing timepoints
#' @param id character, name of pheno data column containing subject identifiers
#' @param color character, name of pheno data column used for coloring sample
#' labels for dendrograms
#' @param assay.type character, assay to be used in case of multiple assays
#' @return None, the function is invoked for its plot-saving side effect.
#'
#' @details If \code{merge} is \code{TRUE} and \code{format} is \code{pdf},
#' then a file containing all the visualizations named \code{prefix.pdf} will 
#' be created. 
#' 
#' The type of visualizations to be saved depends on the type of object.
#' Here is a comprehensive list of the visualizations:
#' \itemize{
#' \item Distribution of quality metrics and flags \code{\link{plot_quality}}
#' \item Boxplots of each sample in injection order 
#' \code{\link{plot_sample_boxplots}}
#' \item PCA scores plot of samples colored by injection order 
#' \code{\link{plot_pca}}
#' \item t-SNE plot of samples colored by injection order 
#' \code{\link{plot_tsne}}
#' \item If the object has over 60 samples, hexbin versions of the PCA and t-
#' SNE plots above
#' \code{\link{plot_pca_hexbin}}, \code{\link{plot_tsne_hexbin}}
#' \item Dendrogram of samples ordered by hierarchical clustering, sample 
#' labels colored by group if present
#' \code{\link{plot_dendrogram}}
#' \item heat map of intersample distances, ordered by hierarchical clustering 
#' \code{\link{plot_sample_heatmap}}
#' \item If the object has QC samples: \itemize{
#' \item Density function of the intersample distances in both QCs and 
#' biological samples \code{\link{plot_dist_density}}
#' \item Histograms of p-values from linear regression of features against 
#' injection order in both QCs and biological samples 
#' \code{\link{plot_p_histogram}}}
#' \item If the object has a group column: \itemize{
#' \item PCA and tSNE plots with points shaped and colored by group 
#' \code{\link{plot_pca}}, \code{\link{plot_tsne}}
#' }
#' \item If the object has a time column: \itemize{
#' \item PCA and tSNE plots with points shaped and colored by time 
#' '\code{\link{plot_pca}}, \code{\link{plot_tsne}}
#' \item Dendrogram of samples ordered by hierarchical clustering, sample 
#' labels colored by time point \code{\link{plot_dendrogram}}
#' }
#' \item If the object has a group column OR a time column: \itemize{
#' \item Boxplots of samples ordered and colored by group and/or time 
#' \code{\link{plot_sample_boxplots}}
#' }
#' \item If the object has a group column AND a time column: \itemize{
#' \item PCA and tSNE plots with points shaped by group and colored by time
#' \code{\link{plot_pca}}, \code{\link{plot_tsne}}
#' }
#' \item If the object has a time column AND a subject column: \itemize{
#' \item PCA and tSNE plots with arrows connecting the samples of each subject 
#' in time point order
#' \code{\link{plot_pca_arrows}}, \code{\link{plot_tsne_arrows}}
#' }
#' }
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' data(toy_notame_set, package = "notame")
#' rp_neg_set <- toy_notame_set[rowData(toy_notame_set)$Split == "RP_neg", ]
#' save_QC_plots(rp_neg_set, prefix="figures/RP_neg", perplexity=5,
#'               group = "Group", color = "Group", time = "Time", 
#'               id = "Subject_ID")
#' \dontshow{setwd(.old_wd)}
#'
#' @export
save_QC_plots <- function(object, prefix, format = "pdf", perplexity = 30,
                          merge = FALSE, remove_singles = FALSE, group = NULL,
                          time = NULL, id = NULL, color = NULL,
                          assay.type = NULL) {
  file_names <- list()
  
  # If not grouped, plot PCA and t-SNE on QC information
  group <- ifelse(is.null(group), "QC", group)
  color <- ifelse(is.null(color), group, color)

  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_QC = TRUE,
    pheno_cols = c(time, id, color), assay.type = from)
  assays(object) <- assays(object)[from]
  if (sum(object$QC == "QC")) {
    file_names <- .save_name(object, prefix, format, fun = plot_dist_density, 
        name = "density_plot", file_names, width = 8, height = 6)
    file_names <- .save_name(object, prefix, format, plot_injection_lm,
               "lm_p_histograms", file_names)
  }
  # Quality metrics
  file_names <- .save_name(object, prefix, format, plot_quality,
      "quality_metrics", file_names)
  # Plots with injection order
  file_names <- .save_name(object, prefix, format, plot_sample_boxplots,
      "boxplots_injection", file_names, order_by = "Injection_order", 
      fill_by = "QC", width = 15)
  file_names <- .save_name(object, prefix, format, plot_pca, "PCA_injection", 
      file_names, color = "Injection_order")
  file_names <- .save_name(object, prefix, format, plot_tsne, "tSNE_injection",
      file_names, perplexity = perplexity, color = "Injection_order")
  # Clustering
  file_names <- .save_name(object, prefix, format, plot_dendrogram, 
      "dendrogram", file_names, width = 15, color = color)
  file_names <- .save_name(object, prefix, format, plot_sample_heatmap, 
      "heatmap_samples",  file_names, width = 15, height = 16, 
             group = group)
  # For large sets, plot hexbin plots
  if (ncol(object) > 60) {
    file_names <- .save_name(object, prefix, format, plot_pca_hexbin,
      "PCA_hexbin", file_names)
    file_names <- .save_name(object, prefix, format, plot_tsne_hexbin, 
      "tSNE_hexbin", file_names, perplexity = perplexity)
  }
  file_names <- .save_name(object, prefix, format, plot_pca, "PCA_group",
      file_names, color = group)
  file_names <- .save_name(object, prefix, format, plot_tsne, "tSNE_group", 
      file_names, perplexity = perplexity, color = group)
  # Time point
  if (!is.null(time)) {
    file_names <- .save_name(object, prefix, format, plot_pca, "PCA_time",
        file_names, color = time)
    file_names <- .save_name(object, prefix, format, plot_tsne, "tSNE_time", 
        file_names, color = time, perplexity = perplexity)
    file_names <- .save_name(object, prefix, format, plot_dendrogram, 
        "dendrogram_time", file_names, color = time, width = 15)
  }
  # Time point OR group
  if (!is.null(group) || !is.null(time)){
    by <- c(group, time)
    file_names <- .save_name(object, prefix, format, plot_sample_boxplots, 
        "boxplots_group", file_names, width = 15, order_by = by, fill_by = by)
  }
  # Time point AND group
  if (!is.null(group) && !is.null(time)) {
    file_names <- .save_name(object, prefix, format, plot_pca, 
        "PCA_group_time", file_names, color = time, shape = group)
    file_names <- .save_name(object, prefix, format, plot_tsne, 
        "tSNE_group_time", file_names, color = time, shape = group,
        perplexity = perplexity)
  }
  # Multiple time points per subject
  if (!is.null(time) && !is.null(id) && sum(object$QC == "QC") == 0) {
    file_names <- .save_name(object, prefix, format, plot_pca_arrows, 
        "PCA_arrows", file_names, color = group, time = time, subject = id)
    file_names <- .save_name(object, prefix, format, plot_tsne_arrows, 
        "tSNE_arrows", file_names, perplexity = perplexity, color = group, 
        time = time, subject = id)
  }
  if (merge && format == "pdf") {
    .merge_to_pdf(prefix, file_names, remove_singles)
  }
}

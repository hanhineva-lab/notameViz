#' @rawNamespace import(ggplot2, except = Position)
#' @importFrom utils citation
#' @import BiocGenerics
#' @import methods
#' @import SummarizedExperiment
#' @importFrom notame drop_flagged drop_qcs flag quality assess_quality
#' combined_data flag merge_notame_sets mark_nas "flag<-" flag_quality log_text 
#' init_log finish_log join_rowData citations
NULL

utils::globalVariables(c('.'))

# Get internal notame functions
.add_citation <- notame:::.add_citation
.get_from_name <- notame:::.get_from_name
.check_object <- notame:::.check_object
.check_feature_data <- notame:::.check_feature_data
.find_mz_rt_cols <- notame:::.find_mz_rt_cols
finite_mean <- notame:::finite_mean
finite_sd <- notame:::finite_sd

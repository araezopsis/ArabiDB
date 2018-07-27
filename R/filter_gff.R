if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' Filtering row by representative model from given data.frame
#' @importFrom magrittr %>%
#' @importFrom dplyr semi_join
#' @importFrom dplyr data_frame
#' @param df data frame
#' @param rep_model character vector
#' @param rep_model_col column name
#' @export
filter_repr_model <-
  function(df, rep_model, rep_model_col = NA_character_){
    if(!is.na(rep_model_col)){
      model_df <- data_frame(rep_model)
      colnames(model_df) <- rep_model_col
      df %>% semi_join(model_df, by = rep_model_col)
    }else{
      stop()
    }
  }

#' Filtering specified feature from gff
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @param gff gff data frame
#' @param feature feature
#' @export
filter_feature <-
  function(gff, feature = "exon"){
    feature <- str_c("^", feature, "$")
    if(length(feature) >= 2) feature <- str_c(feature, collapse = "|")
    gff %>%
      filter(str_detect(.$feature, feature))
  }


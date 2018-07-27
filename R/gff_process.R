
gff3_colnames <-
  c("seqname", "source", "feature",
    "start", "end", "score", "strand",
    "frame", "attribute")

#' extract gff3 attribute
#' @importFrom stringr str_extract_all
#' @importFrom stringr str_remove_all
#' @importFrom magrittr %>%
#' @param gff_attributes attributes. character vector
extract_gff3_attr <-
  function(gff_attributes){
    str_extract_all(gff_attributes, "(^|;).+?=") %>%
      unlist %>%
      str_remove_all("=|;") %>%
      unique
  }

#' Splitting gff3 attribute column
#' @importFrom magrittr %>%
#' @importFrom stringr str_extract
#' @importFrom stringr str_replace_all
#' @importFrom dplyr as_data_frame
#' @param gff_attributes attributes. character vector
split_attr_col <-
  function(gff_attributes){
    unique_attr <- gff_attributes %>% extract_gff3_attr()

    temp <- list()
    for(i in unique_attr){
      cat(paste0("Processing '", i, "' column\n"))
      temp[[i]] <-
        gff_attributes %>%
        {str_extract(., str_c(i, "=.*?(;|$)"))} %>%
        {str_remove_all(., str_c("(", i, "=)|(;|$)"))}
    }
    as_data_frame(temp)
  }

#' Split attribute column
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
#' @importFrom dplyr select
#' @importFrom dplyr bind_cols
#' @param gff gff3 data.frame
#' @export
split_attr_gff <-
  function(gff){
    colnames(gff) <- gff3_colnames

    gff %>%
      select(-9) %>%
      {bind_cols(., split_attr_col(gff$attribute))}
  }

#' Multiply overlaped feature rows
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr do
#' @importFrom dplyr bind_rows
#' @importFrom dplyr arrange
#' @importFrom stringr str_count
#' @importFrom stringr str_extract_all
#' @param gff gff3 data.frame
#' @export
multiply_overlap_feature <-
  function(gff){
    overlap_feature <-
      gff %>%
      mutate(overlap = str_count(.$Parent, "AT.G\\d{5}[.]\\d+")) %>%
      filter(.$overlap >= 2)

    multiply_row <-
      function(df_row){
        parent <- str_extract_all(df_row$Parent, "AT.G\\d{5}[.]\\d+")[[1]]
        df <-
          df_row[rep(1, length(parent)),] %>%
          mutate(Parent = parent)
        df
      }

    category <- (1:nrow(overlap_feature)) %>% cut(breaks = 10)
    result_list <- list()
    for(i in 1:10){
      print(paste("Running ", i, "/10 loop", sep = ""))
      temp <- overlap_feature[category == levels(category)[i],]

      result_list[[as.character(i)]] <-
        temp %>%
        group_by(.$ID) %>%
        do(multiply_row(.))
    }

    gff %>%
      mutate(overlap = str_count(.$Parent, "AT.G\\d{5}[.]\\d+")) %>%
      filter(is.na(.$overlap) | (.$overlap == 1L)) %>%
      {bind_rows(., bind_rows(result_list))} %>%
      arrange(.$seqname, .$start)
  }


#' Summarize GFF3 feature column
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr everything
#' @importFrom purrr map_int
#' @importFrom purrr map_df
#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @param gff gff3 data.frame
feature_table <-
  function(gff){
    feature_name <- gff$feature %>% unique

    mf <-
      function(x){
        map_int(feature_name, ~ sum(str_detect(x, str_c("^", ., "$"))))
      }

    feature_summary <-
      gff %>%
      split(., gff$seqname) %>%
      map_df(~ mf(.x$feature)) %>%
      mutate(feature_name) %>%
      select(feature_name, everything())

    feature_summary %>%
      select(-feature_name) %>%
      {mutate(feature_summary, All = as.integer(rowSums(.)))}
  }

#' Summarize GFF3
#' @param gff gff3 data.frame
#' @export
summarise_gff <-
  function(gff){
    list(
      Source = table(gff$source),
      Score = table(gff$score),
      Strand = table(gff$strand),
      Frame = table(gff$frame),
      Feature = feature_table(gff)
    )
  }

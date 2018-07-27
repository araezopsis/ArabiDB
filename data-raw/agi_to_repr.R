
library(magrittr)
library(tidyverse)

read_lines(mydata_urls$tair_repr_model_path, n_max = 3)
tair_repr_df <-
  read_tsv(
    mydata_urls$tair_repr_model_path,
    skip = 3,
    col_names = "transcript_id"
  ) %>%
  mutate(AGI = str_extract(transcript_id, "AT[12345CM]G\\d{5}"))
tair_repr_df %>% skimr::skim()

left_join(
  tair_repr_df,
  rename(tair10_func_df, transcript_id = name),
  by = "transcript_id"
) %>%
  filter(gene_model_type == "protein_coding")

tair_gff <- read_tsv(mydata_urls$tair_gff3_path, col_names = F)

gff3_colnames <-
  c("seqname", "source", "feature",
    "start", "end", "score", "strand",
    "frame", "attribute")

summarise_gff <-
  function(gff){
    if(ncol(gff) != 9) stop("number of columns is wrong. must be 9.")

    colnames(gff) <- gff3_colnames
    summarise_result <- list()

    summarise_result[["feature_summary"]] <-
      gff %>% group_by(feature) %>% summarise(number = n())

    summarise_result[["attr_splitted_gff"]] <-
      split_attr_gff(gff)

    summarise_result[["sum_feat_len"]] <-
      summarise_result[["attr_splitted_gff"]] %>%
      mutate(feat_len = end - start + 1L) %>%
      split(.$feature) %>%
      keep(~ all(!is.na(.$Parent))) %>%
      map(~ group_by(., Parent)) %>%
      map(~ summarise(., sum_feat_len = sum(feat_len)))

    return(summarise_result)
  }
hoge <- summarise_gff(tair_gff)

#' Mapping AGI code to its representative model
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom stringr str_sub
#' @param filepath path to text file of representative gene models
#' @export
map_agi2repr <-
  function(filepath){
    mapping <-
      read_tsv(filepath, skip = 5, col_names = F) %>%
      rename(rep_model = "X1") %>%
      mutate(AGI = str_sub(.$rep_model, 1, 9))

    if(max(table(mapping$AGI)) == 1){
      return(mapping)
    }else{
      stop()
    }
  }
agi2repr_df <- map_agi2repr(mydata_urls$tair_repr_model_path)
# usethis::use_data(agi2repr_df)

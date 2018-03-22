
#' Mapping AGI code to its representative model
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom stringr str_sub
#' @param filepath path to text file of representative gene models
#' @export
#' 
map_agi2repr <-
  function(filepath = "data-raw/list/TAIR10_representative_gene_models.txt.gz"){
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
agi2repr_df <- map_agi2repr()
usethis::use_data(agi2repr_df)
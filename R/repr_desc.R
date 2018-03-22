if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' Return gene description of representative model
#' @importFrom dplyr left_join
#' @importFrom dplyr rename
#' @export
repr_desc <-
  function(){
    left_join(
      ArabiDB::agi2repr_df,
      rename(ArabiDB::desc_df, rep_model = "name"),
      by = "rep_model"
    )
  }


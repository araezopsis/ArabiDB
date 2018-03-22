
inf <- "data-raw/function/TAIR10_functional_descriptions_20130831.txt.gz"
modified_inf <- "data-raw/function/modified_TAIR10_functional_descriptions_20130831.txt.gz"

# File modification -------------------------------------------------------

### 2 parsing failures
read_tsv(file = inf, quote = "\\")

### searching failures
temp <- readLines("data-raw/function/TAIR10_functional_descriptions_20130831.txt.gz")
str_count(temp, "\t") %>% table
str_count(temp, "\t") %>% {temp[. == 5]} %>% str_view_all("\t")

str_count(temp, "\t") %>% {temp[. == 5]} %>% str_view_all("\tto")
linenum <- str_count(temp, "\t") %>% {which(. == 5)}
temp[linenum] <- temp[linenum] %>% str_replace_all("\tto", " to")
str_count(temp, "\t") %>% table
write_lines(temp, modified_inf)


# loading description -----------------------------------------------------

#' Functional description
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom stringr str_detect
#' @param filepath path to text file of functional description
#' @export
#' 
load_func_desc <-
  function(filepath = "data-raw/function/modified_TAIR10_functional_descriptions_20130831.txt.gz"){
    read_tsv(file = filepath, quote = "\\") %>% 
      filter(str_detect(.$name, "AT.G\\d{5}[.]\\d+")) %>% 
      arrange(name)
  }
desc_df <- load_func_desc()
usethis::use_data(desc_df)


library(magrittr)
library(fs)
library(tidyverse)


# Data URL ----------------------------------------------------------------
url_mydata <- "https://github.com/araezopsis/mydata/raw/master/"
append_url <- function(x) paste0(url_mydata, x)

mydata_urls <-
  list(
    tair_gff3_path =
      append_url("tair/gff/TAIR10_GFF3_genes.gff.gz"),
    tair_gff3_aip_path =
      append_url("aip/TAIR10_GFF3_genes_transposons.AIP.gff.gz"),
    araport_gtf_path =
      append_url("aip/Araport11_GFF3_genes_transposons.201606.gtf.gz"),
    araport_gff3_path =
      append_url("aip/Araport11_GFF3_genes_transposons.201606.gff.gz"),
    tair_func_desc_path =
      append_url("aip/function/TAIR10_functional_descriptions_20130831.txt.gz"),
    tair_repr_model_path =
      append_url("aip/list/TAIR10_representative_gene_models.txt.gz")
  )

# Modification functional descriptions -------------------------------------------------------
### 2 parsing failures
read_tsv(file = mydata_urls$tair_func_desc_path, quote = "\\")

### searching failures
temp <- read_lines(mydata_urls$tair_func_desc_path)

# 普通は\tが4つなのに、2行だけ\tが5つある
str_count(temp, "\t") %>% table

# <Space>toとするべきところを\ttoとしているために、一つ\tが増えている
str_count(temp, "\t") %>% {temp[. == 5]} %>% str_view_all("\t")
str_count(temp, "\t") %>% {temp[. == 5]} %>% str_view_all("\tto")

# 問題の場所を置換
linenum <- str_count(temp, "\t") %>% {which(. == 5)}
temp[linenum] <- temp[linenum] %>% str_replace_all("\tto", " to")
str_count(temp, "\t") %>% table

modified_inf <- "./temp.txt.gz"
write_lines(temp, modified_inf)
tair10_func_df <- read_tsv(modified_inf, quote = "\\")
usethis::use_data(tair10_func_df)


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
# usethis::use_data(desc_df)

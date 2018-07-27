
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
usethis::use_data(mydata_urls)

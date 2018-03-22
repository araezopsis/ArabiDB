library(magrittr)
library(tidyverse)

tair_gff3_path <- "https://github.com/araezopsis/mydata/raw/master/tair/gff/TAIR10_GFF3_genes.gff.gz"
tair_gff3_aip_path <- "https://github.com/araezopsis/mydata/raw/master/aip/TAIR10_GFF3_genes_transposons.AIP.gff.gz"
araport_gtf_path <- "https://github.com/araezopsis/mydata/raw/master/aip/Araport11_GFF3_genes_transposons.201606.gtf.gz"
araport_gff3_path <- "https://github.com/araezopsis/mydata/raw/master/aip/Araport11_GFF3_genes_transposons.201606.gff.gz"


# TAIR10 GFF3 AIP -------------------------------------------------------------
tair_gff3_aip <- 
  read_tsv(tair_gff3_aip_path, col_names = F) %>% 
  split_attr_gff

# modified TAIR10 GFF3 AIP ---------------------------------------------------
tair_gff3_aip_modified <- multiply_overlap_feature(tair_gff3_aip)
# usethis::use_data(tair_gff3_aip, overwrite = TRUE)

# usethis::use_data(tair_gff3_aip_modified, overwrite = TRUE)

# TAIR10 GFF3 -------------------------------------------------------------
tair_gff3 <-
  read_tsv(tair_gff3_path, col_names = F) %>% 
  split_attr_gff
# usethis::use_data(tair_gff3, overwrite = TRUE)

# Araport GTF -------------------------------------------------------------
araport_gtf <- 
  read_tsv(araport_gtf_path, col_names = F)
# glimpse(araport_gtf)

# Araport GFF3 ------------------------------------------------------------
araport_gff3 <- 
  read_tsv(araport_gff3_path,
           col_names = F,
           skip = 1,
           quote = "\"",
           comment = "###"
           ) %>% 
  split_attr_gff
# usethis::use_data(araport_gff3, overwrite = TRUE)

# modified Araport GFF3 ---------------------------------------------------
araport_gff3_modified <- multiply_overlap_feature(araport_gff3)
# usethis::use_data(araport_gff3_modified, overwrite = TRUE)

# Check gff file context --------------------------------------------------
araport_gtf$X1 %>% table
araport_gtf$X2 %>% table
araport_gtf %>% filter(X2 == ".")
araport_gtf$X3 %>% table

feature_name <- araport_gtf$X3 %>% unique
mf <- 
  function(x){
    map_int(feature_name, ~ sum(str_detect(x, str_c("^", ., "$"))))
  }
feature_table <-
  araport_gtf %>%
  split(., .$X1) %>%
  map_df(~ mf(.$X3)) %>%
  mutate(feature_name) %>%
  select(feature_name, everything())
feature_table
rm(feature_table)

araport_gtf$X6 %>% table
araport_gtf$X7 %>% table
araport_gtf$X8 %>% table

araport_gtf <- araport_gtf_sepX9
usethis::use_data(araport_gtf, overwrite = TRUE)

# Filtering by representative gene model ----------------------------------
gff_df_repr <- 
  gff_df_sepV9 %>% 
  semi_join(data_frame(Parent = agi2repr_df$rep_model), by = "Parent")
gff_df_repr
gff_df_repr$X3 %>% table

# Filtering exon feature --------------------------------------------------
gff_df_exon <- 
  gff_df_repr %>% 
  filter(str_detect(X3, "^exon$")) %>% 
  select(1:8, Parent)
# Check
gff_df_exon$Parent %>% 
  unique %>% 
  {identical(intersect(., TAIR10REPR), .)}


# Deviding by coding strand -----------------------------------------------
gff_df_ps <- gff_df_exon %>% filter(X7 == "+")
gff_df_ms <- gff_df_exon %>% filter(X7 == "-")


# Calculating position ----------------------------------------------------
gff_df_ps2 <-
  gff_df_ps %>%
  arrange(X4) %>% # plus strand
  group_by(Parent) %>% 
  mutate(ep = X5-X4) %>% 
  mutate(cum_ep = cumsum(ep)) %>% 
  mutate(cum_sp = cum_ep - ep + 1L) %>% 
  select(1:10, ep, cum_ep, cum_sp) %>% 
  ungroup()

gff_df_ms2 <-
  gff_df_ms %>%
  arrange(desc(X5)) %>% # minus strand
  group_by(Parent) %>% 
  mutate(ep = X5-X4) %>% 
  mutate(cum_ep = cumsum(ep)) %>% 
  mutate(cum_sp = cum_ep - ep + 1L) %>% 
  select(1:10, ep, cum_ep, cum_sp) %>% 
  ungroup()

gff_df_2 <- bind_rows(gff_df_ps2, gff_df_ms2)
gff_df_2

walk(
  chr_name,
  ~ write_csv(filter(gff_df_2, X1 == .),
              str_c(dirname(gff_path), "/modified_", ., "_gff.csv"))
  )
rm(gff_df, gff_df_2, gff_df_exon, gff_df_ms, gff_df_ps, gff_df_ms2, gff_df_ps2)
list.files(dirname(gff_path), pattern = "gff.csv")

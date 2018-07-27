library(magrittr)
library(tidyverse)

# TAIR10 GFF3 -------------------------------------------------------------
tair_gff3 <-
  read_tsv(mydata_urls$tair_gff3_path, col_names = F) %>%
  split_attr_gff
# usethis::use_data(tair_gff3, overwrite = TRUE)

# TAIR10 GFF3 AIP -------------------------------------------------------------
tair_gff3_aip <-
  read_tsv(mydata_urls$tair_gff3_aip_path, col_names = F) %>%
  split_attr_gff
# usethis::use_data(tair_gff3_aip, overwrite = TRUE)

# modified TAIR10 GFF3 AIP ---------------------------------------------------
tair_gff3_aip_modified <- multiply_overlap_feature(tair_gff3_aip)
# usethis::use_data(tair_gff3_aip_modified, overwrite = TRUE)

# Araport GTF -------------------------------------------------------------
araport_gtf <-
  read_tsv(mydata_urls$araport_gtf_path, col_names = F)

# Araport GFF3 ------------------------------------------------------------
araport_gff3 <-
  read_tsv(mydata_urls$araport_gff3_path,
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

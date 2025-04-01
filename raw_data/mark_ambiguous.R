library(tidyverse)
library(here)
path_bed = here("raw_data/igtcr_gene.38.merged.bed4.bed")
df = read_tsv(path_bed, col_names=FALSE)
df

get_genes <- function(x){
  map_chr(strsplit(x, split = ","), \(s) {
    s2 = strsplit(s, split ="\\|");
    paste0(unique(lapply(s2, \(s3){s3[[2]]}), collapse = ","))} 
  )
}

df <- df |>
  #mutate(multi = str_detect(X4, ",")) |>
  mutate(genes = get_genes(X4)) 

df <- df |>
  mutate(genes = if_else(str_detect(genes, ","), "ambiguous", genes))

df |>
  count(genes)

df_bedstyle <- df |> 
  select(1:3, genes)

df_bedstyle |>
  write_tsv(here("src/data/igtcr_gene.38.bed4.bed"), col_names = FALSE)
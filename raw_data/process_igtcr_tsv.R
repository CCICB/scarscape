library(tidyverse)
df <- read_tsv("inst/igtcr_gene.38.tsv");
df_bed6 <- df |> 
  mutate(type = str_replace(gene, "(..).*", "\\1"), subtype = str_replace(gene, "(...).*", "\\1")) |> 
  mutate(name = paste(type, subtype, gene, region, sep="|"), score=".", posStart = posStart-1) |> 
  select(chromosome, posStart, posEnd, name, score, strand)  |> 
  filter(!is.na(posStart))
  
df_bed6 |> 
  write_tsv("inst/igtcr_gene.38.bed6.bed", col_names=F)

df_bed6 |> 
  select(-score, -strand) |> 
  write_tsv("inst/igtcr_gene.38.bed4.bed", col_names=F)


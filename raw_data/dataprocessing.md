
Downloaded data from CIDER resources.

Created initial bed4 and bed6 files using process_igtcr_tsv.R

Sorted and merged overlapping bed regions, collapsing Ids into comma separated lines as follows into 

```
bedtools sort -i igtcr_gene.38.bed4.bed -g hg38.genome | bedtools merge -c 4 -o distinct 
```

Then marked ambiguous samples and simplified names using mark_ambiguous.R
# ScarScape

> \[!WARNING\]  
> This package is in early development and not ready for use

Extract the landscape of tumorigenic scars from cancer omics data.

## Metrics

The metrics we extract from mutation data are designed to reflect different tumorigenic mechanisms.

### Hyperactivity of Endogenous rearrangements

B & T cell derived tumours can be driven by hyperactive activity of endogenous mutagenic enzymes involved in generating receptor diversity (e.g. RAG recombinase for VDJ rearrangement and APOBEC for somatic hypermutation). The following metrics are designed to separate B and T-cell tumours from each other & solid tumours while allowing quantification of diversity-related endogenous mutagenic activity since we expect it to be higher in cancers driven by RAG / APOBEC activity 

1) Counts SVs in VDJ regions of IG (B cell) and TCR (T cell) loci (regions from cider). 
2) Quantify Somatic Hypermutation in IG and TCR loci (only occurs in B cells)
3) SV mutation count


# Installation

See latest release for installation instructions.

# Quick Start

ScarScape extracts features from common mutation file formats. We require a manifest csv with the following column:

1) **sample**: sample identifier. Must be the same identifier used in snv VCFs
2) **snv**: path to a VCF file describing SNVs and INDELs
3) **sv**: path to a VCF file describing structural variants (1 entry per breakend).
4) **cnv**: path to a TSV describing segment copy number. Must contain the columns: chromosome,start (1-based), end (inclusive), copyNumber, minorAlleleCopyNumber

VCFs must be bgzipped and indexed (`tabix -p vcf path/to/vcf.gz`)

You need at least one of these files or the results won't be computed.
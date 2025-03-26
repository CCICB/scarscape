# ScarScape

> \[!WARNING\]  
> This package is in early development and not ready for use

ScarScape extracts the landscape of tumorigenic scars from cancer omics data. For each tumour, it processes common biological data types, classifying each instance of omic damage by its affected loci, sequence context, clustering pattern, and other contextual features. It then reports the frequency of these distinct damage types. By carefully defining how mutations are classified and counted, ScarScape aims to clearly reveal the active mechanisms driving tumorigenesis.

---


## Installation

Refer to the latest release for installation instructions.

---

## Quick Start

Prepare a [manifest](testfiles/manifest.tsv) TSV with the following columns:

- **sample:** Unique sample identifier (must match SNV VCF identifiers)
- **snv:** Path to SNV/INDEL VCF (bgzipped and indexed via `tabix -p vcf path/to/vcf.gz`)
- **sv:** Path to structural variant VCF (one entry per breakend)
- **cnv:** Path to copy number TSV (must include: `chromosome`, `start` (1-based), `end` (inclusive), `copyNumber`, `minorAlleleCopyNumber`)

*Note:* Each sample must have at least one file type specified (snv, sv, or cnv).

```{bash}
scarseek --manifest <path_to_manifest.tsv> --genome hg38
```


---

## Metrics

The metrics captured by ScarScape are carefully chosen to reflect diverse tumorigenic mechanisms, so that downstream analyses stand the best chance of mapping these features back to their aetiological processes. 


### Endogenous Rearrangement Hyperactivity

To characterise B & T cell-derived tumors, which may be driven by endogenous mutagenic enzymes (e.g., RAG for VDJ rearrangement, APOBEC for hypermutation), ScarScape quantifies:

- **SV counts in VDJ regions** to (hopefully) reveal bias towards IG loci in B cell tumours and TCR loci in T cells)
- **Somatic hypermutation** of IG/TCR loci (a B cells-specific process)
- **Total SV mutation count** for normalisation and detection of SV hypermutated samples.

## Output

Each sample produces an `<sample>.svcounts.csv` file with the following columns:

| Column | Description                  |
| ------ | ---------------------------- |
| sample | Sample identifier            |
| total  | Total SV count (pass + fail) |
| pass   | Count of pass SVs            |
| fail   | Count of fail SVs            |
| igtcr  | SVs in IG or TCR regions     |
| igh    | SVs in IGH loci (BCR)        |
| igk    | SVs in IGK loci (BCR)        |
| igl    | SVs in IGL loci (BCR)        |
| tra    | SVs in TRA loci (TCR)        |
| trb    | SVs in TRB loci (TCR)        |
| trd    | SVs in TRD loci (TCR)        |
| trg    | SVs in TRG loci (TCR)        |
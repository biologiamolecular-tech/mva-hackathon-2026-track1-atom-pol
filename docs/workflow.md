# Track 1 analysis workflow

## Scope

This workflow prioritizes candidate small variants from a single-sample
GRCh38 whole-genome VCF using variant annotation, population frequency,
predicted molecular consequence, phenotype associations and Exomiser
phenotype-driven prioritization.

The submitted prediction represents a candidate compound-heterozygous
hypothesis. Read-backed phasing was attempted, but the cis/trans relationship
could not be resolved with the available short-read data.

## Data governance

The analysis used controlled genomic and phenotype data. These inputs are not
distributed in this repository.

Excluded controlled files include:

- FASTQ and BAM files;
- original and annotated VCF files;
- individual-level phenotype documents;
- HPO term lists and phenopackets;
- variant-level private review tables;
- sample identifiers;
- reference databases and software caches.

Only the final challenge CSV, an analysis report, workflow code and aggregate
non-identifying summaries are provided.

## Workflow overview

```mermaid
flowchart TD
    A["Controlled GRCh38 VCF"] --> B["VCF validation and QC"]
    B --> C["Normalization and PASS selection"]
    C --> D["VEP 116 annotation"]
    D --> E["Frequency and consequence filtering"]
    F["Controlled phenotype terms"] --> G["HPO gene and disease matching"]
    E --> H["Core candidate set"]
    G --> H
    H --> I["Exomiser 15.1.0 prioritization"]
    I --> J["Evidence integration and manual review"]
    J --> K["Read-backed phasing attempts"]
    K --> L["Track 1 prediction CSV"]

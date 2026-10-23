# Workflow source code

This directory contains selected scripts used to generate the aggregate
results and final Track 1 prioritization.

The scripts are provided for methodological transparency and provenance. They
preserve the analysis logic used during the project, but require controlled
genomic and phenotype inputs that are not distributed in this repository.

## Execution order

### 1. VEP annotation profiling

Run `bash workflow/vep/run_annotation_profile.sh`.

This script summarizes VEP consequence, impact, population-frequency and
ClinVar fields from the annotated PASS VCF.

### 2. Candidate-tier counting

Run `bash workflow/vep/run_candidate_tier_counts.sh`.

This script counts rare HIGH-impact, damaging MODERATE-impact, splice-related
and ClinVar-supported candidates.

### 3. Candidate extraction

Run `bash workflow/vep/run_candidate_extraction.sh`.

This script extracts the variant-level candidate table and creates an
aggregate candidate summary. The variant-level output is private and is not
included here.

### 4. Direct phenotype intersection

Run `python3 workflow/phenotype/direct_hpo_candidate_intersection.py`.

This script intersects candidate genes with phenotype-to-gene and
phenotype-to-disease associations. Individual phenotype terms and
variant-level outputs remain private.

### 5. Disease-level phenotype coherence

Run `python3 workflow/phenotype/disease_coherence.py`.

This script evaluates whether multiple exact phenotype matches are supported
within the same gene–disease model.

### 6. Exomiser intersection

Run `python3 workflow/exomiser/intersect_core_exomiser.py`.

This script matches phenotype-derived central candidates against the Exomiser
variant-ranking output.

## Required controlled inputs

Depending on the stage, execution requires authorized access to:

- the normalized and VEP-annotated GRCh38 VCF;
- case, personal-history and family-history phenotype term files;
- the HPO `phenotype_to_genes.txt` association file;
- variant-level candidate tables;
- Exomiser variant-ranking output.

These inputs are intentionally excluded from the public repository.

## Software requirements

The scripts require the corresponding tools listed in
`config/software_versions.txt`, including:

- Bash and GNU-compatible command-line utilities;
- bcftools with the `split-vep` plugin;
- Python 3;
- VEP-annotated input generated with Ensembl VEP 116.1.

## Reproducibility boundary

The scripts retain the original relative project structure used during the
analysis. They are not a bundled dataset and are not intended to download or
reconstruct controlled participant data.

Syntax was validated before publication. Reproduction of participant-level
results requires the original authorized inputs and reference resources.

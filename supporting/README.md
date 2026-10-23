# Supporting aggregate results

This directory contains aggregate summaries generated during the Track 1
variant-prioritization workflow.

The files document intermediate filtering, phenotype matching, Exomiser
integration, evidence auditing and final prioritization. They are included to
support transparency and reproducibility of the submitted prediction.

## Contents

- `candidate_extraction_summary.tsv`: candidate counts after VEP-based extraction.
- `core_candidates_summary.tsv`: summary of the central candidate set.
- `direct_hpo_intersection_summary.tsv`: direct phenotype-to-gene intersection counts.
- `disease_coherence_summary.tsv`: phenotype coherence within gene–disease models.
- `core_candidates.exomiser.summary.tsv`: overlap with the genome-wide Exomiser analysis.
- `core25_FULL_intersection.summary.tsv`: targeted Exomiser evaluation of the core candidates.
- `filter_counts.summary.tsv`: Exomiser filter counts.
- `failed_filters_by_tier.summary.tsv`: failed filters grouped by candidate tier.
- `core_25_filter_mapping.summary.tsv`: mapping of failed Exomiser evaluations to core candidates.
- `core_25_evidence_audit.summary.tsv`: aggregate evidence audit.
- `core_25_master_review.summary.tsv`: master-review group counts.
- `core_25_overlap_and_pairs.summary.tsv`: candidate overlap and multi-variant gene counts.
- `core_25_priority_bands.summary.tsv`: final priority-band counts.

## Data governance

Only aggregate counts and non-identifying workflow metrics are included.

This directory does not contain:

- raw sequencing data;
- BAM, VCF or FASTQ files;
- sample identifiers;
- individual HPO terms;
- phenotype documents;
- variant-level private review tables;
- absolute filesystem paths.

The controlled genomic and clinical inputs are intentionally not distributed.

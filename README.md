# MVA Hackathon 2026 — Track 1

Phenotype-guided prioritization of candidate variants from a single-sample
GRCh38 whole-genome VCF.

- Team / display name: `atom_pol`
- Hugging Face username: `PBIOL`
- Affiliation: CEPROCOR
- Genome assembly: GRCh38
- Submission track: Track 1

## Objective

The objective was to identify and prioritize variants potentially associated
with the supplied phenotype using functional annotation, population
frequency, phenotype associations, inheritance models and Exomiser
prioritization.

The analysis is intended for research and hackathon evaluation. It is not a
clinical diagnosis or a substitute for clinical variant interpretation.

## Primary result

The submitted result contains one ranked primary hypothesis involving two
variants in the same gene under a potential compound-heterozygous model.

Both variants were detected as heterozygous. Read-backed phasing with
WhatsHap and HapCUT2 did not determine whether they are in cis or in trans
because no sequencing template directly bridged both sites.

The EPCR value in the submission CSV represents challenge-ranking confidence,
not a clinically calibrated probability.

## Team

**Team name:** `atom_pol`  
**Track 1 submission account:** [`PBIOL`](https://huggingface.co/PBIOL)

| Contributor | Hugging Face | GitHub |
|---|---|---|
| Guillermo Giaj Merlera | [`ggiajmerlera`](https://huggingface.co/ggiajmerlera) | [`guillermogiajmerlera`](https://github.com/guillermogiajmerlera) |
| Pablo Sebastian Velez | [`polface`](https://huggingface.co/polface) | [`poldohko`](https://github.com/poldohko) |

The submission was prepared collaboratively under the team name `atom_pol`. The individual profiles listed above identify the team members, while the Track 1 prediction was submitted through the `PBIOL` Hugging Face account.

## Repository contents

- `results/atom_pol_phenotype-exomiser-v1.csv`: official Track 1 prediction.
- `docs/atom_pol_track1_report.md`: detailed analysis report.
- `docs/workflow.md`: description of the end-to-end analysis workflow.
- `workflow/`: selected analysis scripts retained for provenance.
- `supporting/`: aggregate, non-identifying intermediate summaries.
- `scripts/`: scripts used to construct the public submission and report.
- `config/software_versions.txt`: principal software and reference versions.
- `LICENSE`: CC BY 4.0 license.

## Analysis overview

The workflow comprised:

1. VCF integrity, index and genotype-level quality control.
2. GRCh38 reference validation and variant normalization.
3. PASS variant selection.
4. Offline annotation with Ensembl VEP 116.1.
5. Frequency, consequence, ClinVar and splice-based candidate extraction.
6. Direct HPO-to-gene matching and disease-level phenotype coherence.
7. Exomiser 15.1.0 phenotype-driven prioritization.
8. Evidence integration and candidate review.
9. Read-backed phasing attempts with WhatsHap and HapCUT2.
10. Construction and validation of the Track 1 CSV.

Detailed counts and limitations are provided in the report and in the
aggregate summaries under `supporting/`.

## Reproducibility

The public scripts document the analysis and packaging logic. Exact
participant-level reproduction requires authorized access to the controlled
VCF, BAM and phenotype inputs, which are not included.

The final public artifacts can be rebuilt from authorized private intermediate
tables using the scripts described in `scripts/README.md`.

## Data governance

This repository intentionally excludes:

- FASTQ, BAM, CRAM and VCF files;
- indexes derived from controlled genomic data;
- individual HPO term lists and phenotype documents;
- phenopackets;
- sample identifiers;
- private variant-level review tables;
- software caches and reference databases;
- execution logs containing local paths.

Only the challenge submission, documentation, selected source code and
aggregate non-identifying summaries are public.

## Artificial-intelligence assistance

OpenAI ChatGPT/Codex was used for shell-script drafting, troubleshooting,
workflow organization and documentation support. Scientific decisions and
final result selection were reviewed by the participant.

The ChatGPT setting “Improve the model for everyone” was disabled.

## License

The public code, documentation and aggregate results are released under the
Creative Commons Attribution 4.0 International license. See `LICENSE`.

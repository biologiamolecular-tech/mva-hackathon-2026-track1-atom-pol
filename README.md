# MVA Hackathon 2026 — Track 1

Track 1 submission from **atom_pol** for the Rare Disease, Real Kid:
MVA Hackathon 2026.

- Hugging Face username: `PBIOL`
- Affiliation: CEPROCOR
- Contact: biologiamolecular@ceprocor.uncor.edu
- Genome assembly: GRCh38

## Objective

This project applies a phenotype-guided, genome-wide workflow to prioritize
candidate variants potentially associated with the phenotype supplied by the
hackathon organizers.

The workflow integrated:

- VCF validation and normalization;
- variant quality filtering;
- Ensembl VEP functional annotation;
- population-frequency filtering;
- ClinVar evidence;
- HPO phenotype–gene and disease-model coherence;
- Exomiser phenotype-based prioritization;
- read-backed phasing attempts with WhatsHap and HapCUT2.

## Primary result

The submitted primary hypothesis consists of two heterozygous variants in the
same candidate gene, represented as a single candidate pair in the Track 1 CSV.

Short-read phasing did not determine whether the variants are in cis or in
trans. The pair is therefore reported as a potential recessive candidate model,
not as confirmed compound heterozygosity.

See [`docs/atom_pol_track1_report.md`](docs/atom_pol_track1_report.md) for the
complete methods, evidence and limitations.

## Repository contents

- `results/atom_pol_phenotype-exomiser-v1.csv`: official Track 1 prediction.
- `docs/atom_pol_track1_report.md`: methods and interpretation report.
- `scripts/build_track1_submission.py`: deterministic construction of the
  submission CSV from the private reviewed candidate table.
- `scripts/build_track1_report.py`: deterministic construction of the report.
- `config/software_versions.txt`: software and database versions.

## Reproducing the final packaging step

Authorized users with access to the private reviewed candidate table can run:

```bash
python3 scripts/build_track1_submission.py \
    --input /path/to/P1_primary_review.private.tsv \
    --output results/atom_pol_phenotype-exomiser-v1.csv \
    --epcr 0.90

python3 scripts/build_track1_report.py \
    --input /path/to/P1_primary_review.private.tsv \
    --output docs/atom_pol_track1_report.md

# Scripts

## build_track1_submission.py

Selects `PRIORITY_ORDER=1` and `PRIORITY_ORDER=3` from the private
reviewed candidate table.

The script validates that both candidates are heterozygous and located
on the same chromosome, adds the required `chr` prefix and creates the
official Track 1 submission CSV.

## build_track1_report.py

Reads the same private reviewed candidate table and constructs the
methods and evidence report.

Exact private input tables and controlled patient data are not
distributed in this repository.

## Scope of reproducibility

These scripts reproduce the deterministic final packaging step.

The preceding genome-wide workflow included VCF normalization, VEP
annotation, phenotype association, candidate prioritization, Exomiser
analysis and read-backed phasing. Its parameters, software versions,
results and limitations are documented in the report.

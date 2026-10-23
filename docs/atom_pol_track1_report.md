# MVA Hackathon 2026 — Track 1 Variant Prioritization

**Participant / team:** atom_pol  
**Hugging Face username:** PBIOL  
**Affiliation:** CEPROCOR  
**Contact:** biologiamolecular@ceprocor.uncor.edu  
**Submission:** PBIOL phenotype-guided VEP/Exomiser workflow  
**Date:** 2026-10-23  
**Genome assembly:** GRCh38  

## 1. Scope and clinical boundary

This report describes a research analysis performed for Track 1 of the
Rare Disease, Real Kid: MVA Hackathon 2026. It is not a clinical diagnostic
report and must not be interpreted as medical advice.

The objective was to prioritize variants potentially associated with the
provided phenotype using genome-wide quality control, functional annotation,
population frequency, phenotype–gene associations, disease-model coherence
and Exomiser prioritization.

## 2. Input data and quality control

The analysis started from the controlled-access single-sample WGS VCF supplied
by the organizers.

Verified processing statistics were:

- 5,012,204 records in the original VCF.
- 4,962,046 normalized records on primary reference contigs.
- 4,727,733 normalized PASS records used for genome-wide annotation.
- No remaining multiallelic or symbolic records in the normalized primary VCF.
- No reference-allele mismatches in the reference audit.
- Complete DP, GQ and AD fields in the normalized primary callset.
- GRCh38 reference-contig compatibility was verified before annotation.

Raw FASTQ, BAM, VCF, clinical documents and other controlled-access inputs are
not included in the public submission repository.

## 3. Functional annotation

Variants were annotated offline with Ensembl Variant Effect Predictor
VEP 116.1 using the GRCh38 Ensembl cache version 116.

Annotations included:

- consequence and predicted impact;
- gene and transcript;
- HGVS coding and protein notation;
- MANE and canonical transcript information;
- SIFT and PolyPhen predictions;
- ClinVar annotations available in the VEP cache;
- gnomAD exome and genome allele frequencies;
- maximum observed population frequency.

The annotated VCF contained 4,727,733 records with a populated CSQ field.

## 4. Phenotype-driven prioritization

Five active HPO terms representing the case phenotype were used for primary
phenotype prioritization. Personal and family history terms were retained as
contextual evidence and were not treated as equivalent to manifestations in
the proband.

The prioritization procedure produced:

- 367 rare or potentially relevant candidate rows.
- 281 unique candidate gene symbols.
- 92 candidates with a direct case-HPO gene association.
- 25 candidates supported by at least two exact HPO terms within a coherent
  gene–disease model.
- 21 unique genes among those 25 candidates.

The HPO resource used was the release dated 2026-06-23.

## 5. Exomiser analysis

Exomiser CLI 15.1.0 was run against the GRCh38 genome resource and phenotype
resource release 2602.

The genome-wide PASS-only analysis processed 4,727,733 VCF records and produced:

- 421 reported gene rows.
- 794 reported variant rows.
- 364 genes containing 589 filtered variants.
- Seven exact overlaps between the 25 central candidates and Exomiser
  contributing variants.

A targeted FULL-mode analysis of the 25 central candidates was additionally
used to audit filters and retain candidates requiring inheritance, frequency,
technical or transcript-level review.

Automated Exomiser ACMG classifications were treated as computational evidence,
not as independently curated clinical classifications.

## 6. Primary candidate hypothesis

The highest-priority hypothesis is a pair of heterozygous variants in
**BUB1B**.

### 6.1 Variant identities

| Candidate | Tier | Chromosome | Position | REF | ALT | Consequence | HGVSc | HGVSp | Genotype |
|---|---|---|---:|---|---|---|---|---|---|
| P1_01 | T1A_CLINVAR_PLP | chr15 | 40209701 | T | G | stop_gained | ENST00000287598.11:c.2210T>G | ENSP00000287598.7:p.Leu737Ter | 0/1 |
| P1_03 | T3_MODERATE_DAMAGING_RARE | chr15 | 40220612 | T | G | missense_variant | ENST00000287598.11:c.3006T>G | ENSP00000287598.7:p.Asn1002Lys | 0/1 |

### 6.2 Annotation and phenotype evidence

| Candidate | MAX_AF | ClinVar/VEP significance | SIFT | PolyPhen | MANE | Exomiser automated ACMG | Exact HPO matches |
|---|---|---|---|---|---|---|---:|
| P1_01 | 9.982e-05 | pathogenic&likely_pathogenic | . | . | NM_001211.6 | PATHOGENIC | 2 |
| P1_03 | . | . | deleterious(0.01) | probably_damaging(0.997) | NM_001211.6 | NOT_AVAILABLE;UNCERTAIN_SIGNIFICANCE | 2 |

### 6.3 Exomiser evidence

| Candidate | WGS rank | Exomiser MOI | Phenotype score | Variant score | Automated ACMG evidence |
|---|---:|---|---:|---:|---|
| P1_01 | 1 | AD;AR | 0.5510;0.7750 | 1.0000 | PVS1,PM2_Supporting,PP4_Moderate,PP5_Strong;PVS1,PP5_Strong |
| P1_03 | 1 | AD;AR | 0.5510;0.7750 | 0.9229 | PM2_Supporting,PP4_Moderate,BP1 |

The pair was prioritized because both variants:

- are heterozygous;
- occur in the same candidate gene;
- are rare or absent at appreciable frequency in the population resources
  interrogated;
- are supported by the case phenotype;
- participate in the same recessive Exomiser gene–disease hypothesis;
- passed the basic sequencing-quality review.

P1_01 and P1_03 were therefore submitted together as one primary
compound-heterozygous hypothesis. This pairing represents a ranked research
hypothesis and not confirmed compound heterozygosity.

## 7. Read-backed phase investigation

Phase was investigated using the available short-read BAM and a regional VCF.

### WhatsHap

- WhatsHap version 2.8 completed successfully.
- The regional VCF contained 323 records.
- Both candidate variants were found.
- Neither candidate received a resolved phased genotype shared with the other.
- No shared phase set connected the two variants.

### HapCUT2

- `extractHAIRS` completed successfully and produced 1,374 fragment rows.
- HapCUT2 completed successfully.
- Both candidate variants were present in the output.
- The pair remained outside a common resolved haplotype block.

### Direct read evidence

- 52 templates covered one candidate position.
- 30 templates covered the other candidate position.
- Zero templates bridged both candidate positions directly.

Consequently, the available short-read data did not establish whether the two
variants are in cis or in trans.

The final phase interpretation is:

> **Phase undetermined.**

Parental segregation, long-read sequencing, targeted long-range sequencing or
another validated orthogonal method would be required to establish phase.

## 8. Interpretation

The P1_01–P1_03 pair is the highest-priority causal hypothesis from this
workflow based on phenotype coherence, rarity, functional annotation,
ClinVar evidence where available and Exomiser prioritization.

The evidence supports reporting the pair as a potential recessive candidate
model. It does not establish:

- that the variants are in trans;
- that the missense effect is functionally confirmed;
- that both variants are clinically confirmed;
- or that the pair is definitively causal.

## 9. Limitations

Principal limitations include:

1. Analysis of a single proband.
2. Absence of parental genotypes and segregation data.
3. Absence of direct read templates spanning both candidate sites.
4. Dependence on the database and software versions available during analysis.
5. Use of in-silico pathogenicity predictions.
6. Automated ACMG evidence that was not independently curated under a clinical
   laboratory framework.
7. No orthogonal molecular confirmation.

## 10. Reproducibility and software

Principal software and resources:

| Component | Version or release |
|---|---|
| bcftools | 1.24 |
| samtools | 1.24 |
| Ensembl VEP | 116.1 |
| VEP cache | Homo sapiens 116, GRCh38 |
| Exomiser CLI | 15.1.0 |
| Exomiser data | 2602 phenotype and 2602 hg38 |
| HPO ontology | 2026-06-23 |
| WhatsHap | 2.8 |
| HapCUT2 | Version not independently recorded |

Public scripts contain no raw patient data and expect authorized participants
to obtain controlled inputs from the official challenge dataset.

## 11. Use of artificial intelligence

OpenAI ChatGPT/Codex was used under a ChatGPT Plus plan. In the account Data
Controls, **“Improve the model for everyone” was disabled**.

The assistant was used for shell-script drafting, troubleshooting and workflow
organization. Scientific review, interpretation of the evidence and submission
decisions remained the responsibility of the participant.

## 12. Data governance

The public repository intentionally excludes raw FASTQ, BAM, VCF, index files,
the original phenotype document, private HPO lists, complete Exomiser outputs
and private analysis logs.

Only the final ranked prediction, aggregate methods, selected reproducible
scripts and this report are intended for public release under CC BY 4.0.

## 13. Conclusion

The submitted P1_01–P1_03 pair is the top-ranked result of a genome-wide,
phenotype-guided workflow. The combination is supported as a potential
recessive disease model, but phase remains unresolved and the result requires
orthogonal confirmation before any clinical use.

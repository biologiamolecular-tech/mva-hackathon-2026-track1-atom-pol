#!/usr/bin/env python3

import csv
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

csv.field_size_limit(sys.maxsize)

project = Path.home() / "MVA_Hackathon_2026"

p2g_file = project / "reference/hpo/2026-09-03/phenotype_to_genes.txt"
hpo_file = project / "05_analysis/phenotype/case_hpo.private.txt"
candidate_file = project / "05_analysis/phenotype/candidate_variants.hpo_direct.private.tsv"

private_output = project / "05_analysis/phenotype/candidate_variants.disease_coherence.private.tsv"
summary_output = project / "05_analysis/phenotype/disease_coherence_summary.tsv"


def read_hpo(path):
    with path.open(encoding="utf-8") as handle:
        return {
            line.strip().split()[0]
            for line in handle
            if line.strip() and not line.startswith("#")
        }


def symbols(value):
    if not value or value in {".", "-"}:
        return set()

    return {
        item.strip().upper()
        for item in re.split(r"[,;&]", value)
        if item.strip() and item.strip() not in {".", "-"}
    }


for path in [p2g_file, hpo_file, candidate_file]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty input: {path}")


case_hpo = read_hpo(hpo_file)

with candidate_file.open(encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    original_fields = list(reader.fieldnames or [])
    candidates = list(reader)

if not {"TIER", "SYMBOL"}.issubset(original_fields):
    raise ValueError(f"Unexpected candidate header: {original_fields}")


candidate_symbols = set()

for row in candidates:
    candidate_symbols.update(symbols(row.get("SYMBOL", "")))


# Clave: (gen, enfermedad); valor: HPO del caso asociados a esa enfermedad.
models = defaultdict(set)
public_rows = 0
relevant_rows = 0

with p2g_file.open(encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    required = {"hpo_id", "gene_symbol", "disease_id"}

    if not required.issubset(set(reader.fieldnames or [])):
        raise ValueError(f"Unexpected HPO-gene header: {reader.fieldnames}")

    for row in reader:
        public_rows += 1

        hpo_id = row["hpo_id"].strip()
        symbol = row["gene_symbol"].strip().upper()
        disease = row["disease_id"].strip()

        if (
            hpo_id in case_hpo
            and symbol in candidate_symbols
            and disease
        ):
            models[(symbol, disease)].add(hpo_id)
            relevant_rows += 1


models_by_gene = defaultdict(list)

for (symbol, disease), matched_hpo in models.items():
    models_by_gene[symbol].append((disease, matched_hpo))


extra_fields = [
    "BEST_DISEASE_EXACT_HPO_COUNT",
    "MATCHING_GENE_DISEASE_MODELS",
    "BEST_GENE_DISEASE_MODELS",
    "BEST_DISEASE_EXACT_HPO_IDS",
]

distribution = Counter()
tier_distribution = Counter()

ge1 = 0
ge2 = 0
ge3 = 0
cross_disease_only = 0

genes_ge1 = set()
genes_ge2 = set()
genes_ge3 = set()

with private_output.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=original_fields + extra_fields,
        delimiter="\t",
        extrasaction="ignore",
    )

    writer.writeheader()

    for row in candidates:
        candidate_models = []

        for symbol in symbols(row.get("SYMBOL", "")):
            for disease, matched_hpo in models_by_gene.get(symbol, []):
                candidate_models.append(
                    (len(matched_hpo), symbol, disease, matched_hpo)
                )

        best_count = max(
            (item[0] for item in candidate_models),
            default=0,
        )

        best_models = {
            f"{symbol}|{disease}"
            for count, symbol, disease, matched_hpo in candidate_models
            if count == best_count and count > 0
        }

        best_hpo = set()
        best_symbols = set()

        for count, symbol, disease, matched_hpo in candidate_models:
            if count == best_count and count > 0:
                best_hpo.update(matched_hpo)
                best_symbols.add(symbol)

        row["BEST_DISEASE_EXACT_HPO_COUNT"] = best_count
        row["MATCHING_GENE_DISEASE_MODELS"] = len(candidate_models)
        row["BEST_GENE_DISEASE_MODELS"] = ",".join(sorted(best_models))
        row["BEST_DISEASE_EXACT_HPO_IDS"] = ",".join(sorted(best_hpo))

        writer.writerow(row)

        tier = row["TIER"]
        distribution[best_count] += 1
        tier_distribution[(tier, best_count)] += 1

        if best_count >= 1:
            ge1 += 1
            genes_ge1.update(best_symbols)

        if best_count >= 2:
            ge2 += 1
            genes_ge2.update(best_symbols)

        if best_count >= 3:
            ge3 += 1
            genes_ge3.update(best_symbols)

        direct_count = int(
            row.get("DIRECT_CASE_HPO_MATCH_COUNT", "0") or 0
        )

        if direct_count >= 2 and best_count < 2:
            cross_disease_only += 1


with summary_output.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")

    writer.writerow(["section", "criterion", "count"])

    writer.writerow(["REFERENCE", "public_rows_scanned", public_rows])
    writer.writerow(["REFERENCE", "relevant_association_rows", relevant_rows])

    writer.writerow(["CANDIDATE", "total", len(candidates)])
    writer.writerow(["CANDIDATE", "same_disease_exact_HPO_ge_1", ge1])
    writer.writerow(["CANDIDATE", "same_disease_exact_HPO_ge_2", ge2])
    writer.writerow(["CANDIDATE", "same_disease_exact_HPO_ge_3", ge3])
    writer.writerow([
        "CANDIDATE",
        "multi_HPO_gene_match_but_not_same_disease",
        cross_disease_only,
    ])

    writer.writerow(["GENE", "same_disease_exact_HPO_ge_1", len(genes_ge1)])
    writer.writerow(["GENE", "same_disease_exact_HPO_ge_2", len(genes_ge2)])
    writer.writerow(["GENE", "same_disease_exact_HPO_ge_3", len(genes_ge3)])

    for count in sorted(distribution):
        writer.writerow([
            "BEST_EXACT_HPO_COUNT",
            count,
            distribution[count],
        ])

    for tier, count in sorted(tier_distribution):
        writer.writerow([
            f"TIER_{tier}",
            f"best_exact_HPO_count_{count}",
            tier_distribution[(tier, count)],
        ])


os.chmod(private_output, 0o600)
os.chmod(summary_output, 0o600)

print(f"public_rows_scanned={public_rows}")
print(f"candidate_rows={len(candidates)}")
print(f"same_disease_exact_HPO_ge_1={ge1}")
print(f"same_disease_exact_HPO_ge_2={ge2}")
print(f"same_disease_exact_HPO_ge_3={ge3}")
print(f"summary_exists={'YES' if summary_output.stat().st_size else 'NO'}")
print("final_status=COMPLETED")

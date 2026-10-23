#!/usr/bin/env python3

import csv
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

csv.field_size_limit(sys.maxsize)

project = Path.home() / "MVA_Hackathon_2026"

reference = project / "reference/hpo/2026-09-03/phenotype_to_genes.txt"
phenotype_dir = project / "05_analysis/phenotype"

case_file = phenotype_dir / "case_hpo.private.txt"
personal_file = phenotype_dir / "personal_history_hpo.private.txt"
family_file = phenotype_dir / "family_history_hpo.private.txt"

candidates_file = (
    project /
    "05_analysis/vep_full/candidate_variants.private.tsv"
)

gene_table = phenotype_dir / "case_hpo_gene_panel.private.tsv"
annotated_candidates = phenotype_dir / "candidate_variants.hpo_direct.private.tsv"
matched_candidates = phenotype_dir / "candidate_variants.hpo_direct_matches.private.tsv"
summary_file = phenotype_dir / "direct_hpo_intersection_summary.tsv"


def read_hpo_file(path):
    terms = set()

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            terms.add(line.split()[0])

    return terms


def new_gene_record():
    return {
        "symbol": "",
        "ncbi_ids": set(),
        "case_hpo": set(),
        "personal_hpo": set(),
        "family_hpo": set(),
        "case_disease": set(),
        "personal_disease": set(),
        "family_disease": set(),
    }


def split_symbols(value):
    if not value or value in {".", "-"}:
        return set()

    return {
        item.strip().upper()
        for item in re.split(r"[,;&]", value)
        if item.strip() and item.strip() not in {".", "-"}
    }


for required_file in [
    reference,
    case_file,
    personal_file,
    family_file,
    candidates_file,
]:
    if not required_file.is_file() or required_file.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty input: {required_file}")


case_hpo = read_hpo_file(case_file)
personal_hpo = read_hpo_file(personal_file)
family_hpo = read_hpo_file(family_file)

query_sets = {
    "case": case_hpo,
    "personal": personal_hpo,
    "family": family_hpo,
}

genes = defaultdict(new_gene_record)
association_rows = Counter()
public_rows_scanned = 0

with reference.open(encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    required_columns = {
        "hpo_id",
        "ncbi_gene_id",
        "gene_symbol",
        "disease_id",
    }

    if not required_columns.issubset(set(reader.fieldnames or [])):
        raise ValueError(
            "Unexpected phenotype_to_genes.txt columns: "
            + repr(reader.fieldnames)
        )

    for row in reader:
        public_rows_scanned += 1

        hpo_id = row["hpo_id"].strip()
        symbol_original = row["gene_symbol"].strip()

        if not symbol_original:
            continue

        matched_categories = [
            category
            for category, terms in query_sets.items()
            if hpo_id in terms
        ]

        if not matched_categories:
            continue

        symbol_key = symbol_original.upper()
        record = genes[symbol_key]

        record["symbol"] = symbol_original

        ncbi_id = row["ncbi_gene_id"].strip()
        disease_id = row["disease_id"].strip()

        if ncbi_id:
            record["ncbi_ids"].add(ncbi_id)

        for category in matched_categories:
            record[f"{category}_hpo"].add(hpo_id)

            if disease_id:
                record[f"{category}_disease"].add(disease_id)

            association_rows[category] += 1


with gene_table.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")

    writer.writerow([
        "gene_symbol",
        "ncbi_gene_ids",
        "case_hpo_match_count",
        "personal_history_match_count",
        "family_history_match_count",
        "case_disease_count",
        "personal_history_disease_count",
        "family_history_disease_count",
        "case_hpo_ids",
        "personal_history_hpo_ids",
        "family_history_hpo_ids",
    ])

    ordered_genes = sorted(
        genes.values(),
        key=lambda record: (
            -len(record["case_hpo"]),
            -len(record["personal_hpo"]),
            -len(record["family_hpo"]),
            record["symbol"],
        ),
    )

    for record in ordered_genes:
        writer.writerow([
            record["symbol"],
            ",".join(sorted(record["ncbi_ids"])),
            len(record["case_hpo"]),
            len(record["personal_hpo"]),
            len(record["family_hpo"]),
            len(record["case_disease"]),
            len(record["personal_disease"]),
            len(record["family_disease"]),
            ",".join(sorted(record["case_hpo"])),
            ",".join(sorted(record["personal_hpo"])),
            ",".join(sorted(record["family_hpo"])),
        ])


case_gene_universe = {
    symbol
    for symbol, record in genes.items()
    if record["case_hpo"]
}

personal_gene_universe = {
    symbol
    for symbol, record in genes.items()
    if record["personal_hpo"]
}

family_gene_universe = {
    symbol
    for symbol, record in genes.items()
    if record["family_hpo"]
}

tier_total = Counter()
tier_primary_match = Counter()
tier_context_only = Counter()

candidate_total = 0
primary_match_rows = 0
context_only_rows = 0
unique_candidate_genes = set()
unique_primary_matched_genes = set()

with candidates_file.open(encoding="utf-8", newline="") as source:
    reader = csv.DictReader(source, delimiter="\t")

    if not {"TIER", "SYMBOL"}.issubset(set(reader.fieldnames or [])):
        raise ValueError(
            "Unexpected candidate table columns: "
            + repr(reader.fieldnames)
        )

    extra_fields = [
        "DIRECT_CASE_HPO_MATCH",
        "DIRECT_CASE_HPO_MATCH_COUNT",
        "DIRECT_PERSONAL_HISTORY_MATCH_COUNT",
        "DIRECT_FAMILY_HISTORY_MATCH_COUNT",
        "DIRECT_CASE_HPO_IDS",
        "DIRECT_PERSONAL_HISTORY_HPO_IDS",
        "DIRECT_FAMILY_HISTORY_HPO_IDS",
    ]

    output_fields = list(reader.fieldnames) + extra_fields

    with annotated_candidates.open(
        "w", encoding="utf-8", newline=""
    ) as all_output, matched_candidates.open(
        "w", encoding="utf-8", newline=""
    ) as matched_output:

        all_writer = csv.DictWriter(
            all_output,
            fieldnames=output_fields,
            delimiter="\t",
            extrasaction="ignore",
        )

        matched_writer = csv.DictWriter(
            matched_output,
            fieldnames=output_fields,
            delimiter="\t",
            extrasaction="ignore",
        )

        all_writer.writeheader()
        matched_writer.writeheader()

        for row in reader:
            candidate_total += 1

            tier = row["TIER"]
            tier_total[tier] += 1

            symbols = split_symbols(row.get("SYMBOL", ""))
            unique_candidate_genes.update(symbols)

            case_matches = set()
            personal_matches = set()
            family_matches = set()

            for symbol in symbols:
                record = genes.get(symbol)

                if not record:
                    continue

                case_matches.update(record["case_hpo"])
                personal_matches.update(record["personal_hpo"])
                family_matches.update(record["family_hpo"])

            primary_match = bool(case_matches)
            context_match = bool(personal_matches or family_matches)

            row["DIRECT_CASE_HPO_MATCH"] = (
                "YES" if primary_match else "NO"
            )

            row["DIRECT_CASE_HPO_MATCH_COUNT"] = len(case_matches)
            row["DIRECT_PERSONAL_HISTORY_MATCH_COUNT"] = len(personal_matches)
            row["DIRECT_FAMILY_HISTORY_MATCH_COUNT"] = len(family_matches)

            row["DIRECT_CASE_HPO_IDS"] = ",".join(sorted(case_matches))
            row["DIRECT_PERSONAL_HISTORY_HPO_IDS"] = ",".join(
                sorted(personal_matches)
            )
            row["DIRECT_FAMILY_HISTORY_HPO_IDS"] = ",".join(
                sorted(family_matches)
            )

            all_writer.writerow(row)

            if primary_match:
                primary_match_rows += 1
                tier_primary_match[tier] += 1
                unique_primary_matched_genes.update(
                    symbols & case_gene_universe
                )
                matched_writer.writerow(row)

            elif context_match:
                context_only_rows += 1
                tier_context_only[tier] += 1


case_terms_with_associations = set()
personal_terms_with_associations = set()
family_terms_with_associations = set()

for record in genes.values():
    case_terms_with_associations.update(record["case_hpo"])
    personal_terms_with_associations.update(record["personal_hpo"])
    family_terms_with_associations.update(record["family_hpo"])


with summary_file.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")

    writer.writerow(["section", "criterion", "count"])

    writer.writerow(["REFERENCE", "public_rows_scanned", public_rows_scanned])

    writer.writerow(["PHENOTYPE", "case_terms_input", len(case_hpo)])
    writer.writerow([
        "PHENOTYPE",
        "case_terms_with_direct_gene_associations",
        len(case_terms_with_associations),
    ])
    writer.writerow([
        "PHENOTYPE",
        "personal_terms_with_direct_gene_associations",
        len(personal_terms_with_associations),
    ])
    writer.writerow([
        "PHENOTYPE",
        "family_terms_with_direct_gene_associations",
        len(family_terms_with_associations),
    ])

    writer.writerow([
        "GENE",
        "case_direct_gene_universe",
        len(case_gene_universe),
    ])
    writer.writerow([
        "GENE",
        "personal_history_gene_universe",
        len(personal_gene_universe),
    ])
    writer.writerow([
        "GENE",
        "family_history_gene_universe",
        len(family_gene_universe),
    ])

    writer.writerow(["CANDIDATE", "total", candidate_total])
    writer.writerow([
        "CANDIDATE",
        "unique_gene_symbols",
        len(unique_candidate_genes),
    ])
    writer.writerow([
        "CANDIDATE",
        "direct_case_hpo_gene_match",
        primary_match_rows,
    ])
    writer.writerow([
        "CANDIDATE",
        "personal_or_family_match_only",
        context_only_rows,
    ])
    writer.writerow([
        "CANDIDATE",
        "unique_direct_case_matched_genes",
        len(unique_primary_matched_genes),
    ])

    for tier in sorted(tier_total):
        writer.writerow([
            "TIER_TOTAL",
            tier,
            tier_total[tier],
        ])

        writer.writerow([
            "TIER_DIRECT_CASE_MATCH",
            tier,
            tier_primary_match[tier],
        ])

        writer.writerow([
            "TIER_CONTEXT_ONLY",
            tier,
            tier_context_only[tier],
        ])


for private_output in [
    gene_table,
    annotated_candidates,
    matched_candidates,
]:
    os.chmod(private_output, 0o600)

os.chmod(summary_file, 0o600)

print(f"public_rows_scanned={public_rows_scanned}")
print(f"candidate_rows={candidate_total}")
print(f"direct_case_matches={primary_match_rows}")
print(f"summary_exists={'YES' if summary_file.stat().st_size else 'NO'}")
print("final_status=COMPLETED")

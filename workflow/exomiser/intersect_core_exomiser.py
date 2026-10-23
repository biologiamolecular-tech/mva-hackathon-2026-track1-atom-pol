#!/usr/bin/env python3

import csv
import os
import sys
from collections import Counter, defaultdict

core_path, exomiser_path, output_path, summary_path = sys.argv[1:5]


def normalized_contig(value):
    value = value.strip()

    if value.lower().startswith("chr"):
        value = value[3:]

    if value.upper() in {"M", "MT"}:
        return "MT"

    return value


def normalized_position(value):
    try:
        return str(int(float(value.strip())))
    except (TypeError, ValueError):
        return value.strip()


def variant_key(contig, position, ref, alt):
    return (
        normalized_contig(contig),
        normalized_position(position),
        ref.strip().upper(),
        alt.strip().upper(),
    )


def integer_rank(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 10**12


def numeric_score(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return -1.0


def is_true(value):
    return value.strip().upper() in {"1", "YES", "TRUE", "Y"}


with open(core_path, newline="", encoding="utf-8-sig") as handle:
    core_reader = csv.DictReader(handle, delimiter="\t")
    core_fields = core_reader.fieldnames or []
    core_rows = list(core_reader)

with open(exomiser_path, newline="", encoding="utf-8-sig") as handle:
    exomiser_reader = csv.DictReader(handle, delimiter="\t")
    exomiser_fields = exomiser_reader.fieldnames or []
    exomiser_rows = list(exomiser_reader)

required_core = {"TIER", "CHROM", "POS", "REF", "ALT", "SYMBOL"}
required_exomiser = {
    "#RANK", "GENE_SYMBOL", "MOI", "CONTIG", "START", "REF", "ALT",
    "EXOMISER_GENE_COMBINED_SCORE", "EXOMISER_GENE_PHENO_SCORE",
    "EXOMISER_GENE_VARIANT_SCORE", "EXOMISER_VARIANT_SCORE",
    "CONTRIBUTING_VARIANT", "EXOMISER_ACMG_CLASSIFICATION"
}

missing_core = sorted(required_core.difference(core_fields))
missing_exomiser = sorted(required_exomiser.difference(exomiser_fields))

if missing_core or missing_exomiser:
    print("missing_core_columns=" + ",".join(missing_core), file=sys.stderr)
    print("missing_exomiser_columns=" + ",".join(missing_exomiser), file=sys.stderr)
    sys.exit(1)

exomiser_by_variant = defaultdict(list)

for row in exomiser_rows:
    key = variant_key(
        row["CONTIG"],
        row["START"],
        row["REF"],
        row["ALT"],
    )
    exomiser_by_variant[key].append(row)

additional_fields = [
    "EXOMISER_MATCH",
    "EXOMISER_MATCH_ROWS",
    "EXOMISER_GENE_CONCORDANT",
    "EXOMISER_BEST_RANK",
    "EXOMISER_BEST_GENE_SYMBOL",
    "EXOMISER_BEST_MOI",
    "EXOMISER_ALL_MOI",
    "EXOMISER_GENE_COMBINED_SCORE",
    "EXOMISER_GENE_PHENO_SCORE",
    "EXOMISER_GENE_VARIANT_SCORE",
    "EXOMISER_VARIANT_SCORE",
    "EXOMISER_CONTRIBUTING_VARIANT",
    "EXOMISER_ACMG_CLASSIFICATION",
]

tier_total = Counter()
tier_matched = Counter()
rank_buckets = Counter()
best_moi = Counter()
best_acmg = Counter()

matched_count = 0
gene_concordant_count = 0
multiple_match_count = 0
contributing_count = 0
output_rows = []

for core_row in core_rows:
    key = variant_key(
        core_row["CHROM"],
        core_row["POS"],
        core_row["REF"],
        core_row["ALT"],
    )

    matches = exomiser_by_variant.get(key, [])
    tier = core_row["TIER"]
    tier_total[tier] += 1

    output_row = dict(core_row)

    if not matches:
        output_row.update({
            "EXOMISER_MATCH": "NO",
            "EXOMISER_MATCH_ROWS": "0",
            "EXOMISER_GENE_CONCORDANT": "NO",
        })

        for field in additional_fields[3:]:
            output_row[field] = ""

        output_rows.append(output_row)
        continue

    matched_count += 1
    tier_matched[tier] += 1

    if len(matches) > 1:
        multiple_match_count += 1

    core_symbol = core_row["SYMBOL"].strip().upper()

    concordant_matches = [
        row for row in matches
        if row["GENE_SYMBOL"].strip().upper() == core_symbol
    ]

    gene_concordant = bool(concordant_matches)

    if gene_concordant:
        gene_concordant_count += 1
        selection_pool = concordant_matches
    else:
        selection_pool = matches

    best = min(
        selection_pool,
        key=lambda row: (
            integer_rank(row["#RANK"]),
            -numeric_score(row["EXOMISER_GENE_COMBINED_SCORE"]),
        ),
    )

    best_rank = integer_rank(best["#RANK"])

    if best_rank <= 10:
        rank_bucket = "rank_1_to_10"
    elif best_rank <= 25:
        rank_bucket = "rank_11_to_25"
    elif best_rank <= 50:
        rank_bucket = "rank_26_to_50"
    elif best_rank <= 100:
        rank_bucket = "rank_51_to_100"
    elif best_rank < 10**12:
        rank_bucket = "rank_above_100"
    else:
        rank_bucket = "rank_missing"

    rank_buckets[rank_bucket] += 1

    moi_values = sorted({
        row["MOI"].strip()
        for row in matches
        if row["MOI"].strip()
    })

    if any(is_true(row["CONTRIBUTING_VARIANT"]) for row in matches):
        contributing_count += 1

    best_moi[best["MOI"].strip() or "missing"] += 1
    best_acmg[
        best["EXOMISER_ACMG_CLASSIFICATION"].strip() or "missing"
    ] += 1

    output_row.update({
        "EXOMISER_MATCH": "YES",
        "EXOMISER_MATCH_ROWS": str(len(matches)),
        "EXOMISER_GENE_CONCORDANT": "YES" if gene_concordant else "NO",
        "EXOMISER_BEST_RANK": best["#RANK"],
        "EXOMISER_BEST_GENE_SYMBOL": best["GENE_SYMBOL"],
        "EXOMISER_BEST_MOI": best["MOI"],
        "EXOMISER_ALL_MOI": ";".join(moi_values),
        "EXOMISER_GENE_COMBINED_SCORE":
            best["EXOMISER_GENE_COMBINED_SCORE"],
        "EXOMISER_GENE_PHENO_SCORE":
            best["EXOMISER_GENE_PHENO_SCORE"],
        "EXOMISER_GENE_VARIANT_SCORE":
            best["EXOMISER_GENE_VARIANT_SCORE"],
        "EXOMISER_VARIANT_SCORE":
            best["EXOMISER_VARIANT_SCORE"],
        "EXOMISER_CONTRIBUTING_VARIANT":
            best["CONTRIBUTING_VARIANT"],
        "EXOMISER_ACMG_CLASSIFICATION":
            best["EXOMISER_ACMG_CLASSIFICATION"],
    })

    output_rows.append(output_row)

with open(output_path, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=core_fields + additional_fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(output_rows)

summary_rows = [
    ("REFERENCE", "core_candidates", len(core_rows)),
    ("REFERENCE", "exomiser_variant_rows", len(exomiser_rows)),
    ("MATCH", "exact_variant_match", matched_count),
    ("MATCH", "not_found_in_exomiser", len(core_rows) - matched_count),
    ("MATCH", "gene_symbol_concordant", gene_concordant_count),
    ("MATCH", "multiple_exomiser_rows", multiple_match_count),
    ("MATCH", "contributing_variant", contributing_count),
]

for tier in sorted(tier_total):
    summary_rows.append(("TIER_TOTAL", tier, tier_total[tier]))
    summary_rows.append(("TIER_MATCHED", tier, tier_matched[tier]))
    summary_rows.append((
        "TIER_UNMATCHED",
        tier,
        tier_total[tier] - tier_matched[tier],
    ))

for bucket in [
    "rank_1_to_10",
    "rank_11_to_25",
    "rank_26_to_50",
    "rank_51_to_100",
    "rank_above_100",
    "rank_missing",
]:
    summary_rows.append(("BEST_RANK", bucket, rank_buckets[bucket]))

for moi, count in sorted(best_moi.items()):
    summary_rows.append(("BEST_MOI", moi, count))

for acmg, count in sorted(best_acmg.items()):
    summary_rows.append(("EXOMISER_ACMG", acmg, count))

with open(summary_path, "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["section", "criterion", "count"])
    writer.writerows(summary_rows)

os.chmod(output_path, 0o600)
os.chmod(summary_path, 0o600)

print("core_candidates=" + str(len(core_rows)))
print("exomiser_variant_rows=" + str(len(exomiser_rows)))
print("exact_variant_matches=" + str(matched_count))
print("private_output_exists=YES")
print("summary_exists=YES")
print("final_status=COMPLETED")

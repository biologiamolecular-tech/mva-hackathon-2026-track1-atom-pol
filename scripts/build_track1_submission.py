#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Build the MVA Hackathon Track 1 submission."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--epcr", required=True, type=float)
    args = parser.parse_args()

    if not 0 < args.epcr <= 1:
        raise SystemExit("ERROR: EPCR must be in the interval (0,1].")

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.is_file():
        raise SystemExit(f"ERROR: input file not found: {input_path}")

    selected = {}

    with input_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "PRIORITY_ORDER",
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "GT",
        }

        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                "ERROR: missing input columns: " + ",".join(sorted(missing))
            )

        for row in reader:
            order = row["PRIORITY_ORDER"].strip()
            if order in {"1", "3"}:
                if order in selected:
                    raise SystemExit(
                        f"ERROR: duplicated PRIORITY_ORDER={order}"
                    )
                selected[order] = row

    if set(selected) != {"1", "3"}:
        raise SystemExit(
            "ERROR: could not find exactly PRIORITY_ORDER 1 and 3."
        )

    first = selected["1"]
    second = selected["3"]

    for label, row in (("P1_01", first), ("P1_03", second)):
        for field in ("CHROM", "POS", "REF", "ALT", "GT"):
            if not row[field].strip():
                raise SystemExit(
                    f"ERROR: empty {field} for {label}."
                )

    chrom_1 = first["CHROM"].strip()
    chrom_2 = second["CHROM"].strip()

    if not chrom_1.lower().startswith("chr"):
        chrom_1 = "chr" + chrom_1

    if not chrom_2.lower().startswith("chr"):
        chrom_2 = "chr" + chrom_2

    if chrom_1 != chrom_2:
        raise SystemExit(
            "ERROR: the two candidates are on different chromosomes."
        )

    fieldnames = [
        "proband_id",
        "chrom_1",
        "pos_1",
        "ref_1",
        "alt_1",
        "chrom_2",
        "pos_2",
        "ref_2",
        "alt_2",
        "epcr",
        "finding_type",
        "notes",
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(
            {
                "proband_id": "PROBAND01",
                "chrom_1": chrom_1,
                "pos_1": first["POS"].strip(),
                "ref_1": first["REF"].strip().upper(),
                "alt_1": first["ALT"].strip().upper(),
                "chrom_2": chrom_2,
                "pos_2": second["POS"].strip(),
                "ref_2": second["REF"].strip().upper(),
                "alt_2": second["ALT"].strip().upper(),
                "epcr": f"{args.epcr:.2f}",
                "finding_type": "primary",
                "notes": (
                    "Two heterozygous candidate variants in the same gene; "
                    "phenotype and recessive-model support; "
                    "read-based phase undetermined."
                ),
            }
        )

    print("submission_creation_status=0")
    print("submission_rows=1")
    print("candidate_1_priority_order=1")
    print("candidate_2_priority_order=3")
    print("same_chromosome=YES")
    print(f"epcr={args.epcr:.2f}")
    print(f"output={output_path}")


if __name__ == "__main__":
    main()

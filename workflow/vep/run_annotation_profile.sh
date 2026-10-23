#!/usr/bin/env bash

set -uo pipefail
umask 077

PROJECT_DIR="$HOME/MVA_Hackathon_2026"

ANNOTATED="$PROJECT_DIR/05_analysis/vep_full/case_primary_PASS.vep116.vcf.gz"
CUSTOM_SEVERITY="$PROJECT_DIR/05_analysis/vep_full/split_vep_severity.vep116.txt"

SUMMARY="$PROJECT_DIR/05_analysis/vep_full/annotation_profile.tsv"
ERROR_LOG="$PROJECT_DIR/logs/23_annotation_profile.stderr.log"
TIME_LOG="$PROJECT_DIR/logs/23_annotation_profile.time.log"
STATUS_LOG="$PROJECT_DIR/logs/23_annotation_profile.status.log"

{
    echo "started_utc=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "stage=validating_inputs"
} > "$STATUS_LOG"

if [ ! -s "$ANNOTATED" ]; then
    echo "final_status=FAILED_ANNOTATED_VCF_NOT_FOUND" >> "$STATUS_LOG"
    exit 1
fi

if [ ! -s "$CUSTOM_SEVERITY" ]; then
    echo "final_status=FAILED_SEVERITY_FILE_NOT_FOUND" >> "$STATUS_LOG"
    exit 1
fi

echo "stage=generating_annotation_profile" >> "$STATUS_LOG"

/usr/bin/time -v \
    -o "$TIME_LOG" \
    nice -n 10 \
    bcftools +split-vep \
        -S "$CUSTOM_SEVERITY" \
        -s worst \
        -f '%Consequence\t%IMPACT\t%MAX_AF\t%CLIN_SIG\n' \
        "$ANNOTATED" \
        2> "$ERROR_LOG" |
awk -F '\t' '
BEGIN {
    OFS = "\t"
}
{
    total++

    consequence = $1
    impact = $2
    af = $3
    clin = tolower($4)

    if (consequence == "" || consequence == "." || consequence == "-") {
        consequence_count["MISSING"]++
    } else {
        normalized = consequence
        gsub(/,/, "&", normalized)
        number_terms = split(normalized, terms, "&")

        for (i = 1; i <= number_terms; i++) {
            if (terms[i] != "")
                consequence_count[terms[i]]++
        }
    }

    if (impact == "" || impact == "." || impact == "-")
        impact_count["MISSING"]++
    else
        impact_count[impact]++

    if (af == "" || af == "." || af == "-") {
        af_count["missing"]++
    } else if (af !~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/) {
        af_count["unparseable"]++
    } else if ((af + 0) == 0) {
        af_count["zero"]++
    } else if ((af + 0) <= 0.00001) {
        af_count["0_to_0.00001"]++
    } else if ((af + 0) <= 0.0001) {
        af_count["0.00001_to_0.0001"]++
    } else if ((af + 0) < 0.001) {
        af_count["0.0001_to_0.001"]++
    } else if ((af + 0) <= 0.01) {
        af_count["0.001_to_0.01"]++
    } else {
        af_count["greater_than_0.01"]++
    }

    if (clin == "" || clin == "." || clin == "-") {
        clin_count["missing"]++
    } else if (clin ~ /conflict/) {
        clin_count["conflicting_interpretations"]++
    } else if (clin ~ /pathogenic/) {
        clin_count["pathogenic_or_likely_pathogenic"]++
    } else if (clin ~ /benign/) {
        clin_count["benign_or_likely_benign"]++
    } else if (clin ~ /uncertain/) {
        clin_count["uncertain_significance"]++
    } else if (clin ~ /risk_factor/) {
        clin_count["risk_factor"]++
    } else if (clin ~ /drug_response/) {
        clin_count["drug_response"]++
    } else {
        clin_count["other_clinical_significance"]++
    }
}
END {
    print "section", "category", "count"
    print "TOTAL", "records", total

    for (key in impact_count)
        print "IMPACT", key, impact_count[key]

    for (key in consequence_count)
        print "CONSEQUENCE_TERM", key, consequence_count[key]

    for (key in af_count)
        print "MAX_AF", key, af_count[key]

    for (key in clin_count)
        print "CLIN_SIG", key, clin_count[key]
}
' > "$SUMMARY"

pipeline_status=("${PIPESTATUS[@]}")
split_vep_status="${pipeline_status[0]:-255}"
awk_status="${pipeline_status[1]:-255}"

{
    echo "split_vep_status=$split_vep_status"
    echo "awk_status=$awk_status"
    echo "summary_exists=$([ -s "$SUMMARY" ] && echo YES || echo NO)"
    echo "summary_lines=$(wc -l < "$SUMMARY" 2>/dev/null || echo 0)"
    echo "stderr_bytes=$(wc -c < "$ERROR_LOG" 2>/dev/null || echo 0)"

    if [ -s "$SUMMARY" ]; then
        awk -F '\t' '
            $1=="TOTAL" && $2=="records" {
                print "profiled_records=" $3
            }
        ' "$SUMMARY"
    fi

    if [ "$split_vep_status" -eq 0 ] &&
       [ "$awk_status" -eq 0 ] &&
       [ -s "$SUMMARY" ]; then
        echo "final_status=COMPLETED"
    else
        echo "final_status=FAILED"
    fi

    echo "finished_utc=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
} >> "$STATUS_LOG"

if [ "$split_vep_status" -ne 0 ] || [ "$awk_status" -ne 0 ]; then
    exit 1
fi

exit 0

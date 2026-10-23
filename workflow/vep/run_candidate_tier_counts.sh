#!/usr/bin/env bash

set -uo pipefail
umask 077

PROJECT_DIR="$HOME/MVA_Hackathon_2026"

ANNOTATED="$PROJECT_DIR/05_analysis/vep_full/case_primary_PASS.vep116.vcf.gz"
SEVERITY="$PROJECT_DIR/05_analysis/vep_full/split_vep_severity.vep116.txt"

SUMMARY="$PROJECT_DIR/05_analysis/vep_full/candidate_tier_counts.tsv"
ERROR_LOG="$PROJECT_DIR/logs/24_candidate_tier_counts.stderr.log"
TIME_LOG="$PROJECT_DIR/logs/24_candidate_tier_counts.time.log"
STATUS_LOG="$PROJECT_DIR/logs/24_candidate_tier_counts.status.log"

{
    echo "started_utc=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "stage=validating_inputs"
} > "$STATUS_LOG"

if [ ! -s "$ANNOTATED" ]; then
    echo "final_status=FAILED_ANNOTATED_VCF_NOT_FOUND" >> "$STATUS_LOG"
    exit 1
fi

if [ ! -s "$SEVERITY" ]; then
    echo "final_status=FAILED_SEVERITY_FILE_NOT_FOUND" >> "$STATUS_LOG"
    exit 1
fi

echo "stage=counting_candidate_tiers" >> "$STATUS_LOG"

/usr/bin/time -v \
    -o "$TIME_LOG" \
    nice -n 10 \
    bcftools +split-vep \
        -S "$SEVERITY" \
        -s worst \
        -f '%Consequence\t%IMPACT\t%MAX_AF\t%CLIN_SIG\t%SIFT\t%PolyPhen\n' \
        "$ANNOTATED" \
        2> "$ERROR_LOG" |
awk -F '\t' '
function missing(value) {
    return value == "" || value == "." || value == "-"
}
function numeric(value) {
    return value ~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/
}
function rare_or_missing(value, threshold) {
    return missing(value) ||
           (numeric(value) && (value + 0) <= threshold)
}
BEGIN {
    OFS = "\t"
}
{
    total++

    consequence = tolower($1)
    impact      = tolower($2)
    af          = $3
    clin        = tolower($4)
    sift        = tolower($5)
    polyphen    = tolower($6)

    strict = rare_or_missing(af, 0.001)
    broad  = rare_or_missing(af, 0.01)

    plp = (clin ~ /pathogenic/ && clin !~ /conflict/ && clin !~ /benign/)

    high = (impact == "high")
    moderate = (impact == "moderate")

    damaging = (sift ~ /deleterious/ || polyphen ~ /probably_damaging/ || polyphen ~ /possibly_damaging/)

    splice_related = (consequence ~ /splice_acceptor/ || consequence ~ /splice_donor/ || consequence ~ /splice_region/ || consequence ~ /splice_polypyrimidine/)

    if (missing(af))
        af_missing++
    else if (!numeric(af))
        af_unparseable++
    else {
        if ((af + 0) <= 0.001)
            af_le_0001++
        if ((af + 0) <= 0.01)
            af_le_001++
    }

    if (strict)
        strict_or_missing++
    if (broad)
        broad_or_missing++

    if (plp) {
        clinvar_plp++
        if (strict)
            clinvar_plp_strict++
        if (broad)
            clinvar_plp_broad++
    }

    if (high) {
        high_total++
        if (strict)
            high_strict++
        if (broad)
            high_broad++
    }

    if (moderate) {
        moderate_total++

        if (strict)
            moderate_strict++
        if (broad)
            moderate_broad++

        if (damaging) {
            moderate_damaging++

            if (strict)
                moderate_damaging_strict++
            if (broad)
                moderate_damaging_broad++
        }
    }

    if (splice_related) {
        splice_total++

        if (strict)
            splice_strict++
        if (broad)
            splice_broad++
    }

    strict_candidate = (plp || (strict && high) || (strict && moderate && damaging) || (strict && splice_related))

    broad_candidate = (plp || (broad && high) || (broad && moderate && damaging) || (broad && splice_related))

    if (strict_candidate)
        candidate_union_strict++

    if (broad_candidate)
        candidate_union_broad++

    # Niveles estrictos mutuamente excluyentes
    if (plp)
        tier1++
    else if (strict && high)
        tier2++
    else if (strict && moderate && damaging)
        tier3++
    else if (strict && splice_related)
        tier4++
}
END {
    print "section", "criterion", "count"

    print "TOTAL", "records", total

    print "FREQUENCY", "MAX_AF_missing", af_missing + 0
    print "FREQUENCY", "MAX_AF_unparseable", af_unparseable + 0
    print "FREQUENCY", "MAX_AF_le_0.001", af_le_0001 + 0
    print "FREQUENCY", "MAX_AF_le_0.01", af_le_001 + 0
    print "FREQUENCY", "MAX_AF_le_0.001_or_missing", strict_or_missing + 0
    print "FREQUENCY", "MAX_AF_le_0.01_or_missing", broad_or_missing + 0

    print "CLINVAR", "pathogenic_or_likely_pathogenic", clinvar_plp + 0
    print "CLINVAR", "pathogenic_rare_or_missing_0.001", clinvar_plp_strict + 0
    print "CLINVAR", "pathogenic_rare_or_missing_0.01", clinvar_plp_broad + 0

    print "HIGH", "all", high_total + 0
    print "HIGH", "rare_or_missing_0.001", high_strict + 0
    print "HIGH", "rare_or_missing_0.01", high_broad + 0

    print "MODERATE", "all", moderate_total + 0
    print "MODERATE", "rare_or_missing_0.001", moderate_strict + 0
    print "MODERATE", "rare_or_missing_0.01", moderate_broad + 0
    print "MODERATE", "damaging_all", moderate_damaging + 0
    print "MODERATE", "damaging_rare_or_missing_0.001", moderate_damaging_strict + 0
    print "MODERATE", "damaging_rare_or_missing_0.01", moderate_damaging_broad + 0

    print "SPLICE", "all", splice_total + 0
    print "SPLICE", "rare_or_missing_0.001", splice_strict + 0
    print "SPLICE", "rare_or_missing_0.01", splice_broad + 0

    print "CANDIDATE_UNION", "strict_0.001_or_missing", candidate_union_strict + 0
    print "CANDIDATE_UNION", "broad_0.01_or_missing", candidate_union_broad + 0

    print "EXCLUSIVE_TIER", "T1_ClinVar_PLP", tier1 + 0
    print "EXCLUSIVE_TIER", "T2_HIGH_strict", tier2 + 0
    print "EXCLUSIVE_TIER", "T3_MODERATE_damaging_strict", tier3 + 0
    print "EXCLUSIVE_TIER", "T4_splice_strict", tier4 + 0
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

#!/usr/bin/env bash

set -uo pipefail
umask 077

PROJECT_DIR="$HOME/MVA_Hackathon_2026"

ANNOTATED="$PROJECT_DIR/05_analysis/vep_full/case_primary_PASS.vep116.vcf.gz"
SEVERITY="$PROJECT_DIR/05_analysis/vep_full/split_vep_severity.vep116.txt"

PRIVATE_TABLE="$PROJECT_DIR/05_analysis/vep_full/candidate_variants.private.tsv"
SUMMARY="$PROJECT_DIR/05_analysis/vep_full/candidate_extraction_summary.tsv"

ERROR_LOG="$PROJECT_DIR/logs/25_candidate_extraction.stderr.log"
TIME_LOG="$PROJECT_DIR/logs/25_candidate_extraction.time.log"
STATUS_LOG="$PROJECT_DIR/logs/25_candidate_extraction.status.log"

{
    echo "started_utc=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "stage=validating_inputs"
} > "$STATUS_LOG"

if [ ! -s "$ANNOTATED" ] || [ ! -s "$SEVERITY" ]; then
    echo "final_status=FAILED_INPUT_VALIDATION" >> "$STATUS_LOG"
    exit 1
fi

echo "stage=extracting_private_candidates" >> "$STATUS_LOG"

/usr/bin/time -v \
    -o "$TIME_LOG" \
    nice -n 10 \
    bcftools +split-vep \
        -S "$SEVERITY" \
        -s worst \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%IMPACT\t%SYMBOL\t%Gene\t%Feature\t%BIOTYPE\t%HGVSc\t%HGVSp\t%MANE_SELECT\t%SIFT\t%PolyPhen\t%MAX_AF\t%gnomADe_AF\t%gnomADg_AF\t%CLIN_SIG\t%PHENO\t%PUBMED\t%GENE_PHENO[\t%GT\t%DP\t%GQ\t%AD]\n' \
        "$ANNOTATED" \
        2> "$ERROR_LOG" |
awk -F '\t' '
function missing(value) {
    return value == "" || value == "." || value == "-"
}
function numeric(value) {
    return value ~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/
}
BEGIN {
    OFS = "\t"
    print "TIER","CHROM","POS","REF","ALT","Consequence","IMPACT","SYMBOL","Gene","Feature","BIOTYPE","HGVSc","HGVSp","MANE_SELECT","SIFT","PolyPhen","MAX_AF","gnomADe_AF","gnomADg_AF","CLIN_SIG","PHENO","PUBMED","GENE_PHENO","GT","DP","GQ","AD"
}
{
    consequence = tolower($5)
    impact = tolower($6)
    sift = tolower($14)
    polyphen = tolower($15)
    af = $16
    clin = tolower($19)

    strict = (missing(af) || (numeric(af) && (af + 0) <= 0.001))
    pathogenic = (clin ~ /pathogenic/)
    benign = (clin ~ /benign/)
    conflict = (clin ~ /conflict/)
    high = (impact == "high")
    moderate = (impact == "moderate")
    damaging = (sift ~ /deleterious/ || polyphen ~ /probably_damaging/ || polyphen ~ /possibly_damaging/)
    splice_related = (consequence ~ /splice_acceptor/ || consequence ~ /splice_donor/ || consequence ~ /splice_region/ || consequence ~ /splice_polypyrimidine/)

    tier = ""

    if (conflict)
        tier = "T1C_CLINVAR_CONFLICT"
    else if (pathogenic && benign)
        tier = "T1B_CLINVAR_MIXED"
    else if (pathogenic)
        tier = "T1A_CLINVAR_PLP"
    else if (strict && high)
        tier = "T2_HIGH_RARE"
    else if (strict && moderate && damaging)
        tier = "T3_MODERATE_DAMAGING_RARE"
    else if (strict && splice_related)
        tier = "T4_SPLICE_RARE"

    if (tier != "")
        print tier, $0
}
' > "$PRIVATE_TABLE"

pipeline_status=("${PIPESTATUS[@]}")
split_vep_status="${pipeline_status[0]:-255}"
awk_status="${pipeline_status[1]:-255}"

chmod 600 "$PRIVATE_TABLE" 2>/dev/null || true

if [ "$split_vep_status" -eq 0 ] && [ "$awk_status" -eq 0 ]; then
    awk -F '\t' '
    BEGIN {
        OFS = "\t"
        print "section","tier","count"
    }
    NR > 1 {
        tier_count[$1]++
        total++
    }
    END {
        for (tier in tier_count)
            print "CANDIDATE",tier,tier_count[tier]

        print "TOTAL","private_candidate_rows",total + 0
    }
    ' "$PRIVATE_TABLE" > "$SUMMARY"
fi

candidate_rows=$(awk 'END {print (NR > 0 ? NR - 1 : 0)}' "$PRIVATE_TABLE")

{
    echo "split_vep_status=$split_vep_status"
    echo "awk_status=$awk_status"
    echo "private_table_exists=$([ -s "$PRIVATE_TABLE" ] && echo YES || echo NO)"
    echo "private_candidate_rows=$candidate_rows"
    echo "summary_exists=$([ -s "$SUMMARY" ] && echo YES || echo NO)"
    echo "stderr_bytes=$(wc -c < "$ERROR_LOG" 2>/dev/null || echo 0)"

    if [ "$split_vep_status" -eq 0 ] &&
       [ "$awk_status" -eq 0 ] &&
       [ "$candidate_rows" -gt 0 ]; then
        echo "final_status=COMPLETED"
    else
        echo "final_status=FAILED"
    fi

    echo "finished_utc=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
} >> "$STATUS_LOG"

exit 0

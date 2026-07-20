#!/usr/bin/env bash
###############################################################################
# PIPELINE 2 — VARIANTS DE QUASI-ESPÈCE (bcftools mpileup)
# (à lancer uniquement APRES validation manuelle des résultats du pipeline 1)
#
# BAM triés/indexés (pipeline 1)
#   -> bcftools mpileup (AD, DP par position)
#   -> bcftools call -mv --ploidy 1   (liste des ALT observés + AD/DP)
#   -> bcftools norm -m -both         (1 ligne = 1 variant biallélique)
#   -> extract_variant_table.py       (calcul AF = AD_alt/DP, classification
#                                       majoritaire AF>=0.5 / minoritaire 0.01<=AF<0.5)
#   -> table finale par échantillon + table fusionnée
#
# Utilisation :
#   chmod +x pipeline_2.sh
#   ./pipeline_2.sh -r Consensus2.fasta -b resultat/05_alignements -o resultat_variants
# chmod +x pipeline2v2.sh extractionvariant.py
#./pipeline2v2.sh -r reference.fasta -b resultat/05_alignements -o resultat_variants
###############################################################################

set -euo pipefail
shopt -s nullglob

########################## PARAMÈTRES ########################################

REF=""
BAM_DIR=""
OUT_DIR="resultats_pipeline"
THREADS=4

# --- paramètres bcftools mpileup (sensibilité de détection) ---
MAX_DEPTH=1000000      # profondeur max considérée par mpileup (défaut bcftools = 250, trop bas !)
MIN_BASEQ=10           # qualité de base minimale (-Q)
MIN_MAPQ=20            # qualité de mapping minimale (-q)

# --- seuils de classification quasi-espèce ---
AF_MIN=0.01            # en dessous : considéré comme bruit, écarté
AF_MAJ=0.5             # au-dessus (ou égal) : variant majoritaire
MIN_ALT_DEPTH=5        # nb minimal de reads supportant l'ALT (filtre anti-bruit ONT)

usage() {
    echo "Usage: $0 -r reference.fasta -b dossier_bam_tries [-o dossier_sortie] [-t threads]"
    echo "          [-d max_depth] [-Q min_baseq] [-q min_mapq]"
    echo "          [-m af_min] [-M af_maj] [-a min_alt_depth]"
    exit 1
}

while getopts "r:b:o:t:d:Q:q:m:M:a:" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        b) BAM_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        d) MAX_DEPTH="$OPTARG" ;;
        Q) MIN_BASEQ="$OPTARG" ;;
        q) MIN_MAPQ="$OPTARG" ;;
        m) AF_MIN="$OPTARG" ;;
        M) AF_MAJ="$OPTARG" ;;
        a) MIN_ALT_DEPTH="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "$REF" || -z "$BAM_DIR" ]] && usage
[[ ! -f "$REF" ]] && { echo "Référence introuvable : $REF"; exit 1; }
[[ ! -d "$BAM_DIR" ]] && { echo "Dossier BAM introuvable : $BAM_DIR"; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXTRACT_SCRIPT="$SCRIPT_DIR/extractionvariant.py"
[[ ! -f "$EXTRACT_SCRIPT" ]] && { echo "extractionvariant.py introuvable à côté de ce script"; exit 1; }

REF=$(realpath "$REF")
BAM_DIR=$(realpath "$BAM_DIR")
mkdir -p "$OUT_DIR"
OUT_DIR=$(realpath "$OUT_DIR")

# la référence doit être indexée pour mpileup
[[ ! -f "${REF}.fai" ]] && samtools faidx "$REF"

########################## ARBORESCENCE #######################################

VCF_BRUT="$OUT_DIR/08_vcf_brut"
VCF_BIALLELIC="$OUT_DIR/09_vcf_biallelique"
TABLES="$OUT_DIR/10_tables_variants"
STATS="$OUT_DIR/04_stats"
FINAL="$OUT_DIR/07_resultats_finaux"
LOGS="$OUT_DIR/logs"

mkdir -p "$VCF_BRUT" "$VCF_BIALLELIC" "$TABLES" "$STATS" "$FINAL" "$LOGS"

echo "======================================================================"
echo " PIPELINE 2 — VARIANTS DE QUASI-ESPÈCE (bcftools mpileup)"
echo "======================================================================"
echo ">>> Référence         : $REF"
echo ">>> BAM (pipeline 1)  : $BAM_DIR"
echo ">>> Sortie            : $OUT_DIR"
echo ">>> mpileup           : max-depth=$MAX_DEPTH, min-baseQ=$MIN_BASEQ, min-mapQ=$MIN_MAPQ"
echo ">>> Classification    : majoritaire AF>=$AF_MAJ | minoritaire ${AF_MIN}<=AF<$AF_MAJ | bruit AF<$AF_MIN"
echo ">>> Profondeur ALT min: $MIN_ALT_DEPTH reads"
echo

BAMS=("$BAM_DIR"/*.sorted.bam)
if [[ ${#BAMS[@]} -eq 0 ]]; then
    echo "Aucun fichier *.sorted.bam trouvé dans $BAM_DIR"
    exit 1
fi

########################## 1. mpileup + call (par échantillon) #################

echo "=== [1/3] bcftools mpileup + call (AD/DP par position) ==="
for bam in "${BAMS[@]}"; do
    sample=$(basename "$bam" .sorted.bam)
    echo "  - $sample"

    bcftools mpileup \
        -f "$REF" \
        -a FORMAT/AD,FORMAT/DP \
        -d "$MAX_DEPTH" \
        -Q "$MIN_BASEQ" \
        -q "$MIN_MAPQ" \
        -Ou "$bam" 2> "$LOGS/mpileup_${sample}.log" \
    | bcftools call \
        -mv \
        --ploidy 1 \
        -Ou 2>> "$LOGS/mpileup_${sample}.log" \
    | bcftools norm \
        -f "$REF" \
        -m -both \
        -Oz -o "$VCF_BIALLELIC/${sample}.norm.vcf.gz" \
        2>> "$LOGS/mpileup_${sample}.log"

    bcftools index -t "$VCF_BIALLELIC/${sample}.norm.vcf.gz"

    # conserve aussi le VCF (post-split) pour audit dans un dossier séparé
    cp "$VCF_BIALLELIC/${sample}.norm.vcf.gz" "$VCF_BRUT/${sample}.vcf.gz" 2>/dev/null || true
done

########################## 2. Extraction table + classification ################

echo "=== [2/3] Extraction AD/DP -> AF -> classification majoritaire/minoritaire ==="
for vcf in "$VCF_BIALLELIC"/*.norm.vcf.gz; do
    sample=$(basename "$vcf" .norm.vcf.gz)

    python3 "$EXTRACT_SCRIPT" \
        --vcf "$vcf" \
        --sample "$sample" \
        --af-min "$AF_MIN" \
        --af-maj "$AF_MAJ" \
        --min-alt-depth "$MIN_ALT_DEPTH" \
        -o "$TABLES/${sample}_variants.tsv" \
        2>> "$LOGS/extract_table_${sample}.log"

    # stats bcftools classiques (Ts/Tv, indels/SNP...) pour info QC
    bcftools stats "$vcf" > "$STATS/bcftools_stats_${sample}.txt"
done

########################## 3. Table fusionnée tous échantillons ################

echo "=== [3/3] Fusion des tables (tous échantillons) ==="
MERGED="$TABLES/all_samples_variants.tsv"
first=1
for tsv in "$TABLES"/*_variants.tsv; do
    if [[ $first -eq 1 ]]; then
        cat "$tsv" > "$MERGED"
        first=0
    else
        tail -n +2 "$tsv" >> "$MERGED"
    fi
done

cp "$MERGED" "$FINAL/" 2>/dev/null || true
cp "$TABLES"/*_variants.tsv "$FINAL/" 2>/dev/null || true

echo
echo "======================================================================"
echo " PIPELINE 2 TERMINÉ"
echo "======================================================================"
echo ">>> Table par échantillon : $FINAL/<sample>_variants.tsv"
echo ">>> Table fusionnée       : $FINAL/all_samples_variants.tsv"
echo ">>> VCF biallélique       : $VCF_BIALLELIC/<sample>.norm.vcf.gz"
echo ">>> Stats bcftools (QC)   : $STATS/bcftools_stats_*.txt"
echo ">>> Logs                  : $LOGS"
echo
echo "Colonnes de la table : sample  chrom  position  REF  ALT  DP  AD_ref  AD_alt  AF  classification"
#!/usr/bin/env bash
###############################################################################
# PIPELINE 1 — NETTOYAGE
# Reads bruts -> NanoPlot avant -> Filtrage (NanoFilt) -> NanoPlot après
#            -> Alignement (minimap2) -> Consensus (medaka_consensus)
#
# A la fin de ce pipeline : VALIDER MANUELLEMENT les résultats avant
# de lancer le pipeline 2 (variants) :
#   - comparer 01_nanoplot_avant vs 03_nanoplot_apres (qualité, longueur)
#   - vérifier 04_stats/flagstat_report.txt (taux de reads mappées correct ?)
#   - vérifier 04_stats/stats_avant_filtre.tsv vs stats_apres_filtre.tsv
#
# Utilisation :
# Avant de lancer faire: chmod +x ./pipeline1.sh
#   ./pipeline_1_nettoyage.sh -r reference.fasta -i /chemin/fastq_bruts -o /chemin/sortie
# ./pipeline1.sh -r reference.fasta -i fastq -o resultat
###############################################################################

set -euo pipefail
shopt -s nullglob

########################## PARAMÈTRES ########################################

REF=""
INPUT_DIR=""
OUT_DIR="resultats_pipeline"
THREADS=4
MIN_LEN=350
MAX_LEN=800
MIN_QUAL=15
MEDAKA_MODEL="r1041_e82_400bps_hac_v5.2.0"   # à adapter à votre kit/flowcell/basecaller

usage() {
    echo "Usage: $0 -r reference.fasta -i dossier_fastq_bruts [-o dossier_sortie] [-t threads] [-m modele_medaka] [-l min_len] [-L max_len] [-q min_qual]"
    exit 1
}

while getopts "r:i:o:t:m:l:L:q:" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MEDAKA_MODEL="$OPTARG" ;;
        l) MIN_LEN="$OPTARG" ;;
        L) MAX_LEN="$OPTARG" ;;
        q) MIN_QUAL="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "$REF" || -z "$INPUT_DIR" ]] && usage
[[ ! -f "$REF" ]] && { echo "Référence introuvable : $REF"; exit 1; }
[[ ! -d "$INPUT_DIR" ]] && { echo "Dossier d'entrée introuvable : $INPUT_DIR"; exit 1; }

REF=$(realpath "$REF")
INPUT_DIR=$(realpath "$INPUT_DIR")
mkdir -p "$OUT_DIR"
OUT_DIR=$(realpath "$OUT_DIR")

########################## ARBORESCENCE #######################################

QC_AVANT="$OUT_DIR/01_nanoplot_avant"
FILTRES="$OUT_DIR/02_fastq_filtres"
QC_APRES="$OUT_DIR/03_nanoplot_apres"
STATS="$OUT_DIR/04_stats"
ALIGN="$OUT_DIR/05_alignements"
CONSENSUS="$OUT_DIR/06_consensus"
FINAL="$OUT_DIR/07_resultats_finaux"
LOGS="$OUT_DIR/logs"

mkdir -p "$QC_AVANT" "$FILTRES" "$QC_APRES" "$STATS" "$ALIGN" "$CONSENSUS" "$FINAL" "$LOGS"

echo "======================================================================"
echo " PIPELINE 1 — NETTOYAGE"
echo "======================================================================"
echo ">>> Référence      : $REF"
echo ">>> Fastq bruts     : $INPUT_DIR"
echo ">>> Sortie          : $OUT_DIR"
echo ">>> Threads         : $THREADS"
echo ">>> Modèle medaka   : $MEDAKA_MODEL"
echo ">>> Filtre longueur : $MIN_LEN - $MAX_LEN pb, qualité min $MIN_QUAL"
echo

cd "$INPUT_DIR" 

echo "Fichiers FASTQ détectés :"
ls *.fastq*
echo "Nombre de FASTQ : $(ls *.fastq* | wc -l)"



########################## 1. NANOPLOT AVANT ##################################

echo "=== [1/7] NanoPlot AVANT nettoyage ==="
for f in *.fastq*; do
    [[ -e "$f" ]] || continue
    echo "Traitement de : $f"

    sample=$(basename "$f")
    sample=${sample%.fastq.gz}
    sample=${sample%.fastq}
    NanoPlot --fastq "$f" -o "$QC_AVANT/${sample}" -t "$THREADS" \
        > "$LOGS/nanoplot_avant_${sample}.log" 2>&1
done

########################## 2. NETTOYAGE (qualité + longueur) ##################

echo "=== [2/7] Nettoyage NanoFilt (q>=$MIN_QUAL, ${MIN_LEN}-${MAX_LEN} pb) ==="
for f in *.fastq*; do
    [[ -e "$f" ]] || continue
    sample=$(basename "$f")
    sample=${sample%.fastq.gz}
    sample=${sample%.fastq}

    if [[ "$f" == *.gz ]]; then
        zcat "$f" | NanoFilt -q "$MIN_QUAL" -l "$MIN_LEN" --maxlength "$MAX_LEN" \
            > "$FILTRES/${sample}_filtered.fastq" 2> "$LOGS/nanofilt_${sample}.log"
    else
        NanoFilt -q "$MIN_QUAL" -l "$MIN_LEN" --maxlength "$MAX_LEN" < "$f" \
            > "$FILTRES/${sample}_filtered.fastq" 2> "$LOGS/nanofilt_${sample}.log"
    fi
done

echo "FASTQ filtrés générés :"
ls "$FILTRES"
########################## 3. NANOPLOT APRES ###################################

echo "=== [3/7] NanoPlot APRES nettoyage ==="
cd "$FILTRES"
for f in *_filtered.fastq; do
    [[ -e "$f" ]] || continue
    sample=${f%_filtered.fastq}
    NanoPlot --fastq "$f" -o "$QC_APRES/${sample}" -t "$THREADS" \
        > "$LOGS/nanoplot_apres_${sample}.log" 2>&1
done

########################## 4. STATISTIQUES + ALIGNEMENT ########################

echo "=== [4/7] Statistiques seqkit avant/après ==="
seqkit stats -a "$INPUT_DIR"/*.fastq* -T > "$STATS/stats_avant_filtre.tsv"
seqkit stats -a "$FILTRES"/*_filtered.fastq -T > "$STATS/stats_apres_filtre.tsv"

echo "=== [4/7] Alignement minimap2 + tri + index ==="
for f in *_filtered.fastq; do
    sample=${f%_filtered.fastq}
    minimap2 -ax map-ont -t "$THREADS" "$REF" "$f" 2> "$LOGS/minimap2_${sample}.log" \
        | samtools sort -@ "$THREADS" -o "$ALIGN/${sample}.sorted.bam" -
    samtools index "$ALIGN/${sample}.sorted.bam"
done

echo "=== [4/7] Vérification alignement (flagstat) ==="
{
for bam in "$ALIGN"/*.sorted.bam; do
    echo "==== $(basename "$bam") ===="
    samtools flagstat "$bam"
    echo
done
} > "$STATS/flagstat_report.txt"


echo "=== [5/7] Vérification alignement (profondeur) ==="
{
for bam in "$ALIGN"/*.sorted.bam; do
    echo "==== $(basename "$bam") ===="
    samtools depth -a "$bam" | gawk '
{
    total++;
    depth=$3;
    sum+=depth;

    if(depth>0){
        covered++;
        a[covered]=depth;
    }

    if(depth>max) max=depth;
    if(total==1 || depth<min) min=depth;

    if(depth>=10) cov10++;
    if(depth>=20) cov20++;
}
END{
    if(total==0){
        print "Aucune position";
        exit
    }

    asort(a);

    if(covered%2==1)
        median=a[int((covered+1)/2)];
    else if(covered>0)
        median=(a[covered/2]+a[covered/2+1])/2;
    else
        median=0;

    print "Mean_depth:", sum/total;
    print "Median_depth:", median;
    print "Min_depth:", min;
    print "Max_depth:", max;
    print "Positions >0x:", covered;
    print "Coverage %:", (covered/total)*100;
    print "Positions >=10x:", cov10;
    print "Positions >=20x:", cov20;
}'
    echo
done
} > "$STATS/couverture.tsv"


    
    

########################## 6. CONSENSUS #########################################

echo "=== [6/7] Consensus medaka_consensus par échantillon ==="
for f in *_filtered.fastq; do
    sample=${f%_filtered.fastq}
    medaka_consensus \
        -i "$f" \
        -d "$REF" \
        -o "$CONSENSUS/${sample}" \
        -m "$MEDAKA_MODEL" \
        -t "$THREADS" \
        > "$LOGS/medaka_consensus_${sample}.log" 2>&1

    if [[ -f "$CONSENSUS/${sample}/consensus.fasta" ]]; then
        seqkit replace -p "^(.*)\$" -r "${sample}_\$1" \
            "$CONSENSUS/${sample}/consensus.fasta" \
            > "$FINAL/${sample}_consensus.fasta"
    fi
done

########################## 6. REGROUPEMENT ######################################

echo "=== [7/7] Regroupement des consensus ==="
cat "$FINAL"/*_consensus.fasta > "$FINAL/all_samples_consensus.fasta" 2>/dev/null || true

echo
echo "======================================================================"
echo " PIPELINE 1 TERMINÉ — À VALIDER AVANT DE LANCER LE PIPELINE 2"
echo "======================================================================"
echo ">>> Comparer avant/après : $QC_AVANT  vs  $QC_APRES"
echo ">>> Stats longueur/qualité: $STATS/stats_avant_filtre.tsv / stats_apres_filtre.tsv"
echo ">>> Taux de reads mappées : $STATS/flagstat_report.txt"
echo ">>> Consensus par échantillon : $FINAL/*_consensus.fasta"
echo ">>> BAM triés/indexés (nécessaires au pipeline 2) : $ALIGN"
echo
echo "Une fois validé, lancez :"
echo "  ./pipeline_2_variants.sh -r $REF -b $ALIGN -o $OUT_DIR"
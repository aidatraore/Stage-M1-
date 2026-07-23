#!/usr/bin/env bash
###############################################################################
# PIPELINE 2 — CLIQUESNV (analyse de quasi-espèces, région NS5B, Nanopore)
#
# Prend en entrée les BAM triés/indexés produits par pipeline_1_nettoyage.sh
# (dossier 05_alignements) + la même référence NS5B, et lance CliqueSNV
# en mode "snv-pacbio" (mode recommandé pour reads longs bruités type Nanopore,
# CliqueSNV n'ayant pas de mode nanopore dédié).
#
# Étapes :
#   1) Vérifications (java, samtools, jar CliqueSNV, fichiers d'entrée)
#   2) BAM -> SAM filtré (retrait unmapped/secondary/supplementary)
#   3) CliqueSNV par échantillon (mode snv-pacbio)
#   4) Synthèse : tableau récapitulatif (nb haplotypes, fréquences) + fasta combiné
#
# A LA FIN : VALIDER MANUELLEMENT avant interprétation biologique :
#   - vérifier 04_stats/reads_utilises.tsv (combien de reads sont entrés dans CliqueSNV)
#   - vérifier 02_cliquesnv/<sample>/*.json champ "error" (doit être "none")
#   - regarder 03_synthese/recap_haplotypes.tsv
#
# Utilisation :
#   chmod +x ./pipelinecliquesnv.sh
#   ./pipelinecliquesnv.sh -r reference.fasta -b resultat/05_alignements -j clique-snv.jar -o resultatcliqueSNV
###############################################################################

set -euo pipefail
shopt -s nullglob

########################## PARAMÈTRES ########################################

REF=""
BAM_DIR=""
JAR=""
OUT_DIR="resultats_pipeline"
THREADS=4
XMX="8G"
CSNV_T=10          # -t   : seuil O22 (support statistique du lien entre SNV)
CSNV_TF="0.01,0.02,0.05"  # -f  : liste de seuils tf à tester, séparés par des virgules  
CSNV_EL=700        # -el  : limite d'arêtes dans le graphe de SNV

usage() {
    echo "Usage: $0 -r reference.fasta -b dossier_bam -j clique-snv.jar [-o dossier_sortie] [-p threads] [-x xmx_java] [-t seuil_t] [-f seuil_tf] [-e el]"
    exit 1
}

while getopts "r:b:j:o:p:x:t:f:e:" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        b) BAM_DIR="$OPTARG" ;;
        j) JAR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        p) THREADS="$OPTARG" ;;
        x) XMX="$OPTARG" ;;
        t) CSNV_T="$OPTARG" ;;
        f) CSNV_TF="$OPTARG" ;;
        e) CSNV_EL="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "$REF" || -z "$BAM_DIR" || -z "$JAR" ]] && usage
[[ ! -f "$REF" ]] && { echo "Référence introuvable : $REF"; exit 1; }
[[ ! -d "$BAM_DIR" ]] && { echo "Dossier BAM introuvable : $BAM_DIR"; exit 1; }
[[ ! -f "$JAR" ]] && { echo "Jar CliqueSNV introuvable : $JAR"; exit 1; }

command -v java >/dev/null 2>&1 || { echo "java introuvable dans le PATH"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools introuvable dans le PATH"; exit 1; }

REF=$(realpath "$REF")
BAM_DIR=$(realpath "$BAM_DIR")
JAR=$(realpath "$JAR")
mkdir -p "$OUT_DIR"
OUT_DIR=$(realpath "$OUT_DIR")

########################## ARBORESCENCE #######################################

SAM_DIR="$OUT_DIR/01_sam_filtres"
CSNV_DIR="$OUT_DIR/02_cliquesnv"
STATS="$OUT_DIR/04_stats"
SYNTHESE="$OUT_DIR/03_synthese"
LOGS="$OUT_DIR/logs"

mkdir -p "$SAM_DIR" "$CSNV_DIR" "$STATS" "$SYNTHESE" "$LOGS"

echo "======================================================================"
echo " PIPELINE 2 — CLIQUESNV (NS5B, mode snv-pacbio pour reads Nanopore)"
echo "======================================================================"
echo ">>> Référence       : $REF"
echo ">>> BAM d'entrée    : $BAM_DIR"
echo ">>> Jar CliqueSNV   : $JAR"
echo ">>> Sortie          : $OUT_DIR"
echo ">>> Threads         : $THREADS"
echo ">>> Java -Xmx       : $XMX"
echo ">>> -t / -tf / -el  : $CSNV_T / $CSNV_TF / $CSNV_EL"
echo

BAMS=("$BAM_DIR"/*.sorted.bam)
[[ ${#BAMS[@]} -eq 0 ]] && { echo "Aucun .sorted.bam trouvé dans $BAM_DIR"; exit 1; }

echo "Échantillons détectés :"
for b in "${BAMS[@]}"; do echo "  - $(basename "$b")"; done
echo

########################## 1. BAM -> SAM FILTRÉ ################################

echo "=== [1/3] Conversion BAM -> SAM (retrait unmapped/secondary/supplementary) ==="
> "$STATS/reads_utilises.tsv"
echo -e "sample\treads_total_bam\treads_conserves_sam" >> "$STATS/reads_utilises.tsv"

for bam in "${BAMS[@]}"; do
    sample=$(basename "$bam")
    sample=${sample%.sorted.bam}

    total=$(samtools view -c "$bam")
    # -F 2308 = exclut unmapped (4) + secondary (256) + supplementary (2048)
    samtools view -h -F 2308 "$bam" > "$SAM_DIR/${sample}.sam" \
        2> "$LOGS/samtools_view_${sample}.log"
    conserves=$(samtools view -c "$SAM_DIR/${sample}.sam" 2>/dev/null || echo 0)

    echo -e "${sample}\t${total}\t${conserves}" >> "$STATS/reads_utilises.tsv"
done

echo "Résumé reads conservés :"
column -t "$STATS/reads_utilises.tsv"
echo

########################## 2. CLIQUESNV PAR ÉCHANTILLON ########################

echo "=== [2/3] CliqueSNV (mode snv-pacbio) par échantillon ==="
for sam in "$SAM_DIR"/*.sam; do
    sample=$(basename "$sam" .sam)
    sample_out="$CSNV_DIR/${sample}"
    mkdir -p "$sample_out"

    echo "  -> $sample"
    java -Xmx"$XMX" -jar "$JAR" \
        -m snv-pacbio \
        -in "$sam" \
        -outDir "$sample_out" \
        -t "$CSNV_T" \
        -tf "$CSNV_TF" \
        -el "$CSNV_EL" \
        -threads "$THREADS" \
        -fdf extended3 \
        -rn -log \
        > "$LOGS/cliquesnv_${sample}.log" 2>&1 || \
        echo "     ATTENTION: CliqueSNV a retourné une erreur pour $sample, voir $LOGS/cliquesnv_${sample}.log"
done
echo

########################## 3. SYNTHÈSE #########################################

echo "=== [3/3] Synthèse des résultats ==="

RECAP="$SYNTHESE/recap_haplotypes.tsv"
echo -e "sample\tstatut\tnb_haplotypes\tfreq_max\tfreq_min" > "$RECAP"

: > "$SYNTHESE/all_haplotypes_combined.fasta"

for sample_out in "$CSNV_DIR"/*/; do
    sample=$(basename "$sample_out")
    json_file=$(ls "$sample_out"/*.json 2>/dev/null | head -n1 || true)
    fasta_file=$(ls "$sample_out"/*.fasta 2>/dev/null | head -n1 || true)

    if [[ -z "$json_file" ]]; then
        echo -e "${sample}\tPAS_DE_SORTIE\tNA\tNA\tNA" >> "$RECAP"
        continue
    fi

    # Extraction via python3 (plus fiable que grep/sed sur du JSON)
    python3 - "$json_file" "$sample" >> "$RECAP" <<'PYEOF'
import json, sys
json_path, sample = sys.argv[1], sys.argv[2]
with open(json_path) as f:
    d = json.load(f)
status = d.get("error", "inconnu")
haps = d.get("haplotypes", [])
freqs = [h.get("frequency", 0.0) for h in haps]
nb = len(haps)
fmax = f"{max(freqs):.4f}" if freqs else "NA"
fmin = f"{min(freqs):.4f}" if freqs else "NA"
print(f"{sample}\t{status}\t{nb}\t{fmax}\t{fmin}")
PYEOF

    # Ajout des haplotypes de cet échantillon au fasta combiné, préfixés par le nom d'échantillon
    if [[ -n "$fasta_file" ]]; then
        sed "s/^>/>${sample}_/" "$fasta_file" >> "$SYNTHESE/all_haplotypes_combined.fasta"
    fi
done

echo "Récapitulatif :"
column -t "$RECAP"
echo

echo "======================================================================"
echo " PIPELINE 2 TERMINÉ"
echo "======================================================================"
echo ">>> SAM filtrés               : $SAM_DIR"
echo ">>> Sorties CliqueSNV brutes  : $CSNV_DIR/<sample>/*.json + *.fasta"
echo ">>> Reads conservés par échant.: $STATS/reads_utilises.tsv"
echo ">>> Récapitulatif haplotypes  : $RECAP"
echo ">>> Fasta combiné (multi-éch.): $SYNTHESE/all_haplotypes_combined.fasta"
echo
echo "Prochaine étape possible : alignement multiple + arbre phylogénétique"
echo "sur $SYNTHESE/all_haplotypes_combined.fasta pour comparer les quasi-espèces"
echo "entre échantillons/points de temps."
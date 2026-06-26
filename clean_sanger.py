#####################################################################################
########## Nettoyage du chromatogramme Sanger + Diversité virale du  VHC ############
#####################################################################################
#####################################################################################
#
# PIPELINE : prétraitement d'un chromatogramme Sanger pour analyse intra-familiale du VHC
#   fichier .ab1
#       │
#       ▼
#   trimmed_seq   : séquence brute trimmée (Q20fichier)
#       │
#       ├──→ primary_seq   : base dominante (A/C/G/T UNIQUEMENT — jamais N, jamais IUPAC)
#       │                    min_height=0 → on tranche toujours sur le signal max
#       │                    identique à sangerseqR primarySeq() dans le logiciel R
#       │
#       ├──→ secondary_seq : 2ème signal le plus fort
#       │                    si N dans trimmed et égalité → code IUPAC
#       │                    si ratio < 25% → = primary (pas de variant)
#       │
#       └──→ consensus_seq : pool primary + secondary → code IUPAC
#
#####################################################################################

from Bio import SeqIO
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys, os

#####################################################################################
# 1) ARGUMENTS
#####################################################################################
input_ab1  = sys.argv[1]
prefix     = sys.argv[2]
plotprefix = sys.argv[3]

############################################################
# 2) Lecture du chromatogramme
############################################################
record = SeqIO.read(input_ab1, "abi")

seq  = str(record.seq)
qual = record.letter_annotations["phred_quality"]
abif = record.annotations["abif_raw"]

# Ordre réel des canaux via FWO_1
fwo = abif.get("FWO_1", b"ACGT").decode("utf-8").strip()
if set(fwo) != {"A", "C", "G", "T"}:
    print(f"[WARN] FWO_1 inattendu : '{fwo}' → fallback sur ACGT")
    fwo = "ACGT"
print(f"[INFO] Ordre des canaux (FWO_1) : {fwo}")

data_channels = [abif["DATA9"], abif["DATA10"], abif["DATA11"], abif["DATA12"]]
channel_map   = {base: data_channels[i] for i, base in enumerate(fwo)}
A = channel_map["A"]
C = channel_map["C"]
G = channel_map["G"]
T = channel_map["T"]

# Positions réelles des pics
if "PLOC2" in abif:
    peak_locations = list(abif["PLOC2"])
    print("[INFO] Utilisation de PLOC2")
elif "PLOC1" in abif:
    peak_locations = list(abif["PLOC1"])
    print("[INFO] Fallback sur PLOC1")
else:
    raise ValueError("Aucun champ PLOC trouvé dans le fichier .ab1")

max_channel_len = min(len(A), len(C), len(G), len(T))

############################################################
# 3) Trimming Q20 
############################################################
def find_trim_bounds(qual, threshold=20):
    
    start = 0
    while start < len(qual) and qual[start] < threshold:
        start += 1

    end = len(qual) - 1
    while end > 0 and qual[end] < threshold:
        end -= 1

    return start, end

start, end = find_trim_bounds(qual)

end = min(end, len(peak_locations) - 1)

if start >= end:
    raise ValueError(
        f"Qualité trop faible : aucune région Q20 trouvée. "
        f"start={start}, end={end}. Vérifier le fichier .ab1."
    )

print(f"[INFO] Trim : {start} → {end}  ({end - start + 1} bp)")

############################################################
# 4) Seuil min_height — UNIQUEMENT pour la secondary
# -------------------------------------------------------
# Pour la PRIMARY : min_height = 0
#   → on tranche toujours sur le signal max
#   → jamais de N dus au seuil
#   → identique à sangerseqR primarySeq()
#
# Pour la SECONDARY : seuil dynamique (10e percentile)
#   → on filtre le bruit pour ne pas appeler de faux variants
############################################################
all_peak_max = [
    max(A[peak_locations[i]], C[peak_locations[i]],
        G[peak_locations[i]], T[peak_locations[i]])
    for i in range(start, end + 1)
]
min_height_secondary = int(np.percentile(all_peak_max, 10))
min_height_primary   = 0   # jamais de N dans la primaire

print(f"[INFO] Seuil min_height primary   : {min_height_primary} (toujours une base)")
print(f"[INFO] Seuil min_height secondary : {min_height_secondary} (filtre le bruit)")

############################################################
# 5) Tables IUPAC
############################################################
IUPAC_ENCODE = {
    frozenset(["A"]): "A", frozenset(["C"]): "C",
    frozenset(["G"]): "G", frozenset(["T"]): "T",
    frozenset(["A", "G"]): "R", frozenset(["C", "T"]): "Y",
    frozenset(["G", "C"]): "S", frozenset(["A", "T"]): "W",
    frozenset(["G", "T"]): "K", frozenset(["A", "C"]): "M",
    frozenset(["A", "C", "G"]): "V", frozenset(["A", "C", "T"]): "H",
    frozenset(["A", "G", "T"]): "D", frozenset(["C", "G", "T"]): "B",
    frozenset(["A", "C", "G", "T"]): "N",
}

IUPAC_DECODE = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"},
    "W": {"A", "T"}, "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"}
}

############################################################
# 6) Fonctions
############################################################
def get_primary(a, c, g, t, min_height=0):
    """
    Base PRIMAIRE : signal le plus fort.
    → Toujours A/C/G/T (jamais N, jamais IUPAC)
    → min_height=0 par défaut : on tranche toujours
    → identique à sangerseqR primarySeq()
    """
    signals = {"A": a, "C": c, "G": g, "T": t}
    primary_base, primary_val = max(signals.items(), key=lambda x: x[1])
    if primary_val < min_height:
        return "N"
    return primary_base


def get_secondary(a, c, g, t, min_height=100, minor_ratio=0.25):
    """
    Base SECONDAIRE : 2ème signal le plus fort.
    → Toujours A/C/G/T ou None (jamais IUPAC)
    → None si signal trop faible ou ratio < minor_ratio
    """
    signals = {"A": a, "C": c, "G": g, "T": t}
    sorted_bases = sorted(signals.items(), key=lambda x: x[1], reverse=True)
    primary_val    = sorted_bases[0][1]
    secondary_base = sorted_bases[1][0]
    secondary_val  = sorted_bases[1][1]
    if primary_val < min_height:
        return None
    if secondary_val / primary_val < minor_ratio:
        return None
    return secondary_base


def correct_N_secondary(a, c, g, t, min_height=100, thr=0.20):
    """
    Correction N dans secondary : code IUPAC si plusieurs bases ex-aequo.
    → Codes IUPAC possibles ici (R, Y, M, etc.)
    """
    signals = {"A": a, "C": c, "G": g, "T": t}
    maxsig  = max(signals.values())
    if maxsig < min_height:
        return "N"
    dom = frozenset(b for b, v in signals.items() if v / maxsig >= thr)
    return IUPAC_ENCODE.get(dom, "N")


def pool_consensus(primary, secondary):
    """
    Pool primary + secondary → code IUPAC consensus.
    - primary   : A/C/G/T (jamais N, jamais IUPAC)
    - secondary : A/C/G/T ou code IUPAC ou None
    - primary si secondary == None (pas de variant)
    - IUPAC union sinon
    """
    if secondary is None:
        return primary
    bases_p = IUPAC_DECODE.get(primary, {primary})
    bases_s = IUPAC_DECODE.get(secondary, {secondary})
    bases_union = frozenset(bases_p | bases_s)
    return IUPAC_ENCODE.get(bases_union, "N")

############################################################
# 7) BOUCLE UNIFIÉE
############################################################
trimmed_seq   = ""
primary_seq   = ""
secondary_seq = ""
consensus_seq = ""
qe_ratios_all = []

for idx in range(start, end + 1):

    raw_base = seq[idx]
    peak     = peak_locations[idx]

    # Pic hors bornes → on prend quand même la base brute pour primary
    if peak >= max_channel_len:
        trimmed_seq   += raw_base
        # Pour primary : on garde raw_base si c'est une base valide
        # sinon on ne peut pas faire mieux → N uniquement dans ce cas extrême
        primary_seq   += raw_base if raw_base in "ACGT" else "N"
        secondary_seq += primary_seq[-1]   # secondary = primary (pas de signal)
        consensus_seq += primary_seq[-1]
        qe_ratios_all.append(0.0)
        continue

    a = A[peak]
    c = C[peak]
    g = G[peak]
    t = T[peak]

    # trimmed_seq : base brute
    trimmed_seq += raw_base

    # -----------------------------------------------------------
    # PRIMARY : min_height=0 → toujours A/C/G/T, jamais N
    # Si base brute valide → on garde
    # Si base brute == N   → on corrige avec le signal max (min_height=0)
    # -----------------------------------------------------------
    if raw_base in "ACGT":
        primary = raw_base
    else:
        primary = get_primary(a, c, g, t, min_height=0)

    primary_seq += primary

    # -----------------------------------------------------------
    # SECONDARY : seuil dynamique pour filtrer le bruit
    # Si base brute valide → 2ème signal (si ratio >= 25%)
    # Si base brute == N   → code IUPAC sur les signaux
    # Si pas de vrai variant → = primary
    # -----------------------------------------------------------
    if raw_base in "ACGT":
        sec = get_secondary(a, c, g, t,
                            min_height=min_height_secondary,
                            minor_ratio=0.25)
        secondary = sec if sec is not None else primary
    else:
        secondary = correct_N_secondary(a, c, g, t,
                                        min_height=min_height_secondary,
                                        thr=0.20)

    secondary_seq += secondary

    # CONSENSUS : pool primary + secondary
    consensus_seq += pool_consensus(primary, secondary)

    # Ratio pour le plot
    vals  = sorted([a, c, g, t], reverse=True)
    ratio = vals[1] / vals[0] if vals[0] > 0 else 0.0
    qe_ratios_all.append(ratio)

############################################################
# 8) Vérifications
############################################################
# primary ne doit contenir que A/C/G/T (plus de N !)
bases_autorisees_primary   = set("ACGTN")   # N toléré uniquement si pic hors bornes
bases_autorisees_secondary = set("ACGTNRYSWKMBDHV")

assert set(primary_seq).issubset(bases_autorisees_primary), \
    f"ERREUR inattendu dans primary_seq : {set(primary_seq) - bases_autorisees_primary}"

assert set(secondary_seq).issubset(bases_autorisees_secondary), \
    f"ERREUR : caractères inattendus dans secondary_seq"

assert len(primary_seq) == len(secondary_seq) == len(consensus_seq) == len(trimmed_seq), \
    f"DÉSALIGNEMENT : primary={len(primary_seq)}, secondary={len(secondary_seq)}, " \
    f"consensus={len(consensus_seq)}, trimmed={len(trimmed_seq)}"

trim_len   = len(trimmed_seq)
nb_N_primary = primary_seq.count("N")

print(f"[INFO] Longueur finale : {trim_len} bp")
print(f"[INFO] N résiduels dans primary : {nb_N_primary} "
      f"(uniquement si pic hors bornes)")

############################################################
# 9) Variants quasi-espèces
############################################################
valid_bases    = {"A", "C", "G", "T"}
diff_positions = []
variant_details = []

for i in range(trim_len):
    p = primary_seq[i]
    s = secondary_seq[i]
    if p not in valid_bases or p == s:
        diff_positions.append(0)
    else:
        diff_positions.append(1)
        variant_details.append((i, p, s, qe_ratios_all[i]))

nb_valid = sum(1 for b in primary_seq if b in valid_bases)
diffs    = sum(diff_positions) 
ratio_qe = diffs / nb_valid * 100 if nb_valid > 0 else 0.0

print(f"[INFO] Quasi-espèces : {diffs} / {nb_valid} bp ({ratio_qe:.2f}%)")

############################################################
# 10) Sauvegarde FASTA
############################################################
sample_name = os.path.basename(prefix)

with open(prefix + "_trimmed.fasta", "w") as f:
    f.write(f">{sample_name}_trimmed\n{trimmed_seq}\n")

with open(prefix + "_primary.fasta", "w") as f:
    f.write(f">{sample_name}_primary\n{primary_seq}\n")

with open(prefix + "_secondary.fasta", "w") as f:
    f.write(f">{sample_name}_secondary\n{secondary_seq}\n")

with open(prefix + "_consensus.fasta", "w") as f:
    f.write(f">{sample_name}_consensus\n{consensus_seq}\n")

print(f"[INFO] FASTA sauvegardés (longueur uniforme : {trim_len} bp)")
print(f"[INFO]   trimmed   = séquence brute trimmée")
print(f"[INFO]   primary   = base dominante (A/C/G/T — plus de N)")
print(f"[INFO]   secondary = 2ème signal + IUPAC si N")
print(f"[INFO]   consensus = pool primary + secondary")

############################################################
# 11) Plots
############################################################

# 11a) Profil Qualité
plt.figure(figsize=(10, 3))
plt.plot(qual, color="#2c7bb6", linewidth=0.8)
plt.axvline(x=start, color="green",  linestyle="--", label=f"Début trim (pos {start})")
plt.axvline(x=end,   color="orange", linestyle="--", label=f"Fin trim (pos {end})")
plt.axhline(y=20, color="red",    linestyle=":", linewidth=0.8, label="Q20")
plt.axhline(y=30, color="purple", linestyle=":", linewidth=0.8, label="Q30")
plt.title("Profil Qualité Phred")
plt.xlabel("Position brute")
plt.ylabel("Qualité Phred")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(plotprefix + "_qualite.png", dpi=150)
plt.close()

# 11b) Chromatogramme trimmé
chrom_start = peak_locations[start]
chrom_end   = peak_locations[end]

plt.figure(figsize=(14, 4))
plt.plot(A[chrom_start:chrom_end], color="green", label="A", linewidth=0.7)
plt.plot(C[chrom_start:chrom_end], color="blue",  label="C", linewidth=0.7)
plt.plot(G[chrom_start:chrom_end], color="black", label="G", linewidth=0.7)
plt.plot(T[chrom_start:chrom_end], color="red",   label="T", linewidth=0.7)
plt.legend()
plt.title("Chromatogramme trimmé (positions réelles des pics)")
plt.xlabel("Position chromatographique")
plt.ylabel("Intensité fluorescence")
plt.tight_layout()
plt.savefig(plotprefix + "_chrom_trim.png", dpi=150)
plt.close()

# 11c) Plot quasi-espèces
fig, axes = plt.subplots(2, 1, figsize=(14, 6), sharex=True)

axes[0].bar(range(trim_len), qe_ratios_all,
            color="#f4a261", alpha=0.7, width=1.0, label="Ratio signal secondaire")
axes[0].axhline(y=0.25, color="red", linestyle="--", linewidth=1.2,
                label="Seuil quasi-espèce (25%)")
axes[0].set_ylabel("Ratio secondaire / primaire")
axes[0].set_title("Signal secondaire — Quasi-Espèces VHC")
axes[0].set_ylim(0, 1)
axes[0].legend(fontsize=9)

axes[1].stem(range(trim_len), diff_positions,
             linefmt="r-", markerfmt="ro", basefmt="k-")
for pos, p, s, r in variant_details:
    axes[1].annotate(f"{p}→{s}\n{r:.0%}",
                     xy=(pos, 1), fontsize=6,
                     ha="center", va="bottom", color="darkred")
axes[1].set_ylabel("Variant détecté")
axes[1].set_xlabel("Position sur la séquence trimmée (bp)")
axes[1].set_title("Variants quasi-espèces VHC (primary ≠ secondary)")

plt.tight_layout()
plt.savefig(plotprefix + "_quasiespeces.png", dpi=150)
plt.close()

############################################################
# 12) Statistiques
############################################################
raw_len    = len(seq)
raw_mean   = np.mean(qual)
raw_median = np.median(qual)
raw_q20    = np.sum(np.array(qual) >= 20) / raw_len * 100
raw_q30    = np.sum(np.array(qual) >= 30) / raw_len * 100

trim_frac   = trim_len / raw_len * 100
trim_mean   = np.mean(qual[start:end + 1])
trim_median = np.median(qual[start:end + 1])
trim_q20    = np.sum(np.array(qual[start:end + 1]) >= 20) / trim_len * 100
trim_q30    = np.sum(np.array(qual[start:end + 1]) >= 30) / trim_len * 100

nb_N_primary   = primary_seq.count("N")
nb_N_secondary = secondary_seq.count("N")
nb_N_consensus = consensus_seq.count("N")
pct_N          = nb_N_primary / trim_len * 100

iupac_codes = set("RYSWKMBDHV")
nb_iupac_secondary = sum(1 for b in secondary_seq if b in iupac_codes)
nb_iupac_consensus = sum(1 for b in consensus_seq if b in iupac_codes)

with open(prefix + ".stats.txt", "w") as f:
    f.write("=== STATISTIQUES BRUTES ===\n")
    f.write(f"Longueur brute          : {raw_len}\n")
    f.write(f"Qualité moyenne brute   : {raw_mean:.2f}\n")
    f.write(f"Qualité médiane brute   : {raw_median:.2f}\n")
    f.write(f"% >= Q20                : {raw_q20:.2f}%\n")
    f.write(f"% >= Q30                : {raw_q30:.2f}%\n\n")

    f.write("=== STATISTIQUES APRÈS TRIMMING ===\n")
    f.write(f"Longueur trimmée        : {trim_len}\n")
    f.write(f"Conservé                : {trim_frac:.2f}%\n")
    f.write(f"Qualité moyenne trimmée : {trim_mean:.2f}\n")
    f.write(f"Qualité médiane trimmée : {trim_median:.2f}\n")
    f.write(f"% >= Q20                : {trim_q20:.2f}%\n")
    f.write(f"% >= Q30                : {trim_q30:.2f}%\n\n")

    f.write("=== AMBIGUÏTÉS ===\n")
    f.write(f"N dans primary_seq      : {nb_N_primary} ({pct_N:.2f}%) "
            f"— uniquement si pic hors bornes\n")
    f.write(f"N dans secondary_seq    : {nb_N_secondary}\n")
    f.write(f"N dans consensus_seq    : {nb_N_consensus}\n")
    f.write(f"Codes IUPAC secondary   : {nb_iupac_secondary}\n")
    f.write(f"Codes IUPAC consensus   : {nb_iupac_consensus}\n\n")

    f.write("=== PARAMÈTRES UTILISÉS ===\n")
    f.write(f"Ordre canaux (FWO_1)        : {fwo}\n")
    f.write(f"Seuil trimming              : Q20 \n")
    f.write(f"min_height primary          : 0 (toujours une base)\n")
    f.write(f"min_height secondary        : {min_height_secondary}\n")
    f.write(f"Seuil minor_ratio           : 0.25 (25%)\n")
    f.write(f"Seuil IUPAC thr             : 0.20 (20%)\n")

############################################################
# 13) Rapport quasi-espèces
############################################################
with open(prefix + ".report.txt", "w") as f:
    f.write("===== RAPPORT QUASI-ESPECES VHC =====\n\n")
    f.write(f"Fichier analysé         : {input_ab1}\n")
    f.write(f"Région                  : {sample_name}\n")
    f.write(f"Trim                    : position {start} → {end}\n")
    f.write(f"Longueur analysée       : {nb_valid} bp (bases valides)\n\n")
    f.write(f"Variants détectés       : {diffs}\n")
    f.write(f"% divergence globale    : {ratio_qe:.2f}%\n\n")

    f.write("--- INTERPRETATION ---\n")
    if ratio_qe < 1:
        f.write("Population virale HOMOGENE.\n")
        f.write("Pas de quasi-espèce détectable au Sanger.\n")
        f.write("→ Si suspicion clinique de résistance, envisager NGS.\n")
    elif ratio_qe < 5:
        f.write("FAIBLE diversité quasi-espèce.\n")
        f.write("Proche de la limite de détection Sanger (~20-25%).\n")
        f.write("→ Confirmer par séquençage profond (NGS).\n")
    elif ratio_qe < 15:
        f.write("QUASI-ESPECE MODEREE détectée.\n")
        f.write("Population hétérogène confirmée au Sanger.\n")
        f.write("→ Pertinent pour résistance NS5B/NS3/NS5A.\n")
    else:
        f.write("QUASI-ESPECE MAJEURE détectée.\n")
        f.write("Forte hétérogénéité virale — population mixte probable.\n")
        f.write("→ Analyse de résistance urgente. NGS recommandé.\n")

    f.write("\n--- ESTIMATION QUASI-ESPECES ---\n")
    f.write(f"Positions variants détectées    : {diffs}\n\n")

    if diffs == 0:
        f.write("Population virale HOMOGENE.\n")
        f.write("Nombre de quasi-espèces estimé : 1\n")
        f.write("(aucun double pic détecté au seuil de 25%)\n")
    else:
        # Nombre minimal = 2 (majoritaire + au moins 1 minoritaire)
        f.write(f"Nombre MINIMAL de quasi-espèces : 2\n")
        f.write(f"(séquence majoritaire + ≥1 séquence minoritaire)\n\n")

        # Calcul des ratios aux positions variants uniquement
        ratios_variants  = [r for _, _, _, r in variant_details]
        diversite_mean   = np.mean(ratios_variants)
        diversite_max    = np.max(ratios_variants)
        diversite_min    = np.min(ratios_variants)

        f.write(f"Ratio moyen aux positions variants : {diversite_mean:.1%}\n")
        f.write(f"Ratio min                          : {diversite_min:.1%}\n")
        f.write(f"Ratio max                          : {diversite_max:.1%}\n\n")

        # Interprétation selon les ratios
        f.write("Interprétation des ratios :\n")
        f.write("  Ratio 25-35% → 1 quasi-espèce minoritaire faible\n")
        f.write("  Ratio 35-45% → population minoritaire modérée\n")
        f.write("  Ratio 45-55% → 2 quasi-espèces en proportion égale\n")
        f.write("  Ratios très variables → populations multiples probables\n\n")

        if diversite_mean < 0.35:
            f.write("→ Population minoritaire FAIBLE détectée au Sanger\n")
        elif diversite_mean < 0.45:
            f.write("→ Population minoritaire MODEREE détectée\n")
        else:
            f.write("→ Mélange de populations en proportions comparables\n")

        f.write("\nNOTE : Le Sanger détecte les POSITIONS variants\n")
        f.write("mais ne peut pas reconstruire les haplotypes complets.\n")
        f.write("→ Pour le nombre exact de quasi-espèces : NGS recommandé.\n")
    

    f.write("\n--- VARIANTS DÉTECTÉS ---\n")
    if variant_details:
        f.write(f"{'Position':>10}  {'Primaire':>8}  {'Secondaire':>10}  {'Ratio':>8}\n")
        f.write("-" * 45 + "\n")
        for pos, p, s, r in variant_details:
            f.write(f"{pos + 1:>10}  {p:>8}  {s:>10}  {r:>7.1%}\n")
    else:
        f.write("Aucun variant détecté au seuil de 25%.\n")

############################################################
# 14) Fichier pivot Snakemake
############################################################
with open(prefix + ".done", "w") as f:
    f.write("OK\n")

print(f"[OK] Analyse terminée.")
print(f"     Variants quasi-espèces : {diffs} ({ratio_qe:.2f}%)")
print(f"     Longueur uniforme      : {trim_len} bp")
print(f"     N résiduels primary    : {nb_N_primary} (pic hors bornes uniquement)")
print(f"     Fichiers générés sous  : {prefix}.*")
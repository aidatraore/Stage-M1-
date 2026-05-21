import glob
import os

############################################################
# 1) Détection des échantillons
############################################################
def get_samples(folder):
    return [
        os.path.splitext(os.path.basename(f))[0]
        for f in glob.glob(f"{folder}/*.ab1")
    ]

SAMPLES_FOYER    = get_samples("Sequence_foyer")
SAMPLES_CONTROLE = get_samples("Sequence_controle")

if not SAMPLES_FOYER and not SAMPLES_CONTROLE:
    raise ValueError("Aucun fichier .ab1 trouvé")

print(f"[INFO] Foyer    : {len(SAMPLES_FOYER)} échantillons")
print(f"[INFO] Controle : {len(SAMPLES_CONTROLE)} échantillons")

############################################################
# 2) Règle finale
############################################################
rule all:
    input:
        expand("results_foyer/{sample}.done",    sample=SAMPLES_FOYER),
        expand("results_controle/{sample}.done", sample=SAMPLES_CONTROLE),
        "results_foyer/all_primary.fasta",
        "results_foyer/all_secondary.fasta",
        "results_foyer/all_consensus.fasta",
        "results_controle/all_primary.fasta",
        "results_controle/all_secondary.fasta",
        "results_controle/all_consensus.fasta"  

############################################################
# 3) Règle FOYER
############################################################
rule clean_sanger_foyer:
    input:
        ab1 = "Sequence_foyer/{sample}.ab1"
    output:
        done      = "results_foyer/{sample}.done",
        clean     = "results_foyer/{sample}_primary.fasta",
        secondary = "results_foyer/{sample}_secondary.fasta",
        trimmed   = "results_foyer/{sample}_trimmed.fasta",
        consensus = "results_foyer/{sample}_consensus.fasta",
        stats     = "results_foyer/{sample}.stats.txt",
        report    = "results_foyer/{sample}.report.txt",
        plot_qual = "results_foyer/{sample}_qualite.png",
        plot_chrom= "results_foyer/{sample}_chrom_trim.png",
        plot_qe   = "results_foyer/{sample}_quasiespeces.png",
    log:
        "logs/foyer/{sample}.log"
    params:
        prefix = "results_foyer/{sample}"
    shell:
        """
        mkdir -p results_foyer logs/foyer
        python scripts/clean_sanger.py \
            {input.ab1} \
            {params.prefix} \
            {params.prefix} \
            > {log} 2>&1
        """

############################################################
# 4) Règle CONTROLE
############################################################
rule clean_sanger_controle:
    input:
        ab1 = "Sequence_controle/{sample}.ab1"
    output:
        done      = "results_controle/{sample}.done",
        clean     = "results_controle/{sample}_primary.fasta",
        secondary = "results_controle/{sample}_secondary.fasta",
        trimmed   = "results_controle/{sample}_trimmed.fasta",
        consensus = "results_controle/{sample}_consensus.fasta",
        stats     = "results_controle/{sample}.stats.txt",
        report    = "results_controle/{sample}.report.txt",
        plot_qual = "results_controle/{sample}_qualite.png",
        plot_chrom= "results_controle/{sample}_chrom_trim.png",
        plot_qe   = "results_controle/{sample}_quasiespeces.png",
    log:
        "logs/controle/{sample}.log"
    params:
        prefix = "results_controle/{sample}"
    shell:
        """
        mkdir -p results_controle logs/controle
        python scripts/clean_sanger.py \
            {input.ab1} \
            {params.prefix} \
            {params.prefix} \
            > {log} 2>&1
        """

############################################################
# 5) Règle de concaténation foyer
############################################################
rule concatenate_foyer:
    input:
        prim_foyer    = expand("results_foyer/{sample}_primary.fasta",sample=SAMPLES_FOYER),
        sec_foyer      = expand("results_foyer/{sample}_secondary.fasta",sample=SAMPLES_FOYER),
        cons_foyer     = expand("results_foyer/{sample}_consensus.fasta",sample=SAMPLES_FOYER),
       
    output:
        foyer_primary    = "results_foyer/all_primary.fasta",
        foyer_secondary = "results_foyer/all_secondary.fasta",
        foyer_consensus = "results_foyer/all_consensus.fasta",
        
    shell:
        """
        cat {input.prim_foyer}    > {output.foyer_primary}
        cat {input.sec_foyer}      > {output.foyer_secondary}
        cat {input.cons_foyer}     > {output.foyer_consensus}
        """

############################################################
# 5) Règle de concaténation controle
############################################################
rule concatenate_controle:
    input:
        prim_controle = expand("results_controle/{sample}_primary.fasta", sample=SAMPLES_CONTROLE),
        sec_controle   = expand("results_controle/{sample}_secondary.fasta", sample=SAMPLES_CONTROLE),
        cons_controle  = expand("results_controle/{sample}_consensus.fasta", sample=SAMPLES_CONTROLE),
    output:
        controle_primary     = "results_controle/all_primary.fasta",
        controle_secondary = "results_controle/all_secondary.fasta",
        controle_consensus = "results_controle/all_consensus.fasta",
    shell:
        """
        cat {input.prim_controle} > {output.controle_primary}
        cat {input.sec_controle}   > {output.controle_secondary}
        cat {input.cons_controle}  > {output.controle_consensus}
        """

###################### Les Codes bash #############

####### trier le fichier ##############

for f in *.bam; do
    samtools sort "$f" -o "${f%.bam}.sorted.bam"
done

############ pour indexer le fichier ###########
 
 for f in *.sorted.bam; do
    samtools index "$f" -o "${f%.bam}.sorted.bam"
done



########## extraire le fichier fasta dans fastq #######

for f in *.fastq; do
    seqtk seq -A "$f" > "${f%.fastq}.fasta"
done

seqkit stats fichier.fastq.gz

## creer conda ########

conda create -n medaka_env 

################# activer l'environnement ###########
conda activate medaka_env

####### INSTALLER MEDAKA ##############y

conda install -c conda-forge -c bioconda medaka

####### verification de la version

medaka --version


medaka_consensus \
-i sample.fastq \
-d reference.fasta \
-o medaka_sample \
-t 4 \
-m r1041_e82_400bps_sup_v5.0.0


##### Pour automatiser ########

for f in *.fastq; do
    sample=${f%.fastq}
    
    medaka_consensus \
    -i "$f" \
    -d reference.fasta \
    -o "medaka_${sample}" \
    -t 4
    
done



########### code filtrage par longueur pour plusieur ############
for f in barcode*.fastq; do
    seqkit seq -m 350 -M 400 "$f" > "${f%.fastq}_350_400.fastq"
done


########### stat avant et apres le filtre ############

for f in barcode*.fastq; do
    seqkit stats -a "$f" -T
done > stats_avant_filtre.tsv


for f in *_350_400.fastq; do
    seqkit stats -a "$f" -T
done > stats_apres_filtre.tsv


########### qualité ############


for f in barcode*.fastq; do
    sample=${f%.fastq}
    NanoPlot --fastq "$f" -o "NanoPlot_${sample}"
done



############ alignement avec minimap 2 

minimap2 -ax map-ont Consensus2.fasta barcode01_350_400.fastq > barcode01.sam

for f in barcode*.fastq; do
    sample=${f%.fastq}
    minimap2 -ax map-ont Consensus2.fasta "$f" > "${sample}.sam"
done

########## Puis convertir en BAM  et trier:##############

samtools view -bS barcode01.sam > barcode01.bam

samtools sort barcode01.bam -o barcode01.sorted.bam

## plusieurs
for f in *.sam; do
    sample=${f%.sam}
    
    samtools view -bS "$f" > "${sample}.bam"
    samtools sort "${sample}.bam" -o "${sample}.sorted.bam"
    samtools index "${sample}.sorted.bam"
done

### Vérification de l'alignement pour tous les échantillons
for f in *.sorted.bam; do
    echo "==== $f ===="
    samtools flagstat "$f" | grep mapped
done

############# indexer

samtools index barcode01.sorted.bam


####### Générer la séquence consensus avec Medaka
medaka_consensus \
-i barcode01_350_400.fastq \
-d Consensus2.fasta \
-o medaka_barcode01 \
-m r1041_e82_400bps_hac_v5.2.0 \
-t 4


for f in barcode*.sorted.bam; do

    sample=${f%.sorted.bam}

    echo "Indexation de $f"
    samtools index "$f"

    echo "Consensus Medaka de ${sample}"
    medaka_consensus \
    -i "${sample}_350_400.fastq" \
    -d Consensus2.fasta \
    -o "medaka_${sample}" \
    -m r1041_e82_400bps_hac_v5.2.0 \
    -t 4

done


####### fastQ en fasta

for f in *_350_400.fastq; do
    seqkit fq2fa "$f" > "${f%.fastq}.fasta"
done



##### renommer automatiquement toute les reads 

seqkit replace -p ".*" -r "read_{nr}" barcode01.fasta > barcode01_rename.fasta


for f in barcode*_350_400*.fasta; do
    sample=${f%.fasta}
    seqkit replace -p ".*" -r "read_{nr}" "$f" > "${sample}_rename.fasta"
done


######## combiner plusieurs filtre pour le nettoyage


cat barcode01.fastq | NanoFilt -q 15 -l 300 > barcode01_filtre.fastq


cat barcode01.fastq | NanoFilt -q 15 -l 300 --maxlength 500 > barcode01_filtre.fastq

for file in *.fastq
do
    NanoFilt -q 15 -l 300 --maxlength 500 < "$file" > "${file%.fastq}_filtered.fastq"
done


for file in *.fastq
do
    NanoFilt -q 15 -l 300 --maxlength 500 < "$file" > "filtered_fastq/${file%.fastq}_filtered.fastq"
done

###-q 15 : élimine les lectures de faible qualité.
-l 300 : supprime les lectures trop courtes.
--maxlength 500 : élimine les lectures anormalement longues, souvent dues à des concatémères ou à des artéfacts.

#### Regrouper les consensus
cat medaka_*/consensus.fasta > all_samples_consensus.fasta


######## profondeur 

cd resultat/05_alignements

for bam in *.sorted.bam
do
    sample=${bam%.sorted.bam}
    coverage=$(samtools depth "$bam" | awk '{sum+=$3} END {if(NR>0) print sum/NR; else print 0}')
    echo -e "${sample}\t${coverage}"
done


######## dezipé
for f in *.vcf.gz; do
    gunzip -c "$f" > "${f%.gz}"
done


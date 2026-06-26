######################################################################
# Phylogénie + Distances génétiques intra/inter ménage VHC: stage M1
######################################################################
### chargement des librairies necessaire ######

library(ape)
library(BiocManager)
library(DECIPHER)
library(seqinr)
library(phangorn)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(reshape2)
library(dplyr)
library(Biostrings)
library(phytools)
library(vegan)
library(RColorBrewer)


####### definir le Repertoire de travaille #################
setwd("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/")

#####################################################################
#  STATISTIQUE DESCRIPTIVE DES DEUX GROUPES (FOYER ET CONTROLE)#####
####################################################################

########### importation de la base avec les geno 2 ###########################
foyer <- read_excel("~/MASTER BIOINFORMATIQUE/stageM1/Projets/base_foyermultilpeall.xls")
head(foyer)
################## DESCRIPTION DE FOYER #################
summary(foyer$age)
################### GROUPE AGE ###############
foyer <- foyer %>%
  mutate(
    groupe_age = case_when(
      age < 20  ~ "0-19",
      TRUE     ~ "20+"
    )
  )

table(foyer$groupe_age)
prop.table(table(foyer$groupe_age)) * 100
############ SEXE#####################
table(foyer$sexe)
prop.table(table(foyer$sexe)) * 100

############ relation chef de menage #####################
table(foyer$Relation_chef)
prop.table(table(foyer$Relation_chef)) * 100

#####################################################################
##############  Diversité virale dans le FOYER ######################
####################################################################

########### importation de la base avec les sequence consensus ###########################
consensus_all <- readDNAStringSet("~/MASTER BIOINFORMATIQUE/stageM1/Projets/results_foyer/all_consensus.fasta")
consensus_all
names(consensus_all)<-gsub("_consensus", "",names(consensus_all))
names(consensus_all)<-gsub("-NS5B", "_NS5B",names(consensus_all))
names(consensus_all)

### Filtrer les séquences au genotype 2 #############
base_idsF <- foyer$ID
base_idsF
consensus <- consensus_all[names(consensus_all) %in% base_idsF]
length(consensus_all)
length(consensus)

################# analyse des sites ambigues ####################
consensusF <- read.dna("consensusgeno2.fasta",
                       format = "fasta",
                       as.character = TRUE)
consensusF
class(consensusF)
str(consensusF)


####### definir les code ambigu ##########################
ambigu <- c("r","y","s","w","k","m","b","d","h","v","n")

####### calcul du nombre  ambigus  ##########################
nb_ambigu <- sapply(consensusF,
                    function(x) sum(x %in% ambigu))
nb_ambigu
####### longueur de chaque sequence  ##########################
longueur <- sapply(consensusF, length)

####### pourcentage des sites ambigu  ##########################
pct_ambigu <- nb_ambigu / longueur * 100

####### resumé des valeurs de la diversité intra hote  ##########################
diversite <- data.frame(
  Sequence = names(consensus),
  Sites_ambigus = nb_ambigu,
  longueur=longueur,
  Pourcentage = round(pct_ambigu, 2)
)
diversite$Sequence<-gsub('_NS5B','',diversite$Sequence)
diversite

######################### extraire les infos dont j'ai besoin ###############
foyer1 <- foyer %>%
  select(ID_TERRAIN, menages)
######################### merger  ###############
diversite <- merge(
  diversite,
  foyer1[, c("ID_TERRAIN", "menages")],
  by.x = "Sequence",
  by.y = "ID_TERRAIN",
  all.x = TRUE
)
diversite$Sequence1<-paste0(diversite$menages,"_", diversite$Sequence)

diversite
########################### graphes ################## 

couleurs_menage <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(unique(diversite$menages))),
  sort(unique(diversite$menage))
)


p <- ggplot(diversite,
            aes(x = reorder(Sequence1, Sites_ambigus),
                y = Sites_ambigus,
                fill = menages)) +
  geom_col() +
  geom_text(aes(label = Sites_ambigus),
            hjust = -0.2,
            size = 3) +
  coord_flip() +
  scale_fill_manual(values = couleurs_menage) +
  labs(
    x = "Individus",
    y = "Nombre de sites ambigus",
    fill = "Foyer"
  ) +
  theme_classic()

p

ggsave(
  filename = "~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/resultat_grapheutile/sites_ambigusFoyer.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)
###################### pour le Pourcentage de la diversité ###############
summary(diversite$Pourcentage)
############resume diversité par menage ##########
resume_diversitmenage <- diversite %>%
  group_by(menages) %>%
  summarise(
    n = n(),
    min = min(Pourcentage, na.rm = TRUE),
    Q1 = quantile(Pourcentage, 0.25, na.rm = TRUE),
    mediane = median(Pourcentage, na.rm = TRUE),
    moyenne = mean(Pourcentage, na.rm = TRUE),
    Q3 = quantile(Pourcentage, 0.75, na.rm = TRUE),
    max = max(Pourcentage, na.rm = TRUE),
    ecart_type = sd(Pourcentage, na.rm = TRUE)
  ) %>%
  arrange(mediane)

resume_diversitmenage
############ test de normalité 
qqnorm(diversite$Pourcentage,
       main = "QQ-plot des distances TN93 intra-foyer")

qqline(diversite$Pourcentage,
       col = "red",
       lwd = 2)
shapiro.test(diversite$Pourcentage)
######### arrondir les decimal à 2 chiffre #########
resume_diversitmenage <- resume_diversitmenage %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))



DFM<-ggplot(diversite,
       aes(x = reorder(menages, Pourcentage,median),
           y = Pourcentage,
           fill = menages)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0) +
  scale_fill_manual(values = couleurs_menage) +
  labs(
    x = "Ménage",
    y = "Diversité virale (%)",
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
DFM

ggsave(
  filename = "~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/resultat_grapheutile/diversitviralFoyer.png",
  plot = DFM,
  width = 8,
  height = 6,
  dpi = 300
)


############################################################
# 2) Matrice de distances génétiques de la diversité foyer
# modèle TN93 = recommandé pour VHC
############################################################
# Charger l'alignement
aln_DVF<- read.dna("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/aligned_sequence_diversitefoyer.fasta", format = "fasta")

dist_DVF <- dist.dna(aln_DVF, model = "TN93",
                     pairwise.deletion = TRUE)

# Conversion en data frame long
dist_DVF_df <- as.data.frame(as.table(as.matrix(dist_DVF)))

# Renommer les colonnes
colnames(dist_DVF_df) <- c("Sequence1", "Sequence2", "Distance")

# Supprimer la diagonale
dist_DVF_df <- subset(dist_DVF_df, Sequence1 != Sequence2)

# Garder une seule fois chaque paire
dist_DVF_df <- dist_DVF_df[as.character(dist_DVF_df$Sequence1) <
                             as.character(dist_DVF_df$Sequence2), ]

############ les stats ##############
summary(dist_DVF_df$Distance)
mean(dist_DVF_df$Distance)
sd(dist_DVF_df$Distance)

# Conversion en data frame long
dist_DVF <- as.matrix(dist_DVF)

# Construire toutes les paires
pairesDF <- expand.grid(ind1 = foyer1$ID_TERRAIN,
                       ind2 = foyer1$ID_TERRAIN,
                       stringsAsFactors = FALSE) %>%
  filter(ind1 < ind2) %>%   # éviter doublons et diagonale
  left_join(foyer1, by = c("ind1" = "ID_TERRAIN")) %>%
  rename(menage1 = menages) %>%
  left_join(foyer1, by = c("ind2" = "ID_TERRAIN")) %>%
  rename(menage2 = menages) %>%
  mutate(
    type     = ifelse(menage1 == menage2, "intra", "inter"),
    distance = mapply(function(i, j) dist_DVF[i, j], ind1, ind2)
  )
pairesDF
# Résumé par ménage
resume_DVmenage <- pairesDF %>%
  filter(type == "intra") %>%
  mutate(
    menages = menage1
  )
############## stat distance genetique de la diversité virale intra menage ###
summary(resume_DVmenage$distance)

# Résumé par ménage
resume_DVmenageinter <- pairesDF %>%
  filter(type == "inter") %>%
  mutate(
    menages = menage1
  )
resume_DVmenageinter
############## stat distance genetique de la diversité virale inter menage ###
summary(resume_DVmenageinter$distance)


#################################################################
################ importation de la base controle ################
#################################################################
########### importation de la base avec les geno 2 ###########################
controle <- read_excel("~/MASTER BIOINFORMATIQUE/stageM1/Projets/controle_all.xls")
head(controle)
####################################################################
############### stat descriptives controle ########################
##################################################################
summary(controle$age)
controle <- controle %>%
  mutate(
    groupe_age = case_when(
      age < 20  ~ "0-19",
      TRUE     ~ "20+"
    )
  )

table(controle$groupe_age)
prop.table(table(controle$groupe_age)) * 100

table(controle$sexe)
prop.table(table(controle$sexe)) * 100

table(controle$Relation_chef)
prop.table(table(controle$Relation_chef)) * 100

###############################################################################
########################### Diversité virale dans le groupe controle ##########
#############################################################################
########### importation de la base avec les sequence consensus ###########################
consensus_allC <- readDNAStringSet("~/MASTER BIOINFORMATIQUE/stageM1/Projets/results_controle//all_consensus.fasta")
consensus_allC
names(consensus_allC)<-gsub("_consensus", "",names(consensus_allC))
names(consensus_allC)<-gsub("-NS5B", "_NS5B",names(consensus_allC))
names(consensus_allC)

### Filtrer les séquences au genotype 2 #############
controle$ID<-gsub("_primary","", controle$ID)
controle$ID<-gsub("-NS5B", "_NS5B",controle$ID)
base_ids <- controle$ID
base_ids
consensusC <- consensus_allC[names(consensus_allC) %in% base_ids]
length(consensus_allC)
length(consensusC)

################# analyse des sites ambigues chez les controles ####################
consensusCnt <- read.dna("consensusControlegeno2.fasta",
                         format = "fasta",
                         as.character = TRUE)
consensusCnt

####### calcul du nombre  ambigus  ##########################
nb_ambiguC <- sapply(consensusCnt,
                     function(x) sum(x %in% ambigu))
nb_ambiguC
####### longueur de chaque sequence  ##########################
longueurC <- sapply(consensusCnt, length)

####### pourcentage des sites ambigu  ##########################
pct_ambiguC <- nb_ambiguC / longueurC * 100

####### resumé des valeurs de la diversité intra hote  ##########################
diversiteControle <- data.frame(
  Sequence = names(consensusC),
  Sites_ambigus = nb_ambiguC,
  longueur= longueurC,
  Pourcentage = round(pct_ambiguC, 2)
)
diversiteControle$Sequence<-gsub('_NS5B','',diversiteControle$Sequence)
diversiteControle
############### graphe ############

pc <- ggplot(diversiteControle,
             aes(x = reorder(Sequence, Sites_ambigus),
                 y = Sites_ambigus)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = Sites_ambigus),
            hjust = -0.2,
            size = 3) +
  coord_flip() +
  labs(
    x = "Individus",
    y = "Nombre de sites ambigus"
  ) +
  theme_classic()
pc

ggsave(
  filename = "~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/resultat_grapheutile/sites_ambigusControle.png",
  plot = pc,
  width = 8,
  height = 6,
  dpi = 300
)

names(diversiteControle)
######### STATISTIQUE sur le pourcentage de la diversité chez les cotroles #########
summary(diversiteControle$Pourcentage)

############### graphe ############

pDVC <- ggplot(diversiteControle,
             aes(x = reorder(Sequence, Pourcentage),
                 y = Pourcentage)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = Pourcentage),
            hjust = -0.2,
            size = 3) +
  coord_flip() +
  labs(
    x = "Individus",
    y = "Diversité Virale (%)"
  ) +
  theme_classic()
pDVC

ggsave(
  filename = "~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/resultat_grapheutile/diversiteviraleControle.png",
  plot = pDVC,
  width = 8,
  height = 6,
  dpi = 300
)

########################### COMPARAISON DE LA DIVERSITÉ VIRALE FOYER ET CONTROLE #########
names(diversite)
names(diversiteControle)
############ test wilcoxon 
test_mw<-wilcox.test(diversite$Pourcentage, diversiteControle$Pourcentage)
############# boxplot de comparaison  
diversite$Groupe <- "Foyer"
diversiteControle$Groupe <- "Controle"

################### selectionner les variables utile pour la comparaison#####
diversite1<- diversite %>%
  select(Sequence,Sites_ambigus,longueur,Pourcentage,Groupe)

diversiteControle1<- diversiteControle %>%
  select(Sequence,Sites_ambigus,longueur,Pourcentage,Groupe)


dataComparaisondiversit <- rbind(diversite1, diversiteControle1)

  
PDCompare<-ggplot(dataComparaisondiversit,
                 aes(x = Groupe, y = Pourcentage, fill = Groupe)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Controle" = "#2196F3", "Foyer" = "#F44336")) +
  annotate("text", x = 1.5, y = max(dataComparaisondiversit$Pourcentage, na.rm = TRUE),
           label = paste0("p = ", round(test_mw$p.value, 2)),
           size = 4) +
  theme_classic() +
  labs(
    x     = "Groupe",
    y     = "Diversité intra-hôte ")
PDCompare
ggsave(
  filename = "~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/resultat_grapheutile/diversitefoyercontrole.png",
  plot = PDCompare,
  width = 8,
  height = 6,
  dpi = 300
)

aggregate(Pourcentage ~ Groupe,
          dataComparaisondiversit,
          summary)

############################################################
#  Matrice de distances génétiques de la diversité intra hote
# modèle TN93 = recommandé pour VHC
############################################################
# Charger l'alignement
aln_DVG<- read.dna("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/aligned_diversiteglobale.fasta", format = "fasta")

dist_DVG <- dist.dna(aln_DVG, model = "TN93",
                     pairwise.deletion = TRUE)

dist_DVG_df     <- as.matrix(dist_DVG)

# Conversion en data frame long
dist_DVG_df <- as.data.frame(as.table(as.matrix(dist_DVG)))

# Renommer les colonnes
colnames(dist_DVG_df) <- c("Sequence1", "Sequence2", "Distance")

# Supprimer la diagonale
dist_DVG_df <- subset(dist_DVG_df, Sequence1 != Sequence2)

# Garder une seule fois chaque paire
dist_DVG_df <- dist_DVG_df[as.character(dist_DVG_df$Sequence1) <
                             as.character(dist_DVG_df$Sequence2), ]

############ les stats ##############

distances <- as.vector(dist_DVG)

summary(distances)


##########################################################################
#### distance genetique (consensus)  chez les controle #########################
#########################################################################
# Charger l'alignement
aln_DVC<- read.dna("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/aligned_sequence_diversitecontrole.fasta", format = "fasta")

dist_DVC <- dist.dna(aln_DVC, model = "TN93",
                     pairwise.deletion = TRUE)

# Conversion en data frame long
dist_DVC_df <- as.data.frame(as.table(as.matrix(dist_DVC)))

# Renommer les colonnes
colnames(dist_DVC_df) <- c("Sequence1", "Sequence2", "Distance")

# Supprimer la diagonale
dist_DVC_df <- subset(dist_DVC_df, Sequence1 != Sequence2)

# Garder une seule fois chaque paire
dist_DVC_df <- dist_DVC_df[as.character(dist_DVC_df$Sequence1) <
                             as.character(dist_DVC_df$Sequence2), ]

############ les stats ##############

summary(dist_DVC_df$Distance)
quantile(dist_DVC_df$Distance,
         probs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95),
         na.rm = TRUE)
################ test la normalité 
qqnorm(dist_DVC_df$Distance,
       main = "QQ-plot des distances TN93 intra-foyer")

qqline(dist_DVC_df$Distance,
       col = "red",
       lwd = 2)

shapiro.test(dist_DVC_df$Distance)

ggplot(dist_DVC_df, aes(x = Distance)) +
  geom_histogram(aes(y = ..density..),
                 bins = 30,
                 fill = "steelblue",
                 colour = "black",
                 alpha = 0.6) +
  geom_density(colour = "red",
               linewidth = 1) +
  theme_classic() +
  labs(
    x = "Distance génétique TN93",
    y = "Densité",
  )

summary(dist_DVC_df$Distance)

############# boxplot de comparaison  et test wilcoxon 
dist_DVF_df$Groupe <- "Foyer"
dist_DVC_df$Groupe <- "Controle"
test_DG<- wilcox.test(dist_DVF_df$Distance, dist_DVC_df$Distance)

dataComparaisonDV <- rbind(dist_DVF_df, dist_DVC_df)
dataComparaisonDV

PDV<-ggplot(dataComparaisonDV,
            aes(x = Groupe, y = Distance, fill = Groupe)) +
  geom_boxplot() +
  annotate("text", x = 1.5, y = max(dataComparaisonDV$Distance, na.rm = TRUE),
           label = paste0("p <0.05 "),
           size = 4) +
  theme_classic() +
  labs(
    x     = "Groupe",
    y     = "Distance génétique (TN93)")

PDV
ggsave(
  filename = "~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/resultat_grapheutile/distance_DVfoyercontrole.png",
  plot = PDV,
  width = 8,
  height = 6,
  dpi = 300
)
############ stats
aggregate(Distance ~ Groupe,
          dataComparaisonDV,
          summary)



################ graphe des evenement de la transmission dans les foyers ########
mediane_globale <- median(dist_DVC_df$Distance, na.rm = TRUE) 



Q1 <- quantile(dist_DVC_df$Distance,
               probs = 0.05,
               na.rm = TRUE)

Q3 <- quantile(dist_DVC_df$Distance,
               probs = 0.75,
               na.rm = TRUE)

couleurs_menage1 <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(unique(resume_DVmenage$menages))),
  sort(unique(resume_DVmenage$menages))
)

PDVintra<- ggplot(resume_DVmenage,
                  aes(x = menages,
                      y = distance ,
                      colour = menages)) +
  
  geom_jitter(width = 0.1,
              size = 3) +
  
  geom_hline(yintercept = mediane_globale,
             colour = "red",
             linetype = "dashed",
             linewidth = 1) +
  geom_hline(yintercept = 0.03,
             colour = "darkgreen",
             linetype = "dotted") + 
  geom_hline(yintercept = Q1,
             colour = "darkgreen",
             linetype = "dotted") +
  
  geom_hline(yintercept = Q3,
             colour = "darkgreen",
             linetype = "dotted") +
  
  theme_classic() +
  scale_colour_manual(values = couleurs_menage1) +
  labs(
    x = "Ménage ",
    y = "Distance génétique"
  )

PDVintra

ggsave(
  filename = "~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/resultat_grapheutile/evenementtransmission.png",
  plot = PDVintra,
  width = 8,
  height = 6,
  dpi = 300
)

##########################################################################
############## distance genetique globale dans la region #########################
#########################################################################


# Charger l'alignement
aln_DG<- read.dna("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/aligned_sequence_FC_BON.fasta", format = "fasta")

dist_DVG <- dist.dna(aln_DG, model = "TN93",
                     pairwise.deletion = TRUE)

dist_DVG_df     <- as.matrix(dist_DVG)

# Conversion en data frame long
dist_DVG_df <- as.data.frame(as.table(as.matrix(dist_DVG)))

# Renommer les colonnes
colnames(dist_DVG_df) <- c("Sequence1", "Sequence2", "Distance")

# Supprimer la diagonale
dist_DVG_df <- subset(dist_DVG_df, Sequence1 != Sequence2)

# Garder une seule fois chaque paire
dist_DVG_df <- dist_DVG_df[as.character(dist_DVG_df$Sequence1) <
                             as.character(dist_DVG_df$Sequence2), ]

############ les stats ##############
distances <- as.vector(dist_DVG)
summary(distances)
quantile(distances,
         pobs=c(0.05,0.90,0.95))




###########################################################################
############## arbre phylogenetique du groupe foyer ######################
##########################################################################

############# alignement de la sequequence #######################
## mafft --auto --reorder sequence_foyer.fasta >aligned_sequence_foyer.fasta

alnfoyer <- read.dna("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/aligned_sequence_foyer.fasta", format = "fasta")

### enlever les NS5B des labels
rownames(alnfoyer)<- gsub("_NS5B","",rownames(alnfoyer))


dist_FP <- dist.dna(alnfoyer, model = "TN93", pairwise.deletion = TRUE)

distFP_df     <- as.matrix(dist_FP)





############################################################
# 3) Calcul distances intra vs inter ménage
############################################################
# Construire toutes les paires
paires <- expand.grid(ind1 = foyer1$ID_TERRAIN,
                      ind2 = foyer1$ID_TERRAIN,
                      stringsAsFactors = FALSE) %>%
  filter(ind1 < ind2) %>%   # éviter doublons et diagonale
  left_join(foyer1, by = c("ind1" = "ID_TERRAIN")) %>%
  rename(menage1 = menages) %>%
  left_join(foyer1, by = c("ind2" = "ID_TERRAIN")) %>%
  rename(menage2 = menages) %>%
  mutate(
    type     = ifelse(menage1 == menage2, "intra", "inter"),
    distance = mapply(function(i, j) distFP_df[i, j], ind1, ind2)
  )

summary(paires$distance)

# Vérification

cat("Paires intra :", sum(paires$type == "intra"), "\n")
cat("Paires inter :", sum(paires$type == "inter"), "\n")

############################################################
# 4) Statistiques par ménage
############################################################

# Résumé par ménage
resume_menage <- paires %>%
  filter(type == "intra") %>%
  group_by(menage1) %>%
  summarise(
    n_paires          = n(),
    dist_moyenne      = mean(distance),
    dist_min          = min(distance),
    nb_transmission   = sum(distance ),
    pct_transmission  = mean(distance ) * 100,
    sd_dist   = sd(distance,   na.rm = TRUE)
  ) %>%
  rename(menage = menage1)


print(resume_menage)

# Moyenne intra et inter globale
stats_global <- paires %>%
  group_by(type) %>%
  summarise(
    n        = n(),
    mean_dist = mean(distance, na.rm = TRUE),
    sd_dist   = sd(distance,   na.rm = TRUE),
    min_dist  = min(distance,  na.rm = TRUE),
    max_dist  = max(distance,  na.rm = TRUE),
    median_dist  = median(distance,  na.rm = TRUE),
    Q1_dist  = quantile(distance, probs=0.25, na.rm = TRUE),
    Q3_dist  = quantile(distance, probs=0.75, na.rm = TRUE)
    
  )

print(stats_global)



############################################################
# 5) Test statistique intra vs inter
# Mann-Whitney (non paramétrique, adapté aux petits échantillons)
############################################################

############################################################
# Étape 2 : Tests statistiques
############################################################

# Mann-Whitney : alternative = "less" : on teste si intra < inter
test_mw <- wilcox.test(
  paires$distance[paires$type == "intra"],
  paires$distance[paires$type == "inter"],
  alternative = "less"
)
cat("\n===== TEST MANN-WHITNEY intra < inter =====\n")
cat("Mann-Whitney p-value :", test_mw$p.value, "\n")

if (test_mw$p.value < 0.05) {
  cat("→ Distance intra SIGNIFICATIVEMENT inférieure à inter\n")
  cat("→ Transmission intra-ménage PROBABLE\n")
} else {
  cat("→ Pas de différence significative intra vs inter\n")
}

############################################################
# 7) Heatmap des distances colorée par ménage
############################################################

# Ordonner par ménage
ordre <- foyer1 %>% arrange(menages,ID_TERRAIN) %>% pull(ID_TERRAIN)
ordre

ordre <- ordre[ordre %in% rownames(distFP_df)]
ordre


dist_ordered <- distFP_df[ordre, ordre]
dist_ordered
#melt() transforme une matrice large en format long — elle "fait fondre" la matrice en 3 colonnes pour ggplot.
dist_melt    <- melt(dist_ordered)
colnames(dist_melt) <- c("ind1", "ind2", "distance")
dist_melt$ind1<-as.character(dist_melt$ind1)
dist_melt$ind2<-as.character(dist_melt$ind2)
# Annoter avec les ménages

# Ajouter le ménage à ind1
dist_melt <- dist_melt %>%
  left_join(
    foyer1 %>% select(ID_TERRAIN, menages),
    by = c("ind1" = "ID_TERRAIN")
  ) %>%
  mutate(ind1 = paste0(menages, "-", ind1)) %>%
  rename(menage_ind1 = menages)

# Ajouter le ménage à ind2
dist_melt <- dist_melt %>%
  left_join(
    foyer1 %>% select(ID_TERRAIN, menages),
    by = c("ind2" = "ID_TERRAIN")
  ) %>%
  mutate(ind2 = paste0(menages, "-", ind2)) %>%
  rename(menage_ind2 = menages)

dist_melt

p_heatmap <- ggplot(dist_melt, aes(x = ind1, y = ind2, fill = distance)) +
  geom_tile() +
  
  scale_fill_gradient(
    low = "white",
    high = "darkblue",
    name = "Distance\ngénétique"
  ) +
  
  theme_minimal() +
  
  

  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 16        
    ),
    
    axis.text.y = element_text(
      size = 16        
    ),
    
    plot.title = element_text(
      size = 16,
      face = "bold"
    ),
    
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 12)
  ) +
  
  labs(
    x = "",
    y = ""
  )

p_heatmap
ggsave(
  "~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/resultat_grapheutile/heatmap_distancesfoyerprimary.png",
  p_heatmap,
  width = 18,
  height = 16,
  dpi = 300
)





## Détermination des séparations entre ménages
labels <- data.frame(ID = ordre) %>%
  mutate(menage = sub("-.*", "", ID))

bornes <- labels %>%
  count(menage) %>%
  mutate(fin = cumsum(n))

separations <- bornes$fin + 0.5
separations <- separations[-length(separations)]


##############################################################

dist_melt$type <- ifelse(dist_melt$menage_ind1 == dist_melt$menage_ind2,
                         "Intra",
                         "Inter")

dist_test <- dist_melt %>%
  filter(ind1 != ind2) %>%          # enlève la diagonale
  filter(ind1 < ind2)               # garde uniquement un triangle de la matrice

table(dist_test$type)
nrow(dist_test)
wilcox.test(distance ~ type,
            data = dist_test)

dist_test %>%
  group_by(type) %>%
  summarise(
    n = n(),
    mediane = median(distance),
    moyenne = mean(distance),
    min = min(distance),
    max = max(distance)
  )

dist_test %>%
  filter(type == "Inter", distance == 0)

## Heatmap
p_heatmap1 <- ggplot(dist_melt,
                    aes(x = ind1,
                        y = ind2,
                        fill = distance)) +
  
  geom_tile(color = "grey85", linewidth = 0.2) +
  
  geom_vline(xintercept = separations,
             colour = "black",
             linewidth = 0.8) +
  
  geom_hline(yintercept = separations,
             colour = "black",
             linewidth = 0.8) +
  
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    name = "Distance\ngénétique"
  ) +
  
  coord_equal() +
  
  labs(x = "", y = "") +
  
  theme_classic(base_size = 16) +
  
  theme(
    
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 17
    ),
    
    axis.text.y = element_text(
      size = 17
    ),
    plot.title = element_text(
      size = 17,
      face = "bold"
    ),
    axis.title = element_blank(),
    
    panel.grid = element_blank(),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    legend.text = element_text(size = 12),
    
    plot.margin = margin(10,10,10,10)
  )

p_heatmap1


ggsave(
  "~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/resultat_grapheutile/heatmap_distancesfoyerprimary1.png",
  p_heatmap1,
  width = 18,
  height = 16,
  dpi = 300
)


## Heatmap
p_heatmap2 <- ggplot(dist_test,
                     aes(x = ind1,
                         y = ind2,
                         fill = distance)) +
  
  geom_tile(color = "grey85", linewidth = 0.2) +
  
  geom_vline(xintercept = separations,
             colour = "black",
             linewidth = 0.8) +
  
  geom_hline(yintercept = separations,
             colour = "black",
             linewidth = 0.8) +
  
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    name = "Distance\ngénétique"
  ) +
  
  coord_equal() +
  
  labs(x = "", y = "") +
  
  theme_classic(base_size = 16) +
  
  theme(
    
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 17
    ),
    
    axis.text.y = element_text(
      size = 17
    ),
    plot.title = element_text(
      size = 17,
      face = "bold"
    ),
    axis.title = element_blank(),
    
    panel.grid = element_blank(),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    legend.text = element_text(size = 12),
    
    plot.margin = margin(10,10,10,10)
  )

p_heatmap2


ggsave(
  "~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/resultat_grapheutile/heatmap_distancesfoyerprimary2.png",
  p_heatmap2,
  width = 18,
  height = 16,
  dpi = 300
)



alnfoyer["M04-3312412", ] ==
  alnfoyer["M06-3322292", ]


############################################################
# 8) Boxplot intra vs inter
############################################################
paires
p_boxplot <- ggplot(paires, aes(x = type, y = distance, fill = type)) +
  geom_boxplot() +
  
  scale_fill_manual(values = c("intra" = "#2196F3", "inter" = "#F44336")) +
  annotate("text", x = 1.5, y = max(paires$distance, na.rm = TRUE),
           label = paste0("p = ", round(test_mw$p.value, 4)),
           size = 4) +
  theme_classic() +
  labs(
    x     = "Type de comparaison",
    y     = "Distance génétique (TN93)")

p_boxplot

ggsave("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/resultat_grapheutile/boxplot_intra_inter.png", p_boxplot,
       width = 8, height = 6, dpi = 150)

############################################################
# 6) Phylogénie : methode distance
############################################################
# Vérifier s'il reste des NA
sum(is.na(as.matrix(dist_FP)))
# Construire l'arbre (Neighbor-Joining)

tree1 <- njs(dist_FP)
tree1 <- ladderize(tree1)


# Colorier les tips par ménage
tip_menage <- foyer1$menages[match(tree1$tip.label,foyer$ID_TERRAIN)]

couleurs_menage2 <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(unique(foyer1$menages))),
  unique(foyer1$menages)
)
tip_colors <- couleurs_menage2[tip_menage]

tree1$tip.label <- paste0(
  foyer1$menages[match(tree1$tip.label, foyer1$ID_TERRAIN)],
  "-",
  tree1$tip.label
)
tree1$tip.label
# Plot arbre
plot(tree1)

png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/resultat_grapheutile/phylogenie_VHCfoyerprimDist.png",
    width = 3000,
    height = 2200,
    res = 300)
par(font=2)
plot(tree1,
     tip.color = tip_colors,
     cex       = 1,
     type = "phylogram",
     edge.width = 1.5)

legend("topright",
       legend = names(couleurs_menage2),
       fill   = couleurs_menage2,
       title  = "Ménage",
       cex    = 0.7)
add.scale.bar()


dev.off()

###################### bootstrap de la distance ###########
# Vérifier les séquences trop courtes ou trop gappées
checkAlignment(alnfoyer)



rownames(alnfoyer)

rownames(alnfoyer) <- paste0(
  foyer1$menages[match(rownames(alnfoyer), foyer1$ID_TERRAIN)],
  "-",
  rownames(alnfoyer)
)

################## exporter  ######################################
write.dna(alnfoyer,
          file = "alignedfoyer.fasta",
          format = "fasta")


set.seed(42)
bs <- boot.phylo(
  phy = tree1,
  x = alnfoyer,
  FUN = function(xx)
    nj(dist.dna(xx, model = "TN93",pairwise.deletion = TRUE)),
  B = 100,
  rooted = FALSE
)
bs
plot(tree1,
     tip.color = tip_colors,
     type = "phylogram",
     cex = 1,
     edge.width = 1.5)

nodelabels(ifelse(bs >= 70, bs, ""),
           cex = 0.5,
           frame = "none")


##############################################################################
# 2) phylogenie : arbre avec la vraisemblance et Sélection du meilleur modèle
#############################################################################



# Tester les modèles — peut prendre quelques minutes
cat("[INFO] Sélection du modèle...\n")
modtest <- modelTest(alnfoyer,
                     tree  = tree1,
                     model = c("JC", "HKY", "TN93", "GTR"),
                     G     = TRUE,
                     I     = TRUE)

# Afficher les résultats triés par AIC
modtest_sorted <- modtest[order(modtest$AIC), ]
print(head(modtest_sorted, 10))

# Meilleur modèle
best_model <- modtest_sorted$Model[1]
cat("[INFO] Meilleur modèle (AIC) :", best_model, "\n")

best_model_clean <- gsub("\\+G\\(4\\)|\\+I", "", best_model)
cat("[INFO] Meilleur modèle nettoyer (AIC) :", best_model_clean, "\n")
#################### lire l'arbre issue de Lire l’arbre RAxML ###########
# raxmlHPC \ -s alignedfoyer.phy \ -n VHC \ -m GTRGAMMA \ -p 12345 \-x 12345 \ -# 1000 \  -f a #

tree_raxml <- read.tree("RAxML_bipartitions.VHC")
tree_raxml <- ladderize(tree_raxml)
plot(tree_raxml, cex = 0.8) 

##### arbre en cercle
plot(tree_raxml,
     type = "fan",
     cex = 0.8)

plot(tree_raxml, cex = 0.8)
nodelabels(tree_raxml$node.label, cex = 0.7)

plot(tree_raxml,
     type = "phylogram",
     cex = 0.7,
     edge.width = 1.5)

############################ definir mes couleur 
# Colorier les tips par ménage
tree_raxml$tip.label
foyer2<-foyer1
foyer2$ID_TERRAIN<-paste0(foyer2$menages,"-",foyer2$ID_TERRAIN)
tip_menage2 <- foyer2$menages[match(tree_raxml$tip.label,foyer2$ID_TERRAIN)]
tip_menage2
couleurs_menage3 <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(unique(foyer2$menages))),
  unique(foyer2$menages)
)

couleurs_menage3

tip_colors <- couleurs_menage3[tip_menage2]
tip_colors

tree_raxml$node.label

boot <- as.numeric(tree_raxml$node.label)
boot

png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/resultat_grapheutile/phylogenie_VHCfoyerprimaryRAxml.png",
    width = 3000,
    height = 2200,
    res = 300)
par(font = 2)
plot(tree_raxml,
     tip.color = tip_colors,
     type = "phylogram",
     cex = 1,
     edge.width = 1.5)
nodelabels(
  text = ifelse(boot >= 70, boot, ""),
  frame = "none",
  cex = 0.7,
  adj = c(-0.2, -0.2)
)

legend("topright",
       legend = names(couleurs_menage3),
       fill   = couleurs_menage3,
       cex    = 0.7)

add.scale.bar(x = 0, y = max(nodeHeights(tree_raxml)))

dev.off()
graphics.off()

########## cluster picker #####################
# java -jar ClusterPicker.jar \ -t RAxML_bipartitions.VHC \ -s alignedfoyer.fasta \ -o output_cluster \ -d 0.03 \ -b 90


################### foyer et controle ###################
############### renommer le name seq foyer ###########

# mafft --auto --reorder sequence_FC_BON.fasta > aligned_sequence_FC_BON.fasta

# Charger l'alignement
alnFC <- read.dna("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/aligned_sequence_FC_BON.fasta", format = "fasta")

write.dna(alnFC,
          file = "alignedfoyercontroraxml.fasta",
          format = "fasta")

# Distances et arbre
dist_mat <- dist.dna(alnFC, 
                     model = "TN93", 
                     pairwise.deletion = TRUE)



dist_mat_df     <- as.matrix(dist_mat)

# Conversion en data frame long
dist_mat_df <- as.data.frame(as.table(as.matrix(dist_mat_df)))

# Renommer les colonnes
colnames(dist_mat_df) <- c("Sequence1", "Sequence2", "Distance")

# Supprimer la diagonale
dist_mat_df <- subset(dist_mat_df, Sequence1 != Sequence2)

# Garder une seule fois chaque paire
dist_mat_df <- dist_mat_df[as.character(dist_mat_df$Sequence1) <
                             as.character(dist_mat_df$Sequence2), ]

############ les stats ##############
distances <- as.vector(dist_mat)
summary(distances)
quantile(distances,
         pobs=c(0.05,0.90,0.95))



treeFC     <- njs(dist_mat)
treeFC     <- ladderize(treeFC)



# Extraire le type (FOYER ou CTRL) depuis les noms des tips
tip_labels <- treeFC$tip.label

tip_type   <- ifelse(grepl("^M[0-9]{2}-", tip_labels), "Foyer", "Controle")

# Extraire le ménage pour les foyers (depuis metadata)

tip_menage <- rep(NA, length(tip_labels))
for (i in seq_along(tip_labels)) {
  if (tip_type[i] == "Foyer") {
    individu <- gsub("^M[0-9]{2}-", "", tip_labels[i])
    m <- foyer1$menages [match(individu, foyer1$ID_TERRAIN)]
    tip_menage[i] <- ifelse(length(m) > 0 && !is.na(m), m, "INCONNU")
  }
}

tip_menage

# Couleurs :
# Controle → gris
# Foyer → couleur par ménage
menages_uniques  <- unique(na.omit(tip_menage))
couleurs_menages <- setNames(
  rainbow(length(menages_uniques)),
  menages_uniques
)

tip_colors <- ifelse(
  tip_type == "Controle",
  "grey60",                          # controle = gris
  couleurs_menages[tip_menage]        # foyer = couleur par ménage
)


# Plot
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/resultat_grapheutile/phylogenie_foyer_controleDist.png", width = 1600, height = 1200, res = 150)
par(font = 2)
plot(treeFC,
     tip.color = tip_colors,
     type="phylogram",
     cex=0.6,
     edge.width = 1.5)

# Légende ménages
legend("topright",
       legend = c(names(couleurs_menages), "Controle"),
       fill   = c(couleurs_menages, "grey60"),
       title  = "Ménage",
       cex    = 0.65)


add.scale.bar(x = 0, y = max(nodeHeights(treeFC)))


dev.off()


##############################################################################
# phylogenie : arbre avec la vraisemblance et Sélection du meilleur modèle
#############################################################################

# Tester les modèles — peut prendre quelques minutes
cat("[INFO] Sélection du modèle...\n")
modtestFC <- modelTest(alnFC,
                       tree  =  treeFC,
                       model = c("JC", "HKY", "TN93", "GTR"),
                       G     = TRUE,
                       I     = TRUE)

# Afficher les résultats triés par AIC
modtest_sortedFC <- modtestFC[order(modtestFC$AIC), ]
print(head(modtest_sortedFC, 10))

# Meilleur modèle
best_modelFC <- modtest_sortedFC$Model[1]
cat("[INFO] Meilleur modèle (AIC) :", best_modelFC, "\n")

############################################################
# 3) Arbre Maximum de Vraisemblance
############################################################
############ faire le Raxml ###########
# raxmlHPC  -s alignedfoyercontroraxml.fasta -n VHCFC  -m GTRGAMMAI  -p 12345 -x 12345  -# 1000 -f a #

tree_raxmlFC <- read.tree("RAxML_bipartitions.VHCFC")

tree_raxmlFC <- ladderize(tree_raxmlFC)


############################ definir mes couleur 
# Extraire le type (FOYER ou CTRL) depuis les noms des tips
tip_labels1 <- tree_raxmlFC$tip.label
tip_labels1
tip_type1   <- ifelse(grepl("^M[0-9]{2}-", tip_labels1), "Foyer", "Controle")
tip_type1
# Extraire le ménage pour les foyers (depuis metadata)

tip_menage1 <- rep(NA, length(tip_labels1))
for (i in seq_along(tip_labels1)) {
  if (tip_type1[i] == "Foyer") {
    individu <- gsub("^M[0-9]{2}-", "", tip_labels1[i])
    m <- foyer1$menages [match(individu, foyer1$ID_TERRAIN)]
    tip_menage1[i] <- ifelse(length(m) > 0 && !is.na(m), m, "INCONNU")
  }
}

tip_menage1

# Couleurs :
# Controle → gris
# Foyer → couleur par ménage
menages_uniques1  <- unique(na.omit(tip_menage1))
couleurs_menages1 <- setNames(
  rainbow(length(menages_uniques1)),
  menages_uniques1
)

tip_colors1 <- ifelse(
  tip_type1 == "Controle",
  "grey60",                          # controle = gris
  couleurs_menages1[tip_menage1]        # foyer = couleur par ménage
)

tip_colors1 <- couleurs_menages1[tip_menage1]

boot1 <- as.numeric(tree_raxmlFC$node.label)
boot1

png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/resultat_grapheutile/phylogenie_VHCfoyercontroleRAxml.png",
    width = 3000,
    height = 2200,
    res = 300)
par(font = 2)   # 2 = gras
plot(tree_raxmlFC,
     tip.color = tip_colors1,
     type = "phylogram",
     cex = 0.6,
     edge.width = 1.5)

nodelabels(
  text = ifelse(boot1 >= 70, boot1, ""),
  frame = "none",
  cex = 0.6,
  adj = c(-0.5, -0.3)
)

legend("topright",
       legend = names(couleurs_menages1),
       fill   = couleurs_menages1,
       title  = "Ménage",
       cex    = 0.7)

add.scale.bar(x = 0, y = max(nodeHeights(tree_raxmlFC)))

dev.off()


########## cluster picker #####################
# java -jar ClusterPicker.jar  -t RAxML_bipartitions.VHCFC  -s alignedfoyer.fasta  -o output_cluster  -d 0.03  -b 90
sequence_FC_BON.fasta
RAxML_bipartitions.VHCFC

#########################################################################
################## distance patristique avec les deux arbres ML  #################
dist_pat<-cophenetic.phylo(tree_raxml)
dist_pat <- as.data.frame(as.table(dist_pat))
colnames(dist_pat) <- c("ID1","ID2","Distance_patristique")
str(dist_pat)
######### convertir en caractère
dist_pat$ID1 <- as.character(dist_pat$ID1)
dist_pat$ID2 <- as.character(dist_pat$ID2)

dist_pat <- dist_pat %>%
  filter(ID1 != ID2)
dist_pat <- dist_pat %>%
  rowwise() %>%
  mutate(
    pair = paste(sort(c(ID1, ID2)), collapse = "_")
  ) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)

foyer2

dist_pat_df <- dist_pat %>%
  left_join(foyer2, by = c("ID1" = "ID_TERRAIN")) %>%
  rename(menage1 = menages) %>%
  left_join(foyer2, by = c("ID2" = "ID_TERRAIN")) %>%
  rename(menage2 = menages) %>%
  mutate(type = ifelse(menage1 == menage2,
                       "Intra-foyer",
                       "Inter-foyer"))

wilcox.test(
  Distance_patristique ~ type,
  data = dist_pat_df
)

table(dist_pat_df$type)
intra <- subset(dist_pat_df, type == "Intra-foyer") 
inter <- subset(dist_pat_df, type == "Inter-foyer")
summary(intra$Distance_patristique)
summary(inter$Distance_patristique)




ggplot(dist_pat_df,
       aes(x = type,
           y = Distance_patristique)) +
  geom_boxplot() +
  theme_minimal()


agg <- aggregate(
  Distance_patristique ~ menage1,
  data = intra,
  mean
)

agg[order(agg$Distance_patristique), ]


############################################ FC ######################
tree_raxmlFC

############################ distance patristique #################
dist_patFC<-cophenetic.phylo(tree_raxmlFC)
dist_patFC <- as.data.frame(as.table(dist_patFC))
colnames(dist_patFC) <- c("ID1","ID2","Distance_patristique")
str(dist_patFC)
######### convertir en caractère
dist_patFC$ID1 <- as.character(dist_patFC$ID1)
dist_patFC$ID2 <- as.character(dist_patFC$ID2)

dist_patFC <- dist_patFC %>%
  filter(ID1 != ID2)
dist_patFC <- dist_patFC %>%
  rowwise() %>%
  mutate(
    pair = paste(sort(c(ID1, ID2)), collapse = "_")
  ) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)




wilcox.test(
  Distance_patristique ~ type,
  data = dist_patFC
)

table(dist_patFC$type)
intraFC <- subset(dist_patFC, type == "Intra-foyer") 
interFC <- subset(dist_patFC, type == "Inter-foyer")
summary(intraFC$Distance_patristique)
summary(interFC$Distance_patristique)




ggplot(dist_patFC,
       aes(x = type,
           y = Distance_patristique)) +
  geom_boxplot() +
  theme_minimal()


aggFC <- aggregate(
  Distance_patristique ~ menage1,
  data = intraFC,
  mean
)

aggFC[order(agg$Distance_patristique), ]







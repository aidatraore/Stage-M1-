######################################################################
# Phylogénie + Distances génétiques intra/inter ménage VHC: stage M1
######################################################################
### chargement des librairies necessaire ######
library(rmarkdown)
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


####### definir le Repertoire de travaille #################
setwd("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/")

############################################################
#  Charger l'alignement et les métadonnées
############################################################

# Alignement déjà fait avec MAFFT en bash :
# mafft --auto --reorder Primaryall.fasta >aligned_primary.fasta
aln <- read.dna("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/aligned_primary.fasta", format = "fasta")
### nom des sequences
### enlever les NS5B des labels
rownames(aln)<- gsub("-NS5B","",rownames(aln))
rownames(aln)
### Tableau avec noms + tailles
data.frame(
  Sequence = rownames(aln),
  Length = apply(as.character(aln), 1, function(x) sum(x != "-"))
)
## Longueur de l’alignement total
ncol(aln)
# Fichier metadata :
foyer <- read_excel("~/MASTER BIOINFORMATIQUE/stageM1/Projets/foyergeno2.xlsx")

## Remplacer les   "-"
foyer$ID <- gsub("-NS5B","", foyer$ID)
foyer$ID
############################################################
# 2) Matrice de distances génétiques
# modèle TN93 = recommandé pour VHC
############################################################

dist_matrix <- dist.dna(aln, model = "TN93", pairwise.deletion = TRUE)
dist_df     <- as.matrix(dist_matrix)
############################################################
# 3) Calcul distances intra vs inter ménage
############################################################
names(foyer)

foyer <- foyer %>%
  select(ID, Ménages)

colnames(foyer)

# Renommer au chargement pour éviter les problèmes d'accents
names(foyer) <- trimws(names(foyer))
foyer <- foyer %>%
  rename(menage = Ménages)
colnames(foyer)
# Construire toutes les paires
paires <- expand.grid(ind1 = foyer$ID,
                      ind2 = foyer$ID,
                      stringsAsFactors = FALSE) %>%
  filter(ind1 < ind2) %>%   # éviter doublons et diagonale
  left_join(foyer, by = c("ind1" = "ID")) %>%
  rename(menage1 = menage) %>%
  left_join(foyer, by = c("ind2" = "ID")) %>%
  rename(menage2 = menage) %>%
  mutate(
    type     = ifelse(menage1 == menage2, "intra", "inter"),
    distance = mapply(function(i, j) dist_df[i, j], ind1, ind2)
  )

# Vérification
head(paires)
cat("Paires intra :", sum(paires$type == "intra"), "\n")
cat("Paires inter :", sum(paires$type == "inter"), "\n")

############################################################
# 4) Statistiques par ménage
############################################################

# Moyenne intra et inter globale
stats_global <- paires %>%
  group_by(type) %>%
  summarise(
    n        = n(),
    mean_dist = mean(distance, na.rm = TRUE),
    sd_dist   = sd(distance,   na.rm = TRUE),
    min_dist  = min(distance,  na.rm = TRUE),
    max_dist  = max(distance,  na.rm = TRUE)
  )

print(stats_global)

# Moyenne intra par ménage
stats_menage <- paires %>%
  filter(type == "intra") %>%
  group_by(menage1) %>%
  summarise(
    n         = n(),
    mean_dist = mean(distance, na.rm = TRUE),
    sd_dist   = sd(distance,   na.rm = TRUE)
  ) %>%
  rename(menage = menage1)

print(stats_menage)

############################################################
# 5) Test statistique intra vs inter
# Mann-Whitney (non paramétrique, adapté aux petits échantillons)
############################################################

intra <- paires %>% filter(type == "intra") %>% pull(distance)
inter <- paires %>% filter(type == "inter") %>% pull(distance)

test <- wilcox.test(intra, inter, alternative = "less")
# alternative = "less" : on teste si intra < inter

cat("\n===== TEST MANN-WHITNEY intra < inter =====\n")
cat("p-value :", test$p.value, "\n")
if (test$p.value < 0.05) {
  cat("→ Distance intra SIGNIFICATIVEMENT inférieure à inter\n")
  cat("→ Transmission intra-ménage PROBABLE\n")
} else {
  cat("→ Pas de différence significative intra vs inter\n")
}

############################################################
# 6) Phylogénie : methode distance
############################################################
# Vérifier s'il reste des NA
sum(is.na(as.matrix(dist_matrix)))
# Construire l'arbre (Neighbor-Joining)

tree1 <- njs(dist_matrix)
tree1 <- ladderize(tree1)
tree1$tip.label

# Colorier les tips par ménage
tip_menage <- foyer$menage[match(tree1$tip.label,foyer$ID)]
library(RColorBrewer)
couleurs_menage <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(unique(foyer$menage))),
  unique(foyer$menage)
)
tip_colors <- couleurs_menage[tip_menage]

# Plot arbre


png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_VHCfoyerprimary.png",
    width = 3000,
    height = 2200,
    res = 300)

plot(tree1,
     tip.color = tip_colors,
     cex       = 1,
     no.margin = TRUE)

legend("topright",
       legend = names(couleurs_menage),
       fill   = couleurs_menage,
       title  = "Ménage",
       cex    = 0.7)
add.scale.bar()


dev.off()
graphics.off()

##############################################################################
# 2) phylogenie : arbre avec la vraisemblance et Sélection du meilleur modèle
#############################################################################

# Tester les modèles — peut prendre quelques minutes
cat("[INFO] Sélection du modèle...\n")
modtest <- modelTest(aln,
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
############################################################
# 3) Arbre Maximum de Vraisemblance
############################################################

cat("[INFO] Optimisation ML...\n")

# Initialiser avec le meilleur modèle
aln_phy <- phyDat(aln, type = "DNA")
fit_init <- pml(tree1, aln_phy)

# Optimiser
fit_ml <- optim.pml(
  fit_init,
  model         = best_model_clean,
  optInv        = TRUE,       # proportion de sites invariants
  optGamma      = TRUE,       # distribution gamma
  rearrangement = "stochastic",
  control       = pml.control(trace = 1)
)

cat("[INFO] Log-vraisemblance :", fit_ml$logLik, "\n")
fit_ml$tree
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_VHCfoyerprimaryML.png",
    width = 3000,
    height = 2200,
    res = 300)

plot(fit_ml$tree,
     tip.color = tip_colors,
     cex       = 1,
     no.margin = TRUE)
legend("topright",
       legend = names(couleurs_menage),
       fill   = couleurs_menage,
       cex    = 0.7)

add.scale.bar(x = 0, y = max(nodeHeights(fit_ml$tree)))


dev.off()
graphics.off()



############################################################
# Bootstrap 1000 réplicats
############################################################

cat("[INFO] Bootstrap 1000 réplicats — patience...\n")
set.seed(42)

bs_ml <- bootstrap.pml(
  fit_ml,
  bs      = 1000,
  optNni  = TRUE,
  control = pml.control(trace = 0)
)

# Arbre avec valeurs bootstrap
tree_bs <- plotBS(fit_ml$tree,
                  bs_ml,
                  type = "none")   # récupérer sans plotter
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_VHCfoyerprimaryMLB.png",
    width = 3000,
    height = 2200,
    res = 300)

plot(tree_bs,
     tip.color = tip_colors,
     cex       = 1,
     no.margin = TRUE)
legend("topright",
       legend = names(couleurs_menage),
       fill   = couleurs_menage,
       cex    = 0.7)

add.scale.bar(x = 0, y = max(nodeHeights(tree_bs)))

dev.off()
graphics.off()
############################################################
# 7) Heatmap des distances colorée par ménage
############################################################

# Ordonner par ménage
ordre <- foyer %>% arrange(menage) %>% pull(ID)
ordre <- ordre[ordre %in% rownames(dist_df)]

dist_ordered <- dist_df[ordre, ordre]
#melt() transforme une matrice large en format long — elle "fait fondre" la matrice en 3 colonnes pour ggplot.
dist_melt    <- melt(dist_ordered)
colnames(dist_melt) <- c("ind1", "ind2", "distance")
dist_melt$ind1<-as.character(dist_melt$ind1)
dist_melt$ind2<-as.character(dist_melt$ind2)
# Annoter avec les ménages
dist_melt <- dist_melt %>%
  left_join(foyer, by = c("ind1" = "ID")) %>%
  rename(menage1 = menage)

p_heatmap <- ggplot(dist_melt, aes(x = ind1, y = ind2, fill = distance)) +
  geom_tile() +
  
  scale_fill_gradient(
    low = "white",
    high = "darkred",
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
    legend.text  = element_text(size = 10)
  ) +
  
  labs(
    title = "Matrice de distances génétiques VHC",
    x = "",
    y = ""
  )
ggsave(
  "~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/heatmap_distancesfoyerprimary.png",
  p_heatmap,
  width = 18,
  height = 16,
  dpi = 300
)

############################################################
# 8) Boxplot intra vs inter
############################################################

p_boxplot <- ggplot(paires, aes(x = type, y = distance, fill = type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = c("intra" = "#2196F3", "inter" = "#F44336")) +
  annotate("text", x = 1.5, y = max(paires$distance, na.rm = TRUE),
           label = paste0("p = ", round(test$p.value, 4)),
           size = 4) +
  theme_minimal() +
  labs(title = "Distances génétiques intra vs inter ménage — VHC",
       x     = "Type de comparaison",
       y     = "Distance génétique (TN93)")

p_boxplot

ggsave("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/boxplot_intra_inter.png", p_boxplot,
       width = 8, height = 6, dpi = 150)

############################################################
# 9) foyer primaire et quasi
############################################################

# Alignement déjà fait avec MAFFT en bash :
# mafft --auto all_primary_secondary_combined.fasta > all_combined_aligned.fasta
# Charger l'alignement
aln_combinedfoyer <- read.dna("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/all_combined_aligned.fasta", format = "fasta")
### nom des sequences
rownames(aln_combinedfoyer)

############################################################
# 2) Matrice de distances génétiques
# modèle TN93 = recommandé pour VHC
############################################################
# Distances et arbre
dist_combinedfoyer <- dist.dna(aln_combinedfoyer, model = "TN93",
                          pairwise.deletion = TRUE)

distCF_df     <- as.matrix(dist_combinedfoyer)
distCF_df

rownames(distCF_df)

############################################################
# Phylogénie
############################################################
# Vérifier s'il reste des NA
sum(is.na(as.matrix(dist_combinedfoyer)))

# Construire l'arbre (Neighbor-Joining)

tree_combinedfoyer <- njs(dist_combinedfoyer)
tree_combinedfoyer <- ladderize(tree_combinedfoyer)

# Extraire infos depuis les noms des tips
# ex: "M01_421233-NS5B_primary"
tip_labels  <- tree_combinedfoyer$tip.label
tip_menage  <- sub("_.*", "", tip_labels)            # M01
tip_type    <- sub(".*_", "", tip_labels)            # primary ou secondary

couleurs_menage2 <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(unique(foyer$menage))),
  unique(foyer$menage)
)
tip_colors2 <- couleurs_menage2[tip_menage]


# Forme des points selon primary/secondary
tip_pch <- ifelse(tip_type == "primary", 16, 17)   # rond vs triangle

# Plot
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_combinéfoyer.png", width = 1400, height = 1000, res = 150)
plot(tree_combinedfoyer,
     tip.color = tip_colors2,
     cex       = 1)
# Ajouter des symboles pour distinguer primary/secondary
#tiplabels( col = tip_colors2, cex = 0.8)
legend("topright",
       legend = c(names(couleurs_menage)),
       fill   = c(couleurs_menage, NA, NA),
       title  = "Ménage",
       cex    = 0.7)
add.scale.bar()
dev.off()

##############################################################################
# 2) phylogenie : arbre avec la vraisemblance et Sélection du meilleur modèle
#############################################################################

# Tester les modèles — peut prendre quelques minutes
cat("[INFO] Sélection du modèle...\n")
modtest1 <- modelTest(aln_combinedfoyer,
                     tree  = tree_combinedfoyer,
                     model = c("JC", "HKY", "TN93", "GTR"),
                     G     = TRUE,
                     I     = TRUE)

# Afficher les résultats triés par AIC
modtest_sorted1 <- modtest1[order(modtest1$AIC), ]
print(head(modtest_sorted1, 10))

# Meilleur modèle
best_model <- modtest_sorted1$Model[1]
cat("[INFO] Meilleur modèle (AIC) :", best_model, "\n")

best_model_clean <- gsub("\\+G\\(4\\)|\\+I", "", best_model)
cat("[INFO] Meilleur modèle nettoyer (AIC) :", best_model_clean, "\n")
############################################################
# 3) Arbre Maximum de Vraisemblance
############################################################

cat("[INFO] Optimisation ML...\n")

# Initialiser avec le meilleur modèle
aln_phy <- phyDat(aln_combinedfoyer, type = "DNA")
fit_init1 <- pml(tree_combinedfoyer, aln_phy)

# Optimiser
fit_ml1 <- optim.pml(
  fit_init1,
  model         = best_model_clean,
  optInv        = TRUE,       # proportion de sites invariants
  optGamma      = TRUE,       # distribution gamma
  rearrangement = "stochastic",
  control       = pml.control(trace = 1)
)

cat("[INFO] Log-vraisemblance :", fit_ml1$logLik, "\n")
fit_ml1$tree
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_VHCfoyerPQML.png",
    width = 3000,
    height = 2200,
    res = 300)

plot(fit_ml1$tree,
     tip.color = tip_colors2,
     cex       = 1,
     no.margin = TRUE)
legend("topright",
       legend = names(couleurs_menage),
       fill   = couleurs_menage,
       cex    = 0.7)

add.scale.bar(x = 0, y = max(nodeHeights(fit_ml1$tree)))
dev.off()
graphics.off()



############################################################
# Bootstrap 1000 réplicats
############################################################

cat("[INFO] Bootstrap 1000 réplicats — patience...\n")
set.seed(42)

bs_ml1 <- bootstrap.pml(
  fit_ml1,
  bs      = 1000,
  optNni  = TRUE,
  control = pml.control(trace = 0)
)

# Arbre avec valeurs bootstrap
tree_bs1 <- plotBS(fit_ml1$tree,
                  bs_ml1,
                  type = "none")   # récupérer sans plotter
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_VHCfoyerPQMLB.png",
    width = 3000,
    height = 2200,
    res = 300)

plot(tree_bs1,
     tip.color = tip_colors2,
     cex       = 1,
     no.margin = TRUE)
legend("topright",
       legend = names(couleurs_menage),
       fill   = couleurs_menage,
       cex    = 0.7)
add.scale.bar(x = 0, y = max(nodeHeights(tree_bs1)))
dev.off()
graphics.off()

############################################################
#  Charger l'alignement avec les controles
############################################################

# Alignement déjà fait avec MAFFT en bash :
# mafft --auto --reorder controle_primaryall.fasta >aligned_primarycontrole.fasta

alnPrim_cont <- read.dna("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/aligned_primarycontrole.fasta", format = "fasta")

### nom des sequences
rownames(alnPrim_cont)
### Tableau avec noms + tailles
data.frame(
  Sequence = rownames(alnPrim_cont),
  Length = apply(as.character(alnPrim_cont), 1, function(x) sum(x != "-"))
)

# supprimer les -NS5B_primary 
rownames(alnPrim_cont) <- gsub("-NS5B_primary", "", rownames(alnPrim_cont))

# nombre de sites
#n_sites <- ncol(align)

############################################################
# 2) Matrice de distances génétiques
# modèle TN93 = recommandé pour VHC
############################################################

dist_matrixprim_cont <- dist.dna(alnPrim_cont, model = "TN93", pairwise.deletion = TRUE)
dist_dfcont     <- as.matrix(dist_matrixprim_cont)



############################################################
# 4) Statistiques 
############################################################
dist_values <- dist_dfcont[upper.tri(dist_dfcont)]
mean(dist_values)
sd(dist_values)
median(dist_values)
min(dist_values)
max(dist_values)
summary(dist_values)
############################################################
# 5) Histogramme des distances
############################################################
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/histogram_controle.png",
    width = 3000,
    height = 2200,
    res = 300)
hist(dist_values,
     probability = TRUE,
     main = "Distribution des distances génétiques",
     xlab = "Distance TN93",
     col = "lightblue",
     breaks = 30)

lines(density(dist_values, adjust = 2),
      col = "blue",
      lwd = 2)

dev.off()
graphics.off()

############################################################
# 6) Phylogénie
############################################################
# Vérifier s'il reste des NA
sum(is.na(as.matrix(dist_matrixprim_cont)))


tree_controle <- njs(dist_matrixprim_cont)
tree_controle <- ladderize(tree_controle)


# Plot arbre


png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_VHCControle.png",
    width = 3000,
    height = 2200,
    res = 300)

plot(tree_controle,
     cex       = 1,
     no.margin = TRUE)

add.scale.bar()


dev.off()
graphics.off()


############################################################
# 7) Heatmap des distances colorée par ménage
############################################################
dist_df <- as.data.frame(as.matrix(dist_matrixprim_cont))
dist_df$seq1 <- rownames(dist_df)

dist_long <- melt(dist_df, id.vars = "seq1",
                  variable.name = "seq2",
                  value.name = "dist")

p <- ggplot(dist_long, aes(x = seq1, y = seq2, fill = dist)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    
    plot.title = element_text(
      size = 16,
      face = "bold"
    ),
    
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    
    axis.ticks = element_blank()
  ) +
  
  labs(
    title = "Heatmap des distances génétiques (TN93)",
    x = "",
    y = ""
  )

ggsave("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/heatmap_distances.png",
       plot = p,
       width = 8,
       height = 7,
       dpi = 300)


##############################################################################
# 2) phylogenie : arbre avec la vraisemblance et Sélection du meilleur modèle
#############################################################################

# Tester les modèles — peut prendre quelques minutes
cat("[INFO] Sélection du modèle...\n")
modtest2 <- modelTest(alnPrim_cont,
                      tree  = tree_controle,
                      model = c("JC", "HKY", "TN93", "GTR"),
                      G     = TRUE,
                      I     = TRUE)

# Afficher les résultats triés par AIC
modtest_sorted2 <- modtest2[order(modtest1$AIC), ]
print(head(modtest_sorted2, 10))

# Meilleur modèle
best_model <- modtest_sorted1$Model[1]
cat("[INFO] Meilleur modèle (AIC) :", best_model, "\n")

best_model_clean <- gsub("\\+G\\(4\\)|\\+I", "", best_model)
cat("[INFO] Meilleur modèle nettoyer (AIC) :", best_model_clean, "\n")
############################################################
# 3) Arbre Maximum de Vraisemblance
############################################################

cat("[INFO] Optimisation ML...\n")

# Initialiser avec le meilleur modèle
aln_phy <- phyDat(alnPrim_cont, type = "DNA")
fit_initC <- pml(tree_controle, aln_phy)

# Optimiser
fit_mlC <- optim.pml(
  fit_initC,
  model         = best_model_clean,
  optInv        = TRUE,       # proportion de sites invariants
  optGamma      = TRUE,       # distribution gamma
  rearrangement = "stochastic",
  control       = pml.control(trace = 1)
)

cat("[INFO] Log-vraisemblance :", fit_mlC$logLik, "\n")
fit_mlC$tree
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_VHCControleML.png",
    width = 3000,
    height = 2200,
    res = 300)

plot(fit_mlC$tree,
     cex       = 1,
     no.margin = TRUE)

add.scale.bar(x = 0, y = max(nodeHeights(fit_mlC$tree)))
dev.off()
graphics.off()

  
############################################################
# Bootstrap 1000 réplicats
############################################################

cat("[INFO] Bootstrap 1000 réplicats — patience...\n")
set.seed(42)

bs_ml1C <- bootstrap.pml(
  fit_mlC,
  bs      = 1000,
  optNni  = TRUE,
  control = pml.control(trace = 0)
)

# Arbre avec valeurs bootstrap
tree_bs1C <- plotBS(fit_mlC$tree,
                   bs_ml1C,
                   type = "none")   # récupérer sans plotter
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_VHCControleMLB.png",
    width = 3000,
    height = 2200,
    res = 300)

plot(tree_bs1C,
     cex       = 1,
     no.margin = TRUE)

add.scale.bar(x = 0, y = max(nodeHeights(tree_bs1C)))
dev.off()
graphics.off()

tree_bs1C$edge.length * ncol(alnPrim_cont)

############################################################
#  quasi-espèce controle 
############################################################
aln_combined <- read.dna("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/combined_controle_quasi_aligned.fasta", format = "fasta")
### nom des sequences
rownames(aln_combined)
### Tableau avec noms + tailles
data.frame(
  Sequence = rownames(aln_combined),
  Length = apply(as.character(aln_combined), 1, function(x) sum(x != "-"))
)
## Longueur de l’alignement total
ncol(aln_combined)

############################################################
#  Matrice de distances génétiques
# modèle TN93 = recommandé pour VHC
############################################################
# Distances et arbre
dist_combined <- dist.dna(aln_combined, model = "TN93",
                          pairwise.deletion = TRUE)

dist_df     <- as.matrix(dist_combined)

# Construire l'arbre (Neighbor-Joining)

tree_combined <- njs(dist_combined)
tree_combined <- ladderize(tree_combined)

# Extraire infos depuis les noms des tips
# ex: "M01_421233_primary" ou "CTRL_421233_primary"
# Labels des tips
tip_labels <- tree_combined$tip.label

# Extraire primary / secondary
tip_type <- sub(".*_", "", tip_labels)

# Forme des points
tip_pch <- ifelse(tip_type == "primary", 16, 17)

# Couleurs des points
tip_col <- ifelse(tip_type == "primary",
                  "steelblue",
                  "tomato")

# Export PNG
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_combineControle.png",
    width = 1400,
    height = 1000,
    res = 150)

# Plot arbre
plot(tree_combined,
     cex  = 0.6,
     main = "Phylogénie VHC Contrôle quasi-espèce")

# Ajouter symboles aux tips
tiplabels(pch = tip_pch,
          col = tip_col,
          cex = 0.8)

# Légende
legend("topright",
       legend = c("primary", "secondary"),
       pch    = c(16, 17),
       col    = c("steelblue", "tomato"),
       title  = "Type",
       cex    = 0.7,
       bty    = "n")


add.scale.bar()
dev.off()

##############################################################################
# 2) phylogenie : arbre avec la vraisemblance et Sélection du meilleur modèle
#############################################################################

# Tester les modèles — peut prendre quelques minutes
cat("[INFO] Sélection du modèle...\n")
modtestC <- modelTest(aln_combined,
                      tree  =  tree_combined,
                      model = c("JC", "HKY", "TN93", "GTR"),
                      G     = TRUE,
                      I     = TRUE)

# Afficher les résultats triés par AIC
modtest_sortedC <- modtestC[order(modtestC$AIC), ]
print(head(modtest_sortedC, 10))

# Meilleur modèle
best_model <- modtest_sortedC$Model[1]
cat("[INFO] Meilleur modèle (AIC) :", best_model, "\n")

best_model_clean <- gsub("\\+G\\(4\\)|\\+I", "", best_model)
cat("[INFO] Meilleur modèle nettoyer (AIC) :", best_model_clean, "\n")
############################################################
# 3) Arbre Maximum de Vraisemblance
############################################################

cat("[INFO] Optimisation ML...\n")

# Initialiser avec le meilleur modèle
aln_phy <- phyDat(aln_combined, type = "DNA")
fit_initC <- pml(tree_combined, aln_phy)

# Optimiser
fit_mlC <- optim.pml(
  fit_initC,
  model         = best_model_clean,
  optInv        = TRUE,       # proportion de sites invariants
  optGamma      = TRUE,       # distribution gamma
  rearrangement = "stochastic",
  control       = pml.control(trace = 1)
)

cat("[INFO] Log-vraisemblance :", fit_mlC$logLik, "\n")
fit_mlC$tree
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_VHCControlePQML.png",
    width = 3000,
    height = 2200,
    res = 300)

plot(fit_mlC$tree,
     cex = 0.8,
     no.margin = TRUE)
# Ajouter symboles aux tips
tiplabels(pch = tip_pch,
          col = tip_col,
          cex = 0.6)
# Légende
legend("topright",
       legend = c("primary", "secondary"),
       pch    = c(16, 17),
       col    = c("steelblue", "tomato"),
       title  = "Type",
       cex    = 0.7,
       bty    = "n")
xmax <- max(nodeHeights(fit_mlC$tree))
add.scale.bar(
  x = xmax * 0.85,
  y = 1
)
dev.off()
graphics.off()


############################################################
# Bootstrap 1000 réplicats
############################################################

cat("[INFO] Bootstrap 1000 réplicats — patience...\n")
set.seed(42)

bs_ml1C <- bootstrap.pml(
  fit_mlC,
  bs      = 1000,
  optNni  = TRUE,
  control = pml.control(trace = 0)
)

# Arbre avec valeurs bootstrap
tree_bs1C <- plotBS(fit_mlC$tree,
                    bs_ml1C,
                    type = "none")   # récupérer sans plotter
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_VHCPQMLB.png",
    width = 3000,
    height = 2200,
    res = 300)

plot(tree_bs1C,
     cex = 0.8,
     no.margin = TRUE)
# Ajouter symboles aux tips
tiplabels(pch = tip_pch,
          col = tip_col,
          cex = 0.6)
# Légende
legend("topright",
       legend = c("primary", "secondary"),
       pch    = c(16, 17),
       col    = c("steelblue", "tomato"),
       title  = "Type",
       cex    = 0.7,
       bty    = "n")
xmax <- max(nodeHeights(tree_bs1C))
add.scale.bar(
  x = xmax * 0.85,
  y = 1
)
dev.off()
graphics.off()

#####################################################################
############## Arbre P+Q Controle ##############################
###################################################################

#####################################################################
############## Arbre avec le controle ##############################
###################################################################

# Charger l'alignement
alnFC <- read.dna("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/all_foyer_controle_aligned.fasta", format = "fasta")

# Distances et arbre
dist_mat <- dist.dna(alnFC, model = "TN93", pairwise.deletion = TRUE)
treeFC     <- njs(dist_mat)
treeFC     <- ladderize(treeFC)

# Extraire le type (FOYER ou CTRL) depuis les noms des tips
tip_labels <- treeFC$tip.label
tip_type   <- ifelse(grepl("^FOYER_", tip_labels), "Foyer", "Controle")

# Extraire le ménage pour les foyers (depuis metadata)
view(foyer)
tip_menage <- rep(NA, length(tip_labels))
for (i in seq_along(tip_labels)) {
  if (tip_type[i] == "Foyer") {
    individu <- gsub("FOYER_", "", tip_labels[i])
    m <- foyer$menage [match(individu, foyer$ID)]
    tip_menage[i] <- ifelse(length(m) > 0 && !is.na(m), m, "INCONNU")
  }
}

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

# Forme des points :
# Controle → triangle (17)
# Foyer    → rond (16)
tip_pch <- ifelse(tip_type == "Foyer", 16, 17)

# Plot
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_foyer_controle.png", width = 1600, height = 1200, res = 150)

plot(treeFC,
     tip.color = tip_colors,
     cex       = 0.7,
     label.offset = 0.001,
     main      = "Phylogénie VHC — Foyers vs Contrôles")

# Ajouter symboles sur les tips
tiplabels(pch = tip_pch,
          col = tip_colors,
          cex = 0.8)

# Légende ménages
legend("topright",
       legend = c(names(couleurs_menages), "Controle"),
       fill   = c(couleurs_menages, "grey60"),
       title  = "Ménage",
       cex    = 0.65,
       bty    = "n")

# Légende symboles
legend("bottomright",
       legend = c("Foyer", "Controle"),
       pch    = c(16, 17),
       col    = c("black", "grey60"),
       title  = "Groupe",
       cex    = 0.65,
       bty    = "n")

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

best_model_clean <- gsub("\\+G\\(4\\)|\\+I", "", best_modelFC)
cat("[INFO] Meilleur modèle nettoyer (AIC) :", best_model_clean, "\n")
############################################################
# 3) Arbre Maximum de Vraisemblance
############################################################

cat("[INFO] Optimisation ML...\n")

# Initialiser avec le meilleur modèle
aln_phyFC <- phyDat(alnFC, type = "DNA")
fit_initFC <- pml(treeFC, aln_phyFC)

# Optimiser
fit_mlFC <- optim.pml(
  fit_initFC,
  model         = best_model_clean,
  optInv        = TRUE,       # proportion de sites invariants
  optGamma      = TRUE,       # distribution gamma
  rearrangement = "stochastic",
  control       = pml.control(trace = 1)
)

cat("[INFO] Log-vraisemblance :", fit_mlFC$logLik, "\n")

png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_VHCFoyerControlML.png",
    width = 3000,
    height = 2200,
    res = 300)

plot(fit_mlFC$tree,
     tip.color = tip_colors,
     cex       = 0.7,
     label.offset = 0.001,
     main      = "Phylogénie VHC — Foyers vs Contrôles")
# Ajouter symboles sur les tips
tiplabels(pch = tip_pch,
          col = tip_colors,
          cex = 0.8)

# Légende ménages
legend("topright",
       legend = c(names(couleurs_menages), "Controle"),
       fill   = c(couleurs_menages, "grey60"),
       title  = "Ménage",
       cex    = 0.65,
       bty    = "n")

xmax <- max(nodeHeights(fit_mlFC$tree))
add.scale.bar(
  x = xmax * 0.85,
  y = 1
)
dev.off()
graphics.off()


############################################################
# Bootstrap 1000 réplicats
############################################################

cat("[INFO] Bootstrap 1000 réplicats — patience...\n")
set.seed(42)

bs_mlFC <- bootstrap.pml(
  fit_mlFC,
  bs      = 1000,
  optNni  = TRUE,
  control = pml.control(trace = 0)
)

# Arbre avec valeurs bootstrap
tree_bs1FC <- plotBS(fit_mlFC$tree,
                    bs_mlFC,
                    type = "none")   # récupérer sans plotter
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_VHCFCMLB.png",
    width = 3000,
    height = 2200,
    res = 300)

plot(tree_bs1FC,
     tip.color = tip_colors,
     cex       = 0.7,
     label.offset = 0.001,
     main      = "Phylogénie VHC — Foyers vs Contrôles")
# Ajouter symboles sur les tips
tiplabels(pch = tip_pch,
          col = tip_colors,
          cex = 0.8)

# Légende ménages
legend("topright",
       legend = c(names(couleurs_menages), "Controle"),
       fill   = c(couleurs_menages, "grey60"),
       title  = "Ménage",
       cex    = 0.65,
       bty    = "n")
xmax <- max(nodeHeights(tree_bs1FC))
add.scale.bar(
  x = xmax * 0.85,
  y = 1
)
dev.off()
graphics.off()

########################################################################
### Mettre tout les quasi espece ensemble #############################
#######################################################################

# Charger l'alignement
aln_combined <- read.dna("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/combined_allquasialigned.fasta", format = "fasta")

############################################################
# 2) Matrice de distances génétiques
# modèle TN93 = recommandé pour VHC
############################################################
# Distances et arbre
dist_combined <- dist.dna(aln_combined, model = "TN93",
                          pairwise.deletion = TRUE)

dist_df     <- as.matrix(dist_combined)

# Construire l'arbre (Neighbor-Joining)

tree_combinedQ <- njs(dist_combined)
tree_combinedQ <- ladderize(tree_combinedQ)

# Extraire infos depuis les noms des tips
# ex: "M01_421233-NS5B_primary"
tip_labels  <- tree_combinedQ$tip.label
tip_menage  <- sub("_.*", "", tip_labels)            # M01

tip_type    <- sub(".*_", "", tip_labels)            # primary ou secondary

# Couleurs par ménage
couleurs_menage <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(unique(tip_menage))),
  unique(tip_menage)
)


tip_colors <- couleurs_menage[tip_menage]

# Forme des points selon primary/secondary
tip_pch <- ifelse(tip_type == "primary", 16, 17)   # rond vs triangle

# Plot
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_combinéquasiall.png", width = 1400, height = 1000, res = 150)
plot(tree_combinedQ,
     tip.color = tip_colors,
     cex       = 0.5)
# Ajouter des symboles pour distinguer primary/secondary
tiplabels(pch = tip_pch, col = tip_colors, cex = 0.8)
legend("topright",
       legend = c(names(couleurs_menage), "primary", "secondary"),
       fill   = c(couleurs_menage, NA, NA),
       pch    = c(rep(NA, length(couleurs_menage)), 16, 17),
       title  = "Ménage / Type",
       cex    = 0.7)
add.scale.bar()
dev.off()

##############################################################################
# phylogenie : arbre avec la vraisemblance et Sélection du meilleur modèle
#############################################################################

# Tester les modèles — peut prendre quelques minutes
cat("[INFO] Sélection du modèle...\n")
modtestPQ <- modelTest(aln_combined,
                       tree  =  tree_combinedQ,
                       model = c("JC", "HKY", "TN93", "GTR"),
                       G     = TRUE,
                       I     = TRUE)

# Afficher les résultats triés par AIC
modtest_sortedPQ <- modtestPQ[order(modtestPQ$AIC), ]
print(head(modtest_sortedPQ, 10))

# Meilleur modèle
best_modelPQ <- modtest_sortedPQ$Model[1]
cat("[INFO] Meilleur modèle (AIC) :", best_modelPQ, "\n")

best_model_clean <- gsub("\\+G\\(4\\)|\\+I", "", best_modelPQ)
cat("[INFO] Meilleur modèle nettoyer (AIC) :", best_model_clean, "\n")
############################################################
# 3) Arbre Maximum de Vraisemblance
############################################################

cat("[INFO] Optimisation ML...\n")

# Initialiser avec le meilleur modèle
aln_phyPQ <- phyDat(aln_combined, type = "DNA")
fit_initPQ <- pml(tree_combinedQ, aln_phyPQ)

# Optimiser
fit_mlPQ <- optim.pml(
  fit_initPQ,
  model         = best_model_clean,
  optInv        = TRUE,       # proportion de sites invariants
  optGamma      = TRUE,       # distribution gamma
  rearrangement = "stochastic",
  control       = pml.control(trace = 1)
)

cat("[INFO] Log-vraisemblance :", fit_mlPQ$logLik, "\n")

png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_VHCPQML.png",
    width = 3000,
    height = 2200,
    res = 300)

plot(fit_mlPQ$tree,
     tip.color = tip_colors,
     cex       = 0.5)
# Ajouter des symboles pour distinguer primary/secondary
tiplabels(pch = tip_pch, col = tip_colors, cex = 0.8)
legend("topright",
       legend = c(names(couleurs_menage), "primary", "secondary"),
       fill   = c(couleurs_menage, NA, NA),
       pch    = c(rep(NA, length(couleurs_menage)), 16, 17),
       title  = "Ménage / Type",
       x.intersp = 0.8,
       pt.cex = 1.2,
       bty = "n")

xmax <- max(nodeHeights(fit_mlPQ$tree))
add.scale.bar(
  x = xmax * 0.85,
  y = 1
)
dev.off()
graphics.off()


############################################################
# Bootstrap 1000 réplicats
############################################################

cat("[INFO] Bootstrap 1000 réplicats — patience...\n")
set.seed(42)

bs_mlFC <- bootstrap.pml(
  fit_mlFC,
  bs      = 1000,
  optNni  = TRUE,
  control = pml.control(trace = 0)
)

# Arbre avec valeurs bootstrap
tree_bs1FC <- plotBS(fit_mlFC$tree,
                     bs_mlFC,
                     type = "none")   # récupérer sans plotter
png("~/MASTER BIOINFORMATIQUE/stageM1/Projets/Resultat_R/graphes/phylogenie_VHCPQMLB.png",
    width = 3000,
    height = 2200,
    res = 300)

plot(tree_bs1FC,
     tip.color = tip_colors,
     cex       = 1,
     label.offset = 0.001)
# Ajouter symboles sur les tips
tiplabels(pch = tip_pch,
          col = tip_colors,
          cex = 0.6)

# Légende ménages
legend("topright",
       legend = c(names(couleurs_menages), "Controle"),
       fill   = c(couleurs_menages, "grey60"),
       title  = "Ménage",
       cex    = 0.65,
       bty    = "n")
xmax <- max(nodeHeights(tree_bs1FC))
add.scale.bar(
  x = xmax * 0.85,
  y = 1
)
dev.off()
graphics.off()



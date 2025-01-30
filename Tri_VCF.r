# Charger les bibliothèques
library(vcfR)

# Lire le fichier VCF
vcf <- read.vcfR("~/Téléchargements/P15-1.trimed1000.sv_sniffles.vcf")

# Extraire les colonnes fixes
donnees_fixes <- as.data.frame(vcf@fix)

# Fonction pour extraire les champs de la colonne INFO d'un fichier VCF ou d'un format similaire
extraire_info <- function(info_col) {
  # Sépare chaque ligne de la colonne INFO en fonction du séparateur ";"
  infos <- strsplit(info_col, ";")
  # Identifie les noms de champs uniques dans toutes les informations de la colonne INFO
  champs_possibles <- unique(unlist(lapply(infos, function(x) {
    # Pour chaque chaîne d'information, extrait les noms des champs (partie avant le "=")
    sapply(strsplit(x, "="), function(y) y[1])
  })))
  # Crée une liste où chaque élément représente une ligne de la colonne INFO, avec les valeurs associées aux champs
  infos_list <- lapply(infos, function(x) {
    # Sépare chaque information en paires clé-valeur sur le sép = 
    pairs <- strsplit(x, "=")
    # Extrait les valeurs des paires, en attribuant NA aux clés sans valeur (si une clé n'a pas de valeur après le "=")
    values <- sapply(pairs, function(y) ifelse(length(y) > 1, y[2], y[1]))
    # Associe les noms des champs aux valeurs : 
    names(values) <- sapply(pairs, function(y) y[1])
    # Crée un vecteur avec tous les champs possibles initialisés à NA :
    complet <- setNames(rep(NA, length(champs_possibles)), champs_possibles)
    # Remplir les champs correspondant avec les valeurs extraites
    complet[names(values)] <- values
    # Retourne le vecteur complet pour cette ligne
    return(complet)
  })
  
  # Combine tous les éléments de la liste en un seul DataFrame, chaque ligne correspondant à une ligne INFO traitée
  infos_df <- do.call(rbind, infos_list)
  
  # Convertit la liste en un DataFrame avec des chaînes de caractères (évite la conversion implicite en facteurs)
  infos_df <- as.data.frame(infos_df, stringsAsFactors = FALSE)
  
  # Retourne le DataFrame final avec les informations extraites
  return(infos_df)
}


# Appliquer la fonction sur la colonne INFO
info_colonnes <- extraire_info(donnees_fixes$INFO)

# Fusionner les colonnes fixes et celles extraites
donnees_finales <- cbind(
  donnees_fixes[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")],
  info_colonnes
)
donnees_finales[is.na(donnees_finales)] <- "" # Remplacer NA par des chaînes vides
##########################################################################################
# MONTAGE DE LA µ profondeur et ajout d'une colonne cryptique pour manipuler pr chaque sv
##########################################################################################
# Calculer la moyenne du COVERAGE
calcMoyCoverage <- function(coverage) {
  if (!is.na(coverage) && nzchar(coverage)) { # Vérifie si la chaîne n'est ni vide ni NA
    val_prof <- as.numeric(unlist(strsplit(coverage, ",")))
    if (length(val_prof) > 0) {
      return(mean(val_prof, na.rm = TRUE))
    }
  }
  return(NA)
}

# Ajouter la colonne MU_COVERAGE
info_colonnes$MU_COVERAGE <- sapply(info_colonnes$COVERAGE, calcMoyCoverage)


##########################################################################################
# COSMETIQUE POUR QUE LE CSV SOIT CORRECTEMENT ORGANISEE ET LISIBLE 
##########################################################################################

# Réorganiser les colonnes pour que MU_COVERAGE soit à côté de COVERAGE
colonnes_reorganisees <- c(
  setdiff(names(info_colonnes), c("COVERAGE", "MU_COVERAGE")), # Autres colonnes sauf COVERAGE et MU_COVERAGE
  "COVERAGE", 
  "MU_COVERAGE"
)
info_colonnes <- info_colonnes[, colonnes_reorganisees]

# Fusionner les colonnes fixes et celles extraites
donnees_finales <- cbind(
  donnees_fixes[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")],
  info_colonnes
)


###################################################
# TODO : FONCTION POUR RENDRE INTERACTIF LE TRUC
###################################################
# Demander le nom du fichier de sortie
#fichier <- readline(prompt = "Veuillez saisir un nom de fichier CSV valide (avec l'extension .csv) : ")
#if (tools::file_ext(fichier) != "csv") {
#  fichier <- paste0(fichier, ".csv")
#}


##########################################
# Export dans un CSV PRELABLEMENT (TOUCH) 
###########################################
# Exporter les données finales
write.csv(donnees_finales, "/home/mickael/Projets_GIT/BILL_2025/Traitement_P90_VCF.csv", row.names = FALSE)
#cat("Fichier  :", fichier, "\n")


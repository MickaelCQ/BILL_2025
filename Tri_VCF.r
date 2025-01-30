library(vcfR)  # Bibliothèque pour lire et traiter les fichiers VCF

# Fonction pour lire et extraire les données du VCF
extraire_vcf <- function(vcf_file) {
  vcf <- vcfR::read.vcfR(vcf_file)  # Lire le fichier VCF
  
  # Extraire les colonnes fixes du fichier VCF (comme CHROM, POS, etc.)
  donnees_fixes <- as.data.frame(vcf@fix)  # Extraction des informations conventionnelles du VCF sous forme de data.frame
  
  ##############################################################
  # MISE A L'ECHELLE DES METADONNEES DU CHAMPS INFO 
  #############################################################
  
  extraire_info <- function(info_col) {
    # Séparer les informations de chaque ligne de la colonne INFO en fonction du séparateur ";"
    infos <- strsplit(info_col, ";")  # strsplit() est une fonction de base de R, qui divise les chaînes
    
    # Identifier les noms de champs uniques dans toutes les informations de la colonne INFO
    champs_possibles <- unique(unlist(lapply(infos, function(x) {
      # Extrait les noms des champs (avant le "=") de chaque chaîne d'information
      sapply(strsplit(x, "="), function(y) y[1])
    })))
    
    # Crée une liste où chaque élément représente une ligne de la colonne INFO avec les valeurs associées aux champs
    infos_list <- lapply(infos, function(x) {
      # Séparer chaque information en paires clé-valeur
      pairs <- strsplit(x, "=")
      values <- sapply(pairs, function(y) ifelse(length(y) > 1, y[2], y[1]))  # Extrait la valeur ou la clé
      names(values) <- sapply(pairs, function(y) y[1])  # Associe les noms aux valeurs
      complet <- setNames(rep(NA, length(champs_possibles)), champs_possibles)  # Initialise un vecteur avec des NA
      complet[names(values)] <- values  # Remplir les champs avec les valeurs extraites
      return(complet)
    })
    
    # Combine tous les éléments de la liste en un seul DataFrame
    infos_df <- do.call(rbind, infos_list)  # Combine les listes en un data.frame
    infos_df <- as.data.frame(infos_df, stringsAsFactors = FALSE)  # Assure que les colonnes sont traitées comme des chaînes
    
    # Retourner le DataFrame final avec les informations extraites
    return(infos_df)
  }
  
  # Appliquer la fonction sur la colonne INFO pour extraire les données pertinentes
  info_colonnes <- extraire_info(donnees_fixes$INFO)
  
  ##################################################################################
  #    AJOUT D'UNE COLONNE CRYPTIQUE POUR FAIRE UNE µ DES PROFONDEURS DE REPLICATS 
  #################################################################################
  # Calculer la moyenne du COVERAGE pour chaque variant
  calcMoyCoverage <- function(coverage) {
    if (!is.na(coverage) && nzchar(coverage)) {  # Vérifie que la chaîne de coverage n'est pas vide ni NA
      val_prof <- as.numeric(unlist(strsplit(coverage, ",")))  # Convertit les valeurs en numérique
      if (length(val_prof) > 0) {
        return(mean(val_prof, na.rm = TRUE))  # Calcul de la moyenne, en ignorant les NA
      }
    }
    return(NA)  # Retourne NA si coverage est invalide
  }
  
  # Ajouter la colonne MU_COVERAGE au DataFrame
  if ("COVERAGE" %in% names(info_colonnes)) {  # Vérifie si la colonne COVERAGE existe dans info_colonnes
    info_colonnes$MU_COVERAGE <- sapply(info_colonnes$COVERAGE, calcMoyCoverage)  # Applique la fonction à chaque valeur de COVERAGE
  } else {
    # Si COVERAGE n'existe pas, créer une colonne MU_COVERAGE vide
    info_colonnes$MU_COVERAGE <- NA
  }
  
  # Réorganiser les colonnes pour que MU_COVERAGE soit à côté de COVERAGE
  colonnes_reorganisees <- c(
    setdiff(names(info_colonnes), c("COVERAGE", "MU_COVERAGE")),  # Inclure toutes les autres colonnes sauf COVERAGE et MU_COVERAGE
    "COVERAGE", 
    "MU_COVERAGE"
  )
  info_colonnes <- info_colonnes[, colonnes_reorganisees]  # Réarranger les colonnes
  
  ##################################################################################
  #                               COSMETIQUE DU CSV    
  ##################################################################################
  # Fusionner les colonnes fixes et celles extraites
  donnees_finales <- cbind(
    donnees_fixes[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")],  # Colonnes fixes du VCF
    info_colonnes)  # Colonnes extraites de la colonne INFO
  
  return(donnees_finales)
}

##################################################################################
#                   TRAITEMENT DE PLUSIEURS FICHIERS VCF A LA FOIS 
##################################################################################
# Liste des fichiers VCF dans le répertoire du projet
#vcf_files <- list.files(path = "/home/mickael/Projets_GIT/BILL_2025/", pattern = "*.vcf", full.names = TRUE)
vcf_files <- list.files(path = "~/Téléchargements/", pattern = "*.vcf", full.names = TRUE) # Tests
toutes_les_donnees <- data.frame()

for (vcf_file in vcf_files) {
  # Traiter chaque fichier VCF et récupérer les données
  donnees_vcf <- extraire_vcf(vcf_file)
  toutes_les_donnees <- rbind(toutes_les_donnees, donnees_vcf)
}

# Exporter les données dans un fichier CSV
#write.csv(toutes_les_donnees, "~/Téléchargements/Test.csv", row.names = FALSE)
write.csv(toutes_les_donnees, "/home/mickael/Projets_GIT/BILL_2025/Traitement_P90_VCF.csv", row.names = FALSE)
          

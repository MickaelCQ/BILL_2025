#01/29/25 Structuring and assembly of VCF files listing structural variants of interest
#01/30/25 Programming a function for read and extract columns natives of the VCFs and generate a csv
#01/30/25 Reflection and programming a function for the addition of INFO column metadata
#01/30/25 GIT update, comment translation and test set on generation P15, P30, P50 and P65
#01/31/25 Add a cryptic column on first index with the willingness to identify for each SV the viral generation

library(vcfR) # no std library R for import function for the VCF treatment

################################# 
# Generic Functions
##################################

# For read and extract data in the VCF file
FunExtractVCF <- function(vcf_file) {
  vcf <- vcfR::read.vcfR(vcf_file)  # Read this file
  
  # Extraire les colonnes fixes du fichier VCF (comme CHROM, POS, etc.)
  nativesCol <- as.data.frame(vcf@fix)  # Extract 
  
  ################################# 
  # "INFO" FIELD METADATA SCALING
  ##################################
  
  FunExtractINFO <- function(info_col) {
    # Séparer les informations de chaque ligne de la colonne INFO en fonction du séparateur ";"
    infos <- strsplit(info_col, ";")
    
    # Identifier les noms de champs uniques dans toutes les informations de la colonne INFO
    champs_possibles <- unique(unlist(lapply(infos, function(x) {
      sapply(strsplit(x, "="), function(y) y[1])
    })))
    
    # Crée une liste où chaque élément représente une ligne de la colonne INFO avec les valeurs associées aux champs
    infos_list <- lapply(infos, function(x) {
      pairs <- strsplit(x, "=")
      values <- sapply(pairs, function(y) ifelse(length(y) > 1, y[2], y[1]))  # Extrait la valeur ou la clé
      names(values) <- sapply(pairs, function(y) y[1])  # Associe les noms aux valeurs
      
      # Remplir les champs manquants avec NA
      complet <- setNames(rep(NA, length(champs_possibles)), champs_possibles)
      complet[names(values)] <- values
      return(complet)
    })
    
    # Combine tous les éléments de la liste en un seul DataFrame
    infos_df <- do.call(rbind, infos_list)
    infos_df <- as.data.frame(infos_df, stringsAsFactors = FALSE)  # Assure que les colonnes sont traitées comme des chaînes
    
    return(infos_df)
  }
  
  # Appliquer la fonction sur la colonne INFO pour extraire les données pertinentes
  info_colonnes <- FunExtractINFO(nativesCol$INFO)
  
  ##################################################################################
  # AJOUT D'UNE COLONNE POUR LA MOYENNE DES PROFONDEURS DE RÉPLICATS (COVERAGE)
  ##################################################################################
  calcMoyCoverage <- function(coverage) {
    if (!is.na(coverage) && nzchar(coverage)) {
      val_prof <- as.numeric(unlist(strsplit(coverage, ",")))
      if (length(val_prof) > 0) {
        return(mean(val_prof, na.rm = TRUE))
      }
    }
    return(NA)
  }
  
  if ("COVERAGE" %in% names(info_colonnes)) {
    info_colonnes$MU_COVERAGE <- sapply(info_colonnes$COVERAGE, calcMoyCoverage)
  } else {
    info_colonnes$MU_COVERAGE <- NA
  }
  
  colonnes_reorganisees <- c(
    setdiff(names(info_colonnes), c("COVERAGE", "MU_COVERAGE")),
    "COVERAGE",
    "MU_COVERAGE"
  )
  info_colonnes <- info_colonnes[, colonnes_reorganisees, drop = FALSE]  # Évite les erreurs si colonnes manquent
  
  ##################################################################################
  # COSMÉTIQUE DU CSV
  ##################################################################################
  donnees_finales <- cbind(
    nativesCol[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")],
    info_colonnes
  )
  
  # Extraire la génération du nom du fichier
  generation <- gsub("^P([0-9]+)-.*$", "P\\1", basename(vcf_file))
  donnees_finales$GENERATION <- generation
  
  return(donnees_finales)
}

##################################################################################
# TRAITEMENT DE PLUSIEURS FICHIERS VCF AVEC CHAMPS VARIABLES DANS INFO
##################################################################################
vcf_files <- list.files(path = "/home/mickael/Projets_GIT/BILL_2025/VCF/", pattern = "*.vcf", full.names = TRUE)
toutes_les_donnees <- NULL  # Initialiser comme NULL pour éviter les conflits

for (vcf_file in vcf_files) {
  donnees_vcf <- FunExtractVCF(vcf_file)  # Extraire les données du fichier VCF
  
  if (is.null(toutes_les_donnees)) {
    # Si c'est le premier fichier traité, initialiser avec ses données
    toutes_les_donnees <- donnees_vcf
  } else {
    # Trouver toutes les colonnes existantes entre les données accumulées et les nouvelles données
    colonnes_toutes <- unique(c(names(toutes_les_donnees), names(donnees_vcf)))
    
    # Ajouter des colonnes manquantes à toutes_les_donnees
    for (col in setdiff(colonnes_toutes, names(toutes_les_donnees))) {
      toutes_les_donnees[[col]] <- NA
    }
    
    # Ajouter des colonnes manquantes à donnees_vcf
    for (col in setdiff(colonnes_toutes, names(donnees_vcf))) {
      donnees_vcf[[col]] <- NA
    }
    
    # Réordonner les colonnes et fusionner
    toutes_les_donnees <- rbind(
      toutes_les_donnees[, colonnes_toutes, drop = FALSE],
      donnees_vcf[, colonnes_toutes, drop = FALSE]
    )
  }
}

# Exporter les données dans un fichier CSV
write.csv(toutes_les_donnees, "/home/mickael/Projets_GIT/BILL_2025/Traitement_VCF.csv", row.names = FALSE)

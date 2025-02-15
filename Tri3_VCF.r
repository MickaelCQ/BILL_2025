# Charger la bibliothèque vcfR pour la gestion des fichiers VCF
library(vcfR)
# Définition de la fonction FunExtractVCF qui prend en entrée un fichier VCF
FunExtractVCF <- function(vcf_file) {
  # Lecture du fichier VCF avec la fonction read.vcfR
  vcf <- vcfR::read.vcfR(vcf_file)
  # Définir les colonnes attendues dans le fichier VCF (y compris la colonne INFO)
  expected_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  # Extraire les colonnes fixes du VCF (hors colonne INFO) sous forme de data.frame
  nativesCol <- as.data.frame(vcf@fix)
  # Vérification et ajout des colonnes manquantes avec NA
  for (col in expected_cols) {
    if (!(col %in% names(nativesCol))) {
      nativesCol[[col]] <- NA  # Si la colonne est absente, on l'ajoute avec des valeurs NA
    }
  }
  # Définition de la fonction FunExtractINFO pour extraire et structurer les informations de la colonne INFO
  FunExtractINFO <- function(info_col) {
    infos <- strsplit(info_col, ";")  # Séparer les informations de la colonne INFO par les ";"
    # Identifier tous les champs uniques dans les entrées de la colonne INFO
    allFields <- unique(unlist(lapply(infos, function(x) sapply(strsplit(x, "="), function(y) y[1]))))
    # Créer une liste d'éléments pour chaque ligne de la colonne INFO
    infos_list <- lapply(infos, function(x) {
      pairs <- strsplit(x, "=")  # Séparer les paires clé=valeur
      values <- sapply(pairs, function(y) ifelse(length(y) > 1, y[2], NA))  # Extraire les valeurs
      names(values) <- sapply(pairs, function(y) y[1])  # Donner un nom à chaque valeur (clé)
      # Créer un vecteur avec des NA pour chaque champ, puis remplir avec les valeurs extraites
      full <- setNames(rep(NA, length(allFields)), allFields)
      full[names(values)] <- values  # Remplir les champs avec les valeurs correspondantes
      return(full)
    })
    # Combiner toutes les entrées traitées dans un data.frame
    infos_df <- as.data.frame(do.call(rbind, infos_list), stringsAsFactors = FALSE)
    return(infos_df)  # Retourner le data.frame contenant les informations structurées
  }
  # Appliquer la fonction d'extraction des informations de la colonne INFO
  infoColumns <- FunExtractINFO(nativesCol$INFO)
  # Vérification et ajout de NA aux colonnes manquantes dans infoColumns (incluant "COVERAGE")
  allFields <- unique(c("COVERAGE", colnames(infoColumns)))  # Inclure "COVERAGE" même s'il est absent
  for (col in allFields) {
    if (!(col %in% names(infoColumns))) {
      infoColumns[[col]] <- NA  # Si une colonne manque, ajouter NA
    }
  }
  # Fonction pour calculer la couverture moyenne
  calcMoyCoverage <- function(coverage) {
    if (!is.na(coverage) && nzchar(coverage)) {
      val_prof <- as.numeric(unlist(strsplit(coverage, ",")))  # Convertir les valeurs de couverture en numérique
      if (length(val_prof) > 0) {
        return(mean(val_prof, na.rm = TRUE))  # Calculer la moyenne des valeurs
      }
    }
    return(NA)  # Retourner NA si la couverture est absente ou invalide
  }
  # Ajouter la colonne "MU_COVERAGE" calculée à partir de la couverture
  infoColumns$MU_COVERAGE <- if ("COVERAGE" %in% names(infoColumns)) {
    sapply(infoColumns$COVERAGE, calcMoyCoverage)  # Appliquer la fonction calcMoyCoverage sur chaque ligne
  } else {
    NA  # Si "COVERAGE" n'existe pas, ajouter NA pour "MU_COVERAGE"
  }
  # Réorganiser les colonnes pour mettre "COVERAGE" et "MU_COVERAGE" à la fin
  colonnes_reorganisees <- c(setdiff(names(infoColumns), c("COVERAGE", "MU_COVERAGE")), "COVERAGE", "MU_COVERAGE")
  infoColumns <- infoColumns[, colonnes_reorganisees, drop = FALSE]  # Appliquer la nouvelle organisation des colonnes
  # Fusionner les données natives (chromosome, position, etc.) avec les données extraites de INFO
  finalData <- cbind(nativesCol[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")], infoColumns)
  # Extraction du numéro de semaine et du groupe à partir du nom du fichier
  finalData$WEEK <- gsub("^P([0-9]+)-.*$", "P\\1", basename(vcf_file))  # Extraire la semaine
  print(finalData$WEEK)  # Afficher le numéro de semaine
  finalData$GROUP <- sub("^([^\\.]+\\.[^\\.]+).*", "\\1", basename(vcf_file))  # Extraire le groupe
  # Retourner le data.frame final contenant toutes les informations extraites
  return(finalData)}
# Traitement de tous les fichiers VCF dans le répertoire spécifié
vcf_files <- list.files(path = "/home/mickael/Projets_GIT/BILL_2025/VCF", pattern = "*.vcf", full.names = TRUE)
print(vcf_files)  # Afficher la liste des fichiers VCF trouvés
allData <- NULL  # Initialiser l'objet qui va contenir toutes les données traitées
# Boucle pour traiter chaque fichier VCF un par un
for (vcf_file in vcf_files) {
  print(paste("Processing:", vcf_file))  # Afficher quel fichier est en cours de traitement
  vcfDATA <- FunExtractVCF(vcf_file)  # Appliquer la fonction FunExtractVCF pour traiter le fichier
  print(dim(vcfDATA))  # Afficher les dimensions (nombre de lignes et colonnes) des données extraites
  # Si c'est le premier fichier traité, initialiser allData avec les données extraites
  if (is.null(allData)) {
    allData <- vcfDATA
  } else {
    print(names(vcfDATA))  # Afficher les colonnes du fichier actuel
    print(names(allData))  # Afficher les colonnes des données déjà accumulées
    # Identifier les colonnes uniques entre allData et vcfDATA
    allColumns <- unique(c(names(allData), names(vcfDATA)))
    # Ajouter les colonnes manquantes dans allData avec des valeurs NA
    for (col in setdiff(allColumns, names(allData))) {
      allData[[col]] <- NA
    }
    # Ajouter les colonnes manquantes dans vcfDATA avec des valeurs NA
    for (col in setdiff(allColumns, names(vcfDATA))) {
      vcfDATA[[col]] <- NA
    }
    # Empiler les nouvelles données sur les anciennes
    allData <- rbind(allData[, allColumns, drop = FALSE], vcfDATA[, allColumns, drop = FALSE])
  }
}
# Correction de l'ordre des colonnes avant l'exportation finale
colonnes_finales <- c("WEEK", "GROUP", setdiff(names(allData), c("WEEK", "GROUP")))
allData <- allData[, colonnes_finales, drop = FALSE]  # Réorganiser les colonnes
# Exporter les données résultantes dans un fichier CSV
fileName <- paste0("/home/mickael/Projets_GIT/BILL_2025/Traitement_VCF_", Sys.Date(), ".csv")
write.csv(allData, fileName, row.names = FALSE)  # Enregistrer les données sous forme de fichier CSV sans les noms de lignes

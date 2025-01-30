# Charger le package
library(vcfR)

# Lire le fichier VCF
vcf <- read.vcfR("~/Téléchargements/P15-1.trimed1000.sv_sniffles.vcf")


# Extraire les colonnes fixes
donnees_fixes <- as.data.frame(vcf@fix)

# Séparer les champs de la colonne INFO
extraire_info <- function(info_col) {
  
  # Stratégie pour gérer les sous élements du champs INFO dans le VCF tabulé (passage en data.frame)
  infos <- strsplit(info_col, ";")
  
  # On va parcourir tous nos items du champs info, on extrait les clés (avant le =), et on les séparent de leurs valeurs si 
  # elles sont éxitentes, éliminaation des doubles avec unique. 
  
  champs_possibles <- unique(unlist(lapply(infos,function(x){
    sapply(strsplit(x,"="), function(y) y[1])
  })))
  
  infos_list <- lapply(infos, function(x) {
    
    pairs <- strsplit(x, "=")
    values <- sapply(pairs, function(y) ifelse(length(y) > 1, y[2], y[1]))
    names(values) <- sapply(pairs, function(y) y[1])
  
    complet <- setNames(rep(NA,length(champs_possibles)),champs_possibles)
    complet[names(values)] <- values
    return(complet)
  })
  
  # Convertir en data.frame
  infos_df <- do.call(rbind, infos_list)
  infos_df <- as.data.frame(infos_df, stringsAsFactors = FALSE)
  
  return(infos_df)
}

# Appliquer la fonction sur la colonne INFO 
info_colonnes <- extraire_info(donnees_fixes$INFO)

# Fusionner les colonnes fixes et les colonnes extraites, et mettre une chaine vide si pour un variant pas de métadonnées: 
donnees_finales <- cbind(donnees_fixes[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")], info_colonnes)
donnees_finales[is.na(donnees_finales)] <- ""

# Stockage dans notre CSV avec une vérification postérieur de l'extension (on ne sait jamais):
fichier <- readline(prompt = "Veuillez saisir un nom de CSV valide (avec lextension .csv)");

if (tools::file_ext(fichier) != "csv") {
  fichier <- paste0(fichier, ".csv")
} # Fin verification extension

##############################################################################
# Gestion des colonnes cryptiques dans le CSV généré
##############################################################################
calcMoyCoverage <- function(coverage)
mu_coverage <- "";




##############################################################################

# Exporter vers un fichier CSV
write.csv(donnees_finales, "~/Téléchargements/Test.csv", row.names = FALSE)

print("Exportation en .csv OK ")


# Charger les fichiers
variants_df <- read.csv("/home/mickael/Projets_GIT/BILL_2025/BRUTE_Traitement_VCF_2025.csv", sep="\t", header=TRUE)
regions_df <- read.csv("/home/mickael/Projets_GIT/BILL_2025/Carto_virus.csv", sep="\t", header=TRUE)

# Créer une fonction pour identifier les informations cryptiques
get_cryptic_info <- function(variant_pos) {
  cryptic_info <- c()  # Liste vide pour stocker les informations
  
  # Vérifier si la position du variant tombe dans une des fenêtres des régions
  for (i in 1:nrow(regions_df)) {
    start <- regions_df$Start[i]
    end <- regions_df$End[i]
    
    # Si la position du variant est dans la fenêtre (Start <= POS <= End)
    if (start <= variant_pos && variant_pos <= end) {
      feature <- regions_df$Feature[i]
      attribute <- regions_df$Attribute[i]
      value <- regions_df$Value[i]
      
      # Ajouter l'information au cryptic_info
      cryptic_info <- c(cryptic_info, paste(feature, value, sep=":"))
    }
  }
  
  # Retourner les informations cryptiques ou "NA" si rien n'a été trouvé
  if (length(cryptic_info) > 0) {
    return(paste(cryptic_info, collapse = ", "))
  } else {
    return("NA")
  }
}

# Appliquer la fonction sur chaque variant pour remplir la colonne cryptique
variants_df$cryptique <- sapply(variants_df$POS, get_cryptic_info)

# Sauvegarder le fichier mis à jour avec la nouvelle colonne
write.csv(variants_df, "/home/mickael/Projets_GIT/BILL_2025/BRUTE_Traitement_VCF_2025.csv", sep="\t", row.names=FALSE)


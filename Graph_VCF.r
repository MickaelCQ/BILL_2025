library(ggplot2) # Graphiques

# Fonction pour charger le fichier CSV et générer des graphiques
Graph <- function(fichier_csv) {
  # Charger les données depuis le fichier CSV
  donnees <- read.csv(fichier_csv)
  
  # Créer un histogramme pour visualiser la distribution de la couverture moyenne
  ggplot2::ggplot(donnees, aes(x = MU_COVERAGE)) +
    geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
    labs(title = "Distribution des profondeurs moyenne", x = "Moyenne de la profondeur", y = "Occurences") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "darkblue"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(color = "lightgray", size = 0.5),
      panel.grid.minor = element_line(color = "gray", size = 0.25)
    )
 # ggplot2::ggplot(donnees, aes(x = CHROM, y = QUAL)) +
    # geom_boxplot() +
     #labs(title = "Boxplot de QUAL par CHROM", x = "Chromosome", y = "Qualité") +
     #theme_minimal()
}

Graph("/home/mickael/Projets_GIT/BILL_2025/Traitement_P90_VCF.csv")


library(ggplot2)
library(dplyr)
library(plotly)  # Pour la visualisation 3D interactive

# Lire le fichier CSV (ajustez le chemin vers votre fichier)
data <- Traitement_VCF_2025.02.08  # Remplacer par le chemin réel de votre fichier

# Assurez-vous que les colonnes "MU_COVERAGE" et "AF" sont numériques
data$MU_COVERAGE <- as.numeric(data$MU_COVERAGE)
data$AF <- as.numeric(data$AF)
data$FILTER <- as.factor(data$FILTER)
data$ID <- as.factor(data$ID)
data$WEEK <- as.factor(data$WEEK)

# Appliquer les filtres sur les données
filtered_data <- data %>%
  filter(MU_COVERAGE > 0, AF > 0)

# Créer un graphique en 3D avec une échelle de couleur
plot_3D <- plot_ly(filtered_data, 
                   x = ~MU_COVERAGE, 
                   y = ~AF, 
                   z = ~as.numeric(FILTER == "PASS"),  
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = 5, 
                                 color = ~AF, 
                                 colorscale = "Viridis", 
                                 opacity = 0.7,
                                 colorbar = list(title = "Fréquence Allélique (AF)")  # Ajout de l'échelle
                   ),
                   text = ~paste("<br>PROFONDEUR : ", MU_COVERAGE, 
                                 "<br>AF: ", AF, 
                                 "<br>FILTER: ", FILTER,
                                 "<br>VARIANT: ", ID,
                                 "<br>GENERATION:", WEEK),  
                   hoverinfo = "text") %>%
  layout(title = paste("Critère : FILTER = PASS"),
         scene = list(
           xaxis = list(title = "Profondeur", type = "log"),  
           yaxis = list(title = "Fréquence allélique"),
           zaxis = list(title = "Filtre = PASS")
         ))

# Afficher le graphique interactif
plot_3D


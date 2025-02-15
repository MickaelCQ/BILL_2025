library(ggplot2)
library(dplyr)
library(plotly)  

# Import de mon fichier
data <- Traitement_VCF_2025.02.08  # Remplacer par le chemin réel de votre fichier

# Vérifier la qualité numérique de mes champs d'intéret et passer en facteur les chaines.
data$MU_COVERAGE <- as.numeric(data$MU_COVERAGE)
data$AF <- as.numeric(data$AF)
data$FILTER <- as.factor(data$FILTER)
data$ID <- as.factor(data$ID)
data$WEEK <- as.factor(data$WEEK)

# Appliquer les filtres sur les données , variables pour générer les différents graphes intéractifs
filtered_data <- data %>%
  filter(MU_COVERAGE > 100, AF > 0.25, FILTER == "PASS" )

# Ajout du graph avec l'échelle de couleur, 
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

# Afficher le graphique interactif pour les filtres courants
plot_3D


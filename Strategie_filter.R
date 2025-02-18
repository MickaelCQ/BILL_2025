library(ggplot2)
library(dplyr)
library(plotly) 

data <- Traitement_VCF_2025.02.08  
data$MU_COVERAGE <- as.numeric(data$MU_COVERAGE)
data$AF <- as.numeric(data$AF)
data$FILTER <- as.factor(data$FILTER)
data$ID <- as.factor(data$ID)
data$WEEK <- as.factor(data$WEEK)
data$SUPPORT <- as.numeric(data$SUPPORT)

filtered_data <- data %>%
  filter(MU_COVERAGE > 200, SUPPORT > 200, AF > 0.25, FILTER == "PASS")

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
                                 colorbar = list(title = "AF > 25")  # Ajout de l'échelle
                   ),
                   text = ~paste("<br>PROFONDEUR : ", MU_COVERAGE, 
                                 "<br>AF: ", AF, 
                                 "<br>FILTER: ", FILTER,
                                 "<br>VARIANT: ", ID,
                                 "<br>GENERATION:", WEEK),  
                   hoverinfo = "text") %>%
  layout(title = paste("Critère : FILTER = PASS, Nombre de lignes : ", nrow(filtered_data)),
         scene = list(
           xaxis = list(title = "Profondeur > 200", type = "log"),  
           yaxis = list(title = "Fréquence allélique > 0.25"),
           zaxis = list(title = "Filtre = 'PASS' ")
         ))

plot_3D

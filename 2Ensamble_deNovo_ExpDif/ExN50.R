######################################################
## Tema: RNA-seq: Grafica ExN50                     ##
## Autor: Olga Andrea Hernandez Miranda, Miranda H  ##
## Fecha: 27/08/2021                                ##
## Nota: Sirve para elaborar graficas de lineas     ##
## para los datos de ExN50                          ##
######################################################

# Librerias:
library(ggplot2)
library("RColorBrewer")

directorio <- "C:/Users/andii/OneDrive/Documents/1RespositorioGithub/2Ensamble_deNovo_ExpDif"
setwd(directorio)
data <- read.table("Ex50.csv", sep = ",", header = T)

ggplot(data = data, aes(x = Ex, y = ExN50)) +
  geom_line(size = .25, color = "#0000FF", linetype = "dashed") +
  geom_point(color = "#0000FF", shape = 1, size = 3) +
  xlab ("Porcentaje expresado (%)") + #<-titulo de x
  ylab ("longitud del contig ExN50 (pb)") + #<-titulo de y
  theme_classic()



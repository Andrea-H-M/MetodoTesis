######################################################
## Tema: Redes de co-expresiÛn                      ##
## Autor: Olga Andrea Hernandez Miranda, Miranda H  ##
## Fecha: 27/08/2021                                ##
## Nota: Elaborar redes de co-expresiÛn             ##
######################################################

#Genotipo CH-VI

#--Preprocesamiento--

#Se cargan los paquetes necesarios


library(tidyverse)
library(igraph)
library(plotly)
library(BoolNet)
library(data.table)
library(ggplot2)
library("RColorBrewer")

#Lo primero que se hizo fue generar un archivo de texto plano separado por TAB en excel.
#Una vez que se tiene hecho, se procede a cargar los datos

directorio <- "C:/Users/andii/OneDrive/Documents/CursoED/Denovo/Redesv2"
setwd(directorio)

diff.expr <- fread("Red.csv", header = T)
diff.expr
idGene<-fread("ID.csv", header = T)
ids<-colnames(diff.expr)
dim(diff.expr)
colnames(diff.expr)<-c("V1",idGene[which(idGene$ID %in% ids)]$Genes)
colnames(diff.expr)<-c("V1",idGene[which(idGene$ID %in% ids)]$Genes)
head(diff.expr)

#Como se puede apreciar, hicimos un data frame que contiene en la primera columna los 
#nombres de los estadios de mediciÛn. Luego, como columnas cada uno de los genes, identificando 
#los nombres de las columnas. Debido a ello, ser· necesario darle una pulida a la estructura de datos,
#ya que se necesitan matrices para las pruebas de correlaciÛn

allexpr<-as.matrix(diff.expr[,-1])
row.names(allexpr)<-diff.expr$V1
exprVI<-allexpr[grep("CH-VI",row.names(allexpr)),]
exprI<-allexpr[grep("CH-I",row.names(allexpr)),]

#A continuaciÛn se selecciona una de las posibles matrices (CH-I, CH-VI o toda la tabla),
#con esa matriz se genera la tabla de correlaciÛn.

expr<-allexpr
expr<-allexpr[grep("CH-VI",row.names(allexpr)),]
gene.correlation <- cor(expr)
gene.correlation

#Crar CSV para quitar NA manualmente 
#write.csv(gene.correlation, "correlationCVI.csv", row.names = T, col.names = T)

#Cargar matriz sin NA
data <- read.table("correlationCVI.csv", sep = ",", header = T)
row.names(data) <- data[,1]
mat_data <- data.matrix(data[,1:ncol(data)]) 
gene.correlation2 <- mat_data[,-1]
gene.correlation2

#Se evalua una serie de valores de correlaciÛn para establecer el punto adecuado de correlaciÛn que recupera 
#una mejor construcciÛn de red. En este ejemplo vamos a cribar un espectro de valores entre 0.5 hasta 0.99

thresholds<-seq(0.5,0.99,0.01)
mean.connectivities<-NULL
scale.free.R2<-NULL
## Recorrer todos los valores de umbral
for(i in 1:length(thresholds))
{
  ## Se construye una red que corresponde al umbral de correlaciÛn especifico que se evalua en este paso
  
  ## Matriz de adyacencia
  current.adjacency <- gene.correlation2 > thresholds[i]
  ## Red
  threshold.network <- graph.adjacency(current.adjacency, mode="undirected")
  
  ## C·lculo de grados de los nodos
  node.degrees <- degree(threshold.network)
  
  ## Se guarda la conectividad promedio de esta red
  mean.connectivities[i] <- mean(node.degrees)
  
  ## EvaluaciÛn de la propiedad de escalamiento libre
  h <- hist(node.degrees,main=paste("Histograma de grado de los nodos para el umbral de correlacion",thresholds[i]))
  ## C·lculo de la frecuencia de grados
  degree.frequencies <- table(node.degrees)
  ## DeterminaciÛn por regresiÛn lineal para la transformaciÛn logarÌtmica de las frecuencias de grado
  lm.r <- lm(log10(as.numeric(names(degree.frequencies[-1]))) ~ log10(degree.frequencies[-1]))
  ## ObtenciÛn de R cuadrados como una estimaciÛn del ajuste a la propiedad de escalamiento libre
  s.lm <- summary(lm.r)
  scale.free.R2[i] <- s.lm[["adj.r.squared"]]
}

#El siguiente paso consiste en evaluar el mejor punto de corte para el grado de correlaciÛn que mejor ajusta a 
#la red generada a los modelos de redes observados comunmente en la naturaleza. Es importante tener una idea de 
#las propiedades deseadas en la red, por ejemplo, si se sabe cual es la conectividad promedio observada en redes de 
#regulaciÛn. AsÌ mismo el criterio del ajuste al modelo de red de libre escala, comunmente observado en la naturaleza,
#puede ser un par·metro ideal para determinar la mejor red.

plot(thresholds,mean.connectivities,type="o",col="red",lwd=3,xlab="Umbral de correlacion",ylab="Conectividad promedio")
plot(thresholds,scale.free.R2,type="o",col="blue",lwd=3,xlim=c(0.50,0.99),xlab="Umbral de correlacion",ylab=expression('Ajuste al modelo de escala libre (R  '^2*')'))


#---ConstrucciÛn de la red Ûptima---
#Se selecciona un valor para trabajar con esa red tmax va a ser la variable que contiene nuestro umbral de 
#correlaciÛn. A continuaciÛn se obtiene la red tomando en consideraciÛn una conexiÛn por cada pareja de genes
#con una correlaciÛn mayor al umbral que hemos determinado. Nuestra red quedar· guardada en nuestra variable gene.
#coexpression.network.

tmx<-which(max(scale.free.R2)==scale.free.R2)
tmax<-thresholds[tmx[length(tmx)]]

print(tmax)

adjacency.tmax <- gene.correlation2 > tmax

for(i in 1:ncol(adjacency.tmax))
{
  #  print(i)
  adjacency.tmax[i,i] <- FALSE
}

#adjacency.tmax[40,7] <- TRUE
#adjacency.tmax[7,40] <- TRUE

gene.coexpression.network <- graph.adjacency(adjacency.tmax, mode="undirected")

#Como una alternativa, se puede exportar una red en formato GML, el cual es importable en diferentes tipos
#de software de visualizaciÛn.
#write.graph(gene.coexpression.network,file="gene_coexpression_network2_CHVI.gml",format="gml")

#--VisualizaciÛn--
#El siguiente paso es darle un vistazo r·pido a la red generada

plot(gene.coexpression.network,vertex.size=node.degrees, vertex.label.color="black",
     vertex.label.dist=1, vertex.label.cex=0.25, layout=layout_nicely,
     vertex.color="blue", edge.color="grey")
plot(gene.coexpression.network)

#Y revisamos nuevamente la estructura de la red. Para ello se revisa la distribuciÛn del grado en los nodos
#y revisaremos los ajustes lineales que se obtienen.


network.degrees <- degree(gene.coexpression.network)
degree.histogram <- hist(network.degrees,col="royalblue",
                         xlab="Grado de nodo (K)",
                         ylab="N˙mero de nodos con K enlaces [P(K)]",
                         main = "CH-VI")

#write.table(network.degrees,file = "grados de nodo_CHVI")

degrees.frequencies <- degree.histogram[["counts"]]
node.degrees <- degree.histogram[["mids"]]
log10.degrees.frequencies <- log10(degrees.frequencies[-3])
log10.node.degrees <- log10(node.degrees[-3])
lm.r <- lm(log10.degrees.frequencies[1:6] ~ log10.node.degrees[1:6])
summary(lm.r)

plot(log10.node.degrees,log10.degrees.frequencies, col = "black", pch=16,
     xlab="log10 grado de nodo (K)",
     ylab="log10 n˙mero de nodos con K enlaces [P(K)]",
     abline(lm.r, col = "blue"))


#Prueba de bondad de ajuste para una exponenecial negativa

network.degree.distribution <- degree.distribution(gene.coexpression.network)
fit.scale.free <- power.law.fit(network.degree.distribution)
fit.scale.free[["KS.p"]]

#---Parametros de la red---
#A continuaciÛn revisamos algunos atributos de la red como el score de hub

## La mayor√≠a de los nodos en la redes libre de escala presentan un n√∫mero peque√±o de vecinos. 
## Sin embargo existen unos pocos nodos destacados que tiene un alto n√∫mero de vecinos. Este 
## √∫ltimo tipo de nodos se denominan hubs. La funci√≥n hub.score de igraph que recibe como 
## entrada una red calcula y almacena en el valor vector una puntuaci√≥n entre uno y cero para 
## cada nodo de la red. Cu√°nto mayor sea esta puntuaci√≥n de un nodo m√°s se ajustan sus 
## caracter√≠sticas a las de un hub de la red.

network.hub.scores <- hub.score(gene.coexpression.network)
hub.score.attributes <-network.hub.scores[["vector"]]
#write.table(hub.score.attributes,file = "hub_score2_CHVI")
plot_ly(y = network.hub.scores$vector[which(network.hub.scores$vector>0.1)], text=names(network.hub.scores$vector[which(network.hub.scores$vector>0.1)]),type="bar")
barplot(network.hub.scores$vector[which(network.hub.scores$vector>0.1)],las=2, col = "light green")


#La transitividad y el camino promedio
#Nodos y aristas
gene.coexpression.network
#619 5998

transitivity(gene.coexpression.network)
average.path.length(gene.coexpression.network)

## Para comprobar si el coeficiente de agrupamiento en lo suficientemente alto y 
## la longitud media del camino m√≠nimo entre nodos es lo suficientemente 
## peque√±a como para considerarla de mundo peque√±o es com√∫n generar redes libres de escala del 
## mismo orden y tama√±o de la estudiada para estimar la probabilidad de que por pura 
## aleatoriedad se obtenga una red similar a la estudiada pero con una longitud media del 
## camino m√≠nimo entre nodos inferior. La funci√≥n barabasi.game permite generar redes libres 
## de escala con el n√∫mero de nodos proporcionado en el argumento n.



number.of.added.edges<-10
clustering.coefficients<-NULL
for(i in 1:10000)
{
  if(i%%100==0){
    print(i)
  }
  random.scale.free.graph <- barabasi.game(n=dim(expr)[2],directed=FALSE)
  clustering.coefficients[i] <- transitivity(random.scale.free.graph)
}

sum(clustering.coefficients > transitivity(gene.coexpression.network,type="global")) / 10000

#Es una red de mundo pequeÒo

#---Buesqueda de patrones---

## Para la identificaci√≥n de clusteres o grupos de genes co-expresados en la red necesitamos 
## instalar los siguientes paquetes.


#install.packages("WGCNA")
library("WGCNA")
library("cluster")
allowWGCNAThreads()

## La identificaci√≥n de cl√∫steres o grupos de elementos se basa en una medida de similitud.
## La medida de similitud seleccionada en este estudio es 1 - correlacion
similarity.matrix <- 1 - gene.correlation2

## La funci√≥n hclust usando la matriz de similitudes como distancias y el m√©todo promedio
## para recalcular distancias calcula el clustering jer√°rquico.
hierarchical.clustering <- hclust(as.dist(similarity.matrix),method="average")

## La funci√≥n cutree permite cortar el √°rbol generado en el clustering jer√°rquico a distintas
## alturas para producir distintos n√∫meros de cl√∫steres.
hclust.2 <- cutree(hierarchical.clustering,k=2)
hclust.3 <- cutree(hierarchical.clustering,k=3)
hclust.4 <- cutree(hierarchical.clustering,k=4)
hclust.5 <- cutree(hierarchical.clustering,k=5)
hclust.6 <- cutree(hierarchical.clustering,k=6)
hclust.7 <- cutree(hierarchical.clustering,k=7)
hclust.8 <- cutree(hierarchical.clustering,k=8)
hclust.9 <- cutree(hierarchical.clustering,k=9)
hclust.10 <- cutree(hierarchical.clustering,k=10)

## La funci√≥n pam usa la matriz de similitudes como distancias para determinar cl√∫steres
## seg√∫n el m√©todo de partici√≥n entorno a medoides.
pam.2 <- pam(as.dist(similarity.matrix),k=2,diss=TRUE)
pam.3 <- pam(as.dist(similarity.matrix),k=3,diss=TRUE)
pam.4 <- pam(as.dist(similarity.matrix),k=4,diss=TRUE)
pam.5 <- pam(as.dist(similarity.matrix),k=5,diss=TRUE)
pam.6 <- pam(as.dist(similarity.matrix),k=6,diss=TRUE)
pam.7 <- pam(as.dist(similarity.matrix),k=7,diss=TRUE)
pam.8 <- pam(as.dist(similarity.matrix),k=8,diss=TRUE)
pam.9 <- pam(as.dist(similarity.matrix),k=9,diss=TRUE)
pam.10 <- pam(as.dist(similarity.matrix),k=10,diss=TRUE)

## La funci√≥n silhouette nos permite calcular la silueta de un clustering que sirve de medida
## para la bondad de dicho clustering.
sil2 <- silhouette(hclust.2,dist=similarity.matrix)
sil3 <- silhouette(hclust.3,dist=similarity.matrix)
sil4 <- silhouette(hclust.4,dist=similarity.matrix)
sil5 <- silhouette(hclust.5,dist=similarity.matrix)
sil6 <- silhouette(hclust.6,dist=similarity.matrix)
sil7 <- silhouette(hclust.7,dist=similarity.matrix)
sil8 <- silhouette(hclust.8,dist=similarity.matrix)
sil9 <- silhouette(hclust.9,dist=similarity.matrix)
sil10 <- silhouette(hclust.10,dist=similarity.matrix)

plot(sil3,col="blue")

hclust.sil.values <- c(summary(sil2)[["avg.width"]],summary(sil3)[["avg.width"]],summary(sil4)[["avg.width"]],summary(sil5)[["avg.width"]],summary(sil6)[["avg.width"]],summary(sil7)[["avg.width"]],summary(sil8)[["avg.width"]],summary(sil9)[["avg.width"]],summary(sil10)[["avg.width"]])

sil2 <- silhouette(pam.2)
sil3 <- silhouette(pam.3)
sil4 <- silhouette(pam.4)
sil5 <- silhouette(pam.5)
sil6 <- silhouette(pam.6)
sil7 <- silhouette(pam.7)
sil8 <- silhouette(pam.8)
sil9 <- silhouette(pam.9)
sil10 <- silhouette(pam.10)

pam.sil.values <- c(summary(sil2)[["avg.width"]],summary(sil3)[["avg.width"]],summary(sil4)[["avg.width"]],summary(sil5)[["avg.width"]],summary(sil6)[["avg.width"]],summary(sil7)[["avg.width"]],summary(sil8)[["avg.width"]],summary(sil9)[["avg.width"]],summary(sil10)[["avg.width"]])

## Representamos para los dos m√©todos de clustering, jer√°rquico y pam, y para diferentes n√∫meros
## de cl√∫steres la silueta correspondiente para elegir la mejor combinaci√≥n de m√©todo de 
## clustering y n√∫mero de cl√∫steres.
plot(2:10,pam.sil.values,type="o",col="blue",pch=0,ylim=c(0.3,0.8),xlab="N˙mero de clusters",ylab="Ancho promedio de la silueta",lwd=3)
lines(2:10,hclust.sil.values,type="o",col="red",pch=1,xlab="",ylab="",lwd=3)
legend("topright",legend=c("PAM","HCLUST"),col=c("blue","red"),pch=c(0,1),lwd=3)

## Visualizaci√≥n de cl√∫steres
pam.4[["clustering"]]

clustering.pam.4 <- pam.4[["clustering"]]
#write.table(clustering.pam.4,file="pam_4_CHVI.txt")

#Cargar y utilizar funciÛn IPAK
#ver vÌdeo https://www.youtube.com/watch?v=UjQz9SxG9rk
#ipak <- function(pkg){
# new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#if (length(new.pkg)) 
# install.packages(new.pkg, dependencies = TRUE)
#sapply(pkg, require, character.only = TRUE)
#}

library(fpc)
library(NbClust)
library(cluster)
library(factoextra)
library(tidyr)
library("RColorBrewer")

#calcular la matriz de distancia

library("viridis")
colores <- viridis(256)
#m.distancia <- get_dist(gene.correlation2 , method = "kendall") #el mÈtodo aceptado tambiÈn puede ser: "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" o "kendall"
fviz_dist(m.distancia, gradient = list(low = "blue", mid = "black", high = "yellow"))

#m.distancia <- get_dist(gene.correlation2 , method = "kendall") #el mÈtodo aceptado tambiÈn puede ser: "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" o "kendall"
fviz_dist(m.distancia,
          gradient = list(low = "#440154", mid = "#009c8c",
                          high = "#FDE725"))
#Morado #28024E
#turquoise4
#Viridis claro #009c8c


#probamos algunas visualizaciones

library(scatterplot3d)
scatterplot3d(gene.correlation2[,2:4], pch=20,
              color = c("red","#2E9FDF","green3","purple")[clustering.pam.4])
#Cluster jerarquico
res4 <- hcut(gene.correlation2, k = 4, stand = TRUE, method = "median")
fviz_dend(res4, rect = TRUE, cex = 0.5,
          k_colors = c("red","#2E9FDF","green3","purple"))


#Genotipo CH-I

#--Preprocesamiento--

#Se cargan los paquetes necesarios

library(tidyverse)
library(igraph)
library(plotly)
library(BoolNet)
library(data.table)

#Lo primero que se hizo fue generar un archivo de texto plano separado por TAB en excel.
#Una vez que se tiene hecho, se procede a cargar los datos

diff.expr <- fread("Red.csv", header = T)
head(diff.expr)
idGene<-fread("ID.csv", header = T)
ids<-colnames(diff.expr)
dim(diff.expr)
colnames(diff.expr)<-c("V1",idGene[which(idGene$ID %in% ids)]$Genes)
colnames(diff.expr)<-c("V1",idGene[which(idGene$ID %in% ids)]$Genes)
head(diff.expr)

#Como se puede apreciar, hicimos un data frame que contiene en la primera columna los 
#nombres de los estadios de mediciÛn. Luego, como columnas cada uno de los genes, identificando 
#los nombres de las columnas. Debido a ello, ser· necesario darle una pulida a la estructura de datos,
#ya que se necesitan matrices para las pruebas de correlaciÛn

allexpr<-as.matrix(diff.expr[,-1])
row.names(allexpr)<-diff.expr$V1
exprVI<-allexpr[grep("CH-VI",row.names(allexpr)),]
exprI<-allexpr[grep("CH-I",row.names(allexpr)),]

#A continuaciÛn se selecciona una de las posibles matrices (CH-I, CH-VI o toda la tabla),
#con esa matriz se genera la tabla de correlaciÛn.

expr<-allexpr
expr<-allexpr[grep("CH-I",row.names(allexpr)),]
gene.correlation <- cor(expr)
gene.correlation
#write.csv(gene.correlation, "correlationCI.csv", row.names = T, col.names = T)

#Cargar matriz sin NA
data <- read.table("correlationCI.csv", sep = ",", header = T)
row.names(data) <- data[,1]
mat_data <- data.matrix(data[,1:ncol(data)]) 
gene.correlation2 <- mat_data[,-1]
gene.correlation2 

#Se evalua una serie de valores de correlaciÛn para establecer el punto adecuado de correlaciÛn que recupera 
#una mejor construcciÛn de red. En este ejemplo vamos a cribar un espectro de valores entre 0.5 hasta 0.99

thresholds<-seq(0.5,0.99,0.01)
mean.connectivities<-NULL
scale.free.R2<-NULL
## Recorrer todos los valores de umbral
for(i in 1:length(thresholds))
{
  ## Se construye una red que corresponde al umbral de correlaciÛn especifico que se evalua en este paso
  
  ## Matriz de adyacencia
  current.adjacency <- gene.correlation2 > thresholds[i]
  ## Red
  threshold.network <- graph.adjacency(current.adjacency, mode="undirected")
  
  ## C·lculo de grados de los nodos
  node.degrees <- degree(threshold.network)
  
  ## Se guarda la conectividad promedio de esta red
  mean.connectivities[i] <- mean(node.degrees)
  
  ## EvaluaciÛn de la propiedad de escalamiento libre
  h <- hist(node.degrees,main=paste("Histograma de grado de los nodos para el umbral de correlacion",thresholds[i]))
  ## C·lculo de la frecuencia de grados
  degree.frequencies <- table(node.degrees)
  ## DeterminaciÛn por regresiÛn lineal para la transformaciÛn logarÌtmica de las frecuencias de grado
  lm.r <- lm(log10(as.numeric(names(degree.frequencies[-1]))) ~ log10(degree.frequencies[-1]))
  ## ObtenciÛn de R cuadrados como una estimaciÛn del ajuste a la propiedad de escalamiento libre
  s.lm <- summary(lm.r)
  scale.free.R2[i] <- s.lm[["adj.r.squared"]]
}

#El siguiente paso consiste en evaluar el mejor punto de corte para el grado de correlaciÛn que mejor ajusta a 
#la red generada a los modelos de redes observados comunmente en la naturaleza. Es importante tener una idea de 
#las propiedades deseadas en la red, por ejemplo, si se sabe cual es la conectividad promedio observada en redes de 
#regulaciÛn. AsÌ mismo el criterio del ajuste al modelo de red de libre escala, comunmente observado en la naturaleza,
#puede ser un par·metro ideal para determinar la mejor red.

plot(thresholds,mean.connectivities,type="o",col="red",
     lwd=3,xlab="Umbral de correlacion",ylab="Conectividad promedio")
plot(thresholds,scale.free.R2,type="o",col="blue",
     lwd=3,xlim=c(0.50,0.99),xlab="Umbral de correlacion",ylab=expression('Ajuste al modelo de escala libre (R  '^2*')'))

#---ConstrucciÛn de la red Ûptima---
#Se selecciona un valor para trabajar con esa red tmax va a ser la variable que contiene nuestro umbral de 
#correlaciÛn. A continuaciÛn se obtiene la red tomando en consideraciÛn una conexiÛn por cada pareja de genes
#con una correlaciÛn mayor al umbral que hemos determinado. Nuestra red quedar· guardada en nuestra variable gene.
#coexpression.network.

tmx<-which(max(scale.free.R2)==scale.free.R2)
tmax<-thresholds[tmx[length(tmx)]]

print(tmax)

adjacency.tmax <- gene.correlation2 > tmax

for(i in 1:ncol(adjacency.tmax))
{
  #  print(i)
  adjacency.tmax[i,i] <- FALSE
}

#adjacency.tmax[40,7] <- TRUE
#adjacency.tmax[7,40] <- TRUE

gene.coexpression.network <- graph.adjacency(adjacency.tmax, mode="undirected")

#Como una alternativa, se puede exportar una red en formato GML, el cual es importable en diferentes tipos
#de software de visualizaciÛn.

#write.graph(gene.coexpression.network,file="gene_coexpression_network2_CHI.gml",format="gml")

##--VisualizaciÛn--
#El siguiente paso es darle un vistazo r·pido a la red generada

plot(gene.coexpression.network,vertex.size=node.degrees, vertex.label.color="black",
     vertex.label.dist=1, vertex.label.cex=0.25, layout=layout_nicely,
     vertex.color="blue", edge.color="grey")
plot(gene.coexpression.network)

#Y revisamos nuevamente la estructura de la red. Para ello se revisa la distribuciÛn del grado en los nodos
#y revisaremos los ajustes lineales que se obtienen.

network.degrees <- degree(gene.coexpression.network)
degree.histogram <- hist(network.degrees,col="royalblue",
                         xlab="Grado de nodo (K)",
                         ylab="N˙mero de nodos con K enlaces [P(K)]",
                         main="CH-I")

#write.table(network.degrees,file = "grados de nodo_CHI")


degrees.frequencies <-degree.histogram[["counts"]]
node.degrees <- degree.histogram[["mids"]]
log10.degrees.frequencies <- log10(degrees.frequencies[-3])
log10.node.degrees <- log10(node.degrees[-3])
lm.r <- lm(log10.degrees.frequencies[1:6] ~ log10.node.degrees[1:6])
summary(lm.r)

plot(log10.node.degrees,log10.degrees.frequencies, col = "black", pch=16,
     xlab="log10 grado de nodo (K)",
     ylab="log10 n˙mero de nodos con K enlaces [P(K)]",
     abline(lm.r, col = "blue"))

ggplot(lm.r, aes(x = log10.node.degrees, y = log10.degrees.frequencies))+
  geom_point( color="blue") +
  geom_smooth(method = lm, color= "blue")+
  theme_classic()


#Prueba de bondad de ajuste para una exponenecial negativa

network.degree.distribution <- degree.distribution(gene.coexpression.network)
fit.scale.free <- power.law.fit(network.degree.distribution)
fit.scale.free[["KS.p"]]

#---Parametros de la red---
#A continuaciÛn revisamos algunos atributos de la red como el score de hub

network.hub.scores <- hub.score(gene.coexpression.network)
hub.score.attributes <-network.hub.scores[["vector"]]
#write.table(hub.score.attributes,file = "hub_score2_CHI")
#plot_ly(y = network.hub.scores$vector[which(network.hub.scores$vector>0.1)], text=names(network.hub.scores$vector[which(network.hub.scores$vector>0.1)]),type="bar")
barplot(network.hub.scores$vector[which(network.hub.scores$vector>0.1)],las=2, col = "blue")

#La transitividad y el camino promedio
#Nodos y aristas
gene.coexpression.network

transitivity(gene.coexpression.network)
average.path.length(gene.coexpression.network)

#Y se prueba si el coeficiente de clustering es alto comparado con lo que se espera encontrarlo en una red construida al azar.

number.of.added.edges<-10
clustering.coefficients<-NULL
for(i in 1:10000)
{
  if(i%%100==0){
    print(i)
  }
  random.scale.free.graph <- barabasi.game(n=dim(expr)[2],directed=FALSE)
  clustering.coefficients[i] <- transitivity(random.scale.free.graph)
}

sum(clustering.coefficients > transitivity(gene.coexpression.network,type="global")) / 10000

#---Buesqueda de patrones---

## Para la identificaci√≥n de clusteres o grupos de genes co-expresados en la red necesitamos 
## instalar los siguientes paquetes.


#install.packages("WGCNA")
library("WGCNA")
library("cluster")
allowWGCNAThreads()

## La identificaci√≥n de cl√∫steres o grupos de elementos se basa en una medida de similitud.
## La medida de similitud seleccionada en este estudio es 1 - correlacion
similarity.matrix <- 1 - gene.correlation2

## La funci√≥n hclust usando la matriz de similitudes como distancias y el m√©todo promedio
## para recalcular distancias calcula el clustering jer√°rquico.
hierarchical.clustering <- hclust(as.dist(similarity.matrix),method="average")

## La funci√≥n cutree permite cortar el √°rbol generado en el clustering jer√°rquico a distintas
## alturas para producir distintos n√∫meros de cl√∫steres.
hclust.2 <- cutree(hierarchical.clustering,k=2)
hclust.3 <- cutree(hierarchical.clustering,k=3)
hclust.4 <- cutree(hierarchical.clustering,k=4)
hclust.5 <- cutree(hierarchical.clustering,k=5)
hclust.6 <- cutree(hierarchical.clustering,k=6)
hclust.7 <- cutree(hierarchical.clustering,k=7)
hclust.8 <- cutree(hierarchical.clustering,k=8)
hclust.9 <- cutree(hierarchical.clustering,k=9)
hclust.10 <- cutree(hierarchical.clustering,k=10)

## La funci√≥n pam usa la matriz de similitudes como distancias para determinar cl√∫steres
## seg√∫n el m√©todo de partici√≥n entorno a medoides.
pam.2 <- pam(as.dist(similarity.matrix),k=2,diss=TRUE)
pam.3 <- pam(as.dist(similarity.matrix),k=3,diss=TRUE)
pam.4 <- pam(as.dist(similarity.matrix),k=4,diss=TRUE)
pam.5 <- pam(as.dist(similarity.matrix),k=5,diss=TRUE)
pam.6 <- pam(as.dist(similarity.matrix),k=6,diss=TRUE)
pam.7 <- pam(as.dist(similarity.matrix),k=7,diss=TRUE)
pam.8 <- pam(as.dist(similarity.matrix),k=8,diss=TRUE)
pam.9 <- pam(as.dist(similarity.matrix),k=9,diss=TRUE)
pam.10 <- pam(as.dist(similarity.matrix),k=10,diss=TRUE)

## La funci√≥n silhouette nos permite calcular la silueta de un clustering que sirve de medida
## para la bondad de dicho clustering.
sil2 <- silhouette(hclust.2,dist=similarity.matrix)
sil3 <- silhouette(hclust.3,dist=similarity.matrix)
sil4 <- silhouette(hclust.4,dist=similarity.matrix)
sil5 <- silhouette(hclust.5,dist=similarity.matrix)
sil6 <- silhouette(hclust.6,dist=similarity.matrix)
sil7 <- silhouette(hclust.7,dist=similarity.matrix)
sil8 <- silhouette(hclust.8,dist=similarity.matrix)
sil9 <- silhouette(hclust.9,dist=similarity.matrix)
sil10 <- silhouette(hclust.10,dist=similarity.matrix)

plot(sil3,col="blue")

hclust.sil.values <- c(summary(sil2)[["avg.width"]],summary(sil3)[["avg.width"]],summary(sil4)[["avg.width"]],summary(sil5)[["avg.width"]],summary(sil6)[["avg.width"]],summary(sil7)[["avg.width"]],summary(sil8)[["avg.width"]],summary(sil9)[["avg.width"]],summary(sil10)[["avg.width"]])

sil2 <- silhouette(pam.2)
sil3 <- silhouette(pam.3)
sil4 <- silhouette(pam.4)
sil5 <- silhouette(pam.5)
sil6 <- silhouette(pam.6)
sil7 <- silhouette(pam.7)
sil8 <- silhouette(pam.8)
sil9 <- silhouette(pam.9)
sil10 <- silhouette(pam.10)

pam.sil.values <- c(summary(sil2)[["avg.width"]],summary(sil3)[["avg.width"]],summary(sil4)[["avg.width"]],summary(sil5)[["avg.width"]],summary(sil6)[["avg.width"]],summary(sil7)[["avg.width"]],summary(sil8)[["avg.width"]],summary(sil9)[["avg.width"]],summary(sil10)[["avg.width"]])

## Representamos para los dos m√©todos de clustering, jer√°rquico y pam, y para diferentes n√∫meros
## de cl√∫steres la silueta correspondiente para elegir la mejor combinaci√≥n de m√©todo de 
## clustering y n√∫mero de cl√∫steres.
plot(2:10,pam.sil.values,type="o",col="blue",pch=0,ylim=c(0.3,0.8),xlab="N˙mero de clusters",ylab="Ancho promedio de la silueta",lwd=3)
lines(2:10,hclust.sil.values,type="o",col="red",pch=1,xlab="",ylab="",lwd=3)
legend("topright",legend=c("PAM","HCLUST"),col=c("blue","red"),pch=c(0,1),lwd=3)

## Visualizaci√≥n de cl√∫steres
pam.5[["clustering"]]

clustering.pam.5 <- pam.5[["clustering"]]
#write.table(clustering.pam.5,file="pam_5CHI.txt")

pam.6[["clustering"]]

clustering.pam.6 <- pam.6[["clustering"]]
#write.table(clustering.pam.6,file="pam_6CHI.txt")




library(fpc)
library(NbClust)
library(cluster)
library(factoextra)
library(tidyr)
library("RColorBrewer")

#calcular la matriz de distancia
#m.distancia <- get_dist(gene.correlation2 , method = "kendall") #el mÈtodo aceptado tambiÈn puede ser: "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" o "kendall"
fviz_dist(m.distancia, gradient = list(low = "blue", mid = "black", high = "yellow"))

#probamos algunas visualizaciones
library(scatterplot3d)
scatterplot3d(gene.correlation2[,2:4], pch=20,
              color = c("red","#2E9FDF","green3","purple","orange")[clustering.pam.5])
#Cluster jerarquico
res4 <- hcut(gene.correlation2, k = 5, stand = TRUE, method = "median")
fviz_dend(res4, rect = TRUE, cex = 0.5,
          k_colors = c("red","#2E9FDF","green3","purple","orange" ))



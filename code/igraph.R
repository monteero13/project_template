# Importamos la librería y cargamos los datos del tsv
library(igraph)
genes <- read.table("../results/genes_filtrados.txt", sep = '\t', header = FALSE, stringsAsFactors = FALSE)
grafo <- graph_from_data_frame(genes, directed=F)
# Mostramos la red
par(cex = 0.4)
plot(grafo)

# Comprobamos que no hay ningun nodo sin conexión.
cat("¿Están todos los nodos conectados?: ", is_connected(grafo), "\n")

# Comprobamos que no es dirigido
cat("¿El grafo es dirigido?: ", is_directed(grafo), "\n")

##COEFICIENTE CLUSTERING LOCAL
# Calcular el coeficiente de clustering por nodo
vertex_clustering <- transitivity(grafo, type = "local")
# Visualizar la red con coeficiente de clustering - si los vecinos estan conectadas, se visualiza en rojo!
V(grafo)$color <- ifelse(vertex_clustering > 0, "red", "lightblue")
set.seed(321)
png("../results/resultados_igraph/redClustering.png", width = 800, height = 600)
plot(grafo, vertex.size = 20, main = "Coeficiente de Clustering Local")
dev.off()

# Visualizar la tabla de valores y guardar el resultado
vertex_names <- V(grafo)$name
tabla_clustering <- data.frame(
  Nombre = vertex_names,
  Coef.Clustering = vertex_clustering,
  row.names = NULL
)
print('Los coeficientes de clustering de cada nodo son:')
tabla_clustering
if (!dir.exists("../results/resultados_igraph")) {
  dir.create("../results/resultados_igraph")
}
write.table(tabla_clustering, "../results/resultados_igraph/coefClustering.txt", sep = "\t", row.names = FALSE, quote = FALSE)

##COEFICIENTE CLUSTERING GLOBAL
global_clustering <- transitivity(grafo, type = "global")
cat("Global Clustering Coefficient:", global_clustering, "\n")

##Grado centralidad: numero de conexiones de cada gen
degree_centrality <- degree(grafo)
# Guardar el resultado
tabla_degree <- data.frame(
  Nombre = vertex_names,
  GradoCentralidad = degree_centrality,
  row.names = NULL
)
print('Los grados de centralidad de cada nodo son:')
tabla_degree
write.table(tabla_degree, "../results/resultados_igraph/gradoCentralidad.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Centralidad de cercanía: Mide la distancia promedio entre un nodo y todos los demás nodos
closeness_centrality <- closeness(grafo)
tabla_cercania <- data.frame(
  Nombre = vertex_names,
  DistCercania = closeness_centrality,
  row.names = NULL
)
print('La distancia promedio de cercanía de cada nodo son:')
tabla_cercania

# Conectividad:Mide la fortaleza de la conexión en el grafo.
connectivity <- edge_density(grafo)
cat("Conectividad:", connectivity, "\n")

## Centralidad de la intermediación: cuántas veces un nodo actúa como “puente” o intermediario entre otros nodos
betweenness_centrality <- betweenness(grafo)
# Normalize the betweenness values
betweenness_normalized <- (betweenness_centrality - min(betweenness_centrality)) / 
  (max(betweenness_centrality) - min(betweenness_centrality))
png("../results/resultados_igraph/redBetweenness.png", width = 800, height = 600)
plot(grafo, 
     vertex.size = 20, 
     vertex.color = heat.colors(100, rev=TRUE)[round(betweenness_normalized * 99) + 1],
     main = "Graph colored by Betweenness Centrality")
dev.off()
tabla_betweenness <- data.frame(
  Nombre = vertex_names,
  IndiceIntermediacion = betweenness_normalized,
  row.names = NULL
)
print('Los indices de betweenness de cada nodo son:')
tabla_betweenness
write.table(tabla_betweenness, "../results/resultados_igraph/indiceBetweenness.txt", sep = "\t", row.names = FALSE, quote = FALSE)

##MODULARIDAD
modularity_value <- modularity(cluster_louvain(grafo))
cat("Modularidad de la red:", modularity_value, "\n")

##DENSIDAD Y DISPERSION
density <- edge_density(grafo)
sparsity <- 1 - density
cat("Densidad de la red:", density, "\n")
cat("Dispersión de la red:", sparsity, "\n")

##CALCULO LONGITUDES
average_path_length <- mean_distance(grafo)
cat("Promedio de longitud de camino: ", average_path_length, "\n")
##DIAMETRO
diameter_modular <- diameter(grafo)
cat("Diámetro de la red: ", diameter_modular, "\n")

##Asortatividad de Grado y homofilia
assortativity_coeff <- assortativity_degree(grafo)
cat("Coeficiente asortatividad:", assortativity_coeff, "\n")

## Filtrar nodos que cumplen con los parametros(percentil 0.8 en degree, percentil 0.8 en betweenness, clustering menor a 0.5)
degree_threshold <- quantile(degree_centrality, 0.6)
betweenness_threshold <- quantile(betweenness_normalized, 0.8)
selected_nodes <- which(
  degree_centrality >= degree_threshold &          # Percentil 0.8 en degree
    betweenness_normalized >= betweenness_threshold & # Percentil 0.8 en betweenness
    vertex_clustering < 0.5                          # Clustering menor a 0.5
)
# Mostrar nodos seleccionados
print("Nodos seleccionados:")
print(selected_nodes)

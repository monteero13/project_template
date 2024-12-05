# Importamos la librería y cargamos los datos del tsv
library(igraph)
genes <- read.table("../results/genes_igraph.txt", sep = '\t', header = FALSE, stringsAsFactors = FALSE)
grafo <- graph_from_data_frame(genes, directed=F)
# Mostramos la red
par(cex = 0.4)
plot(grafo)

# Comprobamos que no hay ningun nodo sin conexión.
is_connected(grafo)

# Comprobamos que no es dirigido
is_directed(grafo)

##COEFICIENTE CLUSTERING LOCAL
# Calcular el coeficiente de clustering por nodo
vertex_clustering <- transitivity(grafo, type = "local")
V(grafo)$color <- ifelse(vertex_clustering > 0, "red", "lightblue")
# Visualizar la red con coeficiente de clustering - si los vecinos estan conectadas, se visualiza en rojo!
set.seed(321)
plot(grafo, vertex.size = 5, vertex.label = NA, main = "Coeficiente de Clustering Local")

##COEFICIENTE CLUSTERING GLOBAL
global_clustering <- transitivity(grafo, type = "global")
cat("Global Clustering Coefficient:", global_clustering, "\n")
#es un valor alto, 0.68

# Grado centralidad: numero de conexiones de cada gen
degree_centrality <- degree(grafo)
degree_centrality

# Centralidad de cercanía: Mide la distancia promedio entre un nodo y todos los demás nodos
closeness_centrality <- closeness(grafo)
closeness_centrality

# Conectividad:Mide la fortaleza de la conexión en el grafo.
connectivity <- edge_density(grafo)
connectivity
##Esto esta bajamente conectado, da un 0.16

# Centralidad de la intermediación: cuántas veces un nodo actúa como “puente” o intermediario entre otros nodos
betweenness_centrality <- betweenness(grafo)
betweenness_centrality
# Normalize the betweenness values for better visualization (optional)
betweenness_normalized <- (betweenness_centrality - min(betweenness_centrality)) / 
  (max(betweenness_centrality) - min(betweenness_centrality))
betweenness_normalized
plot(grafo, 
     vertex.size = 20, 
     vertex.color = heat.colors(100, rev=TRUE)[round(betweenness_normalized * 99) + 1], # Color nodes based on centrality # Optionally remove labels
     main = "Graph colored by Betweenness Centrality")
#Principales genes: PAX6, NF1, SEM1, WT1


##MODULARIDAD
modularity_value <- modularity(cluster_louvain(grafo))
cat("Modularidad de la red:", modularity_value, "\n")
#Obtenemos un 0.42, que no es un valor muy alto, no tendrá comunidades bien definidas

##DENSIDAD Y DISPERSION
density <- edge_density(grafo)
sparsity <- 1 - density
cat("Densidad de la red:", density, "\n")
cat("Dispersión de la red:", sparsity, "\n")
#Presentamos una alta dispersion de los datos
#Representa un sistema de señalización específico
#Sugiere que las proteinas siguen una ruta secuencial y específica
#es decir, que las proteinas tienen una funcion concreta, asi se evita que proteinas no implicadas intervengan
#provocando errores

##CALCULO LONGITUDES
average_path_length <- mean_distance(grafo)
average_path_length
#Obtenemos un promedio de camino largo, 2.74
#Se supone que esto indica que estamos ante una red modular y jerárquica, pero antes 
#obtenido también una modularidad baja
#¿?Se puede extraer entonces que nuestra red inicial es ciertamente dispersa y que
#las proteinas o genes intervienen en procesos muy concretos

##DIAMETRO
diameter_modular <- diameter(grafo)
cat("Diámetro de la red: ", diameter_modular, "\n")
#Pero ahora obtenemos un valor de 6, que indica que es una red de "mundo pequeño"

##Asortatividad de Grado y homofilia
#Asortatividad de Grado (Degree Assortativity): Mide la tendencia de los nodos a conectarse con otros nodos de grado similar.
#Homofilia (Homophily): Describe la tendencia de los nodos a conectarse con otros nodos similares en términos de atributos específicos.

assortativity_coeff <- assortativity_degree(grafo)
cat("Assortativity coefficient:", assortativity_coeff, "\n")
#se obtiene 0.48, que indica que las proteinas forman un complejo funcional especifico

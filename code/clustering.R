# Importamos librerías necesarias
library(igraph)
library(linkcomm)

# Cargar datos
file_path <- "../results/subgraph_diamond.tsv"
genes <- read.table(file_path, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# Filtrar filas inválidas
genes <- genes[!(genes$source == "source" & genes$target == "target" & genes$weight == "weight"), ]
genes$weight <- as.numeric(genes$weight)

# Crear grafo no dirigido
net <- graph_from_data_frame(genes, directed = FALSE)

# Verificar conectividad
cat("¿El grafo es conexo?:", is.connected(net), "\n")
cat("¿El grafo es dirigido?:", is.directed(net), "\n")

# Centralidades
cat("Grado de centralidad:\n")
print(degree(net))
cat("Centralidad de cercanía:\n")
print(closeness(net))
cat("Conectividad (densidad de bordes):", edge_density(net), "\n")

# Clustering por edge_betweenness
community <- cluster_edge_betweenness(net)
png(filename = "../results/clustering_edge_betweenness.png", width = 800, height = 600)
plot_dendrogram(community, mode = "hclust")  # Reemplaza dendPlot por plot_dendrogram
dev.off()
png(filename = "../results/clustering_edge_betweenness_graph.png", width = 800, height = 600)
plot(community, net)
dev.off()

# Clustering rápido
cfg <- cluster_fast_greedy(as.undirected(net))
png(filename = "../results/clustering_fast_greedy.png", width = 800, height = 600)
plot(cfg, as.undirected(net))
dev.off()

# Nodos coloreados por comunidad
V(net)$community <- cfg$membership
colrs <- adjustcolor(c("gray50", "tomato", "gold", "yellowgreen", "blue"), alpha = .6)
png(filename = "../results/clustering_coloreado.png", width = 800, height = 600)
plot(net, vertex.color = colrs[V(net)$community])
dev.off()

# Clustering por propagación de etiquetas
clp <- cluster_label_prop(net)
png(filename = "../results/clustering_label_prop.png", width = 800, height = 600)
plot(clp, net)
dev.off()

# Clustering Louvain
community_louvain <- cluster_louvain(net)
png(filename = "../results/clustering_louvain.png", width = 800, height = 600)
plot(community_louvain, net)
dev.off()

# Link communities con `linkcomm`
lc <- getLinkCommunities(genes)
png(filename = "../results/clustering_link_communities.png", width = 800, height = 600)
plot(lc, type = "graph", layout = "spencer.circle")
dev.off()

# Guardar genes del cluster de interés
genes_interes <- c("WNT10B", "WT1", "SEM1", "PAX6", "NF1")
nodos_genes_interes <- V(net)$name %in% genes_interes
cluster_containing_genes <- membership(cfg)[nodos_genes_interes]
genes_cluster <- V(net)$name[membership(cfg) == cluster_containing_genes]
ruta_archivo <- "../results/genes_cluster.txt"
write.table(genes_cluster, file = ruta_archivo, quote = FALSE, col.names = FALSE, row.names = FALSE)
cat("Los genes del cluster se han guardado en el archivo:", ruta_archivo, "\n")

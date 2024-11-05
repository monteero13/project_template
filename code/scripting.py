import requests
import pandas as pd
import io
import igraph as ig
import gseapy as gp
import matplotlib.pyplot as plt

# Lista de genes de interés: PAX6 y genes asociados al fenotipo de aniridia
genes = ["PAX6", "FOXC1", "WT1", "COL4A1", "PITX2"]

# Función para obtener interacciones de STRINGdb
def get_interactions(genes, species=9606, min_score=0.7):
    url = "https://string-db.org/api/tsv/network"
    params = {
        "identifiers": "%0d".join(genes),  # Separador %0d para STRING
        "species": species,                # Código de especie para humanos
        "required_score": int(min_score * 1000),  # Convertir puntaje a escala de 1000
        "network_type": "functional",
        "caller_identity": "aniridia_analysis"
    }
    response = requests.get(url, params=params)
    interactions = pd.read_csv(io.StringIO(response.text), sep="\t")
    return interactions

# Obtener interacciones de genes asociados a PAX6
interactions_df = get_interactions(genes)
print(interactions_df.head())

# Crear grafo desde las interacciones
g = ig.Graph.TupleList(interactions_df[["preferredName_A", "preferredName_B"]].itertuples(index=False), directed=False)

# Configuración de nodos en el grafo
for vertex in g.vs:
    vertex["color"] = "red" if vertex["name"] == "PAX6" else "blue"

# Aplicar clustering en la red para identificar módulos funcionales
clusters = g.community_multilevel()

# Visualizar la red
fig, ax = plt.subplots(figsize=(10, 10))
layout = g.layout("fr")
ig.plot(
    g, target=ax, layout=layout, vertex_label=g.vs["name"],
    vertex_color=g.vs["color"], vertex_size=20, edge_width=0.5,
    bbox=(800, 800), margin=50
)
plt.show()

# Análisis de enriquecimiento funcional para cada cluster usando ENRICHr a través de gseapy
def enrichment_analysis(genes):
    enr = gp.enrichr(gene_list=genes,
                     gene_sets="GO_Biological_Process_2021",
                     outdir=None,  # Eliminar 'description' y outdir=None para evitar salida no deseada
                     cutoff=0.05)  # Solo mostrar resultados con p-valor significativo
    return enr.results[["Term", "P-value", "Combined Score"]]

# Ejecutar análisis de enriquecimiento para cada cluster y mostrar resultados
for i, cluster in enumerate(clusters):
    genes_in_cluster = [g.vs[vertex]["name"] for vertex in cluster]
    print(f"\nCluster {i + 1}: Genes en el cluster -> {genes_in_cluster}")
    
    # Realizar análisis de enriquecimiento
    enrichment_results = enrichment_analysis(genes_in_cluster)
    print(f"Resultados de enriquecimiento para el Cluster {i + 1}")
    print(enrichment_results)

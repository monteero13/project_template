import pandas as pd
import networkx as nx
import requests
import logging
import os

# === CONFIGURACIÓN INICIAL ===
# Configuración de logging para mensajes de progreso y errores
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Definición de parámetros clave
STRING_FILE = "9606.protein.links.full.v12.0.txt.gz"  # Archivo de interacciones de STRINGdb
ALIASES_FILE = "9606.protein.aliases.v12.0.txt.gz"  # Archivo de aliases de STRINGdb
SEED_GENES = ["PAX6", "FOXC1", "WT1", "COL4A1", "PITX2"]  # Genes semilla
STRING_SCORE_THRESHOLD = 700  # Umbral de puntaje de interacción
GO_FILE = "go-basic.obo"  # Archivo OBO para el análisis de GO

# === 1. FUNCIONES DE CARGA Y FILTRADO DE DATOS ===
def load_ppi_data(STRING_FILE, STRING_SCORE_THRESHOLD):
    """
    Cargar las interacciones de STRINGdb desde un archivo CSV comprimido.
    :param STRING_FILE: Ruta del archivo de interacciones.
    :param STRING_SCORE_THRESHOLD: Umbral de puntuación para las interacciones.
    :return: DataFrame filtrado de interacciones.
    """
    logging.info(f"Cargando interacciones desde: {STRING_FILE}")
    
    if not os.path.exists(STRING_FILE):
        logging.error(f"El archivo {STRING_FILE} no se encuentra en la ruta especificada.")
        raise FileNotFoundError(f"El archivo {STRING_FILE} no se encuentra en la ruta especificada.")
    
    try:
        interactions = pd.read_csv(STRING_FILE, sep=" ", compression="gzip", 
                                   names=["protein1", "protein2", "combined_score"], low_memory=False)
        interactions['combined_score'] = pd.to_numeric(interactions['combined_score'], errors='coerce')
        interactions = interactions[interactions['combined_score'] >= STRING_SCORE_THRESHOLD]
        logging.info(f"Interacciones cargadas: {interactions.shape[0]} interacciones.")
        return interactions
    except Exception as e:
        logging.error(f"Error al cargar las interacciones: {e}")
        raise

def load_aliases_data(ALIASES_FILE):
    """
    Cargar los aliases de los genes desde un archivo de texto comprimido.
    :param ALIASES_FILE: Ruta del archivo de aliases de genes.
    :return: Diccionario de aliases de genes.
    """
    logging.info(f"Cargando aliases desde: {ALIASES_FILE}")
    
    if not os.path.exists(ALIASES_FILE):
        logging.error(f"El archivo {ALIASES_FILE} no se encuentra en la ruta especificada.")
        raise FileNotFoundError(f"El archivo {ALIASES_FILE} no se encuentra en la ruta especificada.")
    
    try:
        aliases = pd.read_csv(ALIASES_FILE, sep="\t", compression="gzip", header=None, 
                               names=["protein_id", "alias", "source"])
        alias_dict = aliases.groupby("protein_id")["alias"].apply(list).to_dict()
        logging.info(f"Aliases cargados: {len(alias_dict)} proteinas con aliases.")
        return alias_dict, aliases
    except Exception as e:
        logging.error(f"Error al cargar los aliases: {e}")
        raise

def map_seed_genes_to_protein_ids(SEED_GENES, aliases):
    """
    Mapear los genes semilla a sus respectivos protein_ids usando los aliases.
    :param SEED_GENES: Lista de genes semilla.
    :param aliases: DataFrame de aliases.
    :return: Conjunto de protein_ids correspondientes a los genes semilla.
    """
    seed_aliases = set()
    for gene in SEED_GENES:
        matching_proteins = aliases[aliases["alias"] == gene]["protein_id"].tolist()
        if not matching_proteins:
            logging.warning(f"Gene semilla '{gene}' no encontrado en el archivo de aliases.")
        seed_aliases.update(matching_proteins)
    
    if not seed_aliases:
        logging.error("No se encontraron protein_ids para los genes semilla proporcionados.")
        raise ValueError("No se encontraron protein_ids para los genes semilla proporcionados.")
    
    logging.info(f"Se mapearon {len(seed_aliases)} genes semilla a protein_ids.")
    return seed_aliases

# === 2. FUNCIONES DE DIAMOND ===
def diamond_algorithm(graph, seed_genes, max_added_nodes=100):
    """
    Implementación del algoritmo DIAMOnD para la expansión de un conjunto de genes semilla en una red PPI.
    :param graph: Grafo de interacciones (NetworkX).
    :param seed_genes: Lista de genes semilla (protein_ids).
    :param max_added_nodes: Número máximo de nodos a añadir.
    :return: Lista de nodos añadidos por DIAMOnD.
    """
    nodes_added = []
    candidate_scores = {}

    # Identificar vecinos de los genes semilla
    for seed in seed_genes:
        if seed not in graph:
            logging.warning(f"El protein_id '{seed}' no está presente en la red.")
            continue
        neighbors = graph.neighbors(seed)
        for neighbor in neighbors:
            if neighbor in seed_genes:
                continue
            candidate_scores[neighbor] = candidate_scores.get(neighbor, 0) + 1

    # Expandir la red con nodos adicionales
    for _ in range(max_added_nodes):
        if not candidate_scores:
            break  # Si no hay más candidatos, salir del bucle

        # Seleccionar el nodo con la puntuación más alta
        best_candidate = max(candidate_scores, key=candidate_scores.get)
        nodes_added.append(best_candidate)

        # Actualizar las puntuaciones de los vecinos
        for neighbor in graph.neighbors(best_candidate):
            if neighbor not in candidate_scores and neighbor not in seed_genes:
                candidate_scores[neighbor] = 0
            candidate_scores[neighbor] = candidate_scores.get(neighbor, 0) + 1

        del candidate_scores[best_candidate]  # Eliminar el nodo seleccionado de los candidatos

    logging.info(f"DIAMOnD añadió {len(nodes_added)} nodos al módulo expandido.")
    return nodes_added

# === 3. PROCESO PRINCIPAL ===
def main():
    try:
        # Cargar interacciones y aliases
        interactions = load_ppi_data(STRING_FILE, STRING_SCORE_THRESHOLD)
        alias_dict, aliases = load_aliases_data(ALIASES_FILE)

        # Mapear genes semilla a protein_ids
        seed_aliases = map_seed_genes_to_protein_ids(SEED_GENES, aliases)

        # Filtrar interacciones para genes semilla
        filtered_interactions = interactions[
            interactions["protein1"].isin(seed_aliases) | interactions["protein2"].isin(seed_aliases)
        ]
        
        # Construir la red PPI
        G = nx.Graph()
        for _, row in filtered_interactions.iterrows():
            G.add_edge(row["protein1"], row["protein2"], weight=row["combined_score"])

        # Añadir aliases como atributos de nodos
        for node in G.nodes:
            G.nodes[node]["aliases"] = alias_dict.get(node, ["Unknown"])
            G.nodes[node]["primary_alias"] = alias_dict.get(node, ["Unknown"])[0]

        # Ejecutar DIAMOnD
        expanded_module = diamond_algorithm(G, seed_genes=seed_aliases, max_added_nodes=100)

        # Construir subgrafo del módulo expandido
        logging.info("Construyendo subgrafo del módulo expandido...")
        subgraph = G.subgraph(expanded_module)
        logging.info(f"Número de nodos en el subgrafo: {subgraph.number_of_nodes()}")
        logging.info(f"Número de aristas en el subgrafo: {subgraph.number_of_edges()}")
    
    except Exception as e:
        logging.error(f"Error durante la ejecución del script: {e}")

# Ejecutar el script principal
if __name__ == "__main__":
    main()

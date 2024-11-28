import pandas as pd
import networkx as nx
import requests
import logging
import os
import argparse

# === FUNCIONES DE CONFIGURACIÓN ===
def parse_args():
    """
    Función para analizar los parámetros desde la línea de comandos.
    """
    parser = argparse.ArgumentParser(description="Propagación de redes a partir de genes semilla.")
    parser.add_argument('--protein-links', type=str, required=True, help="Ruta al archivo de interacciones STRINGdb.")
    parser.add_argument('--aliases', type=str, required=True, help="Ruta al archivo de aliases STRINGdb.")
    parser.add_argument('--go', type=str, required=True, help="Ruta al archivo GO (OBO).")
    parser.add_argument('--seed-genes', type=str, nargs='+', required=True, help="Lista de genes semilla.")
    parser.add_argument('--score-threshold', type=int, default=700, help="Umbral de puntaje para las interacciones.")
    return parser.parse_args()

# === FUNCIONES DE CARGA Y FILTRADO DE DATOS ===
def load_ppi_data(STRING_FILE, STRING_SCORE_THRESHOLD):
    """
    Cargar las interacciones de STRINGdb desde un archivo CSV comprimido.
    """
    logging.info(f"Cargando interacciones desde: {STRING_FILE}")
    
    if not os.path.exists(STRING_FILE):
        logging.error(f"El archivo {STRING_FILE} no se encuentra en la ruta especificada.")
        raise FileNotFoundError(f"El archivo {STRING_FILE} no se encuentra en la ruta especificada.")
    
    try:
        # Cargar todas las columnas del archivo
        interactions = pd.read_csv(STRING_FILE, sep=" ", compression="gzip", 
                                   header=None, low_memory=False)
        
        # Renombrar las columnas para que coincidan con el formato esperado
        interactions.columns = [
            "protein1", "protein2", "neighborhood", "neighborhood_transferred", 
            "fusion", "cooccurence", "homology", "coexpression", 
            "coexpression_transferred", "experiments", "experiments_transferred", 
            "database", "database_transferred", "textmining", 
            "textmining_transferred", "combined_score"
        ]
        
        # Filtrar las interacciones por el umbral de combined_score
        interactions['combined_score'] = pd.to_numeric(interactions['combined_score'], errors='coerce')
        interactions = interactions[interactions['combined_score'] >= STRING_SCORE_THRESHOLD]
        
        return interactions
    
    except Exception as e:
        logging.error(f"Error al cargar las interacciones: {e}")
        raise


def load_aliases_data(ALIASES_FILE):
    """
    Cargar los aliases de los genes desde un archivo de texto comprimido.
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


def get_string_id_from_gene(gene, species=9606):
    """
    Obtiene el stringId para un gen a partir de la API de STRINGdb.
    """
    url = f"https://string-db.org/api/json/get_string_ids?identifiers={gene}&species={species}"
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Error en la consulta para {gene}: {response.status_code} - {response.text}")
        return []
    
def map_seed_genes_to_protein_ids(SEED_GENES, graph):
    """
    Mapear los genes semilla a sus respectivos protein_ids usando la API de STRINGdb
    y verificar que estén presentes en la red de interacciones.
    """
    seed_aliases = set()
    not_found_genes = []

    # Mapear los genes semilla a protein_ids usando la API
    for gene in SEED_GENES:
        # Obtener los stringIds de los genes
        string_ids = get_string_id_from_gene(gene)
        
        if not string_ids:
            logging.warning(f"Gene semilla '{gene}' no encontrado en la API de STRING.")
            not_found_genes.append(gene)
        else:
            # Tomamos el primer stringId (en caso de que haya varios resultados)
            seed_aliases.update([string_id['stringId'] for string_id in string_ids])

    # Verificar que los protein_ids estén presentes en la red
    valid_seed_aliases = {pid for pid in seed_aliases if pid in graph.nodes}
    
    if not valid_seed_aliases:
        logging.error("No se encontraron protein_ids válidos para los genes semilla en la red.")
        raise ValueError(f"Los siguientes genes semilla no se encontraron en la red: {', '.join(not_found_genes)}")
    
    logging.info(f"Se mapearon {len(valid_seed_aliases)} genes semilla a protein_ids presentes en la red.")
    return valid_seed_aliases

# === 2. FUNCIONES DE DIAMOND ===
def diamond_algorithm(graph, seed_genes, max_added_nodes=100):
    """
    Implementación del algoritmo DIAMOnD para la expansión de un conjunto de genes semilla en una red PPI.
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

# === PROCESO PRINCIPAL ===
def main():
    try:
        # Analizar argumentos desde la línea de comandos
        args = parse_args()

        # Cargar interacciones y aliases
        interactions = load_ppi_data(args.protein_links, args.score_threshold)
        alias_dict, aliases = load_aliases_data(args.aliases)

        # Construir el grafo de interacciones
        G = nx.Graph()  # Crear un grafo vacío de NetworkX
        for _, row in interactions.iterrows():
            G.add_edge(row["protein1"], row["protein2"], weight=row["combined_score"])

        # Mapear genes semilla a protein_ids usando la API de STRINGdb
        logging.info("Mapeando genes semilla a protein_ids...")
        seed_aliases = map_seed_genes_to_protein_ids(args.seed_genes, G)  # Usando el grafo aquí
        
        # Verificar si los genes semilla están presentes en las interacciones
        logging.info("Verificando interacciones para los genes semilla...")
        seed_interactions = interactions[
            interactions["protein1"].isin(seed_aliases) | interactions["protein2"].isin(seed_aliases)
        ]
        logging.info(f"Interacciones encontradas para genes semilla: {seed_interactions.shape[0]}")

        # Si no hay interacciones, salir
        if seed_interactions.empty:
            logging.warning("No se encontraron interacciones para los genes semilla. Verifica las interacciones o los genes semilla.")
            return
        
        # Añadir aliases como atributos de nodos
        logging.info("Añadiendo aliases a los nodos...")
        for node in G.nodes:
            G.nodes[node]["aliases"] = alias_dict.get(node, ["Unknown"])
            G.nodes[node]["primary_alias"] = alias_dict.get(node, ["Unknown"])[0]

        # Ejecutar DIAMOnD
        logging.info("Aplicando algoritmo DIAMOnD...")
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

import pandas as pd
import networkx as nx
import logging
import os
import argparse

# === CONFIGURACIÓN INICIAL ===
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args():
    parser = argparse.ArgumentParser(description="Script para aplicar el algoritmo DIAMOnD.")
    parser.add_argument("--input", type=str, required=True, help="Ruta al archivo de interacciones (STRING FILE).")
    parser.add_argument("--output", type=str, default="results/subgraph_diamond.tsv", help="Archivo de salida.")
    parser.add_argument("--seed_genes", type=str, nargs='+', default=["WNT10B", "SEM1", "WT1", "PAX6", "NF1"],
                        help="Lista de genes semilla.")
    parser.add_argument("--threshold", type=int, default=400, help="Umbral de puntaje combinado para interacciones.")
    return parser.parse_args()

# === 1. FUNCIONES DE CARGA Y FILTRADO DE DATOS ===
def load_ppi_data(STRING_FILE, STRING_SCORE_THRESHOLD):
    logging.info(f"Cargando interacciones desde: {STRING_FILE}")
    if not os.path.exists(STRING_FILE):
        logging.error(f"El archivo {STRING_FILE} no se encuentra en la ruta especificada.")
        raise FileNotFoundError(f"El archivo {STRING_FILE} no se encuentra en la ruta especificada.")
    
    try:
        interactions = pd.read_csv(STRING_FILE, sep="\t", header=0)
        interactions.columns = [
            "protein1_hugo", "protein2_hugo", "combined_score"
        ]
        interactions['combined_score'] = pd.to_numeric(interactions['combined_score'], errors='coerce')
        interactions = interactions[interactions['combined_score'] >= STRING_SCORE_THRESHOLD]
        return interactions
    except Exception as e:
        logging.error(f"Error al cargar las interacciones: {e}")
        raise

# === 2. FUNCIONES DE DIAMOND ===
def diamond_algorithm(graph, seed_genes, max_added_nodes):
    nodes_added = []
    candidate_scores = {}

    # Verificación de la presencia de genes semilla en el grafo
    logging.info("Verificando si los genes semilla están en el grafo...")
    for seed in seed_genes:
        if seed not in graph:
            logging.warning(f"El gene semilla {seed} no está en la red.")
        else:
            logging.info(f"El gene semilla {seed} está en la red.")

    for seed in seed_genes:
        if seed not in graph:
            continue
        neighbors = graph.neighbors(seed)
        for neighbor in neighbors:
            if neighbor in seed_genes:
                continue
            candidate_scores[neighbor] = candidate_scores.get(neighbor, 0) + 1

    for _ in range(max_added_nodes):
        if not candidate_scores:
            break
        best_candidate = max(candidate_scores, key=candidate_scores.get)
        nodes_added.append(best_candidate)
        
        for neighbor in graph.neighbors(best_candidate):
            if neighbor not in candidate_scores and neighbor not in seed_genes:
                candidate_scores[neighbor] = 0
            candidate_scores[neighbor] = candidate_scores.get(neighbor, 0) + 1
        
        del candidate_scores[best_candidate]

    logging.info(f"DIAMOnD añadió {len(nodes_added)} nodos al módulo expandido.")
    return nodes_added

def export_subgraph(subgraph, output_file):
    edges = nx.to_pandas_edgelist(subgraph)
    edges.to_csv(output_file, sep="\t", index=False)

# === PROCESO PRINCIPAL ===
def main():
    args = parse_args()

    STRING_FILE = args.input
    OUTPUT_FILE = args.output
    SEED_GENES = args.seed_genes
    STRING_SCORE_THRESHOLD = args.threshold

    try:
        interactions = load_ppi_data(STRING_FILE, STRING_SCORE_THRESHOLD)

        G = nx.Graph()
        for _, row in interactions.iterrows():
            G.add_edge(row["protein1_hugo"], row["protein2_hugo"], weight=row["combined_score"])

        logging.info("Verificando interacciones para los genes semilla...")
        seed_interactions = interactions[interactions["protein1_hugo"].isin(SEED_GENES) | interactions["protein2_hugo"].isin(SEED_GENES)]
        logging.info(f"Interacciones encontradas para genes semilla: {seed_interactions.shape[0]}")

        if seed_interactions.empty:
            logging.warning("No se encontraron interacciones para los genes semilla.")
            return

        logging.info("Aplicando algoritmo DIAMOnD...")
        expanded_module = diamond_algorithm(G, seed_genes=SEED_GENES, max_added_nodes=200)

        # Asegurar que los genes semilla estén en el módulo expandido
        expanded_module = list(set(expanded_module) | set(SEED_GENES))

        logging.info("Construyendo subgrafo del módulo expandido...")
        subgraph = G.subgraph(expanded_module)
        logging.info(f"Número de nodos en el subgrafo: {subgraph.number_of_nodes()}")
        logging.info(f"Número de aristas en el subgrafo: {subgraph.number_of_edges()}")

        results_dir = os.path.dirname(OUTPUT_FILE)
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)

        export_subgraph(subgraph, OUTPUT_FILE)
        logging.info(f"Subgrafo exportado a: {OUTPUT_FILE}")
    
    except Exception as e:
        logging.error(f"Error durante la ejecución del script: {e}")

if __name__ == "__main__":
    main()

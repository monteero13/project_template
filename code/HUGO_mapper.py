import argparse
import os
import pandas as pd
import networkx as nx
import logging
import requests
import shutil
import gzip
from tqdm import tqdm
import mygene
import time

# === CONFIGURACIÓN DE LOGGING ===
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# === FUNCIONES PARA DESCARGAR Y DESCOMPRIMIR ===
def download_string_network(url, output_gz, output_txt):
    """Descargar y descomprimir la red STRINGdb."""
    if not os.path.exists("material"):
        os.makedirs("material")
    if not os.path.exists(output_gz):
        logging.info(f"Descargando la red desde STRINGdb: {url}")
        response = requests.get(url, stream=True)
        with open(output_gz, 'wb') as f:
            shutil.copyfileobj(response.raw, f)
        logging.info(f"Archivo descargado: {output_gz}")
    if not os.path.exists(output_txt):
        logging.info(f"Descomprimiendo {output_gz}...")
        with gzip.open(output_gz, 'rb') as f_in:
            with open(output_txt, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        logging.info(f"Archivo descomprimido: {output_txt}")

# === FUNCIONES PARA OBTENER EL MAPEADO DE STRING A HUGO ===
def get_string_to_hugo_mapping(string_ids):
    """
    Get mapping from STRING protein IDs to HUGO symbols using mygene
    
    Args:
        string_ids (set): Set of STRING protein IDs to map
        
    Returns:
        dict: Mapping from STRING protein IDs to HUGO symbols
    """
    logging.info("Getting STRING ID to HUGO mapping using mygene...")
    
    # Initialize mygene client
    mg = mygene.MyGeneInfo()
    
    # Convert STRING IDs to Ensembl protein IDs
    ensembl_ids = set()
    string_to_ensembl = {}
    for string_id in string_ids:
        if string_id.startswith('9606.ENSP'):
            ensembl_id = string_id.replace('9606.ENSP', 'ENSP')
            ensembl_ids.add(ensembl_id)
            string_to_ensembl[string_id] = ensembl_id
    
    # Query mygene in batches to avoid timeout
    batch_size = 1000
    ensembl_list = list(ensembl_ids)
    mapping = {}
    
    logging.info(f"Querying mygene for {len(ensembl_list)} proteins...")
    for i in tqdm(range(0, len(ensembl_list), batch_size)):
        batch = ensembl_list[i:i + batch_size]
        
        try:
            # Query mygene
            results = mg.querymany(batch, 
                                 scopes='ensembl.protein',
                                 fields='symbol',
                                 species='human',
                                 returnall=True)
            
            # Process results
            for hit in results['out']:
                if 'symbol' in hit:
                    ensembl_id = hit['query']
                    hugo = hit['symbol']
                    # Convert back to STRING ID format
                    string_id = [sid for sid, eid in string_to_ensembl.items() if eid == ensembl_id]
                    if string_id:
                        mapping[string_id[0]] = hugo
            
            # Small delay to avoid overwhelming the server
            time.sleep(0.1)
            
        except Exception as e:
            logging.warning(f"Error in batch {i}-{i+batch_size}: {str(e)}")
            continue
    
    logging.info(f"Successfully mapped {len(mapping)} proteins to HUGO symbols")
    return mapping

# === FUNCIONES DE CARGA Y FILTRADO DE DATOS ===
def load_ppi_data(STRING_FILE, STRING_SCORE_THRESHOLD):
    """Cargar interacciones de STRINGdb desde el archivo."""
    logging.info(f"Cargando interacciones desde: {STRING_FILE}")
    try:
        interactions = pd.read_csv(STRING_FILE, sep="\t", header=0)
        interactions.columns = ["protein1_hugo", "protein2_hugo", "combined_score"]
        interactions['combined_score'] = pd.to_numeric(interactions['combined_score'], errors='coerce')
        interactions = interactions.dropna(subset=['combined_score'])
        interactions = interactions[interactions['combined_score'] >= STRING_SCORE_THRESHOLD]
        return interactions
    except Exception as e:
        logging.error(f"Error al cargar las interacciones: {e}")
        raise

# === CONFIGURACIÓN DE `argparse` ===
def parse_args():
    parser = argparse.ArgumentParser(description="Script para descargar y procesar la red STRINGdb para Homo sapiens.")
    
    # Parámetros
    parser.add_argument("--input", type=str, required=True, help="Ruta del archivo STRINGdb descomprimido (por ejemplo, 9606.protein.links.v12.0.txt).")
    parser.add_argument("--output", type=str, help="Archivo de salida para el PPI filtrado.")
    parser.add_argument("--score_threshold", type=int, default=400, help="Umbral de puntaje de interacción (default: 400).")
    parser.add_argument("--url", type=str, default="https://string-db.org/cgi/download?sessionId=bMogNhQqaAZX&species_text=Homo+sapiens&settings_expanded=0&min_download_score=0&filter_redundant_pairs=0&delimiter_type=txt", help="URL para descargar la red STRINGdb.")
    
    return parser.parse_args()

# === PROCESO PRINCIPAL ===
def main():
    args = parse_args()
    
    STRING_FILE = args.input
    OUTPUT_FILE = args.output
    STRING_SCORE_THRESHOLD = args.score_threshold
    STRING_URL = args.url

    try:
        # Descargar y descomprimir la red (si es necesario)
        if not os.path.exists(STRING_FILE):
            download_string_network(STRING_URL, "material/9606.protein.links.v12.0.txt.gz", STRING_FILE)
        
        # Cargar y procesar datos de interacciones
        interactions = load_ppi_data(STRING_FILE, STRING_SCORE_THRESHOLD)

        # Obtener el mapeo de STRING a HUGO
        unique_proteins = set(interactions["protein1_hugo"]).union(set(interactions["protein2_hugo"]))
        string_to_hugo = get_string_to_hugo_mapping(unique_proteins)
        
        # Filtrar y convertir los IDs de STRING a HUGO
        interactions["protein1_hugo"] = interactions["protein1_hugo"].map(string_to_hugo)
        interactions["protein2_hugo"] = interactions["protein2_hugo"].map(string_to_hugo)
        interactions = interactions.dropna(subset=["protein1_hugo", "protein2_hugo"])

        # Construir el grafo
        G = nx.Graph()
        for _, row in interactions.iterrows():
            G.add_edge(row["protein1_hugo"], row["protein2_hugo"], weight=row["combined_score"])

        # Exportar el archivo filtrado
        interactions.to_csv(OUTPUT_FILE, sep="\t", index=False)
        logging.info(f"Archivo de salida generado: {OUTPUT_FILE}")

        # Analizar estadísticas
        analyze_network_statistics(OUTPUT_FILE)

    except Exception as e:
        logging.error(f"Error durante la ejecución del script: {e}")

if __name__ == "__main__":
    main()

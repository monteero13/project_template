import requests
import json
import mygene
import argparse
import os

# Función para obtener los IDs de UniProt a partir de símbolos de genes
def obtener_ids_uniprot(gene_symbols):
    """
    Convierte una lista de símbolos de genes en identificadores de UniProt usando MyGene.info.
    """
    mg = mygene.MyGeneInfo()
    query_result = mg.querymany(gene_symbols, scopes='symbol', fields='uniprot', species='human')

    study_gene_uniprot_ids = []
    for gene in query_result:
        uniprot_data = gene.get('uniprot', {}).get('Swiss-Prot')
        if isinstance(uniprot_data, list):  # Si hay múltiples identificadores, agregarlos todos
            study_gene_uniprot_ids.extend(uniprot_data)
        elif isinstance(uniprot_data, str):  # Si es un único identificador, agregarlo directamente
            study_gene_uniprot_ids.append(uniprot_data)
    return study_gene_uniprot_ids

# Función para realizar el análisis de enriquecimiento en STRINGdb
def realizar_enriquecimiento_string(uniprot_ids, fdr_threshold=0.1, categoria="Process"):
    """
    Realiza el análisis de enriquecimiento en STRINGdb para una categoría específica.
    """
    if not (0 <= fdr_threshold <= 1):
        raise ValueError("El FDR debe estar entre 0 y 1.")

    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "json"
    method = "enrichment"
    request_url = f"{string_api_url}/{output_format}/{method}"

    params = {
        "identifiers": "%0d".join(map(str, uniprot_ids)),  # Asegura que todos los IDs sean cadenas
        "species": 9606,
        "caller_identity": "teresavegamar@uma.es"
    }

    response = requests.post(request_url, data=params)
    if response.status_code != 200:
        raise ConnectionError("Error al conectar con la API de STRINGdb.")

    data = response.json()
    resultados = []
    for row in data:
        term = row["term"]
        preferred_names = ",".join(row["preferredNames"])
        fdr = float(row["fdr"])
        description = row["description"]
        result_category = row["category"]

        # Filtra por la categoría especificada y el umbral de FDR
        if result_category == categoria and fdr < fdr_threshold:
            resultados.append({
                "term": term,
                "genes": preferred_names,
                "fdr": fdr,
                "description": description,
                "category": result_category
            })

    return resultados

# Función para exportar los resultados a un archivo CSV
def exportar_resultados(resultados, filename):
    """
    Exporta los resultados a un archivo de texto.
    """
    with open(filename, "w") as f:
        f.write("Term\tGenes\tFDR\tCategory\tDescription\n")
        for result in resultados:
            linea = f"{result['term']}\t{result['genes']}\t{result['fdr']}\t{result['category']}\t{result['description']}\n"
            f.write(linea)

# Función para leer los genes desde un archivo
def leer_genes_desde_archivo(nombre_archivo):
    """
    Lee los nombres de genes desde un archivo de texto.
    Cada línea debe contener un gen.
    """
    with open(nombre_archivo, "r") as file:
        genes = [line.strip() for line in file if line.strip()]  # Elimina líneas vacías y espacios
    return genes

# Función principal para procesar el archivo de genes y realizar el análisis de enriquecimiento
def main():
    # Usamos argparse para recibir los parámetros desde la línea de comandos
    parser = argparse.ArgumentParser(description="Análisis de enriquecimiento de genes utilizando STRINGdb.")
    parser.add_argument("--input", type=str, required=True, help="Ruta al archivo de genes de entrada (por ejemplo, genes_cluster.txt).")
    parser.add_argument("--output_process", type=str, required=True, help="Ruta para guardar los resultados de la categoría 'Process'.")
    parser.add_argument("--output_hpo", type=str, required=True, help="Ruta para guardar los resultados de la categoría 'HPO'.")
    parser.add_argument("--fdr_threshold", type=float, default=0.1, help="Umbral de FDR (default: 0.1).")
    args = parser.parse_args()

    # Leer los genes del archivo de entrada
    genes = leer_genes_desde_archivo(args.input)

    # Obtener los IDs de UniProt para los genes
    uniprot_ids = obtener_ids_uniprot(genes)

    # Realizar el análisis de enriquecimiento para la categoría "Process"
    resultados_process = realizar_enriquecimiento_string(uniprot_ids, fdr_threshold=args.fdr_threshold, categoria="Process")
    exportar_resultados(resultados_process, args.output_process)
    print(f"Resultados para 'Process' exportados a '{args.output_process}'.")

    # Realizar el análisis de enriquecimiento para la categoría "HPO"
    resultados_hpo = realizar_enriquecimiento_string(uniprot_ids, fdr_threshold=args.fdr_threshold, categoria="HPO")
    exportar_resultados(resultados_hpo, args.output_hpo)
    print(f"Resultados para 'HPO' exportados a '{args.output_hpo}'.")

if __name__ == "__main__":
    main()

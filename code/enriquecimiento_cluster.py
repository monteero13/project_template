import os
import argparse
import mygene
import requests
import csv

# Función para obtener los IDs de UniProt a partir de símbolos de genes
def obtener_ids_uniprot(gene_symbols):
    mg = mygene.MyGeneInfo()
    query_result = mg.querymany(gene_symbols, scopes='symbol', fields='uniprot', species='human')

    study_gene_uniprot_ids = []
    for gene in query_result:
        uniprot_data = gene.get('uniprot', {}).get('Swiss-Prot')
        if isinstance(uniprot_data, list):
            study_gene_uniprot_ids.extend(uniprot_data)
        elif isinstance(uniprot_data, str):
            study_gene_uniprot_ids.append(uniprot_data)
    return study_gene_uniprot_ids

# Función para realizar el análisis de enriquecimiento en STRINGdb
def realizar_enriquecimiento_string(uniprot_ids, fdr_threshold=0.001, categoria="Process"):
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
    with open(filename, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=["term", "genes", "fdr", "category", "description"])
        writer.writeheader()
        for result in resultados:
            writer.writerow(result)

# Función principal para procesar el archivo de genes y realizar el análisis de enriquecimiento
def main():
    # Usamos argparse para recibir los parámetros desde la línea de comandos
    parser = argparse.ArgumentParser(description="Análisis de enriquecimiento de genes utilizando STRINGdb.")
    parser.add_argument("--input", type=str, required=True, help="Ruta al archivo de genes de entrada (por ejemplo, genes_cluster.txt).")
    parser.add_argument("--output", type=str, required=True, help="Ruta para guardar los resultados del enriquecimiento.")
    parser.add_argument("--fdr_threshold", type=float, default=0.001, help="Umbral de FDR (default: 0.001).")
    parser.add_argument("--categoria", type=str, default="Process", choices=["Process", "HPO"], help="Categoría de enriquecimiento (default: Process).")
    args = parser.parse_args()

    # Leer los genes del archivo de entrada
    with open(args.input, "r") as f:
        genes = [line.strip() for line in f if line.strip()]

    # Obtener los IDs de UniProt para los genes
    uniprot_ids = obtener_ids_uniprot(genes)

    # Realizar el análisis de enriquecimiento en STRINGdb
    resultados = realizar_enriquecimiento_string(uniprot_ids, fdr_threshold=args.fdr_threshold, categoria=args.categoria)

    # Exportar los resultados a un archivo CSV
    exportar_resultados(resultados, args.output)
    print(f"Análisis completado. Los resultados se han guardado en {args.output}")

if __name__ == "__main__":
    main()

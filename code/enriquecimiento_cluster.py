import os
import requests
import json
import mygene
import csv

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

def realizar_enriquecimiento_string(uniprot_ids, fdr_threshold=0.001, categoria="Process"):
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

def exportar_resultados(resultados, filename):
    """
    Exporta los resultados a un archivo CSV.
    """
    with open(filename, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=["term", "genes", "fdr", "category", "description"])
        writer.writeheader()
        for result in resultados:
            writer.writerow(result)

def procesar_archivo_cluster(ruta_archivo):
    """
    Procesa un archivo específico y realiza el análisis de enriquecimiento.
    """
    # Leer los genes del archivo
    with open(ruta_archivo, "r") as f:
        genes = [line.strip() for line in f if line.strip()]

    # Obtener IDs de UniProt
    uniprot_ids = obtener_ids_uniprot(genes)

    # Análisis de enriquecimiento para cada categoría
    for categoria in ["Process", "HPO"]:
        resultados = realizar_enriquecimiento_string(uniprot_ids, categoria=categoria)
        nombre_resultados = f"{os.path.splitext(ruta_archivo)[0]}_enriquecimiento_{categoria.lower()}.csv"
        exportar_resultados(resultados, nombre_resultados)
        print(f"Resultados para '{categoria}' exportados a '{nombre_resultados}'.")

# Ejecutar procesamiento del archivo
try:
    archivo_genes = "genes_cluster4.txt"  # Cambia esto por la ruta de tu archivo
    procesar_archivo_cluster(archivo_genes)
except Exception as e:
    print(f"Ha ocurrido un error inesperado: {e}")

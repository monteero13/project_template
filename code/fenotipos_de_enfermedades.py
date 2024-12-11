import requests
import pandas as pd
import argparse
import os

def obtener_fenotipos_asociados_a_enfermedad(disease_id):
    """
    Consulta la API de HPO para obtener los fenotipos asociados a una enfermedad.
    :param disease_id: Identificador de la enfermedad (por ejemplo, 'OMIM:106210' o 'ORPHA:252183').
    :return: Una lista de diccionarios con los fenotipos asociados a la enfermedad.
    """
    url = f'https://ontology.jax.org/api/network/annotation/{disease_id}'
    fenotipos_asociados = []
    
    try:
        response = requests.get(url)
        response.raise_for_status()  # Lanza un error si la respuesta no es exitosa
        data = response.json()
        
        # Extraer categorías y sus fenotipos
        categories = data.get("categories", {})
        for category_name, phenotypes in categories.items():
            for phenotype in phenotypes:
                fenotipos_asociados.append({
                    "disease_id": disease_id,
                    "disease_name": None,  # Esto lo llenaremos más tarde
                    "category": category_name,
                    "phenotype_id": phenotype.get("id"),
                    "phenotype_name": phenotype.get("name"),
                    "frequency": phenotype.get("metadata", {}).get("frequency", "Unknown")
                })
        print(f"Fenotipos asociados obtenidos para la enfermedad {disease_id}")
    except requests.exceptions.RequestException as e:
        print(f"Error al conectarse a la API de HPO para la enfermedad {disease_id}: {e}")
    
    return fenotipos_asociados

def obtener_fenotipos_para_lista_enfermedades(lista_enfermedades):
    """
    Obtiene los fenotipos asociados para una lista de enfermedades.
    :param lista_enfermedades: Lista de identificadores de enfermedades.
    :return: Lista de diccionarios con las enfermedades y sus fenotipos asociados.
    """
    resultados = []
    
    for enfermedad in lista_enfermedades:
        fenotipos = obtener_fenotipos_asociados_a_enfermedad(enfermedad['id'])
        # Asociamos el nombre de la enfermedad
        for fenotipo in fenotipos:
            fenotipo['disease_name'] = enfermedad['name']
        resultados.extend(fenotipos)
    
    return resultados

def leer_enfermedades(archivo):
    """
    Lee el archivo de enfermedades y extrae la información.
    :param archivo: Ruta del archivo con las enfermedades.
    :return: Lista de diccionarios con los IDs y nombres de enfermedades.
    """
    enfermedades = []
    try:
        with open(archivo, 'r') as file:
            for line in file:
                parts = line.split()  # Divide la línea en una lista
                if len(parts) >= 2:
                    id_enfermedad = parts[0]
                    nombre_enfermedad = ' '.join(parts[1:])  # El resto es el nombre de la enfermedad
                    enfermedades.append({"id": id_enfermedad, "name": nombre_enfermedad})
    except FileNotFoundError:
        print(f"Error: El archivo {archivo} no fue encontrado.")
        raise
    return enfermedades

# Función principal
def main():
    # Usamos argparse para recibir los parámetros desde la línea de comandos
    parser = argparse.ArgumentParser(description="Obtiene fenotipos asociados a enfermedades y los guarda en un archivo.")
    
    parser.add_argument("--input", type=str, required=True, help="Ruta al archivo con las enfermedades.")
    parser.add_argument("--output", type=str, required=True, help="Ruta para guardar los resultados de los fenotipos.")
    
    args = parser.parse_args()

    # Verificar que los archivos existen
    if not os.path.exists(args.input):
        print(f"Error: El archivo de entrada {args.input} no existe.")
        return

    # Leer el archivo de enfermedades proporcionado
    enfermedades = leer_enfermedades(args.input)

    # Obtener los fenotipos asociados a las enfermedades
    fenotipos_asociados = obtener_fenotipos_para_lista_enfermedades(enfermedades)

    # Si no se obtuvieron fenotipos, mostrar mensaje y salir
    if not fenotipos_asociados:
        print("No se han obtenido fenotipos para las enfermedades proporcionadas.")
        return

    # Crear un DataFrame para estructurar los resultados
    df_fenotipos = pd.DataFrame(fenotipos_asociados)

    # Exportar resultados a CSV
    df_fenotipos.to_csv(args.output, index=False, encoding='utf-8')

    # Mostrar un resumen de los resultados
    print(f"Fenotipos asociados guardados en: {args.output}")
    print(df_fenotipos.head())

if __name__ == "__main__":
    main()

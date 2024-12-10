import requests
import pandas as pd

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

# Leer el archivo de texto línea por línea
enfermedades = []
with open('../results/diseases_for_HP_0000526', 'r') as file:
    for line in file:
        # Suponiendo que los datos están separados por espacios o tabuladores
        parts = line.split()  # Esto divide la línea en una lista
        if len(parts) >= 2:
            id_enfermedad = parts[0]
            nombre_enfermedad = ' '.join(parts[1:])  # El resto es el nombre de la enfermedad
            enfermedades.append({"id": id_enfermedad, "name": nombre_enfermedad})

# Obtener fenotipos asociados
fenotipos_asociados = obtener_fenotipos_para_lista_enfermedades(enfermedades)

# Crear un DataFrame para estructurar los resultados
df_fenotipos = pd.DataFrame(fenotipos_asociados)

# Exportar resultados a CSV
df_fenotipos.to_csv('../results/fenotipos_de_enfermedades.csv', index=False, encoding='utf-8')

# Mostrar un resumen de los resultados
print(df_fenotipos.head())

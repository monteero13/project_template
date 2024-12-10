import requests
import shutil
import os
import pandas as pd

#Paso 1: obtener los genes y enfermedades asociadas al HPO dado (HP:0000526)
def buscarGenesEnfermedadesHPO(HPO):
    url = f'https://ontology.jax.org/api/network/annotation/{HPO}'
    genes_asociados=[]
    try:
        response = requests.get(url)
        response.raise_for_status()
        data=response.json()
        genes_asociados = data.get('genes', [])
        enfermedades_asociadas = data.get('diseases', [])
        print(f'Se han obtenido los genes y enfermedades asociadas')
    except requests.exceptions.RequestException as e:
        print(f"Error al conectarse a la API de HPO: {e}")
    return genes_asociados, enfermedades_asociadas

def guardarGenesEnArchivo(genes_asociados, HPO):
    no, codigo=HPO.split(':')
    directorio=f'../results/genes_for_HP_{codigo}'
    try:
        with open(directorio, "w", encoding='utf-8') as file:
            # Escribir la cabecera
            file.write("id\tname\n")
            
            # Escribir cada gen en una nueva línea
            for gen in genes_asociados:
                file.write(f"{gen['id']}\t{gen['name']}\n")
        file.close()
        print(f"Genes guardados en el archivo genes_for_HP_{codigo}.tsv")
    except Exception as e:
        print(f"Error al guardar los genes en el archivo: {e}")

def guardarEnfermedadesEnArchivo(enfermedades_asociadas, HPO):
    no, codigo=HPO.split(':')
    directorio=f'../results/diseases_for_HP_{codigo}'
    try:
        with open(directorio, "w", encoding='utf-8') as file:
            # Escribir la cabecera
            file.write("id\tname\n")
            
            # Escribir cada gen en una nueva línea
            for enfermedad in enfermedades_asociadas:
                file.write(f"{enfermedad['id']}\t{enfermedad['name']}\n")
        file.close()
        print(f"Enfermedades guardadas en el archivo diseases_for_HP_{codigo}.tsv")
    except Exception as e:
        print(f"Error al guardar las enfermedades en el archivo: {e}")

#Paso 2: conectarnos a STRINGDB
def descargarRedString(genes, especie, min_score):
    urlTSV = f'https://string-db.org/api/tsv/network'
    urlImagen = f'https://string-db.org/api/image/network'
    
    parametros = {
        # Crear una cadena con los símbolos de genes separados por '%'
        "identifiers" : "%0d".join([gen['name'] for gen in genes]),
        "species" : especie,
        "network_type" : "functional",
        "required_score": int(min_score * 1000)
    }
    
    try:
        response1 = requests.get(urlTSV, params=parametros, stream=True)
        response1.raise_for_status()
        os.makedirs('../results', exist_ok=True)
        with open('../results/red_inicial.tsv', 'w', encoding='utf-8') as file:
            file.write(response1.text)
            print(f'Se ha descargado el archivo red_inicial.tsv')
        
        response2=requests.get(urlImagen, params=parametros, stream=True)
        response2.raise_for_status()
        with open('../results/red_inicial.png', 'wb') as out_file:
            shutil.copyfileobj(response2.raw, out_file)
            print(f'Se ha descargado la imagen red_inicial.png')
        
    except requests.excceptions.RequestException as e:
        print(f"Error al conectarse a la API de STRINGDB: {e}")
        

#Paso 3: Limpiar los datos
def limpiarDatos():
    data = pd.read_csv('../results/red_inicial.tsv', sep='\t')
    # Seleccionar solo las columnas 'preferredName_A' y 'preferredName_B'
    selected_columns = data[['preferredName_A', 'preferredName_B']]
    # Eliminar duplicados
    selected_columns = selected_columns.drop_duplicates()
    # Guardar estas columnas en un nuevo archivo de texto
    selected_columns.to_csv('../results/genes_filtrados.txt', sep='\t', index=False, header=False)
    print("Se ha limpiado el archivo .tsv en genes_filtrados.txt")


fenotipo='Aniridia'
codigoFenotipo='HP:0000526'
especie=9606
min_score=0.4

print(f'Buscando genes asociados a {fenotipo} en HPO...')
genesAsociados, enfermedadesAsociadas=buscarGenesEnfermedadesHPO(codigoFenotipo)
print('Guardando genes y enfermedades...')
guardarGenesEnArchivo(genesAsociados, codigoFenotipo)
guardarEnfermedadesEnArchivo(enfermedadesAsociadas, codigoFenotipo)

if genesAsociados:
    print(f'Descargando red de proteinas de STRINGDB...')
    descargarRedString(genesAsociados, especie, min_score)
    print(f'Limpiando datos para análisis...')
    limpiarDatos()

else:
    print(f'No hay genes asociados a {fenotipo}')

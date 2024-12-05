import requests
import shutil
import os
import pandas as pd

###DEBERIA USAR TRY, MODIFICAR DESPUES

#El primer paso es obtener los genes asociados al HPO dado (HP:0000526)

def buscarGenesHPO(HPO):
    url = f'https://ontology.jax.org/api/network/annotation/{HPO}'
    response = requests.get(url)
    genes_asociados=[]
    if response.status_code==200:
        data=response.json()
        genes_asociados = data.get('genes', [])
    else:
        ###Se podría escribir una excepcion aquí, a revisar###
        print("No se puede obtener los genes asociados")

    return genes_asociados


codigoFenotipo='HP:0000526'
genesAniridia=buscarGenesHPO(codigoFenotipo)

###### FALTA LA FUNCION DE DESCARGAR EL PNG
#with open('../results/red_inicial.png', 'wb') as out_file:
#            shutil.copyfileobj(response.raw, out_file)

#El segundo paso es conectarnos a STRINGDB
def descargarRedString(diccionarioGenes, especie='human', min_score=0.4):
    url = f'https://string-db.org/api/tsv/network'
    
    parametros = {
        # Crear una cadena con los símbolos de genes separados por '%'
        "identifiers" : "%0d".join([gen['name'] for gen in diccionarioGenes]),
        "species" : especie,
        "network_type" : "functional",
        "required_score": int(min_score * 1000)
    }
    
    response = requests.get(url, params=parametros, stream=True)
    
    if response.status_code==200:
        os.makedirs('../results', exist_ok=True)
        
        with open('../results/red_inicial.tsv', 'w', encoding='utf-8') as file:
            file.write(response.text)
    else:
        ###Se podría escribir una excepcion aquí, a revisar###
        print("No se puede obtener los genes asociados")

descargarRedString(genesAniridia)

##PONER COMO FUNCION ESTO
#Paso 3: Limpiar los datos
data = pd.read_csv('../results/red_inicial.tsv', sep='\t')

# Seleccionar solo las columnas 'preferredName_A' y 'preferredName_B'
selected_columns = data[['preferredName_A', 'preferredName_B']]

# Eliminar duplicados
selected_columns.drop_duplicates(inplace=True)

# Guardar estas columnas en un nuevo archivo de texto
selected_columns.to_csv('../results/genes_igraph.txt', sep='\t', index=False, header=False)

print("Se han guardado las columnas 'preferredName_A' y 'preferredName_B' en results/genes_igraph.txt")
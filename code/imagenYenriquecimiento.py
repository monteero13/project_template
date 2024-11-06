import pandas as pd
import requests
import shutil

# Cargar los genes desde un archivo TSV
genes_df = pd.read_csv('./genes.tsv', sep='\t')  # Asegúrate de usar la ruta correcta

# Limpiar los nombres de las columnas, eliminando espacios extras
genes_df.columns = genes_df.columns.str.strip()

# Verificar las primeras filas y las columnas del DataFrame después de limpiar los nombres
print(genes_df.head())  # Muestra las primeras filas
print(genes_df.columns)  # Muestra los nombres de las columnas

# Extraer los símbolos de los genes (en la columna "name")
genes_asociados_aniridia = genes_df['name'].tolist()

# Construir la URL para la solicitud a STRING DB
base_url = 'https://string-db.org/api/image/'

# Los símbolos de genes deben ser separados por '%0d' (salto de línea)
symbols = "%0d".join(genes_asociados_aniridia)

species = "9606"  # Especie humana (ID en STRING DB para humanos)

# Construir la URL completa para la solicitud
url_request = f"{base_url}network?identifiers={symbols}&species={species}&network_type=functional"
print("Requesting to " + url_request)

# Realizar la solicitud HTTP GET a STRING DB
response = requests.get(url_request, stream=True)

# Guardar la imagen de la red de interacción en un archivo
with open('./red_aniridia.png', 'wb') as out_file:
    shutil.copyfileobj(response.raw, out_file)

# Liberar la respuesta de la solicitud
del response

print("Red de interacción guardada como 'red_aniridia.png'")

# Para realizar un análisis de enriquecimiento, podemos hacer una solicitud diferente (por ejemplo, para obtener términos GO)
enrichment_url = f"https://string-db.org/api/enrichment?identifiers={symbols}&species={species}&threshold=0.05&caller=python"
response_enrichment = requests.get(enrichment_url)

# Guardar los resultados de enriquecimiento en un archivo
with open('./enrichment_aniridia.tsv', 'wb') as out_file:
    out_file.write(response_enrichment.content)

# Mostrar el contenido de los resultados de enriquecimiento
enrichment_data = pd.read_csv('./enrichment_aniridia.tsv', sep='\t')
print(enrichment_data)
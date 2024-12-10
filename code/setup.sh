#!/bin/bash

# ============================
# Configuración del entorno
# ============================

set -e  # Termina el script si ocurre un error

PROJECT_DIR="$(cd "$(dirname "$0")"/.. && pwd)"  # Directorio raíz del proyecto (donde se encuentra setup.sh)
SCRIPT_DIR="$PROJECT_DIR/code"  # Directorio donde están los scripts
RESULTS_DIR="$PROJECT_DIR/results"  # Directorio donde se guardan los resultados
MATERIAL_DIR="$SCRIPT_DIR/materials"  # Directorio de materiales

VENV_DIR="$SCRIPT_DIR/venv"  # Ruta del entorno virtual

if [ ! -d "$VENV_DIR" ]; then
    echo "Creando entorno virtual en $VENV_DIR..."
    python3 -m venv "$VENV_DIR"
fi

# Activar el entorno virtual
source "$VENV_DIR/bin/activate"

echo "Instalando dependencias..."
pip install --upgrade pip
pip install -r "$SCRIPT_DIR/requirements.txt"


# ============================
# Paso 1: Propagación de la Red
# ============================

echo "Mapeando la red al formato HUGO"

python "$SCRIPT_DIR"/HUGO_mapper.py --input "$MATERIAL_DIR/9606.protein.links.v12.0.txt" --output "$RESULTS_DIR/homo_sapiens_string_ppi_filtered.tsv" --score_threshold 400

echo "Proceso completado. Los resultados están en la carpeta $RESULTS_DIR."

echo "Iniciando la propagación de redes a partir de los genes semilla (WNT10B, SEM1, WT1, PAX6, NF1)"

INPUT_FILE="$RESULTS_DIR/homo_sapiens_string_ppi_filtered.tsv"
OUTPUT_FILE="$RESULTS_DIR/subgraph_diamond.tsv" 
SEED_GENES="WNT10B SEM1 WT1 PAX6 NF1"
THRESHOLD=400

# Ejecutar el script de propagación de redes
python "$SCRIPT_DIR"/Network_propagation.py \
    --input "$INPUT_FILE" \
    --output "$OUTPUT_FILE" \
    --seed_genes $SEED_GENES \
    --threshold $THRESHOLD

# Verificar si el paso de propagación fue exitoso
if [ $? -ne 0 ]; then
    echo "Error en la propagación de redes. Abortando."
    exit 1
fi

# ============================
# Paso 2: Clusterización de la Red
# ============================

echo "Iniciando la clusterización de la red expandida..."

# Verificar si R está instalado
if ! command -v R &> /dev/null
then
    echo "R no está instalado. Por favor, instálalo antes de continuar."
    exit 1
fi

# Instalar paquetes R necesarios si no están instalados
Rscript -e "if (!require(igraph)) install.packages('igraph', repos='http://cran.rstudio.com/')" 
Rscript -e "if (!require(linkcomm)) install.packages('linkcomm', repos='http://cran.rstudio.com/')"

# Ejecutar el script R
Rscript "$SCRIPT_DIR"/clustering.R

echo "Proceso completado. Los resultados están en la carpeta $RESULTS_DIR."

# ============================
# Paso 3: Representación del Grafo
# ============================

echo "Generando la representación gráfica de la red..."

# Ejecutar el script de representación del grafo en Python
python "$SCRIPT_DIR"/graph_representation.py \
    --input "$RESULTS_DIR/subgraph_diamond.tsv" \
    --output "$RESULTS_DIR"


# Verificar si la representación fue exitosa
if [ $? -ne 0 ]; then
    echo "Error en la representación gráfica de la red. Abortando."
    exit 1
fi

# ============================
# Paso 4: Análisis de Enriquecimiento Funcional
# ============================

echo "Iniciando el análisis de enriquecimiento funcional (GO)..."

# Ejecutar el script de análisis de enriquecimiento funcional
python3 "$SCRIPT_DIR"/enriquecimiento_cluster.py \
    --input "$RESULTS_DIR/genes_cluster.txt" \
    --output "$RESULTS_DIR/enrichment_results.csv" \
    --fdr_threshold 0.001 --categoria "Process"

# Verificar si el análisis de enriquecimiento fue exitoso
if [ $? -ne 0 ]; then
    echo "Error en el análisis de enriquecimiento funcional. Abortando."
    exit 1
fi

# ============================
# Finalización
# ============================

echo "¡El proceso completo se ha ejecutado con éxito!"


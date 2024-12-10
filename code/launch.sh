#!/bin/bash

# ============================
# Configuración del entorno
# ============================

set -e  # Termina el script si ocurre un error

PROJECT_DIR="$(cd "$(dirname "$0")"/.. && pwd)"  # Directorio raíz del proyecto
SCRIPT_DIR="$PROJECT_DIR/code"  # Directorio donde están los scripts
RESULTS_DIR="$PROJECT_DIR/results"  # Directorio donde se guardan los resultados
MATERIAL_DIR="$SCRIPT_DIR/material"  # Directorio de materiales

# ============================
# Paso 1: Ejecutar el setup.sh
# ============================

echo "Ejecutando setup.sh para instalar dependencias..."
# Ejecutar setup.sh para instalar dependencias
bash "$SCRIPT_DIR/setup.sh"

echo "Dependencias instaladas con éxito."

# ============================
# Paso 2: Propagación de la Red
# ============================

echo "Mapeando la red al formato HUGO..."

python3 "$SCRIPT_DIR/HUGO_mapper.py" --input "$STRING_TXT" --output "$RESULTS_DIR/homo_sapiens_string_ppi_filtered.tsv" --score_threshold 400

echo "Proceso de mapeo completado. Los resultados están en la carpeta $RESULTS_DIR."

# ============================
# Paso 3: Propagación de Redes a partir de los genes semilla
# ============================

echo "Iniciando la propagación de redes a partir de los genes semilla (WNT10B, SEM1, WT1, PAX6, NF1)..."

INPUT_FILE="$RESULTS_DIR/homo_sapiens_string_ppi_filtered.tsv"
OUTPUT_FILE="$RESULTS_DIR/subgraph_diamond.tsv" 
SEED_GENES="WNT10B SEM1 WT1 PAX6 NF1"
THRESHOLD=400

# Ejecutar el script de propagación de redes
python3 "$SCRIPT_DIR/Network_propagation.py" \
    --input "$INPUT_FILE" \
    --output "$OUTPUT_FILE" \
    --seed_genes $SEED_GENES \
    --threshold $THRESHOLD

# Verificar si la propagación fue exitosa
if [ $? -ne 0 ]; then
    echo "Error en la propagación de redes. Abortando."
    exit 1
fi

# ============================
# Paso 4: Clusterización de la Red
# ============================

echo "Iniciando la clusterización de la red expandida..."

# Ejecutar el script R para la clusterización
Rscript "$SCRIPT_DIR/clustering.R"

echo "Proceso de clusterización completado. Los resultados están en la carpeta $RESULTS_DIR."

# ============================
# Paso 5: Representación del Grafo
# ============================

echo "Generando la representación gráfica de la red..."

# Ejecutar el script de representación del grafo en Python
python3 "$SCRIPT_DIR/graph_representation.py" --input "$RESULTS_DIR/subgraph_diamond.tsv" --output "$RESULTS_DIR"

# Verificar si la representación fue exitosa
if [ $? -ne 0 ]; then
    echo "Error en la representación gráfica de la red. Abortando."
    exit 1
fi

# ============================
# Paso 6: Análisis de Enriquecimiento Funcional
# ============================

echo "Iniciando el análisis de enriquecimiento funcional (GO)..."

# Ejecutar el script de análisis de enriquecimiento funcional
python3 "$SCRIPT_DIR/enriquecimiento_cluster.py" \
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

#!/bin/bash

# ============================
# Configuración del entorno
# ============================

# Definir variables
PROTEIN_LINKS="9606.protein.links.full.v12.0.txt.gz"
ALIASES_FILE="9606.protein.aliases.v12.0.txt.gz"
GO_FILE="go-basic.obo"
SEED_GENES=("PAX6" "FOXC1" "WT1" "COL4A1" "PITX2")  # Genes semilla

# Directorio donde están los scripts Python
SCRIPT_DIR="."

# ============================
# Paso 1: Propagación de la Red
# ============================

echo "Iniciando la propagación de redes a partir de los genes semilla..."

# Ejecutar el script de propagación de redes
python3 "$SCRIPT_DIR"/propagation_network.py \
    --protein-links "$PROTEIN_LINKS" \
    --aliases "$ALIASES_FILE" \
    --go "$GO_FILE" \
    --seed-genes "${SEED_GENES[@]}"  # Pasar los genes semilla al script

# Verificar si el paso de propagación fue exitoso
if [ $? -ne 0 ]; then
    echo "Error en la propagación de redes. Abortando."
    exit 1
fi

# ============================
# Paso 2: Clusterización de la Red
# ============================

echo "Iniciando la clusterización de la red expandida..."

# Ejecutar el script de clusterización (puedes usar k-means o community detection)
python3 "$SCRIPT_DIR"/network_clustering.py \
    --input-network expanded_network.gml  # Ajusta si usas otro formato de archivo

# Verificar si la clusterización fue exitosa
if [ $? -ne 0 ]; then
    echo "Error en la clusterización de la red. Abortando."
    exit 1
fi

# ============================
# Paso 3: Análisis de Enriquecimiento Funcional
# ============================

echo "Iniciando el análisis de enriquecimiento funcional (GO)..."

# Ejecutar el script de análisis de enriquecimiento funcional
python3 "$SCRIPT_DIR"/enrichment_analysis.py \
    --network expanded_network.gml \
    --go-file "$GO_FILE"

# Verificar si el análisis de enriquecimiento fue exitoso
if [ $? -ne 0 ]; then
    echo "Error en el análisis de enriquecimiento funcional. Abortando."
    exit 1
fi

# ============================
# Finalización
# ============================

echo "¡El proceso completo se ha ejecutado con éxito!"

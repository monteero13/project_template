#!/bin/bash

# ============================
# Instalación de dependencias de Python
# ============================

set -e  # Termina el script si ocurre un error

PROJECT_DIR="$(cd "$(dirname "$0")"/.. && pwd)"  # Directorio raíz del proyecto
SCRIPT_DIR="$PROJECT_DIR/code"  # Directorio donde están los scripts

echo "Iniciando la instalación de dependencias..."

# Verificar si Python está instalado
if ! command -v python3 &>/dev/null; then
  echo "Python no está instalado. Por favor, instala Python antes de continuar."
  exit 1
fi

# Actualizar pip
echo "Actualizando pip..."
pip install --upgrade pip

# Instalar dependencias desde requirements.txt
echo "Instalando dependencias de Python..."
pip install -r "$SCRIPT_DIR/requirements.txt"

echo "Dependencias instaladas con éxito."

# Verificar si R está instalado
if ! command -v R &> /dev/null
then
    echo "R no está instalado. Por favor, instálalo antes de continuar."
    exit 1
fi

# Instalar paquetes R necesarios si no están instalados
echo "Instalando paquetes R necesarios..."
Rscript -e "if (!require(igraph)) install.packages('igraph', repos='http://cran.rstudio.com/')" 
Rscript -e "if (!require(linkcomm)) install.packages('linkcomm', repos='http://cran.rstudio.com/')"

echo "Paquetes R instalados con éxito."

# Limpieza de caché de Python (opcional)
# python -m pip cache purge

# Imprimir mensaje de finalización
echo "Instalación de dependencias completada con éxito."


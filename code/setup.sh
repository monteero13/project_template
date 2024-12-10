#!/bin/bash

# ============================
# Instalación de depedencias de python
# ============================

set -e  # Termina el script si ocurre un error

echo "Instalando dependencias..."
pip install --upgrade pip
pip install -r "$SCRIPT_DIR/requirements.txt"


# Verificar si Python está instalado
if ! command -v python &>/dev/null; then
  echo "Python no está instalado. Por favor, instala Python antes de continuar."
  exit 1
fi


# Verificar si R está instalado
if ! command -v R &> /dev/null
then
    echo "R no está instalado. Por favor, instálalo antes de continuar."
    exit 1
fi


# Instalar paquetes R necesarios si no están instalados
Rscript -e "if (!require(igraph)) install.packages('igraph', repos='http://cran.rstudio.com/')" 
Rscript -e "if (!require(linkcomm)) install.packages('linkcomm', repos='http://cran.rstudio.com/')"


# Limpieza de caché de Python
#python -m pip cache purge


# Imprimir mensaje de finalización
echo "Instalación de dependencias completada con éxito."
setup.sh

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

# Verificar si pip está instalado
if ! command -v pip3 &>/dev/null; then
    echo "pip no está instalado. Instalando pip..."
    # Para sistemas basados en Debian/Ubuntu
    if [ -f /etc/debian_version ]; then
        sudo apt-get update
        sudo apt-get install -y python3-pip
    # Para sistemas basados en RedHat (CentOS, Fedora, etc.)
    elif [ -f /etc/redhat-release ]; then
        sudo yum install -y python3-pip
    else
        echo "No se pudo detectar el gestor de paquetes adecuado. Instale pip manualmente."
        exit 1
    fi
fi

# Actualizar pip
echo "Actualizando pip..."
pip3 install --upgrade pip

# Instalar dependencias desde requirements.txt
echo "Instalando dependencias de Python..."
pip3 install -r "$SCRIPT_DIR/requirements.txt"

# Verificar si curl está instalado, si no, instalarlo
if ! command -v curl &>/dev/null; then
    echo "curl no está instalado. Instalando curl..."
    # Detectar el gestor de paquetes y usar el adecuado
    if [ -f /etc/debian_version ]; then
        # Para sistemas basados en Debian (Ubuntu, etc.)
        sudo apt-get update && sudo apt-get install -y curl
    elif [ -f /etc/redhat-release ]; then
        # Para sistemas basados en RedHat (CentOS, Fedora, etc.)
        sudo yum install -y curl
    else
        echo "No se pudo detectar un gestor de paquetes adecuado. Instale curl manualmente."
        exit 1
    fi
fi

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


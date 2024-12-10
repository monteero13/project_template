# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 21:33:24 2024

@author: usuario
"""

import pandas as pd

# Cargar los datos de los archivos
fenotipos_enfermedades = pd.read_csv('../results/fenotipos_de_enfermedades.csv')
genes_enriquecimiento = pd.read_csv('../results/genes_cluster4_enriquecimiento_hpo.csv')

# Realizar el cruce entre los dos archivos
# Nos enfocamos en el "phenotype_id" en fenotipos_enfermedades y "term" en genes_enriquecimiento
fenotipos_comunes = fenotipos_enfermedades.merge(
    genes_enriquecimiento, 
    left_on='phenotype_id', 
    right_on='term', 
    how='inner'
)

# Seleccionar las columnas relevantes para el archivo de salida
resultado = fenotipos_comunes[['phenotype_id', 'phenotype_name', 'disease_id', 'disease_name']]

# Guardar el resultado en un archivo CSV
output_path = 'fenotipos_comunes.csv'
resultado.to_csv(output_path, index=False)

print(f"Archivo creado con éxito en: {output_path}")

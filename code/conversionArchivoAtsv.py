# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 00:09:49 2024

@author: alexa
"""

import pandas as pd

# Cargar el archivo de texto (usando '\t' como delimitador de tabuladores)
data = pd.read_csv('C:./genes_for_HP_0000526', delimiter='\t')

# Guardar el archivo como .tsv (separado por tabuladores)
data.to_csv('genes.tsv', sep='\t', index=False)

print("Archivo convertido exitosamente a .tsv")
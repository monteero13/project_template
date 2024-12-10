import pandas as pd
import argparse

# Función para realizar el cruce entre los archivos de fenotipos y enriquecimiento de genes
def cruzar_fenotipos(archivo_fenotipos, archivo_enriquecimiento, output_path):
    # Cargar los datos de los archivos
    fenotipos_enfermedades = pd.read_csv(archivo_fenotipos)
    genes_enriquecimiento = pd.read_csv(archivo_enriquecimiento)

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
    resultado.to_csv(output_path, index=False)

    print(f"Archivo creado con éxito en: {output_path}")

def main():
    # Usamos argparse para recibir los parámetros desde la línea de comandos
    parser = argparse.ArgumentParser(description="Cruce de fenotipos y enriquecimiento de genes para crear un archivo de resultados.")
    
    parser.add_argument("--fenotipos", type=str, required=True, help="Ruta al archivo de fenotipos y enfermedades (por ejemplo, fenotipos_de_enfermedades.csv).")
    parser.add_argument("--genes_enriquecimiento", type=str, required=True, help="Ruta al archivo de enriquecimiento de genes (por ejemplo, genes_cluster_enriquecimiento_hpo.csv).")
    parser.add_argument("--output", type=str, required=True, help="Ruta para guardar el archivo de salida (por ejemplo, fenotipos_comunes.csv).")

    args = parser.parse_args()

    # Llamamos a la función para hacer el cruce de los archivos
    cruzar_fenotipos(args.fenotipos, args.genes_enriquecimiento, args.output)

if __name__ == "__main__":
    main()

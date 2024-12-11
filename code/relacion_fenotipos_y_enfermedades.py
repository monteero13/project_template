import pandas as pd
import argparse

def main():
    # Usamos argparse para recibir los parámetros desde la línea de comandos
    parser = argparse.ArgumentParser(description="Cruce de fenotipos y enriquecimiento de genes para crear un archivo de resultados.")
    
    parser.add_argument("--fenotipos", type=str, required=True, help="Ruta al archivo de fenotipos y enfermedades (por ejemplo, fenotipos_de_enfermedades.csv).")
    parser.add_argument("--genes_enriquecimiento", type=str, required=True, help="Ruta al archivo de enriquecimiento de genes (por ejemplo, genes_cluster4_enriquecimiento_hpo.csv).")
    parser.add_argument("--output", type=str, required=True, help="Ruta para guardar el archivo de salida (por ejemplo, fenotipos_comunes.csv).")

    args = parser.parse_args()

    # Leer los datos de los archivos con manejo de líneas problemáticas
    try:
        fenotipos_enfermedades = pd.read_csv(args.fenotipos)  # Ignorar líneas problemáticas
        genes_enriquecimiento = pd.read_csv(args.genes_enriquecimiento, sep='\t')  # Ignorar líneas problemáticas
    except Exception as e:
        print(f"Error al leer los archivos: {e}")
        return

    # Realizar el cruce entre los dos archivos
    try:
        fenotipos_comunes = fenotipos_enfermedades.merge(
            genes_enriquecimiento, 
            left_on='phenotype_id', 
            right_on='Term', 
            how='inner'
        )
    except Exception as e:
        print(f"Error al realizar el cruce: {e}")
        return

    # Seleccionar las columnas relevantes para el archivo de salida
    resultado = fenotipos_comunes[['phenotype_id', 'phenotype_name', 'disease_id', 'disease_name']]

    # Guardar el resultado en un archivo CSV
    try:
        resultado.to_csv(args.output, index=False)
        print(f"Archivo creado con éxito en: {args.output}")
    except Exception as e:
        print(f"Error al guardar el archivo de salida: {e}")

if __name__ == "__main__":
    main()

import os
import argparse
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Función para crear la representación gráfica
def crear_representacion_grafo(file_path, output_dir):
    # Cargar los datos
    data = pd.read_csv(file_path, sep="\t")
    print(data.head())  # Para verificar que los datos se cargaron correctamente

    # Crear un grafo no dirigido
    G = nx.Graph()

    # Añadir aristas al grafo (source, target, weight)
    for _, row in data.iterrows():
        G.add_edge(row['source'], row['target'], weight=row['weight'])

    print(f"Número de nodos: {G.number_of_nodes()}")
    print(f"Número de aristas: {G.number_of_edges()}")

    # Dibuja la red con pesos
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G)  # Disposición de los nodos
    nx.draw(G, pos, node_size=10, edge_color="gray", with_labels=False)

    # Mostrar pesos en las aristas (opcional)
    weights = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=weights)

    # Guardar la imagen
    output_file = os.path.join(output_dir, "graph_with_weights.png")
    plt.title("Red de interacción proteína-proteína")
    plt.savefig(output_file)
    plt.close()  # Cerrar la figura para liberar memoria
    print(f"Gráfico guardado en {output_file}")

    # Dibuja la red sin mostrar los pesos
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G)  # Calcula la posición de los nodos

    # Dibuja nodos y aristas sin etiquetas de pesos
    nx.draw(
        G,
        pos,
        node_size=50,          # Tamaño de los nodos
        node_color="skyblue",  # Color de los nodos
        edge_color="gray",     # Color de las aristas
        with_labels=False      # Sin etiquetas para nodos
    )

    # Guardar la imagen
    output_file_no_weights = os.path.join(output_dir, "graph_no_weights.png")
    plt.title("Red de interacción proteína-proteína (sin pesos)", fontsize=16)
    plt.savefig(output_file_no_weights)
    plt.close()
    print(f"Gráfico sin pesos guardado en {output_file_no_weights}")

# Función principal
def main():
    # Usar argparse para gestionar los parámetros
    parser = argparse.ArgumentParser(description="Generar la representación gráfica de un grafo de interacciones proteína-proteína.")
    parser.add_argument("--input", type=str, required=True, help="Ruta del archivo de entrada con las interacciones (subgraph_diamond.tsv).")
    parser.add_argument("--output", type=str, required=True, help="Directorio donde se guardarán los gráficos generados.")
    args = parser.parse_args()

    # Crear la carpeta de salida si no existe
    os.makedirs(args.output, exist_ok=True)

    # Llamar a la función para crear la representación gráfica
    crear_representacion_grafo(args.input, args.output)

if __name__ == "__main__":
    main()

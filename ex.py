import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import load_iris

# Cargar un conjunto de datos de ejemplo (Iris)
data = load_iris()
X = data.data  # Datos originales (características)
y = data.target  # Etiquetas (clases)

# Normalizar los datos para que tengan media cero y varianza unitaria
scaler = StandardScaler()
X_normalized = scaler.fit_transform(X)

# Calcular la matriz de covarianza
cov_matrix = np.cov(X_normalized, rowvar=False)

# Calcular autovectores y autovalores
eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)

# Ordenar autovectores y autovalores en orden descendente
eigen_pairs = [(eigenvalues[i], eigenvectors[:, i]) for i in range(len(eigenvalues))]
eigen_pairs.sort(reverse=True, key=lambda x: x[0])

# Obtener las componentes principales (autovectores) para reducir la dimensionalidad
num_components = 2  # Número de componentes principales que deseamos mantener
top_eigen_pairs = eigen_pairs[:num_components]
top_eigenvalues = [pair[0] for pair in top_eigen_pairs]
top_eigenvectors = [pair[1] for pair in top_eigen_pairs]

# Transformar los datos originales en el nuevo espacio de características
transformed_data = np.dot(X_normalized, np.array(top_eigenvectors).T)

# Imprimir las componentes principales y sus autovalores
print("Componentes principales:")
for i, eigenvector in enumerate(top_eigenvectors):
    print(f"Componente {i + 1}: {eigenvector} (Autovalor: {top_eigenvalues[i]})")

# Imprimir los datos transformados en el nuevo espacio de características
print("\nDatos transformados:")
print(transformed_data)

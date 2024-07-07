import cv2
import numpy as np
import os
import shutil

def crear_mascara_roja(image):
    hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
    lower_red1 = np.array([0, 50, 50])
    upper_red1 = np.array([10, 255, 255])
    lower_red2 = np.array([170, 50, 50])
    upper_red2 = np.array([180, 255, 255])
    mask1 = cv2.inRange(hsv, lower_red1, upper_red1)
    mask2 = cv2.inRange(hsv, lower_red2, upper_red2)
    mask = cv2.bitwise_or(mask1, mask2)
    return mask

def extraer_keypoints_y_descriptores(image, mask):
    orb = cv2.ORB_create()
    keypoints, descriptores = orb.detectAndCompute(image, mask)
    return keypoints, descriptores

def comparar_keypoints(desc1, desc2):
    bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=True)
    matches = bf.match(desc1, desc2)
    matches = sorted(matches, key=lambda x: x.distance)
    return matches

def cargar_imagenes(directorio):
    imagenes = []
    for nombre_archivo in os.listdir(directorio):
        ruta = os.path.join(directorio, nombre_archivo)
        if os.path.isfile(ruta):
            imagen = cv2.imread(ruta)
            imagenes.append((nombre_archivo, imagen))
    return imagenes

def detectar_circulos(image, mask):
    masked_image = cv2.bitwise_and(image, image, mask=mask)
    gray = cv2.cvtColor(masked_image, cv2.COLOR_BGR2GRAY)
    gray = cv2.medianBlur(gray, 5)
    circles = cv2.HoughCircles(gray, cv2.HOUGH_GRADIENT, dp=1.2, minDist=50, param1=50, param2=30, minRadius=0, maxRadius=0)
    if circles is not None:
        circles = np.uint16(np.around(circles))
        return len(circles[0, :])
    return 0

# Ruta a la imagen de referencia
ruta_imagen_referencia = '/home/alex/StochasticSlipGenerator/bin/python_codes/coupling_ia.png'
imagen_referencia = cv2.imread(ruta_imagen_referencia)
mascara_referencia = crear_mascara_roja(imagen_referencia)
keypoints_referencia, descriptores_referencia = extraer_keypoints_y_descriptores(imagen_referencia, mascara_referencia)

# Procesar las carpetas de 8.1 a 9.3
for i in range(81, 94):
    version = i / 10.0
    directorio_imagenes = f'/home/alex/StochasticSlipGenerator/Output_data/Simulation_{version}_1000_coupling/img/'

    # Cargar las imágenes a comparar
    imagenes = cargar_imagenes(directorio_imagenes)

    # Comparar las imágenes y almacenar los resultados
    resultados = []
    for nombre, imagen in imagenes:
        mascara = crear_mascara_roja(imagen)
        keypoints, descriptores = extraer_keypoints_y_descriptores(imagen, mascara)
        if descriptores is not None and len(descriptores) > 0:
            matches = comparar_keypoints(descriptores_referencia, descriptores)
            similitud = len(matches)
            num_circulos = detectar_circulos(imagen, mascara)
            resultados.append((nombre, similitud, num_circulos))

    # Ordenar los resultados primero por similitud y luego por el número de círculos (ambos de mayor a menor)
    resultados = sorted(resultados, key=lambda x: (x[1], x[2]), reverse=True)

    # Crear la carpeta img_filtered si no existe
    directorio_filtrado = f'/home/alex/StochasticSlipGenerator/Output_data/Simulation_{version}_1000_coupling/img_filtered/'
    if not os.path.exists(directorio_filtrado):
        os.makedirs(directorio_filtrado)

    # Imprimir los primeros 10 resultados y copiar las imágenes a la nueva carpeta
    for nombre, similitud, num_circulos in resultados[:10]:
        print(f'Imagen: {nombre}, Similitud: {similitud}, Círculos: {num_circulos}')
        ruta_origen = os.path.join(directorio_imagenes, nombre)
        ruta_destino = os.path.join(directorio_filtrado, nombre)
        shutil.copy(ruta_origen, ruta_destino)

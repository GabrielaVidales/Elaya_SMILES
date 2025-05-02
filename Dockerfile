# Usa una imagen base con soporte para ciencia de datos y RDKit
FROM python:3.10-slim

# Instala dependencias del sistema necesarias para RDKit, OpenBabel y Auto3D
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    libgl1 \
    libglib2.0-0 \
    libxrender1 \
    libxext6 \
    libsm6 \
    libboost-all-dev \
    && apt-get clean

# Instala OpenBabel manualmente
RUN pip install openbabel-wheel

# Instala Auto3D y RDKit desde wheels
RUN pip install --upgrade pip \
 && pip install \
    rdkit \
    auto3d \
    flask \
    flask-cors \
    dscribe \
    ase \
    networkx \
    py3Dmol

# Crea directorio de trabajo
WORKDIR /app

# Copia todos los archivos del proyecto
COPY . .

# Expone el puerto de Flask
EXPOSE 5000

# Comando para iniciar la app
CMD ["python", "app.py"]

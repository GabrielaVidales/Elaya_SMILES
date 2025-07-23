# Imagen base con Python
FROM python:3.10-slim

# Instala dependencias del sistema necesarias
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

# Establecer directorio de trabajo
WORKDIR /app

# Copiar archivos del proyecto
COPY . .

# Instalar pip actualizado
RUN pip install --upgrade pip

# Instalar todas las dependencias necesarias
RUN pip install \
    flask \
    flask-cors \
    rdkit \
    openbabel-wheel \
    auto3d \
    torch \
    pandas \
    py3Dmol \
    dscribe \
    ase \
    networkx \
    ipython

EXPOSE 5000

# Comando de inicio de la aplicación
CMD ["python", "app.py"]

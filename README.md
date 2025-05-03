# 🧪 ELAYA-Smiles: Herramienta de Conversión y Análisis Molecular

![Logo Elaya](Logo_Elaya.jpg)

**ELAYA** (_Energy-minimized Linear-to-structure Atom Yielding Algorithm_) es una aplicación web para la conversión de representaciones moleculares lineales (SMILES) a estructuras tridimensionales optimizadas, permitiendo su visualización y análisis comparativo. Desarrollada como una herramienta de apoyo en química computacional, esta plataforma integra múltiples motores de conversión y métodos de análisis estructural.

---

## 🚀 Funcionalidades principales

- ✅ Conversión de cadenas SMILES a estructuras 3D con **RDKit**, **OpenBabel**, **NetworkX** o **Auto3D**.
- ✅ Visualización interactiva con **Py3Dmol**.
- ✅ Análisis de similitud molecular mediante **Tanimoto**, **SOAP** y **Valle-Oganov**.
- ✅ Carga individual o por archivo.
- ✅ Exportación de estructuras en formato `.xyz`.
- ✅ Interfaz moderna y amigable.

---

## 🧰 Tecnologías utilizadas

| Categoría             | Herramientas                                       |
|-----------------------|----------------------------------------------------|
| Backend               | Flask, Flask-CORS                                  |
| Química computacional | RDKit, OpenBabel, Auto3D, ASE, dscribe             |
| Visualización         | Py3Dmol, Chart.js                                  |
| Frontend              | HTML5, CSS3, JavaScript                            |
| Infraestructura       | Docker, Render                                     |
| Machine Learning      | PyTorch (Auto3D dependency)                        |

---

## 🗂️ Estructura del proyecto

```
ElayaSmiles/
├── app.py                # Servidor Flask principal
├── elaya_smiles.py       # Lógica de conversión molecular
├── index.html            # Página web principal
├── styles.css            # Estilos personalizados
├── app.js                # Lógica del frontend
├── Dockerfile            # Imagen para despliegue
├── render.yaml           # Configuración Render (modo Docker)
├── requirements.txt      # Dependencias del proyecto
└── assets/               # Archivos estáticos (logo, íconos)
```

---

## ⚙️ Instalación local (modo desarrollador)

```bash
git clone https://github.com/TU-USUARIO/ElayaSmiles.git
cd ElayaSmiles
pip install -r requirements.txt
python app.py
```

> Asegúrate de tener instalados RDKit, OpenBabel, Auto3D y sus dependencias (preferiblemente mediante Docker o Conda).

---

## 🐳 Despliegue en Render (Docker)

1. Sube el repositorio a GitHub con `Dockerfile` y `render.yaml`.
2. Crea un Web Service en [Render](https://render.com).
3. Selecciona “Deploy with Docker”.
4. Verifica que el servidor Flask esté configurado como:  
```python
app.run(host="0.0.0.0", port=5000)
```

---

## 🧪 Ejemplo de uso

### Conversión SMILES a 3D

1. Ingresa una cadena SMILES como: `C1=CC=CC=C1`.
2. Elige el método de conversión (RDKit, OpenBabel, etc.).
3. Visualiza y descarga el archivo `.xyz` generado.

### Análisis de Similitud

1. Selecciona dos moléculas convertidas.
2. Elige el método de análisis.
3. Visualiza el resultado numérico y gráfico.

---

## 👨‍🔬 Créditos

Proyecto desarrollado por **Gabriela Vidales** en colaboración con **Theochem Mérida - CINVESTAV**, como una propuesta para facilitar el análisis molecular estructural desde el entorno web.

---

## 📄 Licencia

Este proyecto está bajo licencia **MIT**. Consulta el archivo `LICENSE` para más información.

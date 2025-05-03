# 🧪 ELAYA-Smiles: Herramienta de Conversión y Análisis Molecular

![Logo Elaya](Logo_Elaya.jpg)

**ELAYA** (_Energy-minimized Linear-to-structure Atom Yielding Algorithm_) es una aplicación web para la conversión de representaciones moleculares lineales (SMILES) a estructuras tridimensionales optimizadas, permitiendo su visualización y análisis comparativo. Desarrollada como una herramienta de apoyo en química computacional, esta plataforma integra múltiples motores de conversión y métodos de análisis estructural.

---

## 🌐 Propósito y relevancia

La representación lineal SMILES ha sido ampliamente utilizada por su simplicidad sintáctica y compatibilidad con bases de datos químicas. Sin embargo, su utilidad para la simulación, visualización tridimensional y predicción molecular depende críticamente de su transformación en geometrías 3D realistas.

**ELAYA-Smiles automatiza y democratiza este proceso** al integrar distintos motores de conversión y algoritmos de optimización estructural, permitiendo además el análisis de similitud molecular desde múltiples enfoques.

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

## ⚗️ Características científicas destacadas

- **Conversión multi-método SMILES → 3D**: soporta RDKit, OpenBabel, Auto3D y NetworkX.
- **Optimización geométrica por energía** (Auto3D) con soporte para desactivación de GPU.
- **Análisis estructural de similitud molecular** mediante:
  - Índices topológicos (Tanimoto sobre fingerprints)
  - Descriptores atómicos y SOAP (Smooth Overlap of Atomic Positions)
  - Distancias en el espacio de configuraciones (Valle-Oganov)
- **Visualización molecular interactiva** mediante Py3Dmol.
- **Exportación en formato XYZ para simulaciones posteriores.**

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

## 👩‍🔬 Desarrollo

Este sistema ha sido desarrollado por **Gabriela Vidales** como parte de un esfuerzo por integrar herramientas de código abierto con flujos de trabajo reproducibles para la representación y análisis molecular tridimensional.

El proyecto toma inspiración tanto de la necesidad académica en el aula como de las exigencias del laboratorio computacional.

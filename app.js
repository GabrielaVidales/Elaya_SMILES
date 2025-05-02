// Configuración global
const API_BASE_URL = window.location.origin + '/api';
let molecules = [];
let currentVisualizations = [];
let currentXYZ = '';
let currentMoleculeIndex = 0;
let rotationSpeed = 0.8;
let isRotationPaused = false;
let rotationAnimationId = null;
let selectedMolecules = new Set();

// Elementos del DOM
let visualizationContainer, viewerContainer, moleculeTitle, moleculeControls;
let prevMoleculeBtn, nextMoleculeBtn, rotationSpeedInput;

// Event Listeners
document.addEventListener('DOMContentLoaded', function() {
    console.log("DOM fully loaded and parsed");

    // Inicializar elementos del DOM
    visualizationContainer = document.querySelector('.visualization-container');
    viewerContainer = document.getElementById('viewer-container');
    moleculeTitle = document.getElementById('molecule-title');
    moleculeControls = document.getElementById('molecule-controls');
    prevMoleculeBtn = document.getElementById('prev-molecule');
    nextMoleculeBtn = document.getElementById('next-molecule');
    rotationSpeedInput = document.getElementById('rotation-speed');

    // Configurar rango de velocidad (0-2)
    rotationSpeedInput.min = 0;
    rotationSpeedInput.max = 3;
    rotationSpeedInput.step = 0.3;
    rotationSpeedInput.value = 0.8;
    
    // Configuración inicial
    checkBackendConnection();
    toggleInputMethod();
    toggleForceField();

    // Event listeners
    document.querySelectorAll('input[name="input-method"]').forEach(radio => {
        radio.addEventListener('change', toggleInputMethod);
    });

    document.getElementById('conversion-method').addEventListener('change', toggleForceField);
    document.getElementById('load-smiles').addEventListener('click', loadSmiles);
    document.getElementById('convert-molecules').addEventListener('click', convertMolecules);
    document.getElementById('download-xyz').addEventListener('click', downloadXYZ);
    document.getElementById('download-image').addEventListener('click', downloadImage);
    document.getElementById('view-file').addEventListener('click', viewFileContent);

    // Nuevos event listeners
    if (prevMoleculeBtn) prevMoleculeBtn.addEventListener('click', () => navigateMolecules(-1));
    if (nextMoleculeBtn) nextMoleculeBtn.addEventListener('click', () => navigateMolecules(1));
    if (rotationSpeedInput) rotationSpeedInput.addEventListener('input', (e) => changeRotationSpeed(e.target.value));

    // Inicializar 3Dmol después de que se cargue la librería
    init3DMol();
});

function init3DMol() {
    // Verificar si 3Dmol está cargado
    if (window.$3Dmol) {
        console.log("3Dmol está cargado correctamente");
    } else {
        console.error("3Dmol no se cargó correctamente");
        // Podrías recargar la librería aquí si es necesario
    }
}

// Funciones principales
async function checkBackendConnection() {
    try {
        console.log("Verificando conexión con el backend...");
        const response = await fetch(API_BASE_URL);
        if (!response.ok) {
            console.error('Error en el backend:', await response.text());
            showError('Error en el servidor backend');
        } else {
            const data = await response.json();
            console.log('Backend conectado:', data);
            showSuccess('Conexión con el backend establecida');
        }
    } catch (error) {
        console.error('Conexión fallida con el backend:', error);
        showError('No se pudo conectar al servidor. ¿Está ejecutándose el backend?');
    }
}

async function convertMolecules() {
    try {
        if (!molecules.length) {
            throw new Error('No hay moléculas cargadas para convertir');
        }

        // Si hay múltiples moléculas, mostrar el modal de selección
        if (molecules.length > 1) {
            const selected = await showMoleculeSelectionModal();
            if (!selected || selected.size === 0) return;
            
            // Convertir las moléculas seleccionadas
            for (const index of selected) {
                currentMoleculeIndex = index;
                await convertSingleMolecule();
            }
            
            // Mostrar la primera molécula seleccionada
            const firstSelected = [...selected][0];
            currentMoleculeIndex = firstSelected;
            await convertSingleMolecule();
        } else {
            // Solo una molécula
            currentMoleculeIndex = 0;
            await convertSingleMolecule();
        }
    } catch (error) {
        console.error("Error en convertMolecules:", error);
        showError(`Error al convertir molécula: ${error.message}`);
    }
}

async function convertSingleMolecule() {
    const method = document.getElementById('conversion-method').value;
    const forceField = document.getElementById('force-field').value;
    const molecule = molecules[currentMoleculeIndex];

    const response = await fetch(`${API_BASE_URL}/convert`, {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify({
            smiles: molecule.smiles,
            identifier: molecule.id,
            method: method,
            force_field: forceField
        })
    });

    if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.error || 'Error en la conversión');
    }

    const data = await response.json();
    if (!data.xyz) {
        throw new Error('No se recibieron datos XYZ del servidor');
    }

    currentXYZ = data.xyz;
    updateXYZDisplay();
    render3DMolecule(data.xyz, molecule.smiles);
    showSuccess(`Molécula ${molecule.id} convertida exitosamente`);
}


function updateXYZDisplay() {
    const xyzDisplay = document.getElementById('xyz-display');
    if (xyzDisplay) {
        xyzDisplay.textContent = currentXYZ;
        xyzDisplay.style.display = 'block';
    }
}

async function selectMoleculeToConvert() {
    return new Promise(resolve => {
        const modal = document.createElement('div');
        modal.style.position = 'fixed';
        modal.style.top = '0';
        modal.style.left = '0';
        modal.style.width = '100%';
        modal.style.height = '100%';
        modal.style.backgroundColor = 'rgba(0,0,0,0.7)';
        modal.style.display = 'flex';
        modal.style.justifyContent = 'center';
        modal.style.alignItems = 'center';
        modal.style.zIndex = '1000';
        
        const content = document.createElement('div');
        content.style.backgroundColor = '#12122B';
        content.style.padding = '20px';
        content.style.borderRadius = '10px';
        content.style.maxWidth = '500px';
        content.style.width = '80%';
        
        content.innerHTML = `
            <h3 style="color: #98ebc4; margin-bottom: 15px;">Seleccione la molécula a convertir</h3>
            <div style="max-height: 300px; overflow-y: auto; margin-bottom: 15px;">
                ${molecules.map((mol, index) => `
                    <div style="padding: 8px; border-bottom: 1px solid #0f6975; cursor: pointer;"
                         onclick="document.getElementById('modal-selected-index').value = ${index};">
                        <strong>${mol.id}:</strong> ${mol.smiles}
                    </div>
                `).join('')}
            </div>
            <input type="hidden" id="modal-selected-index" value="-1">
            <div style="display: flex; justify-content: space-between;">
                <button onclick="this.parentNode.parentNode.parentNode.removeChild(this.parentNode.parentNode); resolve(null);"
                        style="background: #FF6B6B; padding: 8px 15px; border: none; border-radius: 5px; cursor: pointer;">
                    Cancelar
                </button>
                <button onclick="const selectedIndex = parseInt(document.getElementById('modal-selected-index').value); 
                                 this.parentNode.parentNode.parentNode.removeChild(this.parentNode.parentNode); 
                                 resolve(selectedIndex >= 0 ? selectedIndex : null);"
                        style="background: #98ebc4; padding: 8px 15px; border: none; border-radius: 5px; cursor: pointer;">
                    Convertir
                </button>
            </div>
        `;
        
        modal.appendChild(content);
        document.body.appendChild(modal);
    });
}

function render3DMolecule(xyzContent, smiles = '') {
    // Limpiar visualización anterior
    if (rotationAnimationId) {
        cancelAnimationFrame(rotationAnimationId);
        rotationAnimationId = null;
    }
    
    viewerContainer.innerHTML = '';
    moleculeTitle.textContent = `SMILES: ${smiles}`;
    moleculeTitle.style.display = 'block';
    
    // Mostrar controles si hay múltiples moléculas
    moleculeControls.style.display = 'flex';

    if (window.$3Dmol && window.$3Dmol.createViewer) {
        const viewer = window.$3Dmol.createViewer(viewerContainer, {
            width: viewerContainer.clientWidth,
            height: viewerContainer.clientHeight,
            backgroundColor: 'black'
        });

        viewer.addModel(xyzContent, "xyz", {
            keepH: true,
            noComputeSecondaryStructure: true
        });

        viewer.setStyle({}, {
            stick: {radius: 0.15},
            sphere: {scale: 0.25}
        });
        
        viewer.zoomTo();
        viewer.render();
        
        // Animación de rotación mejorada
        const animate = function() {
            if (!isRotationPaused) {
                viewer.rotate(rotationSpeed, {x:0, y:1, z:0});
                viewer.render();
            }
            rotationAnimationId = requestAnimationFrame(animate);
        };
        animate();
        
        currentVisualizations = [viewer];
    } else {
        viewerContainer.innerHTML = `
            <div class="placeholder error">
                Error: La biblioteca 3Dmol no está disponible
            </div>
        `;
    }
}

function navigateMolecules(direction) {
    const newIndex = currentMoleculeIndex + direction;
    if (newIndex >= 0 && newIndex < molecules.length) {
        currentMoleculeIndex = newIndex;
        convertMolecules();
    }
}

function toggleRotation() {
    isRotationPaused = !isRotationPaused;
    const btn = document.getElementById('pause-btn');
    if (btn) {
        btn.textContent = isRotationPaused ? 'Reanudar' : 'Pausar';
    }
}

function changeRotationSpeed(speed) {
    rotationSpeed = parseFloat(speed);
}

function downloadXYZ() {
    if (!currentXYZ) {
        showError('No hay datos XYZ para descargar');
        return;
    }

    try {
        const blob = new Blob([currentXYZ], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `molecule_${molecules[currentMoleculeIndex]?.id || 'unknown'}.xyz`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
        showSuccess('Archivo XYZ descargado');
    } catch (error) {
        console.error("Error al descargar XYZ:", error);
        showError('Error al descargar el archivo XYZ');
    }
}

function downloadImage() {
    if (!currentVisualizations.length) {
        showError('No hay visualización para descargar');
        return;
    }

    try {
        const viewer = currentVisualizations[0];
        const imgData = viewer.pngURI();
        const a = document.createElement('a');
        a.href = imgData;
        a.download = `molecule_${molecules[currentMoleculeIndex]?.id || 'unknown'}.png`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        showSuccess('Imagen descargada');
    } catch (error) {
        console.error("Error al descargar imagen:", error);
        showError('Error al descargar la imagen');
    }
}

function toggleInputMethod() {
    try {
        const method = document.querySelector('input[name="input-method"]:checked').value;
        
        document.getElementById('single-input').style.display = method === 'single' ? 'block' : 'none';
        document.getElementById('single-id').style.display = method === 'single' ? 'block' : 'none';
        document.getElementById('file-input').style.display = method === 'file' ? 'block' : 'none';
        
        // Mover el input de archivo a la izquierda
        const fileInputDiv = document.getElementById('file-input');
        if (fileInputDiv) {
            fileInputDiv.style.marginLeft = '0';
            fileInputDiv.style.marginRight = 'auto';
        }
        
        document.getElementById('view-file').style.display = method === 'file' ? 'inline-block' : 'none';
        document.getElementById('file-format-example').style.display = method === 'file' ? 'block' : 'none';
    } catch (error) {
        console.error("Error en toggleInputMethod:", error);
    }
}

async function loadSmiles() {
    try {
        console.log("Iniciando carga de SMILES...");
        const method = document.querySelector('input[name="input-method"]:checked').value;
        console.log("Método seleccionado:", method);

        if (method === 'single') {
            const smiles = document.getElementById('smiles').value.trim();
            const molId = document.getElementById('mol-id').value.trim();
            console.log("SMILES ingresado:", smiles, "ID:", molId);

            if (!smiles) {
                throw new Error('Por favor ingrese un string SMILES');
            }

            const finalId = molId || `mol_${Date.now()}`;
            molecules = [{ smiles, id: finalId }];
        } else {
            const fileInput = document.getElementById('smiles-file');
            console.log("Archivo seleccionado:", fileInput?.files[0]);

            if (!fileInput?.files?.length) {
                throw new Error('Por favor seleccione un archivo');
            }

            const text = await fileInput.files[0].text();
            molecules = [];
            console.log("Contenido del archivo:", text);

            text.split('\n').forEach(line => {
                const trimmedLine = line.trim();
                if (trimmedLine) {
                    const parts = trimmedLine.split(/\s+/);
                    if (parts.length >= 1) {
                        const smiles = parts[0];
                        const id = parts.length >= 2 ? parts[1] : `mol_${molecules.length + 1}`;
                        molecules.push({ smiles, id });
                    }
                }
            });

            if (molecules.length === 0) {
                throw new Error('El archivo no contiene moléculas válidas');
            }
        }

        console.log("Moléculas cargadas:", molecules);
        const loadConfirmation = document.getElementById('load-confirmation');
        if (loadConfirmation) {
            const moleculeInfo = molecules.map(mol => 
                `<div><strong>SMILES:</strong> ${mol.smiles}</div><div><strong>ID:</strong> ${mol.id}</div>`
            ).join('<hr>');
            
            loadConfirmation.innerHTML = `
                <div class="alert success">
                    ${molecules.length === 1 ? 'Molécula cargada exitosamente:' : `${molecules.length} moléculas cargadas exitosamente:`}
                    ${moleculeInfo}
                </div>
            `;
            loadConfirmation.style.display = 'block';
            setTimeout(() => {
                loadConfirmation.style.display = 'none';
            }, 5000);
        }

        const convertBtn = document.getElementById('convert-molecules');
        if (convertBtn) {
            convertBtn.disabled = false;
            console.log("Botón Convert habilitado");
        }
    } catch (error) {
        console.error("Error en loadSmiles:", error);
        showError(`Error al cargar SMILES: ${error.message}`);
        throw error;
    }
}

function updateLoadConfirmation() {
    if (!loadConfirmation) {
        loadConfirmation = document.getElementById('load-confirmation');
    }

    if (loadConfirmation) {
        const moleculeInfo = molecules.map(mol => 
            `<div><strong>SMILES:</strong> ${mol.smiles}</div><div><strong>ID:</strong> ${mol.id}</div>`
        ).join('<hr>');
        
        loadConfirmation.innerHTML = `
            <div class="alert success">
                ${molecules.length === 1 ? 'Molécula cargada exitosamente:' : `${molecules.length} moléculas cargadas exitosamente:`}
                ${moleculeInfo}
            </div>
        `;
        loadConfirmation.style.display = 'block';
        setTimeout(() => {
            loadConfirmation.style.display = 'none';
        }, 5000);
    }
}

function enableConversionControls() {
    const convertBtn = document.getElementById('convert-molecules');
    const mol1Select = document.getElementById('mol1');
    const mol2Select = document.getElementById('mol2');
    const compareMolSelect = document.getElementById('compare-molecule');

    if (convertBtn) convertBtn.disabled = false;
    
    if (mol1Select && mol2Select && compareMolSelect) {
        // Limpiar y llenar selects
        mol1Select.innerHTML = '';
        mol2Select.innerHTML = '';
        compareMolSelect.innerHTML = '';
        
        molecules.forEach(mol => {
            const option = document.createElement('option');
            option.value = mol.id;
            option.textContent = `${mol.id} (${mol.smiles})`;
            
            mol1Select.appendChild(option.cloneNode(true));
            mol2Select.appendChild(option.cloneNode(true));
            compareMolSelect.appendChild(option);
        });
        
        mol1Select.disabled = false;
        mol2Select.disabled = false;
        compareMolSelect.disabled = false;
        document.getElementById('calculate-similarity').disabled = false;
        document.getElementById('compare-methods').disabled = false;
    }
}

function toggleForceField() {
    const method = document.getElementById('conversion-method').value;
    const forceFieldDiv = document.getElementById('force-field-options');
    
    if (method === 'networkx' || method === 'auto3d') {
        forceFieldDiv.style.display = 'none';
    } else {
        forceFieldDiv.style.display = 'block';
    }
}

function viewFileContent() {
    const fileInput = document.getElementById('smiles-file');
    if (!fileInput?.files?.length) return;

    const file = fileInput.files[0];
    const reader = new FileReader();
    
    reader.onload = function(e) {
        alert(`Contenido del archivo:\n\n${e.target.result}`);
    };
    
    reader.readAsText(file);
}

function retryVisualization() {
    if (currentXYZ && molecules[currentMoleculeIndex]) {
        render3DMolecule(currentXYZ, molecules[currentMoleculeIndex].smiles);
    }
}

// Funciones de ayuda
function showSuccess(message) {
    console.log("Mostrando mensaje de éxito:", message);
    const alertDiv = document.createElement('div');
    alertDiv.className = 'alert success';
    alertDiv.innerHTML = `✓ ${message}`;
    document.body.appendChild(alertDiv);
    setTimeout(() => {
        alertDiv.remove();
    }, 3000);
}

function showError(message) {
    console.error("Mostrando mensaje de error:", message);
    const alertDiv = document.createElement('div');
    alertDiv.className = 'alert error';
    alertDiv.innerHTML = `⚠ ${message}`;
    document.body.appendChild(alertDiv);
    setTimeout(() => {
        alertDiv.remove();
    }, 5000);
}

async function showMoleculeSelectionModal() {
    return new Promise((resolve) => {
        const modal = document.createElement('div');
        modal.className = 'modal-overlay';
        
        modal.innerHTML = `
            <div class="modal-content">
                <h3>Seleccione los SMILES a convertir</h3>
                <div id="molecule-list" class="molecule-list"></div>
                <div class="modal-buttons">
                    <button id="modal-cancel">Cancelar</button>
                    <button id="modal-confirm">Convertir seleccionadas</button>
                </div>
            </div>
        `;
        
        const moleculeList = modal.querySelector('#molecule-list');
        molecules.forEach((mol, index) => {
            const item = document.createElement('div');
            item.className = 'molecule-item';
            item.textContent = `${mol.id}: ${mol.smiles}`;
            item.dataset.index = index;
            moleculeList.appendChild(item);
        });
        
        // Manejar clics en elementos
        modal.addEventListener('click', (e) => {
            const item = e.target.closest('.molecule-item');
            if (item) {
                item.classList.toggle('selected');
                const index = parseInt(item.dataset.index);
                if (item.classList.contains('selected')) {
                    selectedMolecules.add(index);
                } else {
                    selectedMolecules.delete(index);
                }
            }
        });
        
        // Cerrar al hacer clic fuera del contenido
        modal.addEventListener('click', (e) => {
            if (e.target === modal) {
                document.body.removeChild(modal);
                document.removeEventListener('keydown', handleKeyDown);
                resolve(null);
            }
        });

        // Botones del modal
        modal.querySelector('#modal-cancel').addEventListener('click', () => {
            document.body.removeChild(modal);
            resolve(null);
        });
        
        modal.querySelector('#modal-confirm').addEventListener('click', () => {
            document.body.removeChild(modal);
            resolve(new Set(selectedMolecules));
            selectedMolecules.clear(); // Limpiar selección después de usar
        });
        
        document.body.appendChild(modal);
    });
}

// Hacer funciones accesibles globalmente para botones HTML
window.retryVisualization = retryVisualization;
window.navigateMolecules = navigateMolecules;
window.toggleRotation = toggleRotation;
window.changeRotationSpeed = changeRotationSpeed;
window.navigateMolecules = navigateMolecules;
window.toggleRotation = toggleRotation;
window.changeRotationSpeed = changeRotationSpeed;
window.retryVisualization = retryVisualization;

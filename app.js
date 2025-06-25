console.log(">>> app.js está cargado");

// Configuración global
const API_BASE_URL = window.location.origin + '/api';
let molecules = [];
let selected = new Set();
let currentVisualizations = [];
let currentXYZ = '';
let currentMoleculeIndex = 0;
let rotationSpeed = 0.8;
let isRotationPaused = false;
let rotationAnimationId = null;
let selectedMolecules = new Set();
let visualizationContainer, viewerContainer, moleculeTitle, moleculeControls;
let prevMoleculeBtn, nextMoleculeBtn, rotationSpeedInput;

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
rotationSpeedInput.max = 5;
rotationSpeedInput.step = 0.5;
rotationSpeedInput.value = 0.8;
    
document.addEventListener("DOMContentLoaded", function () {
    console.log("DOM fully loaded and parsed");

    checkBackendConnection();
    toggleInputMethod();
    toggleForceField();
    loadSmilesHistory();
    setupConversionHandler();

    const homeLink = document.querySelector('a[href="#home"]');
    const simulatorLink = document.querySelector('a[href="#simulator"]');
    const homeSection = document.getElementById('home-info');
    const header = document.querySelector('header');
    const inputSection = document.querySelector('.input-section');
    const conversionSection = document.querySelector('.conversion-section');

    // Botón "Home"
    if (homeLink) {
        homeLink.addEventListener('click', function (e) {
            e.preventDefault();
            // Oculta todo
            document.querySelectorAll('main > section').forEach(section => section.style.display = 'none');
            // Muestra la sección home
            if (homeSection) {
                homeSection.style.display = 'block';
                setTimeout(() => {
                    const topOffset = homeSection.offsetTop - 40; // Puedes ajustar este valor
                    window.scrollTo({ top: topOffset, behavior: 'smooth' });
                }, 100);
            }
            if (header) header.style.display = 'none'; 
        });
    }

    // Botón "Simulator"
    if (simulatorLink) {
        simulatorLink.addEventListener('click', function (e) {
            e.preventDefault();
        // Ocultar todas las secciones
        document.querySelectorAll('main > section').forEach(section => section.style.display = 'none');

        // Mostrar solo las secciones necesarias del simuladorr
            if (inputSection) inputSection.style.display = 'block';
            if (conversionSection) conversionSection.style.display = 'flex';
            if (header) header.style.display = 'block'; 
            window.scrollTo({ top: 0, behavior: 'smooth' });
        });
    }

    const backgroundLink = document.querySelector('a[href="#background"]');
    const backgroundSection = document.getElementById('background-info');

    if (backgroundLink) {
    backgroundLink.addEventListener('click', function (e) {
        e.preventDefault();
        document.querySelectorAll('main > section').forEach(section => section.style.display = 'none');
        if (backgroundSection) {
        backgroundSection.style.display = 'block';
        setTimeout(() => {
            window.scrollTo({ top: backgroundSection.offsetTop - 40, behavior: 'smooth' });
        }, 100);
        }
        if (header) header.style.display = 'none';
    });
    }
});

    // Event listeners
document.querySelectorAll('input[name="input-method"]').forEach(radio => {
    radio.addEventListener('change', toggleInputMethod);
});

document.getElementById('conversion-method').addEventListener('change', toggleForceField);
document.getElementById('load-smiles').addEventListener('click', loadSmiles);
document.getElementById('convert-molecules').addEventListener('click', convertMolecules);
document.getElementById('download-xyz').addEventListener('click', downloadXYZ);
document.getElementById('download-image').addEventListener('click', downloadImage);

document.getElementById('smiles-file').addEventListener('change', function () {
    this.classList.add('touched');
    const methodRadio = document.querySelector('input[name="input-method"][value="file"]');
    if (methodRadio) methodRadio.checked = true;
    toggleInputMethod();
});

    // Nuevos event listeners
if (prevMoleculeBtn) prevMoleculeBtn.addEventListener('click', () => navigateMolecules(-1));
if (nextMoleculeBtn) nextMoleculeBtn.addEventListener('click', () => navigateMolecules(1));
if (rotationSpeedInput) rotationSpeedInput.addEventListener('input', (e) => changeRotationSpeed(e.target.value));

    // Inicializar 3Dmol después de que se cargue la librería
init3DMol();

document.querySelectorAll('input, select').forEach(el => {
    const type = el.getAttribute('type');
    
    // Para select y file usamos 'change', para text también usamos 'input'
    const eventType = (el.tagName === 'SELECT' || type === 'file') ? 'change' : 'input';
    
    el.addEventListener(eventType, () => {
        el.classList.add('touched');
    });
});

document.getElementById('load-smiles-file').addEventListener('click', () => {
    loadSmiles();
});

loadSmilesHistory();

// Mostrar automáticamente la oxitocina al iniciar por primera vez en la sesión
if (!sessionStorage.getItem('oxytocinShown')) {
    const oxytocinSmiles = "CC(C)C[C@H](NC(=O)[C@H](CC(C)C)NC(=O)CN)C(=O)N";
    document.getElementById('smiles').value = oxytocinSmiles;
    molecules = [{ smiles: oxytocinSmiles, id: oxytocinSmiles }];
    updateSmilesHistory(oxytocinSmiles);
        
    // Llama a convertSingleMolecule sin await
    convertSingleMolecule().catch(err => {
        console.error("Error al convertir oxitocina automáticamente:", err);
    });

    sessionStorage.setItem('oxytocinShown', 'true');
}

        /* Cargar por defecto la oxitocina al iniciar
    if (!molecules.length) {
        const oxytocinSmiles = "CC(C)C[C@H](NC(=O)[C@H](CC(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@H](CSSC[C@H](C(=O)N[C@H](C)C(=O)N[C@H](CCC(=O)O)C(=O)N[C@H](CO)C(=O)N1CCC[C@H]1C(=O)N[C@H](CCCCN)C(=O)N2CCCC2)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](N)C(C)C)NC(=O)[C@H](N)CO)NC(=O)[C@H](N)CC(=O)O)C(=O)O";
        document.getElementById('smiles').value = oxytocinSmiles;
        molecules = [{ smiles: oxytocinSmiles, id: "Oxitocina" }];
        updateSmilesHistory(oxytocinSmiles);  // opcional si quieres agregarla al historial
        convertSingleMolecule(); // Lanza la conversión automática
        }*/

// Efecto de rayo dorado al hacer clic en el título
const goldFlashLink = document.querySelector('.gold-flash');
const modal = document.createElement('div');
modal.className = 'modal-overlay';

if (goldFlashLink) {
    goldFlashLink.addEventListener('click', function (e) {
    e.preventDefault();
    this.classList.remove('clicked'); // reinicia si estaba activo
    void this.offsetWidth; // fuerza reflow para reiniciar animación
    this.classList.add('clicked');

    // quitar clase después de que termine la animación (~1s)
    setTimeout(() => {
            this.classList.remove('clicked');
        }, 1000); // debe coincidir con la duración del @keyframes
    });
}

console.log("Conexión con backend establecida");

let is3DEnabled = true;
let is2DEnabled = false;

// === Funciones de grabación y captura de imagen ===
const recordBtn = document.getElementById('record-btn');
const screenshotBtn = document.getElementById('screenshot-btn');
const countdownEl = document.getElementById('countdown');
const viewer = document.getElementById('viewer-container');

let mediaRecorder;
let recordedChunks = [];
let isRecording = false;

if (recordBtn && screenshotBtn) {
    recordBtn.addEventListener('click', () => {
        if (!isRecording) {
            // Mostrar cuenta regresiva y luego iniciar grabación
            countdown(3, () => {
                startRecording();
                toggleRecordingUI(true);
                isRecording = true;
            });
        } else {
            stopRecording();
            toggleRecordingUI(false);
            isRecording = false;
        }
    });

    screenshotBtn.addEventListener('click', downloadImage);
}

function toggleRecordingUI(recording) {
    const svg = recordBtn.querySelector('svg');
    svg.innerHTML = recording
        ? `<rect x="8" y="8" width="8" height="8" />`
        : `<circle cx="12" cy="12" r="8" />`;

    recordBtn.classList.toggle('recording-blink', recording);
    recordBtn.title = recording ? "Stop Recording" : "Start Recording";
}

function countdown(seconds, callback) {
    countdownEl.style.display = 'block';
    let current = seconds;
    const timer = setInterval(() => {
        countdownEl.textContent = current;
        current--;
        if (current < 0) {
            clearInterval(timer);
            countdownEl.style.display = 'none';
            callback();
        }
    }, 1000);
}

let timerInterval;
let startTime;

function startRecording() {
    const canvas = document.querySelector('#viewer-container canvas');
    if (!canvas) {
        console.error("No se encontró canvas para capturar.");
        return;
    }
    const stream = canvas.captureStream(30);
    console.log("Stream capturado desde canvas:", stream);
    
    mediaRecorder = new MediaRecorder(stream, { mimeType: 'video/webm' });

    mediaRecorder.ondataavailable = function (e) {
        if (e.data.size > 0) {
            recordedChunks.push(e.data);
        }
    };

    mediaRecorder.onstop = function () {
        clearInterval(timerInterval);
        document.getElementById('countdown').style.display = 'none';

        const blob = new Blob(recordedChunks, { type: 'video/webm' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');

        const smilesId = window.molecules?.[window.currentMoleculeIndex]?.smiles?.replace(/[^a-zA-Z0-9]/g, "_") || `video_${Date.now()}`;
        const filename = `molecule_video_${smilesId}.webm`;

        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        URL.revokeObjectURL(url);
        recordedChunks = [];
    };

    mediaRecorder.start();
    startTimer();
}

function startTimer() {
    startTime = Date.now();
    const countdownEl = document.getElementById('countdown');
    countdownEl.style.display = 'block';
    timerInterval = setInterval(() => {
        const elapsed = Math.floor((Date.now() - startTime) / 1000);
        countdownEl.textContent = formatTime(elapsed);
    }, 1000);
}

function formatTime(seconds) {
    const mins = Math.floor(seconds / 60).toString().padStart(2, '0');
    const secs = (seconds % 60).toString().padStart(2, '0');
    return `${mins}:${secs}`;
}


function stopRecording() {
    if (mediaRecorder && mediaRecorder.state !== 'inactive') {
        mediaRecorder.stop();
    }
}

document.querySelectorAll('.info-icon').forEach(icon => {
    let timer;
    const tooltipId = icon.getAttribute('data-tooltip-id');
    const tooltip = document.getElementById(tooltipId);
    const url = icon.getAttribute('data-url');

    // Manejar hover para tooltip
    icon.addEventListener('mouseenter', () => {
        timer = setTimeout(() => {
            if (tooltip) tooltip.classList.add('visible');
        }, 300);
    });

    icon.addEventListener('mouseleave', () => {
        clearTimeout(timer);
        if (tooltip) tooltip.classList.remove('visible');
    });

    // Manejar click para abrir URL
    icon.addEventListener('click', (e) => {
        e.stopPropagation(); // Evitar que el evento se propague
        if (url) window.open(url, '_blank');
    });
});


const btn3D = document.getElementById("toggle-3d");
const btn2D = document.getElementById("toggle-2d");
    
btn3D.addEventListener("click", () => toggleView("3d"));
btn2D.addEventListener("click", () => toggleView("2d"));

function toggleView(mode) {
    if (mode === "3d") {
        if (is3DEnabled && !is2DEnabled) return; // no permitir que ambos se apaguen
        is3DEnabled = !is3DEnabled;
        btn3D.classList.toggle("active", is3DEnabled);
    } else if (mode === "2d") {
        if (is2DEnabled && !is3DEnabled) return; // no permitir que ambos se apaguen
        is2DEnabled = !is2DEnabled;
        btn2D.classList.toggle("active", is2DEnabled);
    }

    updateViewerContainerLayout();
}

function updateViewerContainerLayout() {
        const viewer = document.querySelector('.visualization-container');
        const speedControl = document.querySelector(".speed-control");
        const mol = molecules[currentMoleculeIndex];
        if (!mol) {
            console.warn("Molécula no definida al cambiar vista.");
            return;
        }

        // Limpiar el contenedor
        viewer.innerHTML = '';

        // 3D + 2D
        if (is3DEnabled && is2DEnabled) {
            const viewer3D = document.createElement('div');
            viewer3D.id = 'viewer-3d';
            viewer3D.style.width = '50%';
            viewer3D.style.height = '100%';
            viewer3D.style.display = 'inline-block';
            viewer3D.style.verticalAlign = 'top';
            viewer3D.style.position = 'relative';

            const divider = document.createElement('div');
            divider.className = 'viewer-divider';
            divider.style.position = 'absolute';
            divider.style.right = '0';
            divider.style.top = '0';
            divider.style.width = '1px';
            divider.style.height = '100%';
            divider.style.backgroundColor = '#666';

            viewer3D.appendChild(divider);

            const viewer2D = document.createElement('div');
            viewer2D.id = 'viewer-2d';
            viewer2D.style.width = '50%';
            viewer2D.style.height = '100%';
            viewer2D.style.display = 'inline-block';
            viewer2D.style.verticalAlign = 'top';
            viewer2D.style.color = 'white';
            viewer2D.style.fontFamily = 'monospace';
            viewer2D.style.overflowY = 'auto';
            viewer2D.style.padding = '10px';
            const mol = molecules[currentMoleculeIndex];
            viewer2D.innerHTML = ''; // Limpia el contenedor
            const chemViewer = new Kekule.ChemWidget.Viewer(viewer2D);
            chemViewer.setPredefinedSetting('fullFunc'); // incluye zoom, hover, etc.
            setTimeout(() => {
                const interval = setInterval(() => {
                    const reader = Kekule.IO.ChemDataReaderManager.getReaderByFormat('smi');
                    if (reader) {
                        clearInterval(interval);
                        try {
                            const chemObj = reader.readData(mol.smiles, 'smi');
                            if (chemObj) {
                                chemViewer.setChemObj(chemObj);
                            } else {
                                viewer2D.textContent = "Error loading 2D structure.";
                            }
                        } catch (err) {
                            console.error("Error al cargar el SMILES para Kekule:", err);
                            viewer2D.textContent = "Error loading 2D structure.";
                        }
                    } else {
                        console.warn("Esperando que Kekule cargue el lector SMILES...");
                    }
                }, 200);
            }, 100); 
            viewer.appendChild(viewer3D);
            viewer.appendChild(viewer2D);
            const format = mol.molBlock ? 'sdf' : 'xyz';
            const modelData = mol.molBlock || mol.xyz;
            render3DMolecule(modelData, mol.smiles, format, viewer3D);

            if (speedControl) speedControl.style.display = 'flex';

        } else if (is3DEnabled) {
            const viewer3D = document.createElement('div');
            viewer3D.id = 'viewer-3d';
            viewer3D.style.width = '100%';
            viewer3D.style.height = '100%';
            viewer3D.style.position = 'relative';

            viewer.appendChild(viewer3D);

            const mol = molecules[currentMoleculeIndex];
            const format = mol.molBlock ? 'sdf' : 'xyz';
            const modelData = mol.molBlock || mol.xyz;
            render3DMolecule(modelData, mol.smiles, format, viewer3D);

            if (speedControl) speedControl.style.display = 'flex';

        } else if (is2DEnabled) {
            const viewer2D = document.createElement('div');
            viewer2D.id = 'viewer-2d';
            viewer2D.style.width = '100%';
            viewer2D.style.height = '100%';
            viewer2D.style.color = 'white';
            viewer2D.style.fontFamily = 'monospace';
            viewer2D.style.overflowY = 'auto';
            viewer2D.style.padding = '10px';
            const mol = molecules[currentMoleculeIndex];
                fetch('/api/draw2d', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ smiles: mol.smiles })
                })
                .then(response => {
                    if (!response.ok) throw new Error("Error al obtener imagen 2D");
                    return response.text();
                })
                .then(svg => {
                    viewer2D.innerHTML = svg;
                })
                .catch(err => {
                    console.error("Error cargando visualización 2D:", err);
                    viewer2D.textContent = "Error al mostrar estructura 2D.";
                });
            viewer.appendChild(viewer2D);

            if (speedControl) speedControl.style.display = 'none';
        }
        if (moleculeControls) {
            moleculeControls.style.display = 'flex';

            const prevBtn = document.getElementById('prev-molecule');
            const nextBtn = document.getElementById('next-molecule');
            const speedControl = document.querySelector(".speed-control");

            if (prevBtn) prevBtn.style.display = molecules.length > 1 ? 'inline-block' : 'none';
            if (nextBtn) nextBtn.style.display = molecules.length > 1 ? 'inline-block' : 'none';
            if (speedControl) speedControl.style.display = is3DEnabled ? 'flex' : 'none';
        }
    updateNavigationButtons();
}

function createBubbles(container, count = 15) {
    container.innerHTML = '';  // limpia burbujas anteriores
    for (let i = 0; i < count; i++) {
        const bubble = document.createElement('div');
        bubble.classList.add('bubble');
        bubble.style.left = `${Math.random() * 100}%`;
        const size = `${6 + Math.random() * 6}px`;
        bubble.style.width = size;
        bubble.style.height = size;
        bubble.style.animationDuration = `${3 + Math.random() * 2}s`;
        bubble.style.animationDelay = `${Math.random() * 5}s`;
        container.appendChild(bubble);
    }
}

function init3DMol() {
    // Verificar si 3Dmol está cargado
    if (window.$3Dmol) {
        console.log("3Dmol está cargado correctamente");
    } else {
        console.error("3Dmol no se cargó correctamente");
        // Podrías recargar la librería aquí si es necesario
    }
}

document.querySelectorAll('.mol-entry').forEach((entry, idx) => {
    entry.addEventListener('click', () => {
        selected.clear();
        selected.add(idx);

        // Remover estilo 'selected' de todas las entradas
        document.querySelectorAll('.mol-entry').forEach(e => e.classList.remove('selected'));
        entry.classList.add('selected');

        // Actualizar el índice y cambiar la visualización a 2D
        currentMoleculeIndex = idx;
        is2DEnabled = true;
        is3DEnabled = false;
        btn2D.classList.add("active");
        btn3D.classList.remove("active");

        // Actualizar visualización
        updateViewerContainerLayout();
    });
});

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
    const convertBtn = document.getElementById('convert-molecules');
    const loader = convertBtn.querySelector('.water-loader');
    const label = convertBtn.querySelector('.button-label');

    loader.style.display = 'block';
    label.style.opacity = '0.3';
    convertBtn.disabled = true;
    createBubbles(loader);
    
    try {
        if (!molecules.length) {
            throw new Error('No hay moléculas cargadas para convertir');
        }

        const prevBtn = document.getElementById('prev-molecule');
        const nextBtn = document.getElementById('next-molecule');
        const speedControl = document.querySelector('.speed-control');

        let selectedIndices = null;

        if (molecules.length > 1) {
            const selected = await showMoleculeSelectionModal();
            if (!selected || selected.size === 0) return;

            selectedIndices = [...selected];

            // Convertir y almacenar XYZ de las seleccionadas
            const selectedMolecules = [];
            for (const index of selectedIndices) {
                const mol = molecules[index];
                const method = document.getElementById('conversion-method').value;
                const forceField = document.getElementById('force-field').value;

                const response = await fetch(`${API_BASE_URL}/convert`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        smiles: mol.smiles,
                        identifier: mol.id,
                        method,
                        force_field: forceField
                    })
                });

                if (!response.ok) {
                    const errorData = await response.json();
                    throw new Error(errorData.error || 'Error en la conversión');
                }

                const data = await response.json();
                selectedMolecules.push({ ...mol, xyz: data.xyz });
                showSuccess(`Molécula convertida exitosamente`);
            }

            molecules = selectedMolecules;
            currentMoleculeIndex = 0;
            currentXYZ = molecules[0].xyz;

            render3DMolecule(currentXYZ, molecules[0].smiles);
            updateXYZDisplay();

        } else {
            // Solo una molécula, convertirla
            const mol = molecules[0];
            const method = document.getElementById('conversion-method').value;
            const forceField = document.getElementById('force-field').value;

            const response = await fetch(`${API_BASE_URL}/convert`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    smiles: mol.smiles,
                    identifier: mol.id,
                    method,
                    force_field: forceField
                })
            });

            if (!response.ok) {
                const errorData = await response.json();
                throw new Error(errorData.error || 'Error en la conversión');
            }

            const data = await response.json();
            mol.xyz = typeof data.xyz === 'string' ? data.xyz : data.xyz?.xyz;
            mol.molBlock = data.mol || null;

            showSuccess(`Molécula convertida exitosamente`);

            currentMoleculeIndex = 0;
            currentXYZ = typeof data.xyz === 'string' ? data.xyz : data.xyz?.xyz || '';

            const format = mol.molBlock ? 'sdf' : 'xyz';
            const modelData = mol.molBlock || mol.xyz;

            render3DMolecule(modelData, mol.smiles, format);
            updateXYZDisplay();
        }

        // Mostrar controles si hay más de una molécula
        if (moleculeControls) {
            const multiple = molecules.length > 1;
            if (prevBtn) prevBtn.style.display = multiple ? 'inline-block' : 'none';
            if (nextBtn) nextBtn.style.display = multiple ? 'inline-block' : 'none';
            if (speedControl) speedControl.style.display = 'flex';
            moleculeControls.style.display = 'flex';
        }

    loader.style.display = 'none';
    label.style.visibility = 'visible';
    convertBtn.disabled = false;

    } catch (error) {
        console.error("Error en convertMolecules:", error);
        loader.style.display = 'none';
        label.style.opacity = '1';
        convertBtn.disabled = false;
        loader.innerHTML = ''; // eliminar burbujas
        showError(`Error al convertir molécula: ${error.message}`);
    }

    loader.style.display = 'none';
    label.style.opacity = '1';
    convertBtn.disabled = false;
    loader.innerHTML = '';
    updateNavigationButtons();
}

// Mostrar controles de navegación si hay múltiples moléculas
if (moleculeControls) {
    const prevBtn = document.getElementById('prev-molecule');
    const nextBtn = document.getElementById('next-molecule');
    const speedControl = document.querySelector('.speed-control');

    const moleculesToShow = selected ? [...selected] : molecules.map((_, idx) => idx);
    const showControls = moleculesToShow.length > 1;

    const selectedEntry = document.querySelector(".mol-entry.selected");
    if (selectedEntry) {
        selectedEntry.classList.remove("selected");
    }
    
    if (prevBtn) prevBtn.style.display = showControls ? 'inline-block' : 'none';
    if (nextBtn) nextBtn.style.display = showControls ? 'inline-block' : 'none';
    if (speedControl) speedControl.style.display = 'flex';

    // Si se seleccionaron moléculas, filtrar `molecules` para evitar navegación a vacíos
    if (selected && selected.size > 0) {
        molecules = moleculesToShow.map(i => molecules[i]);
    }
}


function insertTitleLine(xyz, title) {
    const lines = xyz.trim().split('\n');
    if (lines.length < 2) return xyz;
    lines[1] = title;
    return lines.join('\n');
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

    currentXYZ = typeof data.xyz === 'string' ? data.xyz : data.xyz?.xyz || '';
    updateXYZDisplay();
    const titleLine = `0  0   ${molecule.smiles}`;
    const xyzWithHeader = insertTitleLine(data.xyz, titleLine);
    currentXYZ = xyzWithHeader;
    updateXYZDisplay();

    const format = data.mol ? 'mol' : 'xyz';
    const modelData = data.mol || xyzWithHeader;

    render3DMolecule(modelData, molecule.smiles, format);
    showSuccess(`Molécula convertida exitosamente`);
}

function updateXYZDisplay() {
    const xyzDisplay = document.getElementById('xyz-display');
    if (xyzDisplay) {
        const lines = currentXYZ.trim().split('\n');
        if (lines.length >= 3) {
            const header = lines.slice(0, 2);
            const atoms = lines.slice(2).map(line => {
                const parts = line.trim().split(/\s+/);
                return parts.length === 4
                ? `${parts[0].padEnd(2)} ${parts[1].padStart(8)} ${parts[2].padStart(8)} ${parts[3].padStart(8)}`
                : line;
            });
            xyzDisplay.textContent = [...header, ...atoms].join('\n');
        } else {
            xyzDisplay.textContent = currentXYZ; // fallback
        }
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
            <h3 style="color: #98ebc4; margin-bottom: 15px;">Select the molecule to convert</h3>
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
                    Cancel
                </button>
                <button onclick="const selectedIndex = parseInt(document.getElementById('modal-selected-index').value); 
                                 this.parentNode.parentNode.parentNode.removeChild(this.parentNode.parentNode); 
                                 resolve(selectedIndex >= 0 ? selectedIndex : null);"
                        style="background: #98ebc4; padding: 8px 15px; border: none; border-radius: 5px; cursor: pointer;">
                    Convert
                </button>
            </div>
        `;
        
        modal.appendChild(content);
        document.body.appendChild(modal);
    });
}

function render3DMolecule(modelData, smiles = '', format = 'xyz', targetContainer = viewerContainer) {
    // Limpiar visualización anterior
    if (rotationAnimationId) {
        cancelAnimationFrame(rotationAnimationId);
        rotationAnimationId = null;
    }

    targetContainer.innerHTML = '';
    moleculeTitle.textContent = smiles;
    moleculeTitle.style.display = 'block';

    if (window.$3Dmol && window.$3Dmol.createViewer) {
        const viewer = window.$3Dmol.createViewer(viewerContainer, {
            width: '100%',
            height: '100%',
            backgroundColor: 'black'
        });

        viewer.addModel(modelData, format, {
            keepH: true,
            noComputeSecondaryStructure: true
        });

        viewer.setStyle({}, {
            stick: { radius: 0.15 },
            sphere: { scale: 0.25 }
        });

        viewer.setStyle({}, { stick: {}, sphere: { scale: 0.25 } });  // mixto

        // Animación de rotación
        const animate = function () {
            if (!isRotationPaused) {
                viewer.rotate(rotationSpeed, { x: 0, y: 1, z: 0 });
            }
            rotationAnimationId = requestAnimationFrame(animate);
        };
        animate();

        currentVisualizations = [viewer];

        // Ocultar placeholder y mostrar el contenedor si todo va bien
        const placeholder = document.querySelector('.placeholder');
        if (placeholder) placeholder.style.display = 'none';

        if (targetContainer) targetContainer.style.display = 'block';

    } else {
        viewerContainer.innerHTML = `
            <div class="placeholder error">
                Error: La biblioteca 3Dmol no está disponible
            </div>
        `;
    }
    if (!targetContainer) return;

    if (rotationAnimationId) {
        cancelAnimationFrame(rotationAnimationId);
        rotationAnimationId = null;
    }

    targetContainer.innerHTML = '';
    if (window.$3Dmol && window.$3Dmol.createViewer) {
        const viewer = window.$3Dmol.createViewer(targetContainer, {
            width: targetContainer.offsetWidth || targetContainer.clientWidth || 500,
            height: targetContainer.offsetHeight || targetContainer.clientHeight || 400,
            backgroundColor: 'black'
        });

        viewer.addModel(modelData, format, {
            keepH: true,
            noComputeSecondaryStructure: true
        });

        viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.25 } });
        viewer.zoomTo();
        viewer.zoom(1.5);
        viewer.render();
        viewer.zoomTo();

        const animate = function () {
            if (!isRotationPaused) {
                viewer.rotate(rotationSpeed, { x: 0, y: 1, z: 0 });
                viewer.render();
            }
            rotationAnimationId = requestAnimationFrame(animate);
        };
        animate();

        currentVisualizations = [viewer];
    } else {
        targetContainer.innerHTML = `<div class="placeholder error">Error: La biblioteca 3Dmol no está disponible</div>`;
    }

    setTimeout(() => {
        viewer.resize();
        viewer.zoomTo();
        viewer.render();
    }, 100);
}

function navigateMolecules(direction) {
    const newIndex = currentMoleculeIndex + direction;
    if (newIndex >= 0 && newIndex < molecules.length) {
        currentMoleculeIndex = newIndex;
        currentXYZ = molecules[newIndex].xyz;
        const mol = molecules[newIndex];
        const format = mol.molBlock ? 'sdf' : 'xyz';
        const modelData = mol.molBlock || mol.xyz;
        render3DMolecule(modelData, mol.smiles, format);
        updateXYZDisplay();
        updateNavigationButtons();
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
        const formattedXYZ = formatXYZWithAlignment(currentXYZ);
        const blob = new Blob([formattedXYZ], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        const smilesName = molecules[currentMoleculeIndex]?.smiles || 'unknown';
        const safeName = smilesName.replace(/[\\/:*?"<>|]/g, '_'); 
        a.download = `${safeName}.xyz`;                 
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

        const singleInput = document.getElementById('single-input');
        const fileInput = document.getElementById('file-input');
        const loadButton = document.getElementById('load-smiles');
        const fileFormatExample = document.getElementById('file-format-example');

        singleInput.style.display = method === 'single' ? 'block' : 'none';
        fileInput.style.display = method === 'file' ? 'block' : 'none';
        fileFormatExample.style.display = method === 'file' ? 'block' : 'none';

        // Asegurar la visibilidad del botón
        if (method === 'single') {
            loadButton.style.display = 'inline-block';
            loadButton.disabled = false;
        } else {
            loadButton.style.display = 'none';
        }

    } catch (error) {
        console.error("Error en toggleInputMethod:", error);
    }
}


async function loadSmiles(showConfirmation = true) {
    try {
        console.log("Iniciando carga de SMILES...");
        if (moleculeControls) {
            const prevBtn = document.getElementById('prev-molecule');
            const nextBtn = document.getElementById('next-molecule');
            const speedControl = document.querySelector('.speed-control');
        
            if (molecules.length > 1) {
                moleculeControls.style.display = 'flex';
                if (prevBtn) prevBtn.style.display = 'inline-block';
                if (nextBtn) nextBtn.style.display = 'inline-block';
            } else {
                moleculeControls.style.display = 'flex'; // Mantener visible la barra
                if (prevBtn) prevBtn.style.display = 'none';
                if (nextBtn) nextBtn.style.display = 'none';
            }
        
            if (speedControl) speedControl.style.display = 'flex'; // Asegurar visibilidad
        }       
        const method = document.querySelector('input[name="input-method"]:checked').value;
        console.log("Método seleccionado:", method);

        if (method === 'single') {
            const smiles = document.getElementById('smiles').value.trim();
            console.log("SMILES ingresado:", smiles);
            updateSmilesHistory(smiles);

            if (!smiles) {
                throw new Error('Por favor ingrese un string SMILES');
            }
            molecules = [{ smiles, id: smiles }];

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

        if (loadConfirmation && showConfirmation) {
            const moleculeInfo = molecules.map(mol => 
                `<div><strong>SMILE:</strong> ${mol.smiles}</div>`
            ).join('<hr>');
            loadConfirmation.innerHTML = `
                <div class="alert success" style="text-align: left; padding: 10px 15px;">
                    <strong>${molecules.length === 1 ? 'Molécula cargada exitosamente:' : `${molecules.length} moléculas cargadas exitosamente:`}</strong><br>
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
    const method = document.getElementById('conversion-method')?.value;
    const forceFieldDiv = document.getElementById('force-field-options');

    if (!method || !forceFieldDiv) {
        console.warn("Elementos necesarios para toggleForceField no están disponibles aún.");
        return;
    }

    if (method === 'networkx' || method === 'auto3d') {
        forceFieldDiv.style.display = 'none';
    } else {
        forceFieldDiv.style.display = 'block';
    }
}

function retryVisualization() {
    if (currentXYZ && molecules[currentMoleculeIndex]) {
        render3DMolecule(currentXYZ, molecules[currentMoleculeIndex].smiles);
    }
}

function formatXYZWithAlignment(xyzText) {
    const lines = xyzText.trim().split('\n');
    if (lines.length < 3) return xyzText;

    const header = lines.slice(0, 2);
    const atomLines = lines.slice(2).map(line => {
        const parts = line.trim().split(/\s+/);
        return parts.length === 4
            ? `${parts[0].padEnd(2)} ${parts[1].padStart(8)} ${parts[2].padStart(8)} ${parts[3].padStart(8)}`
            : line;
    });

    return [...header, ...atomLines].join('\n');
}


// Funciones de ayuda
function showSuccess(message) {
    console.log("Mostrando mensaje de éxito:", message);
    const alertDiv = document.createElement('div');
    alertDiv.className = 'alert success';
    const [first, ...rest] = message.split(':');
    const formattedMessage = rest.length
        ? `<strong>${first}:</strong><br>${rest.join(':').trim()}`  //✓
        : `<strong>${message}</strong>`;
    alertDiv.innerHTML = formattedMessage;
    document.body.appendChild(alertDiv);
    setTimeout(() => {
        alertDiv.remove();
    }, 3000);
}


function showError(message) {
    console.error("Mostrando mensaje de error:", message);
    const alertDiv = document.createElement('div');
    alertDiv.className = 'alert error';
    const [first, ...rest] = message.split(':');
    const formattedMessage = rest.length
        ? `⚠ <strong>${first}:</strong><br>${rest.join(':').trim()}`
        : `⚠ ${message}`;

    alertDiv.innerHTML = formattedMessage;
    document.body.appendChild(alertDiv);
    setTimeout(() => {
        alertDiv.remove();
    }, 5000);
}

async function showMoleculeSelectionModal() {
    return new Promise((resolve) => {
        selectedMolecules = new Set(); // ← resetear selección

        const modal = document.createElement('div');
        modal.className = 'modal-overlay';

        // Variables de navegación por teclado
        let focusedIndex = 0;
        let shiftAnchor = null;

        modal.innerHTML = `
            <div class="modal-content" tabindex="0">
                <h3>Seleccione los SMILES a convertir</h3>
                <div id="molecule-list" class="molecule-list" style="max-height: 300px; overflow-y: auto;"></div>
                <div class="modal-buttons">
                    <button id="modal-cancel">Cancelar</button>
                    <button id="modal-confirm">Convertir seleccionadas</button>
                </div>
            </div>
        `;

        document.body.appendChild(modal);
        const moleculeList = modal.querySelector('#molecule-list');

        // Renderizar lista de moléculas
        molecules.forEach((mol, index) => {
            const item = document.createElement('div');
            item.className = 'molecule-item';
            item.textContent = `${mol.id}: ${mol.smiles}`;
            item.dataset.index = index;
            moleculeList.appendChild(item);
        });

        const items = moleculeList.querySelectorAll('.molecule-item');

        // Aplicar foco visual
        const updateFocus = () => {
            items.forEach((item, i) => {
                item.classList.toggle('focused', i === focusedIndex);
            });
            if (focusedIndex >= 0) {
                items[focusedIndex].scrollIntoView({ block: "nearest", behavior: "smooth" });
            }
        };

        // Manejar selección con teclado
        const handleKeyDown = (e) => {
            const maxIndex = items.length - 1;

            if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
                e.preventDefault();
                const oldIndex = focusedIndex;

                focusedIndex += (e.key === 'ArrowDown') ? 1 : -1;
                focusedIndex = Math.max(0, Math.min(focusedIndex, maxIndex));

                if (e.shiftKey) {
                    if (shiftAnchor === null) shiftAnchor = oldIndex;
                    const [start, end] = [shiftAnchor, focusedIndex].sort((a, b) => a - b);
                    selectedMolecules.clear();
                    items.forEach((item, i) => {
                        if (i >= start && i <= end) {
                            item.classList.add('selected');
                            selectedMolecules.add(i);
                        } else {
                            item.classList.remove('selected');
                        }
                    });
                } else {
                    shiftAnchor = null;
                }

                updateFocus();
            }

            if (e.key === 'Enter') {
                cleanup(true);
            }

            if (e.key === 'Escape') {
                cleanup(false);
            }
        };

        // Manejo con mouse
        items.forEach((item, index) => {
            item.addEventListener('click', () => {
                item.classList.toggle('selected');
                focusedIndex = index;
                shiftAnchor = index;

                if (item.classList.contains('selected')) {
                    selectedMolecules.add(index);
                } else {
                    selectedMolecules.delete(index);
                }

                updateFocus();
            });
        });

        // Confirmar o cancelar con botones
        modal.querySelector('#modal-cancel').addEventListener('click', () => cleanup(false));
        modal.querySelector('#modal-confirm').addEventListener('click', () => cleanup(true));

        // Cerrar haciendo clic fuera del modal
        modal.addEventListener('click', (event) => {
            if (event.target === modal) cleanup(false);
        });

        // Función para limpiar y cerrar
        function cleanup(confirm) {
            document.removeEventListener('keydown', handleKeyDown);
            document.body.removeChild(modal);
            resolve(confirm ? selectedMolecules : null);
        }

        // Escuchar teclado global
        document.addEventListener('keydown', handleKeyDown);
        updateFocus();
        modal.querySelector('.modal-content').focus();
    });
}

// Cerrar el modal al hacer clic fuera del contenido
modal.addEventListener('click', (e) => {
    const modalContent = modal.querySelector('.modal-content');
    if (!modalContent.contains(e.target)) {
        document.removeEventListener('keydown', handleKeyDown);
        modal.remove();
        resolve(null);
    }
});


function updateSmilesHistory(newSmiles) {
    const key = 'smilesHistory';
    let history = JSON.parse(localStorage.getItem(key)) || [];

    // Elimina duplicados y agrega nuevo al principio
    history = [newSmiles, ...history.filter(s => s !== newSmiles)];

    // Limita a 5 entradas
    history = history.slice(0, 5);

    // Guarda en localStorage
    localStorage.setItem(key, JSON.stringify(history));

    // Actualiza el datalist
    const datalist = document.getElementById('smiles-history');
    datalist.innerHTML = '';
    history.forEach(smiles => {
        const option = document.createElement('option');
        option.value = smiles;
        datalist.appendChild(option);
    });
}

function loadSmilesHistory() {
    const key = 'smilesHistory';
    const history = JSON.parse(localStorage.getItem(key)) || [];
    const datalist = document.getElementById('smiles-history');
    if (datalist) {
        datalist.innerHTML = '';
        history.forEach(smiles => {
            const option = document.createElement('option');
            option.value = smiles;
            datalist.appendChild(option);
        });
    }
}

function setupConversionHandler() {
    const convertButton = document.getElementById("convert-molecules");
    if (convertButton) {
        convertButton.addEventListener("click", async () => {
            try {
                await convertMolecules();
            } catch (error) {
                showError(error.message || "Conversion failed");
            } finally {
                showLoader(false);
            }
        });
    }
}

function showLoader(show) {
    const loader = document.querySelector('.water-loader');
    if (loader) {
        loader.style.display = show ? 'block' : 'none';
    }
}

function updateNavigationButtons() {
    if (prevMoleculeBtn) {
        prevMoleculeBtn.style.display = currentMoleculeIndex > 0 ? 'inline-block' : 'none';
    }
    if (nextMoleculeBtn) {
        nextMoleculeBtn.style.display = currentMoleculeIndex < molecules.length - 1 ? 'inline-block' : 'none';
    }
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
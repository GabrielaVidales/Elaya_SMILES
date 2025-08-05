# -*- coding: utf-8 -*-
"""                                ELAYA SMILES
                 Molecular Conversion and Analysis Tool
        "Energy-minimized Linear-to-structure Atom Yielding Algorithm"
Elaborado por Gabriela Vidales, Luis Gonzalez, Filiberto Ortiz, and Gabriel Merino.

Funcionalidades principales
- Conversión SMILES a 3D (RDKit, OpenBabel, NetworkX, Auto3D)
- Visualización 3D (py3Dmol)
- Carga individual o por archivo.
- Exportación de estructuras en formato .xyz.
- Interfaz moderna y amigable.
"""

import os
import shutil
import contextlib
import numpy as np
from scipy.spatial.distance import pdist, squareform
import networkx as nx

# RDKit and cheminformatics
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops, Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import MolToMolBlock
from rdkit.Geometry import Point3D

# 3D visualization
import py3Dmol

# Molecular dynamics and descriptors
from ase.io import read, write
from dscribe.descriptors import SOAP, ValleOganov
from dscribe.kernels import AverageKernel

# OpenBabel
from openbabel import openbabel, pybel
from openbabel.pybel import readstring
from openbabel import openbabel as ob

# Auto3D for AI-based optimization
import Auto3D
from Auto3D.auto3D import options, main

class MolecularTools:
    def __init__(self):
        """Initialize with ASCII art and setup"""
        os.environ['BABEL_DATADIR'] = os.path.abspath('openbabel_data')
        self.print_banner()
        self.setup_directories()

    def print_banner(self):
        """Display program banner"""
        print("""
                              ELAYA SMILES
                 Molecular Conversion and Analysis Tool
        "Energy-minimized Linear-to-structure Atom Yielding Algorithm"
        """)

    def setup_directories(self):
        """Create necessary directories"""
        self.dirs = {
            'smi': 'output_smi_individuales',
            'xyz_rdkit': 'output_xyz_rdkit',
            'xyz_openbabel': 'output_xyz_openbabel',
            'xyz_networkx': 'output_xyz_networkx',
            'xyz_auto3d': 'output_xyz_auto3d',
            'similarity': 'similarity_analysis'
        }

        for dir_name, dir_path in self.dirs.items():
            os.makedirs(dir_path, exist_ok=True)
            self.dirs[dir_name] = os.path.abspath(dir_path)
            print(f"Directorio {dir_name}: {self.dirs[dir_name]}")

    # SMILES Processing Methods
    def load_smiles(self, input_type=1, file_path=None):
        """Load SMILES either from single input or file"""
        if input_type == 1:
            self.smiles_list = [input("Enter SMILES: ")]
            self.identifiers = ["mol1"]
        else:
            if not file_path:
                raise ValueError("No file path provided")
                
            self.smiles_list = []
            self.identifiers = []
            with open(file_path, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        self.smiles_list.append(parts[0])
                        self.identifiers.append(parts[1])

    # 3D Conversion Methods
    def rdkit_conversion(self, smiles, identifier, force_field='uff', optimize=True):
        """Convert SMILES to 3D using RDKit"""
        try:
            print(f"Intentando convertir SMILES: {smiles}")
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"RDKit no pudo parsear el SMILES: {smiles}")
            
            mol = Chem.AddHs(mol)

            # Generate 3D coordinates
            params = AllChem.ETKDGv3()
            params.randomSeed = 665
            status = AllChem.EmbedMolecule(mol, params=params)
            if status == -1:
                raise ValueError("Embedding failed - could not generate 3D coordinates")
            print("Coordenadas 3D generadas")

            # Optimization
            if optimize:
                if force_field.lower() == 'uff':
                    print("Optimizando con UFF...")
                    status = AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
                    if status == -1:
                        print(f"Warning: UFF optimization failed for {identifier}")
                else:  # MMFF94
                    print("Optimizando con MMFF94...")
                    mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
                    if mmff_props is None:
                        print(f"Warning: MMFF94 not available for {identifier}, using UFF instead")
                        AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
                    else:
                        ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props)
                        ff.Minimize(maxIts=2000)

            # Save XYZ file
            xyz_content = f"{mol.GetNumAtoms()}\n0  0   {smiles}\n"
            conf = mol.GetConformer()
            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                xyz_content += f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"

            safe_name = smiles.translate(str.maketrans({
                '\\': '_', '/': '_', ':': '_', '*': '_', '?': '_', '"': '_', '<': '_', '>': '_', '|': '_'
            }))
            output_path = os.path.join(self.dirs['xyz_rdkit'], f"{safe_name}.xyz")
            with open(output_path, 'w') as f:
                f.write(xyz_content)
            
            mol_block = MolToMolBlock(mol)

            print(f"Conversión exitosa para {identifier}")
            return {"xyz": xyz_content, "mol": mol_block}

        except Exception as e:
            print(f"Error converting {smiles}: {str(e)}")
            raise


    def openbabel_conversion(self, smiles, identifier, force_field='uff'):
        """Convert SMILES to 3D using OpenBabel with enhanced bond detection and warm-up"""
        try:
            print(f"Convirtiendo {smiles} con OpenBabel...")

            # --- Warm-up opcional: evita fallos internos en primeras llamadas a make3D ---
            try:
                warm_conv = openbabel.OBConversion()
                warm_conv.SetInAndOutFormats("smi", "mol")
                warm_mol = openbabel.OBMol()
                warm_conv.ReadString(warm_mol, '[O]')
                warm_pybel = pybel.readstring("mol", warm_conv.WriteString(warm_mol))
                warm_pybel.make3D(forcefield=force_field, steps=1)
                print("Warm-up de OpenBabel completado")
            except Exception as warm_error:
                print("Warm-up fallido o innecesario:", warm_error)

            # --- Conversión principal ---
            conv = openbabel.OBConversion()
            conv.SetInAndOutFormats("smi", "mol")

            mol = openbabel.OBMol()
            if not conv.ReadString(mol, smiles):
                raise ValueError(f"OpenBabel no pudo interpretar el SMILES: {smiles}")

            mol.AddHydrogens()
            mol.PerceiveBondOrders()
            self.add_lone_pairs_openbabel(mol)

            # Generar coordenadas y optimizar
            builder = openbabel.OBBuilder()
            builder.Build(mol)

            ff = openbabel.OBForceField.FindForceField(force_field.upper())
            if ff is None:
                raise ValueError(f"Force field '{force_field}' no encontrado en OpenBabel")

            ff.Setup(mol)
            ff.ConjugateGradients(2000)
            ff.GetCoordinates(mol)

            # Generar archivo MOL
            conv_mol = openbabel.OBConversion()
            conv_mol.SetOutFormat("mol")
            mol_block = conv_mol.WriteString(mol)

            # Generar archivo XYZ
            xyz_content = f"{mol.NumAtoms()}\n0  0   {smiles}\n"
            for atom in ob.OBMolAtomIter(mol):
                x, y, z = atom.GetX(), atom.GetY(), atom.GetZ()
                xyz_content += f"{atom.GetType()[0]} {x:.4f} {y:.4f} {z:.4f}\n"

            # Guardar XYZ
            safe_name = smiles.translate(str.maketrans('\\/:*?"<>|', '_________'))
            output_path = os.path.join(self.dirs['xyz_openbabel'], f"{safe_name}.xyz")
            with open(output_path, 'w') as f:
                f.write(xyz_content)

            print(f"Conversión exitosa para {identifier}")
            return {"xyz": xyz_content, "mol": mol_block}

        except Exception as e:
            print(f"Error en OpenBabel: {str(e)}")
            raise

    def networkx_conversion(self, smiles, identifier):
        """Generar coordenadas aproximadas en 3D usando NetworkX a partir del SMILES"""
        try:
            print(f"Generando coordenadas aproximadas para {smiles}")

            # Convertir SMILES a Mol y añadir hidrógenos
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"SMILES inválido: {smiles}")
            mol = Chem.AddHs(mol)

            # Construcción de nodos y enlaces para el grafo molecular
            self.atoms = [(atom.GetIdx(), atom.GetSymbol()) for atom in mol.GetAtoms()]
            self.bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]

            # Crear grafo y calcular layout 3D aproximado
            G = nx.Graph()
            G.add_nodes_from(self.atoms)
            G.add_edges_from(self.bonds)
            coord_xyz = nx.spring_layout(G, dim=3, k=5, iterations=800, scale=3.4, seed=586)
            
            # Generar contenido XYZ
            self.xyz_networkx = []
            xyz_content = f"{len(self.atoms)}\n0  0   {smiles}\n"
            for atom in self.atoms:
                idx = atom[0]
                symbol = atom[1]
                x, y, z = coord_xyz[idx]
                self.xyz_networkx.append([round(x, 4), round(y, 4), round(z, 4)])
                xyz_content += f"{symbol} {x:.4f} {y:.4f} {z:.4f}\n"

            # Guardar archivo XYZ
            safe_name = smiles.translate(str.maketrans('\\/:*?"<>|', '_________'))
            output_path = os.path.join(self.dirs['xyz_openbabel'], f"{safe_name}.xyz")
            with open(output_path, 'w') as f:
                f.write(xyz_content)

            print(f"Conversión con NetworkX completada para {identifier}")
            return xyz_content

        except Exception as e:
            print(f"Error en NetworkX para {identifier}: {str(e)}")
            raise

    def auto3d_conversion(self, smiles, identifier):
        """Use Auto3D for AI-based 3D structure generation"""
        try:
            print(f"Usando Auto3D para {smiles}")
            # Create temporary SMILES file
            temp_file = os.path.join(self.dirs['smi'], f"temp_{identifier}.smi")
            with open(temp_file, 'w') as f:
                f.write(f"{smiles} {identifier}\n")

            # Run Auto3D
            args = options(temp_file, k=1, use_gpu=False)
            out_path = main(args)

            # Move output to our directory
            final_path = os.path.join(self.dirs['xyz_auto3d'], f"{identifier}.xyz")
            if os.path.isfile(out_path):
                shutil.move(out_path, final_path)
            else:
                raise FileNotFoundError(f"No se encontró el archivo de salida de Auto3D: {out_path}")

            # Read and return the content
            with open(final_path, 'r') as f:
                return f.read()
        except Exception as e:
            print(f"Error en Auto3D: {str(e)}")
            raise

    """
    def add_lone_pairs(self, mol):
        rw_mol = Chem.RWMol(mol)
        original_num_atoms = rw_mol.GetNumAtoms()

        # Creamos lista de posiciones originales
        conf_old = mol.GetConformer()
        positions = [conf_old.GetAtomPosition(i) for i in range(original_num_atoms)]

        # Añadir átomos virtuales y calcular sus posiciones
        new_positions = []
        for atom_idx in range(original_num_atoms):
            atom = mol.GetAtomWithIdx(atom_idx)
            symbol = atom.GetSymbol()
            if symbol in ['O', 'N', 'F']:
                pos = positions[atom_idx]
                for dx, dy, dz in [(0.3, 0.3, 0.0), (-0.3, -0.3, 0.0)]:
                    dummy = Chem.Atom(0)  # átomo dummy
                    idx = rw_mol.AddAtom(dummy)
                    new_positions.append((idx, Point3D(pos.x + dx, pos.y + dy, pos.z + dz)))

        # Crear un nuevo conformador del tamaño correcto
        total_atoms = rw_mol.GetNumAtoms()
        conf_new = Chem.Conformer(total_atoms)

        # Asignar posiciones originales
        for i in range(original_num_atoms):
            conf_new.SetAtomPosition(i, positions[i])

        # Asignar posiciones de los nuevos átomos
        for idx, pos in new_positions:
            conf_new.SetAtomPosition(idx, pos)

        # Asignar el conformador al mol
        rw_mol.RemoveAllConformers()
        rw_mol.AddConformer(conf_new)

        return rw_mol.GetMol()
    """
    
        # Visualization Methods
    def visualize_3d(self, xyz_content, width=400, height=400):
        """Visualize molecule from XYZ content"""
        try:
            view = py3Dmol.view(width=width, height=height)
            view.addModel(xyz_content, "xyz")
            view.addBonds()  # <<<<<< AÑADIDO: intenta inferir enlaces desde posiciones
            view.setStyle({'sphere': {'scale': 0.3}, 'stick': {'radius': 0.2}})
            view.zoomTo()
            return view
        except Exception as e:
            print(f"Error en visualización: {str(e)}")
            raise
    
    def add_lone_pairs_openbabel(self, mol):
        """Agrega átomos virtuales tipo He para simular pares libres (visualización)"""
        from openbabel import OBAtom

        atom_count = mol.NumAtoms()
        for i in range(1, atom_count + 1):
            atom = mol.GetAtom(i)
            symbol = atom.GetType()
            if symbol in ['N', 'O', 'F']:
                x, y, z = atom.GetX(), atom.GetY(), atom.GetZ()
                
                # Crear átomo dummy tipo He como par libre (mejor si no interfiere con enlaces)
                for dx, dy, dz in [(0.3, 0.3, 0.0), (-0.3, -0.3, 0.0)]:
                    dummy = OBAtom()
                    dummy.SetAtomicNum(2)  # Helio (Z=2)
                    dummy.SetVector(x + dx, y + dy, z + dz)
                    mol.AddAtom(dummy)

    # Similarity Analysis Methods
    def tanimoto_similarity(self, xyz_file1, xyz_file2):
        """Calculate Tanimoto similarity between two XYZ structures"""
        try:
            def parse_xyz(xyz_file):
                coords = set()
                with open(xyz_file, 'r') as f:
                    lines = f.readlines()[2:]  # Skip first two lines
                    for line in lines:
                        parts = line.split()
                        if len(parts) >= 4:
                            element = parts[0]
                            x, y, z = map(lambda v: round(float(v), 2), parts[1:4])
                            coords.add((element, x, y, z))
                return coords

            set1 = parse_xyz(xyz_file1)
            set2 = parse_xyz(xyz_file2)

            intersection = len(set1.intersection(set2))
            union = len(set1.union(set2))
            return intersection / union if union != 0 else 0
        except Exception as e:
            print(f"Error en similitud Tanimoto: {str(e)}")
            raise

    def soap_similarity(self, xyz_files, species=None, r_cut=5.0):
        """Calculate SOAP similarity matrix for multiple XYZ files"""
        try:
            if not species:
                # Auto-detect species from first file
                with open(xyz_files[0], 'r') as f:
                    elements = set()
                    for line in f.readlines()[2:]:
                        parts = line.split()
                        if parts:
                            elements.add(parts[0])
                    species = list(elements)

            # Create SOAP descriptor
            soap = SOAP(
                species=species,
                r_cut=r_cut,
                n_max=4,
                l_max=4,
                periodic=False
            )

            # Process all molecules
            features = []
            for xyz_file in xyz_files:
                mol = read(xyz_file, format='xyz')
                features.append(soap.create(mol))

            # Calculate similarity matrix
            kernel = AverageKernel(metric="linear")
            return kernel.create(features)
        except Exception as e:
            print(f"Error en SOAP: {str(e)}")
            raise

    def valle_oganov_similarity(self, xyz_files, species=None):
        """Calculate Valle-Oganov similarity matrix"""
        try:
            if not species:
                # Auto-detect species from first file
                with open(xyz_files[0], 'r') as f:
                    elements = set()
                    for line in f.readlines()[2:]:
                        parts = line.split()
                        if parts:
                            elements.add(parts[0])
                    species = list(elements)

            # Create Valle-Oganov descriptor
            vo = ValleOganov(
                species=species,
                function='distance',
                n=100,
                sigma=1E-5,
                r_cut=10
            )

            # Process all molecules
            features = []
            for xyz_file in xyz_files:
                mol = read(xyz_file, format='xyz')
                features.append(vo.create(mol))

            # Calculate similarity matrix
            sim_matrix = np.zeros((len(features), len(features)))
            for i in range(len(features)):
                for j in range(i, len(features)):
                    norm_i = np.linalg.norm(features[i])
                    norm_j = np.linalg.norm(features[j])
                    dot_product = np.dot(features[i], features[j])
                    similarity = dot_product / (norm_i * norm_j)
                    sim_matrix[i, j] = similarity
                    sim_matrix[j, i] = similarity

            return sim_matrix
        except Exception as e:
            print(f"Error en Valle-Oganov: {str(e)}")
            raise

    # Connectivity Matrix Methods
    def get_connectivity_matrix(self, xyz_file, threshold=2.0):
        """Generate connectivity matrix from XYZ file"""
        try:
            with open(xyz_file, 'r') as f:
                num_atoms = int(f.readline().strip())
                _ = f.readline()  # Skip comment line
                coords = []
                elements = []
                for _ in range(num_atoms):
                    parts = f.readline().split()
                    elements.append(parts[0])
                    coords.append(list(map(float, parts[1:4])))

            coords = np.array(coords)
            distances = squareform(pdist(coords))
            connectivity = (distances < threshold).astype(int)
            np.fill_diagonal(connectivity, 0)

            return elements, connectivity
        except Exception as e:
            print(f"Error en matriz de conectividad: {str(e)}")
            raise

    # Batch Processing Methods
    def process_all_smiles(self, method='rdkit', force_field='uff'):
        """Process all loaded SMILES with specified method"""
        results = {}
        for smi, ident in zip(self.smiles_list, self.identifiers):
            try:
                if method == 'rdkit':
                    result = tool.rdkit_conversion(smiles, identifier, force_field)
                elif method == 'openbabel':
                    result = tool.openbabel_conversion(smiles, identifier)
                elif method == 'networkx':
                    result = {"xyz": tool.networkx_conversion(smiles, identifier), "mol": None}
                else:
                    return jsonify({"error": "Método no soportado"}), 400

                results[ident] = {
                    'xyz': xyz,
                    'visualization': self.visualize_3d(xyz)
                }
            except Exception as e:
                print(f"Error procesando {ident}: {str(e)}")
                continue
                
        return results

    def compare_all_methods(self, identifier):
        """Compare results from all conversion methods for one molecule"""
        try:
            # Find the SMILES for this identifier
            idx = self.identifiers.index(identifier)
            smi = self.smiles_list[idx]

            # Generate with all methods
            rdkit_xyz = self.rdkit_conversion(smi, f"{identifier}_rdkit")
            obabel_xyz = self.openbabel_conversion(smi, f"{identifier}_obabel")
            networkx_xyz = self.networkx_conversion(smi, f"{identifier}_networkx")
            auto3d_xyz = self.auto3d_conversion(smi, f"{identifier}_auto3d")

            return {
                'rdkit': rdkit_xyz,
                'openbabel': obabel_xyz,
                'networkx': networkx_xyz,
                'auto3d': auto3d_xyz
            }
        except Exception as e:
            print(f"Error comparando métodos: {str(e)}")
            raise

if __name__ == "__main__":
    tool = MolecularTools()
    print("MolecularTools inicializado correctamente")

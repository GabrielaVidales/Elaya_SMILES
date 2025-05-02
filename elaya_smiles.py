# -*- coding: utf-8 -*-
"""Elaya_smiles

Tool completo para conversión y análisis molecular
"""

import os
import shutil
import numpy as np
from scipy.spatial.distance import pdist, squareform
import networkx as nx

# RDKit and cheminformatics
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops, Draw
from rdkit.Chem.Draw import IPythonConsole

# 3D visualization
import py3Dmol

# Molecular dynamics and descriptors
from ase.io import read, write
from dscribe.descriptors import SOAP, ValleOganov
from dscribe.kernels import AverageKernel

# OpenBabel
from openbabel import openbabel, pybel
from openbabel.pybel import readstring

# Auto3D for AI-based optimization
import Auto3D
from Auto3D.auto3D import options, main

class MolecularTools:
    def __init__(self):
        """Initialize with ASCII art and setup"""
        self.print_banner()
        self.setup_directories()

    def print_banner(self):
        """Display program banner"""
        print("""
        CENTRO DE INVESTIGACIÓN Y ESTUDIOS AVANZADOS DEL IPN (CINVESTAV)
                      THEOCHEM MÉRIDA YUCATÁN
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
            print(f"Molécula con hidrógenos añadidos: {mol.GetNumAtoms()} átomos")

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
            xyz_content = f"{mol.GetNumAtoms()}\n{identifier}\n"
            conf = mol.GetConformer()
            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                xyz_content += f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"

            output_path = os.path.join(self.dirs['xyz_rdkit'], f"{identifier}.xyz")
            with open(output_path, 'w') as f:
                f.write(xyz_content)

            print(f"Conversión exitosa para {identifier}")
            return xyz_content

        except Exception as e:
            print(f"Error converting {smiles}: {str(e)}")
            raise

    def openbabel_conversion(self, smiles, identifier, force_field='uff'):
        """Convert SMILES to 3D using OpenBabel"""
        try:
            print(f"Convirtiendo {smiles} con OpenBabel")
            mol = pybel.readstring("smi", smiles)
            mol.addh()
            mol.make3D(forcefield=force_field, steps=2000)

            # Save XYZ file
            xyz_content = f"{len(mol.atoms)}\n{identifier}\n"
            for atom in mol.atoms:
                x, y, z = atom.coords
                xyz_content += f"{atom.atomicnum} {x:.4f} {y:.4f} {z:.4f}\n"

            output_path = os.path.join(self.dirs['xyz_openbabel'], f"{identifier}.xyz")
            with open(output_path, 'w') as f:
                f.write(xyz_content)

            return xyz_content
        except Exception as e:
            print(f"Error en OpenBabel: {str(e)}")
            raise

    def networkx_conversion(self, smiles, identifier):
        """Generate approximate 3D coordinates using NetworkX"""
        try:
            print(f"Generando coordenadas aproximadas para {smiles}")
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Invalid SMILES: {smiles}")
                
            mol = Chem.AddHs(mol)

            # Create molecular graph
            G = nx.Graph()
            atoms = [(atom.GetIdx(), {"symbol": atom.GetSymbol()}) for atom in mol.GetAtoms()]
            bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]
            G.add_nodes_from(atoms)
            G.add_edges_from(bonds)

            # Generate 3D layout
            pos = nx.spring_layout(G, dim=3, seed=42)

            # Save XYZ file
            xyz_content = f"{len(atoms)}\n{identifier}\n"
            for atom in atoms:
                idx = atom[0]
                x, y, z = pos[idx]
                xyz_content += f"{atom[1]['symbol']} {x:.4f} {y:.4f} {z:.4f}\n"

            output_path = os.path.join(self.dirs['xyz_networkx'], f"{identifier}.xyz")
            with open(output_path, 'w') as f:
                f.write(xyz_content)

            return xyz_content
        except Exception as e:
            print(f"Error en NetworkX: {str(e)}")
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
            shutil.move(out_path, final_path)

            # Read and return the content
            with open(final_path, 'r') as f:
                return f.read()
        except Exception as e:
            print(f"Error en Auto3D: {str(e)}")
            raise

    # Visualization Methods
    def visualize_3d(self, xyz_content, width=400, height=400):
        """Visualize molecule from XYZ content"""
        try:
            view = py3Dmol.view(width=width, height=height)
            view.addModel(xyz_content, "xyz")
            view.setStyle({'sphere': {'scale': 0.3}, 'stick': {'radius': 0.2}})
            view.zoomTo()
            return view
        except Exception as e:
            print(f"Error en visualización: {str(e)}")
            raise

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
                if method.lower() == 'rdkit':
                    xyz = self.rdkit_conversion(smi, ident, force_field)
                elif method.lower() == 'openbabel':
                    xyz = self.openbabel_conversion(smi, ident, force_field)
                elif method.lower() == 'networkx':
                    xyz = self.networkx_conversion(smi, ident)
                elif method.lower() == 'auto3d':
                    xyz = self.auto3d_conversion(smi, ident)
                else:
                    raise ValueError("Invalid method specified")

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
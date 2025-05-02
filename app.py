from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
from elaya_smiles import MolecularTools
import os

app = Flask(__name__, static_folder='.', static_url_path='')
CORS(app)

tool = MolecularTools()

@app.route('/')
def index():
    return send_from_directory('.', 'index.html')

@app.route('/<path:path>')
def static_file(path):
    return send_from_directory('.', path)

@app.route('/api')
def api_base():
    return jsonify({
        "status": "API ready",
        "endpoints": {
            "convert": "/api/convert [POST]",
            "compare": "/api/compare [POST]",
            "similarity": "/api/similarity [POST]",
            "connectivity": "/api/connectivity [POST]"
        }
    })

@app.route('/api/convert', methods=['POST'])
def convert():
    try:
        data = request.json
        print("Datos recibidos:", data)
        
        smiles = data['smiles']
        identifier = data['identifier']
        method = data.get('method', 'rdkit')
        force_field = data.get('force_field', 'uff')
        
        print(f"Convirtiendo {smiles} con {method}")
        
        if method == 'rdkit':
            xyz = tool.rdkit_conversion(smiles, identifier, force_field)
        elif method == 'openbabel':
            xyz = tool.openbabel_conversion(smiles, identifier, force_field)
        elif method == 'networkx':
            xyz = tool.networkx_conversion(smiles, identifier)
        else:
            return jsonify({"error": "Método no soportado"}), 400
            
        return jsonify({"xyz": xyz})
    except Exception as e:
        print("Error en la conversión:", str(e))
        return jsonify({"error": str(e)}), 500

@app.route('/api/compare', methods=['POST'])
def compare():
    try:
        data = request.json
        smiles = data['smiles']
        identifier = data['identifier']
        
        rdkit_xyz = tool.rdkit_conversion(smiles, f"{identifier}_rdkit")
        obabel_xyz = tool.openbabel_conversion(smiles, f"{identifier}_obabel")
        networkx_xyz = tool.networkx_conversion(smiles, f"{identifier}_networkx")
        
        return jsonify({
            "rdkit": rdkit_xyz,
            "openbabel": obabel_xyz,
            "networkx": networkx_xyz
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    app.run(port=5000, debug=True)
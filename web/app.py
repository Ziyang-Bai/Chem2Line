import base64
from io import BytesIO
from flask import Flask, render_template, request, jsonify
from rdkit import Chem
from rdkit.Chem import AllChem
from ctlcore import formula_to_bondline, load_smiles_database, get_smiles_options  # 使用 ctlcore.py 中的逻辑

app = Flask(__name__)

# 加载 SMILES 数据库
smiles_dict = load_smiles_database("lib/db/default_database.xml")

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/molecule", methods=["POST"])
def molecule():
    data = request.json
    input_text = data.get("smiles", "C")
    smiles_list = get_smiles_options(input_text, smiles_dict)
    if not smiles_list:
        return jsonify({"error": "No matching SMILES found"}), 400

    selected_smiles = data.get("selected_smiles")
    if selected_smiles and selected_smiles in smiles_list:
        smiles = selected_smiles
    else:
        smiles = smiles_list[0]

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES string"}), 400

    # Generate 3D coordinates
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)

    # Extract atom positions and bonds
    conf = mol.GetConformer()
    atoms = [{"symbol": atom.GetSymbol(), "x": conf.GetAtomPosition(i).x, 
              "y": conf.GetAtomPosition(i).y, "z": conf.GetAtomPosition(i).z}
             for i, atom in enumerate(mol.GetAtoms())]
    bonds = []
    for bond in mol.GetBonds():
        start_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        start_pos = conf.GetAtomPosition(start_idx)
        end_pos = conf.GetAtomPosition(end_idx)

        # 计算键的方向向量
        dx = end_pos.x - start_pos.x
        dy = end_pos.y - start_pos.y
        dz = end_pos.z - start_pos.z

        bonds.append({
            "start": {"x": start_pos.x, "y": start_pos.y, "z": start_pos.z},
            "end": {"x": end_pos.x, "y": end_pos.y, "z": end_pos.z},
            "direction": {"x": dx, "y": dy, "z": dz},
            "type": int(bond.GetBondTypeAsDouble())
        })

    return jsonify({"atoms": atoms, "bonds": bonds, "smiles_options": smiles_list, "used_smiles": smiles})

@app.route("/smiles_svg", methods=["POST"])
def smiles_svg():
    data = request.json
    input_text = data.get("smiles", "C")
    smiles_list = get_smiles_options(input_text, smiles_dict)
    if not smiles_list:
        return jsonify({"error": "No matching SMILES found"}), 400

    selected_smiles = data.get("selected_smiles")
    if selected_smiles and selected_smiles in smiles_list:
        smiles = selected_smiles
    else:
        smiles = smiles_list[0]
    try:
        svg = formula_to_bondline(smiles, image_type="SVG")
        return jsonify({"svg": svg, "smiles_options": smiles_list, "used_smiles": smiles})
    except Exception as e:
        return jsonify({"error": f"Failed to generate SVG: {str(e)}"}), 500

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)

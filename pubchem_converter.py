import warnings
from rdkit import Chem
from rdkit.Chem import AllChem
import xml.etree.ElementTree as ET

def process_sdf_and_generate_smiles(sdf_file, output_xml, batch_size=1000):
    # 禁用警告
    warnings.filterwarnings("ignore", category=UserWarning, module="rdkit")
    
    # 创建 XML 根元素
    root = ET.Element("smiles_database")

    # 读取 SDF 文件
    supplier = Chem.SDMolSupplier(sdf_file)
    
    # 处理分子的计数器
    molecule_count = 0

    # 遍历所有分子
    for i, mol in enumerate(supplier):
        if mol is None:
            continue  # 跳过无效分子
        
        # 输出当前正在处理的分子
        molecule_name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"Compound_{i+1}"
        mol_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        cas_number = mol.GetProp('CAS') if mol.HasProp('CAS') else "N/A"
        
        print(f"Processing molecule: {molecule_name} (Molecular Formula: {mol_formula}, CAS: {cas_number})")
        
        # 检查是否为有机分子（例如，检查是否含有C元素）
        if 'C' not in mol_formula:
            print(f"Skipping non-organic molecule: {molecule_name}")
            continue
        
        # 生成 SMILES
        try:
            smiles = Chem.MolToSmiles(mol)
            compound_element = ET.SubElement(root, "compound")
            
            # 添加信息到 XML
            name_element = ET.SubElement(compound_element, "name")
            smiles_element = ET.SubElement(compound_element, "smiles")
            formula_element = ET.SubElement(compound_element, "formula")
            cas_element = ET.SubElement(compound_element, "cas")
            
            name_element.text = molecule_name
            smiles_element.text = smiles
            formula_element.text = mol_formula
            cas_element.text = cas_number
        except Exception as e:
            print(f"Error processing molecule {molecule_name}: {e}")
            continue
        
        molecule_count += 1
        
        # 每处理1000个分子，写入一次文件
        if molecule_count % batch_size == 0:
            print(f"Processed {molecule_count} molecules, saving progress to file...")
            tree = ET.ElementTree(root)
            tree.write(output_xml)
    
    # 处理完所有分子后，写入文件
    print(f"Processing completed. Saving final data to {output_xml}")
    tree = ET.ElementTree(root)
    tree.write(output_xml)

# 使用示例
process_sdf_and_generate_smiles('pubchem_data.sdf', 'smiles_database.xml')

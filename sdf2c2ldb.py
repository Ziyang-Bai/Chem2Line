from rdkit import Chem
from rdkit.Chem import rdMolDescriptors  # 添加这一行
import xml.etree.ElementTree as ET

def sdf_to_xml(sdf_file, xml_output_file):
    # 创建XML根节点
    root = ET.Element("smiles_database")
    name = ET.SubElement(root, "name")
    name.text = "PubChem物质大全"
    publisher = ET.SubElement(root, "publisher")
    publisher.text = "Ziyang-Bai"
    description = ET.SubElement(root, "description")
    description.text = "这个数据库包含了PubChem截止2024年12月份的所有物质。"
    
    # 以二进制模式打开SDF文件
    with open(sdf_file, 'rb') as f:  # 这里改成'rb'模式
        suppl = Chem.ForwardSDMolSupplier(f)  # 使用流式读取方式加载SDF
        count = 0
        for mol in suppl:
            if mol is None:
                continue  # 如果分子无效，则跳过
            
            try:
                # 获取分子信息
                smiles = Chem.MolToSmiles(mol)  # 获取SMILES
                name = mol.GetProp("PUBCHEM_IUPAC_OPENEYE_NAME") if mol.HasProp("PUBCHEM_IUPAC_OPENEYE_NAME") else "Unknown Name"
                formula = rdMolDescriptors.CalcMolFormula(mol)  # 获取分子式，改为使用rdMolDescriptors.CalcMolFormula
                cas = mol.GetProp("PUBCHEM_IUPAC_CAS_NAME") if mol.HasProp("PUBCHEM_IUPAC_CAS_NAME") else "Unknown CAS"
                print(f"OK: {name}-{smiles}-{formula}-{cas}")
            except Exception as e:
                print(f"Error processing molecule: {e}")
                continue  # 跳过无法处理的分子
            
            # 创建XML节点
            compound = ET.SubElement(root, "compound")
            
            # 添加分子信息到XML节点
            compound_name = ET.SubElement(compound, "name")
            compound_name.text = name
            compound_smiles = ET.SubElement(compound, "smiles")
            compound_smiles.text = smiles
            compound_formula = ET.SubElement(compound, "formula")
            compound_formula.text = formula
            compound_cas = ET.SubElement(compound, "cas")
            compound_cas.text = cas

            count += 1
            if count % 1000 == 0:
                # 每处理1000个分子追加保存到文件
                tree = ET.ElementTree(root)
                with open(xml_output_file, 'ab') as xml_file:
                    tree.write(xml_file, encoding="utf-8", xml_declaration=False)
                print(f"File Save: {xml_output_file}")

    # 创建最终的XML树并写入文件
    tree = ET.ElementTree(root)
    tree.write(xml_output_file, encoding="utf-8", xml_declaration=True)

    print(f"File Saved: {xml_output_file}")

# 调用函数，指定SDF文件路径和XML输出路径
sdf_file = 'pubchemfull.sdf'  # 请替换为你的SDF文件路径
xml_output_file = 'output.xml'    # 输出XML文件路径
sdf_to_xml(sdf_file, xml_output_file)

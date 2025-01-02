import warnings
import tqdm
from rdkit import Chem
from rdkit.Chem import PandasTools
import pandas as pd

# 禁用警告
warnings.filterwarnings("ignore", category=UserWarning, module="rdkit")

def load_sdf(file_path):
    """
    逐分子加载SDF文件并返回一个生成器
    """
    sdf_supplier = Chem.SDMolSupplier(file_path)
    return sdf_supplier  # 返回一个分子生成器，而不是整个DataFrame

def search_by_name(file_path, name):
    # 按名称查找分子
    results = []
    # 获取分子总数以设置进度条
    sdf_supplier = Chem.SDMolSupplier(file_path)
    total_molecules = sum(1 for _ in sdf_supplier)  # 先遍历一次文件计算总数
    sdf_supplier = Chem.SDMolSupplier(file_path)  # 重新加载SDF文件
    # 显示进度条
    with tqdm.tqdm(total=total_molecules, desc="正在查询...", unit="分子") as pbar:
        for mol in sdf_supplier:
            if mol and 'Name' in mol.GetPropsAsDict() and name.lower() in mol.GetProp('Name').lower():
                results.append(mol)
            pbar.update(1)
    return results

def search_by_formula(file_path, formula):
    # 按化学式查找分子
    results = []
    # 获取分子总数以设置进度条
    sdf_supplier = Chem.SDMolSupplier(file_path)
    total_molecules = sum(1 for _ in sdf_supplier)  # 先遍历一次文件计算总数
    sdf_supplier = Chem.SDMolSupplier(file_path)  # 重新加载SDF文件
    # 显示进度条
    with tqdm.tqdm(total=total_molecules, desc="正在查询...", unit="分子") as pbar:
        for mol in sdf_supplier:
            if mol and formula.lower() in Chem.rdMolDescriptors.CalcMolFormula(mol).lower():
                results.append(mol)
            pbar.update(1)
    return results

def search_by_cas(file_path, cas):
    # 按CAS号查找分子
    results = []
    # 获取分子总数以设置进度条
    sdf_supplier = Chem.SDMolSupplier(file_path)
    total_molecules = sum(1 for _ in sdf_supplier)  # 先遍历一次文件计算总数
    sdf_supplier = Chem.SDMolSupplier(file_path)  # 重新加载SDF文件
    # 显示进度条
    with tqdm.tqdm(total=total_molecules, desc="正在查询...", unit="分子") as pbar:
        for mol in sdf_supplier:
            if mol and 'CAS' in mol.GetPropsAsDict() and cas in mol.GetProp('CAS'):
                results.append(mol)
            pbar.update(1)
    return results

def display_results(results):
    if not results:
        print("没有找到匹配的分子")
    else:
        for mol in results:
            print(f"名称: {mol.GetProp('Name')}")
            print(f"分子式: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
            print(f"CAS号: {mol.GetProp('CAS')}")
            print(f"SMILES: {Chem.MolToSmiles(mol)}")
            print(f"分子量: {Chem.rdMolDescriptors.CalcExactMolWt(mol)}")
            print("-" * 40)

def main():
    # 加载SDF文件
    file_path = "pubchem_data.sdf"
    print("正在加载SDF文件，请稍候...")
    
    # 定义菜单选项
    options = ["名称", "化学式", "CAS号"]
    
    while True:
        print("\n请选择查询方式：")
        for idx, option in enumerate(options):
            print(f"{idx + 1}. {option}")

        choice = input("请输入选项(1/2/3): ")

        if choice == '1':
            query = "请输入分子名称: "
            query_func = search_by_name
        elif choice == '2':
            query = "请输入化学式: "
            query_func = search_by_formula
        elif choice == '3':
            query = "请输入CAS号: "
            query_func = search_by_cas
        else:
            print("无效的选项，请重新选择。")
            continue

        user_input = input(query).strip()

        print("正在查询，请稍候...")
        # 执行查询并显示结果
        results = query_func(file_path, user_input)
        display_results(results)

if __name__ == "__main__":
    main()

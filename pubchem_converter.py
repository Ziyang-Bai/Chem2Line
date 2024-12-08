import xml.etree.ElementTree as ET
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from tqdm import tqdm

# 定义用于处理单个 XML 片段的函数
def process_chunk(chunk):
    """解析XML片段并提取所需数据"""
    try:
        # 为确保每个块是有效的XML，添加一个根元素
        chunk = f'<root>{chunk}</root>'
        
        # 解析每个小片段
        root = ET.fromstring(chunk)
        compounds = []
        
        # 查找所有 <PC-Compound> 元素
        for elem in root.findall('.//PC-Compound'):
            # 提取 Name 和 SMILES 数据
            name_element = elem.find('.//PC-InfoData[PC-Urn/PC-Urn_label="Name"]/PC-InfoData_value/PC-InfoData_value_sval')
            smiles_element = elem.find('.//PC-InfoData[PC-Urn/PC-Urn_label="SMILES"]/PC-InfoData_value/PC-InfoData_value_sval')

            # 如果找到元素，提取并存储信息
            if name_element is not None and smiles_element is not None:
                compounds.append({
                    'name': name_element.text,
                    'smiles': smiles_element.text
                })
        
        print(f"Processed chunk with {len(compounds)} compounds.")  # 调试信息
        return compounds
    except Exception as e:
        print(f"Error processing chunk: {e}")
        return []

# 将大文件按行读取并分割为块
def split_file(file_path, chunk_size=100000):
    """按指定大小分割大文件，每个chunk为指定字符数"""
    with open(file_path, 'r', encoding='utf-8') as f:
        chunk = ""
        for line in f:
            chunk += line
            
            # 确保每个块的结束处是一个完整的 XML 元素
            if len(chunk) > chunk_size and chunk.endswith('</root>'):
                print(f"Yielding chunk of size {len(chunk)}")  # 输出调试信息
                yield chunk
                chunk = ""
        if chunk:
            print(f"Yielding final chunk of size {len(chunk)}")  # 输出调试信息
            yield chunk

def main():
    # 输入文件路径和输出文件路径
    input_file = 'Compound_000000001_000500000.xml'
    output_file = 'formatted_smiles_database.xml'
    
    # 使用多进程处理
    with ProcessPoolExecutor() as executor:
        futures = []
        chunk_generator = split_file(input_file)
        
        # 向池中提交任务，逐块处理文件
        for chunk in chunk_generator:
            futures.append(executor.submit(process_chunk, chunk))
        
        results = []
        with tqdm(total=len(futures)) as pbar:
            for future in as_completed(futures):
                try:
                    result = future.result()
                    results.extend(result)
                except Exception as e:
                    print(f"Error processing a chunk: {e}")
                pbar.update(1)

    # 将结果写入输出文件
    output_root = ET.Element('smiles_database')
    for compound in results:
        compound_element = ET.SubElement(output_root, 'compound')
        name = ET.SubElement(compound_element, 'name')
        name.text = compound['name']
        smiles = ET.SubElement(compound_element, 'smiles')
        smiles.text = compound['smiles']

    # 写入新的 XML 文件
    tree = ET.ElementTree(output_root)
    tree.write(output_file, encoding='utf-8', xml_declaration=True)

    print(f"转换完成，结果已保存到 '{output_file}'")

# 确保代码在主程序中运行
if __name__ == '__main__':
    main()

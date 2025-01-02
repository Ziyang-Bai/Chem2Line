import xml.dom.minidom

def format_xml(input_file, output_file):
	# 读取 XML 文件
	print("reading")
	with open(input_file, 'r', encoding='utf-8') as f:
		xml_content = f.read()
	print("thinking")
	# 使用 minidom 解析 XML
	dom = xml.dom.minidom.parseString(xml_content)
	print("formatting")
	# 格式化 XML 内容
	formatted_xml = dom.toprettyxml(indent="  ")  # 使用两个空格缩进
	print("writing")
	# 将格式化后的 XML 写入输出文件
	with open(output_file, 'w', encoding='utf-8') as f:
		f.write(formatted_xml)
	print(f"XML 文件已格式化并保存为 {output_file}")

# 示例用法
input_xml_file = 'smiles_database.xml'  # 输入的 XML 文件路径
output_xml_file = 'formatted_output.xml'  # 输出的格式化 XML 文件路径

format_xml(input_xml_file, output_xml_file)

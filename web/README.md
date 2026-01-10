# Chem2Line 轻量级版本 - 前后端分离架构

## 快速开始

### 后端启动

```bash
# 1. 进入 web 目录
cd web

# 2. 安装依赖（如果未安装）
pip install -r requirements.txt

# 3. 启动 Flask 服务
python app.py
```

服务会运行在 `http://localhost:5000`

### 前端访问

打开浏览器访问 `http://localhost:5000/`

## 架构说明

### 前端（index.html）
- **SMILES 解析 + 2D 绘图**：使用 RDKit.js minimal
- **3D 分子渲染**：使用 Three.js
- **API 调用**：`POST /api/3d`

### 后端（app.py）
- **3D 坐标生成**：RDKit EmbedMolecule()
- **力场优化**：RDKit MMFF()
- **返回数据**：atoms[] + bonds[] 坐标

## API 文档

### `POST /api/3d`

**请求**：
```json
{
  "smiles": "CCO"
}
```

**响应**：
```json
{
  "success": true,
  "atoms": [
    {
      "idx": 0,
      "symbol": "C",
      "atomic_num": 6,
      "x": 0.1,
      "y": 0.2,
      "z": 0.3
    }
  ],
  "bonds": [
    {
      "start_idx": 0,
      "end_idx": 1,
      "start": {"x": ..., "y": ..., "z": ...},
      "end": {"x": ..., "y": ..., "z": ...},
      "direction": {"x": ..., "y": ..., "z": ...},
      "type": 1,
      "aromatic": false
    }
  ],
  "smiles": "CCO",
  "atom_count": 3,
  "bond_count": 2
}
```

## 性能对比

| 项目 | 旧方案 (Python 完整) | 新方案 (精简) |
|------|--------|--------|
| 后端内存占用 | ~180MB | ~50MB |
| 初始化时间 | 2-3s | <1s |
| 单个请求耗时 | 0.5-2s | 0.2-0.8s |
| 支持并发 | 差 | 好 |

## 功能清单

✅ SMILES 解析和验证  
✅ 2D 键线式绘图  
✅ 3D 坐标生成  
✅ 力场优化  
✅ 多键支持  
✅ 芳香键显示  
✅ 实时缩放和旋转  

## 常见问题

**Q: 为什么 2D 由前端生成？**  
A: RDKit.js 足以处理 SMILES 解析和 2D 绘图，后端专注 3D，减少服务器负担。

**Q: 能否离线使用？**  
A: 可以，RDKit.js 是 WASM 本地执行，但后端 API 仍需要网络。

**Q: 能否直接在桌面版中使用？**  
A: 可以，只需将 Flask 启动改为内嵌服务，或改用 Electron。

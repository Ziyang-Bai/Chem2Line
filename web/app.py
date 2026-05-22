from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import inchi
import requests
import logging
import urllib.parse
from collections import deque
from datetime import datetime

app = Flask(__name__, static_folder='.', static_url_path='')
CORS(app)

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

class InMemoryLogHandler(logging.Handler):
    def __init__(self, max_lines=100):
        super().__init__()
        self.logs = deque(maxlen=max_lines)

    def emit(self, record):
        try:
            msg = self.format(record)
            timestamp = datetime.now().strftime('%H:%M:%S')
            self.logs.append(f"[{timestamp}] {record.levelname}: {msg}")
        except Exception:
            pass

log_handler = InMemoryLogHandler(max_lines=100)
log_handler.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(log_handler)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
PUBCHEM_USER_AGENT = "Chem2Line/1.0 (+https://github.com/Ziyang-Bai/Chem2Line)"
DEBUG_MODE = False


def set_debug_mode(enabled: bool):
    global DEBUG_MODE
    DEBUG_MODE = bool(enabled)
    target_level = logging.DEBUG if DEBUG_MODE else logging.WARNING
    logger.setLevel(target_level)
    app.logger.setLevel(target_level)
    logging.getLogger().setLevel(target_level)
    logger.warning(f"调试模式已{'开启' if DEBUG_MODE else '关闭'}，日志级别={logging.getLevelName(target_level)}")

@app.route("/")
def serve_index():
    return send_from_directory('.', 'index.html')

@app.route("/api/3d", methods=["POST"])
def generate_3d():
    try:
        data = request.get_json(silent=True) or {}
        smiles = data.get("smiles", "").strip()
        if DEBUG_MODE:
            logger.debug(f"[/api/3d] payload_keys={list((data or {}).keys()) if isinstance(data, dict) else 'n/a'}")
        
        if not smiles:
            return jsonify({"error": "SMILES 不能为空", "success": False}), 400
        
        logger.info(f"处理 SMILES: {smiles}")
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({"error": f"无效的 SMILES: {smiles}", "success": False}), 400
        
        logger.info(f"SMILES 解析成功，原子数: {mol.GetNumAtoms()}")
        
        mol = Chem.AddHs(mol)
        
        try:
            AllChem.EmbedMolecule(mol, randomSeed=42)
        except Exception as e:
            logger.warning(f"EmbedMolecule 失败: {e}，使用 2D 坐标")
            AllChem.Compute2DCoords(mol)
        
        try:
            AllChem.MMFFOptimizeMolecule(mol)
            logger.info("力场优化完成")
        except Exception as e:
            logger.warning(f"力场优化失败: {e}")
        
        conf = mol.GetConformer()
        atoms = []
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            atoms.append({
                "idx": i,
                "symbol": atom.GetSymbol(),
                "atomic_num": atom.GetAtomicNum(),
                "x": float(pos.x),
                "y": float(pos.y),
                "z": float(pos.z)
            })
        
        bonds = []
        for bond in mol.GetBonds():
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            start_pos = conf.GetAtomPosition(start_idx)
            end_pos = conf.GetAtomPosition(end_idx)
            
            bonds.append({
                "start_idx": start_idx,
                "end_idx": end_idx,
                "start": {
                    "x": float(start_pos.x),
                    "y": float(start_pos.y),
                    "z": float(start_pos.z)
                },
                "end": {
                    "x": float(end_pos.x),
                    "y": float(end_pos.y),
                    "z": float(end_pos.z)
                },
                "direction": {
                    "x": float(end_pos.x - start_pos.x),
                    "y": float(end_pos.y - start_pos.y),
                    "z": float(end_pos.z - start_pos.z)
                },
                "type": float(bond.GetBondTypeAsDouble()),
                "aromatic": bond.GetIsAromatic()
            })
        
        logger.info(f"返回结果: {len(atoms)} 个原子，{len(bonds)} 个键")
        if DEBUG_MODE:
            logger.debug(f"[/api/3d] smiles={smiles} atom_count={len(atoms)} bond_count={len(bonds)}")
        
        return jsonify({
            "success": True,
            "atoms": atoms,
            "bonds": bonds,
            "smiles": smiles,
            "atom_count": len(atoms),
            "bond_count": len(bonds)
        })
    
    except Exception as e:
        logger.error(f"错误: {str(e)}")
        return jsonify({"error": str(e), "success": False}), 500


@app.route("/api/name", methods=["POST"])
def generate_name():
    try:
        data = request.get_json(silent=True) or {}
        smiles = (data.get("smiles") or "").strip()
        if DEBUG_MODE:
            logger.debug(f"[/api/name] smiles={smiles}")
        if not smiles:
            return jsonify({"error": "SMILES 不能为空", "success": False}), 400

        # 解析 SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({"error": f"无效的 SMILES: {smiles}", "success": False}), 400

        # 尝试生成 InChI / InChIKey（若 RDKit 支持）
        inchistr = None
        inchikey = None
        try:
            inchistr = inchi.MolToInchi(mol)
            inchikey = inchi.MolToInchiKey(mol)
            if DEBUG_MODE:
                logger.debug(f"[/api/name] inchikey={inchikey}")
        except Exception as e:
            logger.warning(f"生成 InChI/Key 失败: {e}")

        result = {"success": False, "source": None, "iupac_name": None, "synonyms": [], "inchI": inchistr, "inchIKey": inchikey}

        # 首选：用 InChIKey 查询 PubChem
        session = requests.Session()
        session.headers.update({"User-Agent": "Chem2Line/1.0 (+https://github.com/Ziyang-Bai/Chem2Line)"})
        try:
            if inchikey:
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/IUPACName/JSON"
                r = session.get(url, timeout=8)
                if r.ok:
                    j = r.json()
                    name = j.get('PropertyTable', {}).get('Properties', [{}])[0].get('IUPACName')
                    if name:
                        result.update({"success": True, "source": "pubchem_inchikey", "iupac_name": name})
                        if DEBUG_MODE:
                            logger.debug("[/api/name] 命中 pubchem_inchikey")
                # synonyms
                if not result['success']:
                    url2 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/synonyms/JSON"
                    r2 = session.get(url2, timeout=8)
                    if r2.ok:
                        j2 = r2.json()
                        syns = j2.get('InformationList', {}).get('Information', [{}])[0].get('Synonym', [])
                        if syns:
                            result.update({"success": True, "source": "pubchem_inchikey_synonyms", "synonyms": syns, "iupac_name": syns[0]})
                            if DEBUG_MODE:
                                logger.debug("[/api/name] 命中 pubchem_inchikey_synonyms")

            # 回退：用 SMILES 查询 PubChem
            if not result['success']:
                import urllib.parse
                esmi = urllib.parse.quote(smiles, safe='')
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{esmi}/property/IUPACName/JSON"
                r = session.get(url, timeout=8)
                if r.ok:
                    j = r.json()
                    name = j.get('PropertyTable', {}).get('Properties', [{}])[0].get('IUPACName')
                    if name:
                        result.update({"success": True, "source": "pubchem_smiles", "iupac_name": name})
                        if DEBUG_MODE:
                            logger.debug("[/api/name] 命中 pubchem_smiles")
                if not result['success']:
                    url2 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{esmi}/synonyms/JSON"
                    r2 = session.get(url2, timeout=8)
                    if r2.ok:
                        j2 = r2.json()
                        syns = j2.get('InformationList', {}).get('Information', [{}])[0].get('Synonym', [])
                        if syns:
                            result.update({"success": True, "source": "pubchem_smiles_synonyms", "synonyms": syns, "iupac_name": syns[0]})
                            if DEBUG_MODE:
                                logger.debug("[/api/name] 命中 pubchem_smiles_synonyms")

        except Exception as e:
            logger.warning(f"调用 PubChem 失败: {e}")

        # 最后回退：返回 InChI 或 InChIKey 作为标准标识（若没有名字）
        if not result['success']:
            if inchikey:
                result.update({"success": True, "source": "inchikey_fallback", "iupac_name": inchikey})
            elif inchistr:
                result.update({"success": True, "source": "inchi_fallback", "iupac_name": inchistr})
            else:
                # 作为最后回退，返回规范 SMILES
                try:
                    can = Chem.MolToSmiles(mol, isomericSmiles=True)
                except Exception:
                    can = smiles
                result.update({"success": True, "source": "smiles_fallback", "iupac_name": can})

        result["smiles"] = smiles
        if DEBUG_MODE:
            logger.debug(f"[/api/name] result_source={result.get('source')} success={result.get('success')}")
        return jsonify(result)

    except Exception as e:
        logger.error(f"/api/name 错误: {e}")
        return jsonify({"error": str(e), "success": False}), 500


@app.route("/api/debug-mode", methods=["GET", "POST"])
def debug_mode():
    """Toggle backend debug logging at runtime."""
    if request.method == "GET":
        return jsonify({"success": True, "debug_mode": DEBUG_MODE, "log_level": logging.getLevelName(logger.level)})

    payload = request.get_json(silent=True) or {}
    enabled = payload.get("enabled", False)
    set_debug_mode(bool(enabled))
    return jsonify({"success": True, "debug_mode": DEBUG_MODE, "log_level": logging.getLevelName(logger.level)})


@app.route("/api/get-logs", methods=["GET"])
def get_logs():
    """Get recent backend logs for frontend display."""
    logs = list(log_handler.logs)
    return jsonify({"success": True, "logs": logs})


@app.route("/api/pubchem/proxy", methods=["POST"])
def pubchem_proxy():
    """PubChem fallback proxy for browser-side direct fetch failures.

    Body example:
    {
      "path": "compound/name/ethanol/property/IsomericSMILES,SMILES,ConnectivitySMILES/JSON",
      "query": {"MaxRecords": 20},
      "timeout": 12
    }
    """
    try:
        payload = request.get_json(silent=True) or {}
        path = (payload.get("path") or "").strip().lstrip("/")
        query = payload.get("query") or {}
        timeout = payload.get("timeout", 12)

        if not path:
            return jsonify({"error": "path 不能为空", "success": False}), 400
        if path.startswith("http://") or path.startswith("https://"):
            return jsonify({"error": "path 不能是完整 URL", "success": False}), 400
        if not path.endswith("/JSON"):
            return jsonify({"error": "仅允许请求 JSON 接口（path 必须以 /JSON 结尾）", "success": False}), 400

        if not isinstance(query, dict):
            return jsonify({"error": "query 必须是对象", "success": False}), 400

        try:
            timeout = int(timeout)
        except Exception:
            timeout = 12
        timeout = min(max(timeout, 3), 25)

        target_url = urllib.parse.urljoin(PUBCHEM_BASE, path)
        parsed = urllib.parse.urlparse(target_url)
        if parsed.scheme != "https" or parsed.netloc != "pubchem.ncbi.nlm.nih.gov":
            return jsonify({"error": "非法目标地址", "success": False}), 400

        if DEBUG_MODE:
            logger.debug(f"[pubchem_proxy] path={path} query={query} timeout={timeout}")

        headers = {"User-Agent": PUBCHEM_USER_AGENT}
        resp = requests.get(target_url, params=query, headers=headers, timeout=timeout)

        if DEBUG_MODE:
            logger.debug(f"[pubchem_proxy] target_url={resp.url} status={resp.status_code}")

        try:
            data = resp.json()
            return jsonify({
                "success": resp.ok,
                "status": resp.status_code,
                "data": data,
                "url": target_url
            }), resp.status_code
        except ValueError:
            return jsonify({
                "success": False,
                "status": resp.status_code,
                "error": "PubChem 返回了非 JSON 响应",
                "text": (resp.text or "")[:500],
                "url": target_url
            }), 502

    except requests.RequestException as e:
        logger.warning(f"PubChem 代理请求失败: {e}")
        return jsonify({"error": f"PubChem 代理请求失败: {e}", "success": False}), 502
    except Exception as e:
        logger.error(f"/api/pubchem/proxy 错误: {e}")
        return jsonify({"error": str(e), "success": False}), 500

@app.route("/health", methods=["GET"])
def health():
    return jsonify({"status": "ok", "message": "Chem2Line 3D 后端正常运行"})

@app.errorhandler(404)
def not_found(error):
    return jsonify({"error": "端点不存在", "success": False}), 404

@app.errorhandler(500)
def internal_error(error):
    return jsonify({"error": "服务器错误", "success": False}), 500

if __name__ == "__main__":
    logger.info("启动 Chem2Line 3D 后端服务...")
    app.run(debug=False, host="0.0.0.0", port=5000, threaded=True)

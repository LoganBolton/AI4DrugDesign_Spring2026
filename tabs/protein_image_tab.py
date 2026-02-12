import io
import re
import requests
from PIL import Image

try:
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

UNIPROT_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", re.IGNORECASE)
PDB_RE = re.compile(r"^[0-9][A-Z0-9]{3}$", re.IGNORECASE)


def detect_id_type(pid):
    pid = pid.strip().upper()
    if UNIPROT_RE.match(pid):
        return "uniprot"
    if PDB_RE.match(pid):
        return "pdb"
    return "search"


def fetch_uniprot(uid):
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}.json", timeout=15)
    r.raise_for_status()
    d = r.json()
    name = d.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", uid)
    organism = d.get("organism", {}).get("scientificName", "Unknown")
    fn = " ".join(
        t.get("value", "") for c in d.get("comments", [])
        if c.get("commentType") == "FUNCTION"
        for t in c.get("texts", [])
    )
    pubchem_cids = [x.get("id") for x in d.get("dbCrossReferences", []) if x.get("database") == "PubChem"]
    return {"name": name, "organism": organism, "function": fn or "N/A", "uniprot_id": uid.upper(), "pubchem_cids": pubchem_cids, "chem_ids": []}


def search_uniprot(query):
    r = requests.get("https://rest.uniprot.org/uniprotkb/search", params={"query": query, "format": "json", "size": 1}, timeout=15)
    r.raise_for_status()
    results = r.json().get("results", [])
    if not results:
        raise ValueError(f"No UniProt entry found for '{query}'")
    return fetch_uniprot(results[0]["primaryAccession"])


def fetch_pdb(pid):
    r = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{pid}", timeout=15)
    r.raise_for_status()
    d = r.json()
    name = d.get("struct", {}).get("title", pid)
    lr = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{pid}/components", timeout=10)
    skip = {"HOH", "SO4", "PO4", "GOL", "EDO", "ACT", None}
    chem_ids = [c.get("chem_comp", {}).get("id") for c in (lr.json() if lr.ok else []) if c.get("chem_comp", {}).get("id") not in skip]
    return {"name": name, "organism": "N/A", "function": f"PDB structure {pid}", "pdb_id": pid, "chem_ids": chem_ids[:5], "pubchem_cids": []}


def smiles_from_pubchem(cid):
    try:
        r = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES/JSON", timeout=10)
        props = r.json().get("PropertyTable", {}).get("Properties", [])
        return props[0].get("IsomericSMILES") if props else None
    except Exception:
        return None


def smiles_from_pdb_ligand(chem_id):
    try:
        r = requests.get(f"https://data.rcsb.org/rest/v1/core/chemcomp/{chem_id}", timeout=10)
        d = r.json().get("rcsb_chem_comp_descriptor", {})
        return d.get("smiles_stereo") or d.get("smiles")
    except Exception:
        return None


def smiles_from_chembl(uniprot_id):
    try:
        tr = requests.get("https://www.ebi.ac.uk/chembl/api/data/target.json", params={"target_components__accession": uniprot_id, "limit": 1}, timeout=10)
        targets = tr.json().get("targets", [])
        if not targets:
            return None, None
        mr = requests.get("https://www.ebi.ac.uk/chembl/api/data/molecule.json", params={"molecule_type": "Small molecule", "limit": 5}, timeout=10)
        for mol in mr.json().get("molecules", []):
            smi = mol.get("molecule_structures", {}).get("canonical_smiles")
            name = mol.get("pref_name") or mol.get("molecule_chembl_id", "Unknown")
            if smi:
                return smi, name
    except Exception:
        pass
    return None, None


def render_image(smiles):
    if RDKIT_AVAILABLE:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                AllChem.Compute2DCoords(mol)
                return Draw.MolToImage(mol, size=(600, 400))
        except Exception:
            pass
    try:
        r = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/PNG", params={"smiles": smiles}, timeout=15)
        if r.ok:
            return Image.open(io.BytesIO(r.content))
    except Exception:
        pass
    return None


def lookup(protein_id):
    protein_id = protein_id.strip()
    if not protein_id:
        return "Enter a Protein ID.", None

    id_type = detect_id_type(protein_id)
    try:
        if id_type == "uniprot":
            info = fetch_uniprot(protein_id)
        elif id_type == "pdb":
            info = fetch_pdb(protein_id)
        else:
            info = search_uniprot(protein_id)
    except Exception as e:
        return f"Error fetching data for '{protein_id}': {e}", None

    smiles, compound_name = None, "Unknown"

    if info.get("pubchem_cids"):
        cid = info["pubchem_cids"][0]
        smiles = smiles_from_pubchem(cid)
        compound_name = f"PubChem CID {cid}"

    if not smiles and info.get("uniprot_id"):
        smiles, compound_name = smiles_from_chembl(info["uniprot_id"])

    if not smiles and info.get("chem_ids"):
        for cid in info["chem_ids"]:
            smiles = smiles_from_pdb_ligand(cid)
            if smiles:
                compound_name = cid
                break

    image = render_image(smiles) if smiles else None

    lines = [
        f"Protein: {info.get('name', protein_id)}",
        f"Organism: {info.get('organism', 'Unknown')}",
    ]
    if info.get("uniprot_id"):
        lines.append(f"UniProt: {info['uniprot_id']}")
    if info.get("pdb_id"):
        lines.append(f"PDB: {info['pdb_id']}")
    lines += ["", f"Function: {info.get('function', 'N/A')}", ""]
    if smiles:
        lines += [f"Compound: {compound_name}", f"SMILES: {smiles[:80]}{'...' if len(smiles) > 80 else ''}"]
    else:
        lines.append("No compound found. Try P00533 (EGFR) or 1IEP (Bcr-Abl).")

    return "\n".join(lines), image


def build_tab():
    import gradio as gr

    with gr.Column():
        gr.Markdown("## Protein → Compound Structure Viewer")
        with gr.Row():
            id_input = gr.Textbox(label="Protein ID", placeholder="e.g. P00533  or  1IEP  or  EGFR", scale=3)
            submit_btn = gr.Button("Look Up", variant="primary", scale=1)
        with gr.Row():
            info_output = gr.Textbox(label="Info", lines=12, interactive=False, scale=1)
            img_output = gr.Image(label="2D Compound Structure", type="pil", scale=1, height=420)
        gr.Examples(examples=[["P00533"], ["1IEP"], ["P0DTD1"], ["EGFR"]], inputs=id_input)

        submit_btn.click(fn=lookup, inputs=id_input, outputs=[info_output, img_output])
        id_input.submit(fn=lookup, inputs=id_input, outputs=[info_output, img_output])

import os
import gradio as gr
import requests


try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


OUTPUT_DIR = "docking_files"


def _fetch_protein_pdb(pdb_id: str) -> tuple[str, str]:
    """Download PDB file from RCSB and save locally. Returns (filepath, message)."""
    pdb_id = pdb_id.strip().upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url, timeout=20)

    if response.status_code == 404:
        raise ValueError(f"PDB ID '{pdb_id}' not found in the RCSB database.")
    response.raise_for_status()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    filepath = os.path.join(OUTPUT_DIR, f"{pdb_id}_protein.pdb")
    with open(filepath, "w") as f:
        f.write(response.text)

    # Count ATOM/HETATM lines as a basic sanity check
    atom_count = sum(1 for line in response.text.splitlines() if line.startswith(("ATOM", "HETATM")))
    return filepath, atom_count


def _smiles_to_pdb(smiles: str, pdb_id: str) -> tuple[str, str]:
    """Convert SMILES to a 3D PDB file using RDKit. Returns (filepath, inchi_key)."""
    if not RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is not installed. Cannot convert SMILES to 3D structure.")

    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise ValueError(f"Invalid SMILES string: '{smiles}'")

    # Add hydrogens and generate 3D coordinates
    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result == -1:
        raise RuntimeError("3D coordinate generation failed. Please check your SMILES string.")

    AllChem.MMFFOptimizeMolecule(mol)

    inchi_key = Chem.InchiInfo(Chem.MolToInchi(mol)).GetInchiKey() if hasattr(Chem, "InchiInfo") else "LIGAND"

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    filepath = os.path.join(OUTPUT_DIR, f"{pdb_id}_ligand.pdb")
    writer = Chem.PDBWriter(filepath)
    writer.write(mol)
    writer.close()

    return filepath, inchi_key


def _smiles_to_pdb_fallback(smiles: str, pdb_id: str) -> tuple[str, str]:
    """Simpler fallback that writes SMILES-derived PDB without InchiKey dependency."""
    mol = Chem.MolFromSmiles(smiles.strip())
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    filepath = os.path.join(OUTPUT_DIR, f"{pdb_id}_ligand.pdb")
    Chem.MolToPDBFile(mol, filepath)
    return filepath


def _prepare_docking(pdb_id: str, smiles: str) -> str:
    """
    Main handler: downloads protein PDB and converts SMILES to 3D ligand PDB.
    Returns a formatted results string for display in the UI.
    """
    # --- Input validation ---
    pdb_id = pdb_id.strip() if pdb_id else ""
    smiles = smiles.strip() if smiles else ""

    if not pdb_id:
        return "❌ Error: Please enter a PDB ID."
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        return f"❌ Error: '{pdb_id}' does not look like a valid PDB ID (should be 4 alphanumeric characters, e.g. 6LU7)."
    if not smiles:
        return "❌ Error: Please enter a SMILES string."

    lines = []
    protein_path = None
    ligand_path = None
    protein_atoms = 0

    # --- Step 1: Fetch protein ---
    lines.append("⏳ Fetching protein structure from RCSB...")
    try:
        protein_path, protein_atoms = _fetch_protein_pdb(pdb_id)
        lines.append(f"✅ Protein PDB downloaded successfully.")
        lines.append(f"   File: {protein_path}")
        lines.append(f"   Atom/HETATM records: {protein_atoms}")
    except Exception as e:
        lines.append(f"❌ Protein fetch failed: {e}")
        return "\n".join(lines)

    # --- Step 2: Convert SMILES to 3D ligand PDB ---
    lines.append("")
    lines.append("⏳ Converting SMILES to 3D ligand structure...")
    try:
        ligand_path = _smiles_to_pdb_fallback(smiles, pdb_id.upper())
        ligand_size = os.path.getsize(ligand_path)
        lines.append(f"✅ Ligand PDB generated successfully.")
        lines.append(f"   File: {ligand_path}")
        lines.append(f"   File size: {ligand_size} bytes")
    except Exception as e:
        lines.append(f"❌ Ligand preparation failed: {e}")
        return "\n".join(lines)

    # --- Step 3: Summary and next steps ---
    lines.append("")
    lines.append("=" * 50)
    lines.append("✅ PREPARATION COMPLETE")
    lines.append("=" * 50)
    lines.append(f"  Protein:  {protein_path}")
    lines.append(f"  Ligand:   {ligand_path}")
    lines.append("")
    lines.append("📋 NEXT STEPS FOR AUTODOCK VINA:")
    lines.append("")
    lines.append("1. Convert files to PDBQT format using MGLTools:")
    lines.append(f"   python prepare_receptor4.py -r {protein_path} -o {pdb_id.upper()}_receptor.pdbqt -A hydrogens")
    lines.append(f"   python prepare_ligand4.py   -l {ligand_path} -o {pdb_id.upper()}_ligand.pdbqt")
    lines.append("")
    lines.append("   (Alternatively use the 'Molecular Docking' tab which handles this automatically.)")
    lines.append("")
    lines.append("2. Define the docking search box in your Vina config:")
    lines.append("   center_x = <X>   center_y = <Y>   center_z = <Z>")
    lines.append("   size_x = 20      size_y = 20      size_z = 20")
    lines.append("")
    lines.append("3. Run AutoDock Vina:")
    lines.append(f"   vina --receptor {pdb_id.upper()}_receptor.pdbqt \\")
    lines.append(f"        --ligand    {pdb_id.upper()}_ligand.pdbqt \\")
    lines.append("        --config    config.txt \\")
    lines.append("        --exhaustiveness 8 \\")
    lines.append(f"        --out       {pdb_id.upper()}_docked.pdbqt")
    lines.append("")
    lines.append("4. Visualize results in PyMOL or UCSF Chimera.")

    return "\n".join(lines)


def create_tab():
    with gr.Tab("Docking Preparation"):
        with gr.Row():
            with gr.Column(scale=1):
                dock_pdb_input = gr.Textbox(
                    label="Target Protein PDB ID",
                    placeholder="e.g., 6LU7",
                )
                dock_smiles_input = gr.Textbox(
                    label="Compound SMILES",
                    placeholder="e.g., CC(=O)OC1=CC=CC=C1C(=O)O",
                    lines=3,
                )
                prepare_btn = gr.Button(
                    "Prepare Files for Vina",
                    size="lg",
                    elem_classes=["analyze-btn"],
                )

            with gr.Column(scale=1):
                dock_output = gr.Textbox(
                    label="Preparation Results",
                    interactive=False,
                    lines=12,
                    placeholder=(
                        "a. Confirmation that Preparation was Successful\n"
                        "b. Protein PDB Saved / Ligand PDB Saved\n"
                        "c. Steps on Proceeding with AutoDock Vina"
                    ),
                )

        prepare_btn.click(
            _prepare_docking,
            inputs=[dock_pdb_input, dock_smiles_input],
            outputs=[dock_output],
        )

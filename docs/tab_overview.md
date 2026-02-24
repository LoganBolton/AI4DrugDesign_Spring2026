# Tab Overview

Quick reference for demoing each tab.

---

## Protein Analysis - Jase

- Takes a **PDB ID** (e.g. `1AZ5`), fetches protein metadata (title, organism, chains, resolution, ligands, binding sites) from the **RCSB PDB** and **PDBe** APIs.
- Sends that data to **OpenAI (GPT-5-nano)** to generate a detailed structural analysis: binding sites, key residues, interaction types, drug design considerations.
- Output: text report with raw protein info + AI-generated analysis.

## Compound Optimization - Sathvik

- Takes a **SMILES string** and a free-text **optimization goal** (e.g. "Increase solubility").
- Uses **RDKit** locally to calculate molecular properties (MW, LogP, H-bond donors/acceptors) and checks **Lipinski's Rule of Five**.
- Generates goal-specific optimization suggestions (no external API calls — all rule-based logic).
- Output: text report with property analysis + suggestions, and a 2D molecule image.

## Compound Visualization - Alex

- Takes a **SMILES string** (e.g. Aspirin, Caffeine).
- Uses **RDKit** locally to render a 2D structure image. No API calls.
- Output: 2D molecule image.

## Molecular Docking - Logan

- Takes a **PDB ID** (protein target), a **SMILES string** (ligand), and optional **binding box coordinates** (X, Y, Z).
- Downloads the PDB file from **RCSB**, converts the ligand to 3D with **RDKit**, prepares receptor/ligand files with **Meeko**, then runs **AutoDock Vina** to simulate docking.
- "Auto-Calculate Center" finds the co-crystallized ligand in the PDB file and uses its geometric center.
- Output: table of binding affinities (kcal/mol) for each pose + downloadable docked PDBQT file.

## Compound Structure Viewer - Dhvani

- Takes a **protein identifier** — UniProt ID (`P00533`), PDB ID (`1IEP`), or gene name (`EGFR`).
- Looks up the protein via **UniProt** or **RCSB PDB**, then searches for a known compound/ligand through **PubChem**, **ChEMBL**, or the PDB's ligand data.
- Renders the first compound it finds using **RDKit** (falls back to PubChem's image API).
- Output: protein info text + 2D image of an associated compound.

## Chat

- Free-form conversational interface powered by **OpenAI (GPT-5-nano)**.
- Maintains multi-turn conversation history.
- Output: AI-generated text responses.

---

## Data Sources & Tools

### External APIs

- **RCSB PDB** (https://data.rcsb.org) — The main public database of 3D protein structures. We fetch protein metadata (title, organism, resolution, ligands) and download `.pdb` structure files for docking. Used by: Protein Analysis, Molecular Docking, Compound Structure Viewer.

- **PDBe** (https://www.ebi.ac.uk/pdbe) — The European mirror of PDB, run by EMBL-EBI. We use it specifically to fetch binding site annotations. Used by: Protein Analysis.

- **UniProt** (https://rest.uniprot.org) — The main protein sequence and function database. We look up protein names, organisms, functions, and cross-references to compound databases. Used by: Compound Structure Viewer.

- **PubChem** (https://pubchem.ncbi.nlm.nih.gov) — NIH's database of chemical molecules. We fetch SMILES strings and fallback compound images from it. Used by: Compound Structure Viewer.

- **ChEMBL** (https://www.ebi.ac.uk/chembl) — A database of bioactive drug-like molecules with assay data. We search it for known compounds that bind to a given protein target. Used by: Compound Structure Viewer.


### Local Libraries

- **RDKit** — Open-source cheminformatics toolkit. Parses SMILES strings, calculates molecular properties (MW, LogP, etc.), generates 2D/3D molecular coordinates, and renders 2D structure images. Used by: Compound Visualization, Compound Optimization, Compound Structure Viewer, Molecular Docking.

- **AutoDock Vina** — Molecular docking engine that simulates how a small molecule (ligand) binds to a protein. Outputs binding affinity scores in kcal/mol. Used by: Molecular Docking.



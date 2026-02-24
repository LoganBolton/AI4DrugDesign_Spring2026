# Molecular Docking Tab

## Overview

The Molecular Docking tab lets you dock a small molecule (ligand) into a protein target using [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina). It takes a PDB ID and a SMILES string as input, and returns predicted binding affinities plus a downloadable file of docked poses.

## Prerequisites

### AutoDock Vina CLI

The `vina` binary must be installed and available on your `PATH`.

1. Download from [GitHub Releases](https://github.com/ccsb-scripps/AutoDock-Vina/releases) (v1.2.5+ for macOS ARM64)
2. Extract and move the `vina` executable somewhere on your PATH:
   ```bash
   # Example: move to /usr/local/bin
   sudo cp vina /usr/local/bin/vina
   sudo chmod +x /usr/local/bin/vina

   # Verify
   vina --version
   ```

### Python Dependencies

Installed automatically via `uv sync`:
- `rdkit` — SMILES parsing and 3D conformer generation
- `meeko` — Ligand and receptor PDBQT file preparation
- `numpy` — Coordinate math
- `scipy`, `gemmi` — Required by meeko

## How It Works

### Step 1: Enter a Target PDB ID

Enter a 4-character PDB ID from the [RCSB Protein Data Bank](https://www.rcsb.org/). This is the protein you want to dock against.

**Good examples:**
| PDB ID | Protein | Notes |
|--------|---------|-------|
| `6LU7` | SARS-CoV-2 main protease | Has co-crystallized inhibitor (PJE) |
| `1HSG` | HIV-1 protease | Has co-crystallized inhibitor (MK1) |
| `3HTB` | CDK2 kinase | Has co-crystallized ligand |

### Step 2: Enter a Ligand SMILES

Enter a SMILES string for the molecule you want to dock.

**Examples:**
| SMILES | Molecule |
|--------|----------|
| `CC(=O)Oc1ccccc1C(=O)O` | Aspirin |
| `CC(C)Cc1ccc(cc1)C(C)C(=O)O` | Ibuprofen |
| `c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34` | Pyrene |

### Step 3: Auto-Calculate Binding Box Center

Click **"Auto-Calculate Center from Co-Ligand"** to automatically determine where to focus the docking search.

**What it does:**
1. Fetches the PDB file from RCSB
2. Scans all HETATM records for non-solvent molecules (excludes water, ions, buffer molecules, etc.)
3. Selects the ligand with the most atoms (most likely the pharmacological ligand)
4. Computes the geometric center (mean X, Y, Z) of that ligand
5. Populates the Center X, Y, Z fields

**Expected output for `6LU7`:**
```
Center X: -11.575
Center Y:  14.612
Center Z:  65.164
```

You can also enter coordinates manually if you know the binding site location.

### Step 4: Run AutoDock Vina

Click **"Run AutoDock Vina"** to start the docking calculation.

**What happens behind the scenes:**
1. **Receptor preparation** — Downloads PDB, converts to PDBQT format using `mk_prepare_receptor.py` (meeko)
2. **Ligand preparation** — Converts SMILES to 3D coordinates (RDKit EmbedMolecule + MMFF optimization), then to PDBQT format (meeko)
3. **Vina docking** — Runs the Vina CLI with a 25x25x25 Angstrom search box centered on the coordinates from Step 3, exhaustiveness=8
4. **Results** — Parses the affinity table from Vina output

**Vina CLI command executed:**
```
vina --receptor receptor.pdbqt --ligand ligand.pdbqt \
     --center_x X --center_y Y --center_z Z \
     --size_x 25 --size_y 25 --size_z 25 \
     --exhaustiveness 8 --out docked_output.pdbqt
```

## Expected Output

### Binding Affinities

A table of docking modes ranked by predicted binding energy:

```
mode |   affinity   | dist from best mode
     | (kcal/mol)   | rmsd l.b.| rmsd u.b.
-----+----------+----------+----------
   1       -5.2          0.0        0.0
   2       -5.0          1.832      2.417
   3       -4.9          1.654      2.103
   4       -4.8          2.391      3.876
   5       -4.7          3.214      5.102
```

- **mode** — Pose rank (1 = best)
- **affinity (kcal/mol)** — Predicted binding energy. More negative = stronger predicted binding.
- **rmsd l.b.** — Lower-bound RMSD from the best mode
- **rmsd u.b.** — Upper-bound RMSD from the best mode

### Docked PDBQT File

A downloadable `.pdbqt` file containing all docked poses. This can be visualized in tools like PyMOL, UCSF Chimera, or AutoDockTools.

## Error Messages

| Error | Cause | Fix |
|-------|-------|-----|
| "Invalid PDB ID" | PDB ID is not 4 alphanumeric characters | Check the ID on [rcsb.org](https://www.rcsb.org/) |
| "PDB ID not found on RCSB" | PDB doesn't exist | Verify the PDB ID is correct |
| "No co-crystallized ligand found" | The PDB has no HETATM records besides solvent/ions | Use a structure with a bound ligand, or enter coordinates manually |
| "Invalid SMILES string" | RDKit can't parse the SMILES | Validate your SMILES (e.g. on [molview.org](https://molview.org/)) |
| "Failed to generate 3D coordinates" | RDKit can't embed the molecule in 3D | Try a simpler or smaller molecule |
| "Receptor preparation failed" | `mk_prepare_receptor.py` error | Check that the PDB file is a valid protein structure |
| "AutoDock Vina binary not found" | `vina` is not installed or not on PATH | See [Prerequisites](#prerequisites) |
| "Vina docking failed" | Vina returned an error | Check coordinates are reasonable; check stderr output |
| "Vina timed out" | Docking took longer than 300 seconds | Try a smaller molecule or lower exhaustiveness |

## Full Example Walkthrough

1. Start the app: `uv run python app.py`
2. Go to the **Molecular Docking** tab
3. Enter PDB ID: **`6LU7`**
4. Enter SMILES: **`CC(=O)Oc1ccccc1C(=O)O`** (aspirin)
5. Click **"Auto-Calculate Center from Co-Ligand"** → X/Y/Z fields fill in
6. Click **"Run AutoDock Vina"** → Affinities table appears + PDBQT file available for download

# AI4DrugDesign_Spring2026

A simple Gradio-based application for AI4DrugDesign course.

## Installation

1. **Install uv** (if not already installed):

   **macOS/Linux:**
   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

   **Windows (PowerShell):**
   ```powershell
   powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
   ```

   After installation, restart your terminal or run:
   - macOS/Linux: `source $HOME/.local/bin/env` or `source $HOME/.cargo/env`
   - Windows: Restart PowerShell

   Verify installation:
   ```bash
   uv --version
   ```


2. **Install system build dependencies** (required to compile AutoDock Vina):

   **macOS (Homebrew):**
   ```bash
   brew install swig boost
   ```

   **Ubuntu/Debian:**
   ```bash
   sudo apt install swig libboost-all-dev
   ```

3. **Install dependencies**:
   ```bash
   uv sync
   ```
   This will:
   - Create a virtual environment in `.venv/`
   - Install Python 3.12 (specified in `.python-version`)
   - Install all dependencies from `pyproject.toml` (including AutoDock Vina)

## Running the App

```bash
uv run python app.py
```

To activate the virtual environment manually:
- **macOS/Linux:** `source .venv/bin/activate`
- **Windows (PowerShell):** `.venv\Scripts\Activate.ps1`
- **Windows (CMD):** `.venv\Scripts\activate.bat`



## Compound Optimization (Sathvik)

### Overview
The Compound Optimization tab provides medicinal chemistry guidance for improving small-molecule drug candidates.

This module:
- Accepts a valid SMILES string
- Calculates molecular descriptors
- Evaluates drug-likeness using Lipinski’s Rule of Five
- Interprets user-defined optimization goals
- Generates structural modification suggestions
- Displays a 2D molecular visualization

### How to Use
1. Navigate to the **Compound Optimization** tab
2. Enter a valid SMILES string
3. Enter optimization goals
4. Click **Optimize Compound**

### Example Inputs

**Example 1 — Aspirin**  
SMILES:
```text
CC(=O)Oc1ccccc1C(=O)O
```
Optimization Goals:
```text
Increase solubility, reduce lipophilicity
```

**Example 2 — Acetaminophen**  
SMILES:
```text
CC(=O)NC1=CC=C(C=C1)O
```
Optimization Goals:
```text
Improve metabolic stability, reduce toxicity
```

**Example 3 — Caffeine**  
SMILES:
```text
Cn1cnc2n(C)c(=O)n(C)c(=O)c12
```
Optimization Goals:
```text
Improve BBB permeability, optimize logP
```

**Example 4 — Ibuprofen**  
SMILES:
```text
CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O
```
Optimization Goals:
```text
Reduce lipophilicity, improve aqueous solubility
```

**Example 5 — Benzene**  
SMILES:
```text
c1ccccc1
```
Optimization Goals:
```text
Increase polarity, improve solubility
```

Output Includes:
   - Molecular weight
   - LogP
   - Hydrogen bond donors
   - Hydrogen bond acceptors
   - Lipinski PASS/FAIL status
   - Goal-based structural suggestions
   - 2D molecule visualization
   - Valid Input Notes
   - SMILES must represent a valid small organic molecule.
   - Invalid syntax will return an error.
   - Optimization suggestions are rule-based medicinal chemistry guidance intended for early-stage design decisions.
# Compound Optimization (Sathvik)

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
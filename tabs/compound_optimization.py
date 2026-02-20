import gradio as gr

from rdkit import Chem
from rdkit.Chem import Descriptors, Draw


def _optimize_compound(smiles: str, goals: str):
    """
    Compound Optimization Tab Logic

    Inputs:
      - smiles: SMILES string for the compound
      - goals: free-text optimization goals (e.g., "Increase solubility, reduce lipophilicity")

    Outputs:
      - formatted text summary of properties + goal-based suggestions + Lipinski assessment
      - 2D molecule image rendered from SMILES
    """
    smiles = (smiles or "").strip()
    goals = (goals or "").strip()

    if not smiles:
        return "Please enter a valid SMILES string.", None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Invalid SMILES string. Please check your input.", None

    # --- Basic properties ---
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)

    # --- Lipinski rule-of-five violations ---
    violations = 0
    violations += 1 if mw > 500 else 0
    violations += 1 if logp > 5 else 0
    violations += 1 if hbd > 5 else 0
    violations += 1 if hba > 10 else 0

    # A common quick heuristic: <= 1 violation is typically acceptable in early discovery
    lipinski_status = "Pass" if violations <= 1 else "Fail"

    # --- Goal-based suggestions ---
    goal_text = goals.lower()
    suggestions = []
    goal_mapping = []

    if "solubility" in goal_text:
        suggestions.append(
            "Increase polarity by adding or substituting polar functional groups such as hydroxyl (-OH), amines (-NH2), or heteroatoms in rings."
        )
        suggestions.append(
            "Reduce crystallinity by introducing small branching or disrupting symmetry, if feasible for the scaffold."
        )
        goal_mapping.append("Solubility → increase polarity, reduce crystallinity, balance LogP")

    if "lipophil" in goal_text or "logp" in goal_text:
        suggestions.append(
            "Reduce lipophilicity by shortening alkyl chains, reducing aromatic ring count, or adding polar substituents while preserving key pharmacophore features."
        )
        suggestions.append(
            "Avoid adding large hydrophobic groups that increase LogP unless BBB penetration is a primary requirement."
        )
        goal_mapping.append("Lipophilicity → reduce hydrophobic surface area and aromatic burden")

    if "bbb" in goal_text or "blood brain" in goal_text or "brain" in goal_text:
        suggestions.append(
            "For BBB permeability, aim for moderate lipophilicity. A common early heuristic is LogP roughly in the 1–3 range while keeping polarity controlled."
        )
        suggestions.append(
            "Minimize strongly ionizable groups when possible and avoid excessive hydrogen-bonding that can reduce CNS exposure."
        )
        goal_mapping.append("BBB permeability → moderate LogP, controlled polarity, fewer strong ionizable groups")

    if "metabolism" in goal_text or "stability" in goal_text or "half-life" in goal_text:
        suggestions.append(
            "To improve metabolic stability, consider reducing easily oxidized motifs (e.g., exposed aromatic methyls) and consider steric shielding near labile positions."
        )
        suggestions.append(
            "Avoid reactive functional groups and consider bioisosteres for metabolically labile linkers."
        )
        goal_mapping.append("Metabolic stability → reduce labile motifs, consider bioisosteres")

    if "tox" in goal_text or "toxicity" in goal_text or "off-target" in goal_text:
        suggestions.append(
            "To reduce off-target risk, avoid high lipophilicity and consider removing structural alerts if present. Keep the molecule as simple as possible while maintaining activity."
        )
        goal_mapping.append("Safety/off-target → avoid excessive lipophilicity, reduce structural alerts")

    if not suggestions:
        suggestions.append(
            "Add specific optimization goals such as 'increase solubility', 'reduce lipophilicity', or 'improve BBB permeability' to receive tailored recommendations."
        )
        goal_mapping.append("No recognized goals → provide more specific goals for targeted suggestions")

    suggestions_text = "\n".join([f"- {s}" for s in suggestions])
    goal_mapping_text = "\n".join([f"- {m}" for m in goal_mapping])

    # --- Drug-likeness notes (simple, practical) ---
    drug_likeness_notes = [
        "Lipinski Rule of Five is a quick screen, not a guarantee of success.",
        "Optimization involves tradeoffs. Improving BBB permeability can reduce solubility.",
        "Use docking and ADME tools together to guide iteration rather than relying on one metric.",
    ]
    drug_likeness_text = "\n".join([f"- {n}" for n in drug_likeness_notes])

    # --- Build response ---
    response = f"""Initial Properties:
- Molecular Weight: {mw:.2f}
- LogP: {logp:.2f}
- H-Bond Donors: {hbd}
- H-Bond Acceptors: {hba}

Optimization Goals:
{goals if goals else "(none provided)"}

Suggested Structural Modifications:
{suggestions_text}

Optimization Goal Mapping:
{goal_mapping_text}

Drug-Likeness Considerations:
- Lipinski Violations: {violations}
- Overall Lipinski Status: {lipinski_status}
{drug_likeness_text}
"""

    # --- 2D image ---
    # Using RDKit's default depiction; this returns a PIL.Image.Image which Gradio accepts.
    img = Draw.MolToImage(mol, size=(350, 350))

    return response.strip(), img


def create_tab():
    with gr.Tab("Compound Optimization"):
        with gr.Row():
            with gr.Column(scale=1):
                optim_smiles_input = gr.Textbox(
                    label="Compound SMILES",
                    placeholder="e.g., CC(=O)Oc1ccccc1C(=O)O",
                    lines=3,
                )
                optim_goals_input = gr.Textbox(
                    label="Optimization Goals",
                    placeholder="e.g., Increase solubility, reduce lipophilicity",
                    lines=3,
                )
                optimize_btn = gr.Button(
                    "Optimize Compound",
                    size="lg",
                    elem_classes=["analyze-btn"],
                )

            with gr.Column(scale=1):
                optim_output = gr.Textbox(
                    label="Optimization Suggestions",
                    interactive=False,
                    lines=14,
                    placeholder=(
                        "a. Initial Properties of the Compound\n"
                        "b. Suggested Structural Modifications\n"
                        "c. Optimization Goal Mapping\n"
                        "d. Drug-Likeness Considerations"
                    ),
                )
                optim_image = gr.Image(
                    label="Molecule Visualization (2D)",
                    type="pil",
                    height=350,
                )

        optimize_btn.click(
            _optimize_compound,
            inputs=[optim_smiles_input, optim_goals_input],
            outputs=[optim_output, optim_image],
        )

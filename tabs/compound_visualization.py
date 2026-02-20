import gradio as gr
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import io


def _visualize_compound(smiles):
    if not smiles:
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        Chem.rdDepictor.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=(400, 400))

        return img

    except Exception:
        return None


def create_tab():
    with gr.Tab("Compound Visualization"):
        with gr.Row():
            with gr.Column(scale=1):
                viz_smiles_input = gr.Textbox(
                    label="Compound SMILES",
                    placeholder="e.g., CC(=O)OC1=CC=CC=C1C(=O)O",
                    lines=3,
                )
                visualize_btn = gr.Button(
                    "Visualize Compound",
                    size="lg",
                    elem_classes=["analyze-btn"],
                )

            with gr.Column(scale=1):
                viz_output = gr.Image(
                    label="2D Compound Representation",
                )

        visualize_btn.click(
            _visualize_compound,
            inputs=[viz_smiles_input],
            outputs=[viz_output],
        )
import gradio as gr


def _prepare_docking(pdb_id, smiles):
    """Placeholder – returns empty until logic is wired up."""
    return ""


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

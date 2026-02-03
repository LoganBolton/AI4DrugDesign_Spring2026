import gradio as gr


def _evaluate_compound(pdb_id, smiles):
    """Placeholder – returns empty until logic is wired up."""
    return ""


def create_tab():
    with gr.Tab("Compound Evaluation"):
        with gr.Row():
            with gr.Column(scale=1):
                eval_pdb_input = gr.Textbox(
                    label="Target Protein PDB ID",
                    placeholder="e.g., 6LU7",
                )
                eval_smiles_input = gr.Textbox(
                    label="Compound SMILES",
                    placeholder="e.g., CC(=O)OC1=CC=CC=C1C(=O)O",
                    lines=3,
                )
                evaluate_btn = gr.Button(
                    "Evaluate Compound",
                    size="lg",
                    elem_classes=["analyze-btn"],
                )

            with gr.Column(scale=1):
                eval_output = gr.Textbox(
                    label="Evaluation Results",
                    interactive=False,
                    lines=12,
                    placeholder=(
                        "a. Likelihood of Target Binding\n"
                        "b. Potential Activity Against the Target\n"
                        "c. Pharmacokinetic Considerations\n"
                        "d. Structural Improvements for Better Target Activity\n"
                        "e. Potential Off-Target Concerns\n"
                        "f. Overall Suitability for Further Development"
                    ),
                )

        evaluate_btn.click(
            _evaluate_compound,
            inputs=[eval_pdb_input, eval_smiles_input],
            outputs=[eval_output],
        )

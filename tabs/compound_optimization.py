import gradio as gr


def _optimize_compound(smiles, goals):
    """Placeholder – returns empty until logic is wired up."""
    return ""


def create_tab():
    with gr.Tab("Compound Optimization"):
        with gr.Row():
            with gr.Column(scale=1):
                optim_smiles_input = gr.Textbox(
                    label="Compound SMILES",
                    placeholder="e.g., CC(=O)OC1=CC=CC=C1C(=O)O",
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
                    lines=12,
                    placeholder=(
                        "a. Initial Properties of the Compound\n"
                        "b. Suggested Structural Modifications\n"
                        "c. Optimization Goal Mapping\n"
                        "d. Drug-Likeness Considerations"
                    ),
                )

        optimize_btn.click(
            _optimize_compound,
            inputs=[optim_smiles_input, optim_goals_input],
            outputs=[optim_output],
        )

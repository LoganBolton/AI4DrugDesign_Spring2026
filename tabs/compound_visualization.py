import gradio as gr


def _visualize_compound(smiles):
    """Placeholder – returns None until 2D rendering is wired up."""
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

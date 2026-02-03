import gradio as gr


def _analyze_protein(pdb_id):
    """Placeholder – returns empty until logic is wired up."""
    return ""


def create_tab():
    with gr.Tab("Protein Analysis"):
        with gr.Row():
            with gr.Column(scale=1):
                pdb_input = gr.Textbox(
                    label="Enter PDB ID",
                    placeholder="e.g., 1AZ5",
                    lines=3,
                )
                analyze_btn = gr.Button(
                    "Analyze Protein",
                    size="lg",
                    elem_classes=["analyze-btn"],
                )

            with gr.Column(scale=1):
                summary_output = gr.Textbox(
                    label="Protein Summary",
                    interactive=False,
                    lines=8,
                )

        analyze_btn.click(_analyze_protein, inputs=[pdb_input], outputs=[summary_output])

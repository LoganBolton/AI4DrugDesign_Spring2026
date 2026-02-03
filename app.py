import gradio as gr
import os
from openai import OpenAI
from dotenv import load_dotenv

load_dotenv()
client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

# Dark theme with orange active-tab accent
theme = gr.themes.Base(
    primary_hue=gr.themes.colors.orange,
    neutral_hue=gr.themes.colors.gray,
).set(
    body_background_fill="rgb(30, 30, 30)",
    background_fill_secondary="rgb(40, 40, 42)",
    block_background_fill="rgb(42, 42, 44)",
    block_border_color="rgb(60, 60, 62)",
    block_label_text_color="rgb(200, 200, 200)",
    input_background_fill="rgb(50, 50, 52)",
    input_border_color="rgb(80, 80, 82)",
    body_text_color="rgb(230, 230, 230)",
    button_primary_background_fill="rgb(90, 90, 92)",
    button_primary_background_fill_hover="rgb(110, 110, 112)",
    button_primary_text_color="rgb(255, 255, 255)",
    button_primary_border_color="rgb(90, 90, 92)",
)

css = """
/* title bar */
.platform-title {
    background-color: rgb(35, 35, 37);
    padding: 12px 20px;
    border-radius: 8px 8px 0 0;
    margin-bottom: 0 !important;
}
.platform-title p {
    color: #fff;
    font-size: 1.05rem;
    font-weight: 600;
    margin: 0 !important;
}

/* active tab → orange text + underline */
.tabs .tab-nav button.selected {
    color: rgb(232, 120, 42) !important;
    border-bottom-color: rgb(232, 120, 42) !important;
}

/* "Analyze Protein" button full-width */
.analyze-btn button {
    width: 100%;
}
"""


def analyze_protein(pdb_id):
    """Placeholder – returns empty until logic is wired up."""
    return ""


with gr.Blocks(title="AI Drug Design Platform") as demo:

    # Title bar
    gr.Markdown("AI Drug Design Platform", elem_classes=["platform-title"])

    with gr.Tabs():
        # ── Protein Analysis ────────────────────────────────────
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

            analyze_btn.click(analyze_protein, inputs=[pdb_input], outputs=[summary_output])

        # ── Placeholder tabs ────────────────────────────────────
        with gr.Tab("Compound Evaluation"):
            gr.Markdown("*Coming soon*")

        with gr.Tab("Compound Optimization"):
            gr.Markdown("*Coming soon*")

        with gr.Tab("Compound Visualization"):
            gr.Markdown("*Coming soon*")

        with gr.Tab("Docking Preparation"):
            gr.Markdown("*Coming soon*")

        with gr.Tab("Molecular Docking"):
            gr.Markdown("*Coming soon*")

demo.launch(theme=theme, css=css)

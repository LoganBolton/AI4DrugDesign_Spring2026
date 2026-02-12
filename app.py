import gradio as gr
from dotenv import load_dotenv
from tabs.protein_image_tab import build_tab as protein_image_tab

load_dotenv()  # must run before tab imports so OpenAI keys are available

from tabs.protein_analysis import create_tab as protein_analysis_tab
from tabs.compound_evaluation import create_tab as compound_evaluation_tab
from tabs.compound_optimization import create_tab as compound_optimization_tab
from tabs.compound_visualization import create_tab as compound_visualization_tab
from tabs.docking_preparation import create_tab as docking_preparation_tab
from tabs.molecular_docking import create_tab as molecular_docking_tab
from tabs.chat import create_tab as chat_tab

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

/* full-width action buttons */
.analyze-btn button {
    width: 100%;
}
"""

with gr.Blocks(title="AI Drug Design Platform") as demo:
    gr.Markdown("AI Drug Design Platform", elem_classes=["platform-title"])

    with gr.Tabs():
        protein_analysis_tab()
        compound_evaluation_tab()
        compound_optimization_tab()
        compound_visualization_tab()
        docking_preparation_tab()
        molecular_docking_tab()
        chat_tab()
        with gr.Tab("Structure Viewer"):
            protein_image_tab()

demo.launch(theme=theme, css=css)

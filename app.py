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


def evaluate_compound(pdb_id, smiles):
    """Placeholder – returns empty until logic is wired up."""
    return ""


def optimize_compound(smiles, goals):
    """Placeholder – returns empty until logic is wired up."""
    return ""


def prepare_docking(pdb_id, smiles):
    """Placeholder – returns empty until logic is wired up."""
    return ""


def auto_calculate_center(pdb_id, smiles):
    """Placeholder – returns default coordinates until logic is wired up."""
    return 0.0, 0.0, 0.0


def run_docking(pdb_id, smiles, center_x, center_y, center_z):
    """Placeholder – returns empty until logic is wired up."""
    return "", None


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
        # ── Compound Evaluation ─────────────────────────────────
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
                evaluate_compound,
                inputs=[eval_pdb_input, eval_smiles_input],
                outputs=[eval_output],
            )

        # ── Compound Optimization ───────────────────────────────
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
                optimize_compound,
                inputs=[optim_smiles_input, optim_goals_input],
                outputs=[optim_output],
            )

        with gr.Tab("Compound Visualization"):
            gr.Markdown("*Coming soon*")

        # ── Docking Preparation ─────────────────────────────────
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
                prepare_docking,
                inputs=[dock_pdb_input, dock_smiles_input],
                outputs=[dock_output],
            )

        # ── Molecular Docking ────────────────────────────────────
        with gr.Tab("Molecular Docking"):
            with gr.Row():
                with gr.Column(scale=1):
                    vina_pdb_input = gr.Textbox(
                        label="Target PDB ID",
                        placeholder="e.g., 6LU7",
                    )
                    vina_smiles_input = gr.Textbox(
                        label="Ligand SMILES",
                        placeholder="e.g., CC(=O)OC1=CC=CC=C1C(=O)O",
                        lines=3,
                    )
                    auto_calc_btn = gr.Button(
                        "Auto-Calculate Center from Co-Ligand",
                        size="sm",
                        elem_classes=["analyze-btn"],
                    )
                    with gr.Row():
                        center_x = gr.Number(label="Center X", value=0.0, precision=3)
                        center_y = gr.Number(label="Center Y", value=0.0, precision=3)
                        center_z = gr.Number(label="Center Z", value=0.0, precision=3)
                    run_vina_btn = gr.Button(
                        "Run AutoDock Vina",
                        size="lg",
                        elem_classes=["analyze-btn"],
                    )

                with gr.Column(scale=1):
                    vina_affinities = gr.Textbox(
                        label="Vina Binding Affinities",
                        interactive=False,
                        lines=8,
                    )
                    vina_output_file = gr.File(
                        label="Docked PDBQT File (Download)",
                    )

            auto_calc_btn.click(
                auto_calculate_center,
                inputs=[vina_pdb_input, vina_smiles_input],
                outputs=[center_x, center_y, center_z],
            )
            run_vina_btn.click(
                run_docking,
                inputs=[vina_pdb_input, vina_smiles_input, center_x, center_y, center_z],
                outputs=[vina_affinities, vina_output_file],
            )

demo.launch(theme=theme, css=css)

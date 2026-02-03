import gradio as gr


def _auto_calculate_center(pdb_id, smiles):
    """Placeholder – returns default coordinates until logic is wired up."""
    return 0.0, 0.0, 0.0


def _run_docking(pdb_id, smiles, center_x, center_y, center_z):
    """Placeholder – returns empty until logic is wired up."""
    return "", None


def create_tab():
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
            _auto_calculate_center,
            inputs=[vina_pdb_input, vina_smiles_input],
            outputs=[center_x, center_y, center_z],
        )
        run_vina_btn.click(
            _run_docking,
            inputs=[vina_pdb_input, vina_smiles_input, center_x, center_y, center_z],
            outputs=[vina_affinities, vina_output_file],
        )

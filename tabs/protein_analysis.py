import gradio as gr
import requests
import os
from openai import OpenAI


def _fetch_pdb_info(pdb_id: str) -> dict:
    """Fetch structured info from RCSB PDB REST API."""
    pdb_id = pdb_id.strip().upper()
    base = "https://data.rcsb.org/rest/v1/core"

    # Core entry
    r = requests.get(f"{base}/entry/{pdb_id}", timeout=15)
    r.raise_for_status()
    entry = r.json()

    struct = entry.get("struct", {})
    title = struct.get("title", "N/A")
    keywords = entry.get("struct_keywords", {}).get("pdbx_keywords", "N/A")
    text_keywords = entry.get("struct_keywords", {}).get("text", "N/A")

    # Organism / source
    source_organisms = []
    try:
        poly_r = requests.get(f"{base}/entry/{pdb_id}/polymer_entities", timeout=15)
        if poly_r.ok:
            for entity in poly_r.json():
                src = entity.get("rcsb_entity_source_organism", [])
                for s in src:
                    org = s.get("ncbi_scientific_name")
                    if org and org not in source_organisms:
                        source_organisms.append(org)
    except Exception:
        pass
    organism = ", ".join(source_organisms) if source_organisms else "N/A"

    # Chains
    chains = []
    try:
        chains_r = requests.get(f"{base}/entry/{pdb_id}/polymer_entity_instances", timeout=15)
        if chains_r.ok:
            for inst in chains_r.json():
                chain_id = inst.get("rcsb_polymer_entity_instance_container_identifiers", {}).get("auth_asym_id")
                if chain_id and chain_id not in chains:
                    chains.append(chain_id)
    except Exception:
        pass

    # Biological assemblies count
    bio_assemblies = entry.get("rcsb_entry_info", {}).get("assembly_count", "N/A")

    # Resolution
    resolution = entry.get("refine", [{}])[0].get("ls_d_res_high", None) if entry.get("refine") else None

    # Non-polymer / ligand entity IDs
    non_poly_ids = entry.get("rcsb_entry_container_identifiers", {}).get("non_polymer_entity_ids") or []
    skip = {"HOH", "SO4", "PO4", "GOL", "EDO", "ACT", "FMT", "IOD", "CL", "MG", "ZN", "CA", "NA", "PEG", "DMS"}
    ligands = []
    for eid in non_poly_ids[:10]:
        try:
            er = requests.get(f"{base}/nonpolymer_entity/{pdb_id}/{eid}", timeout=10)
            if er.ok:
                comp_id = er.json().get("pdbx_entity_nonpoly", {}).get("comp_id", "")
                name = er.json().get("pdbx_entity_nonpoly", {}).get("name", comp_id)
                if comp_id and comp_id not in skip:
                    ligands.append(f"{comp_id} ({name})")
        except Exception:
            continue

    # Binding sites
    binding_sites = []
    try:
        site_r = requests.get(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/{pdb_id.lower()}", timeout=15)
        if site_r.ok:
            site_data = site_r.json().get(pdb_id.lower(), {})
            for site_id, site_info in site_data.items():
                details = site_info[0] if isinstance(site_info, list) and site_info else site_info
                site_desc = details.get("site_description", site_id)
                ligand_id = details.get("ligand_entity_id", "")
                residues = details.get("site_residues", [])
                res_list = ", ".join(
                    f"{res.get('residue_name','')}{res.get('residue_number','')}{res.get('chain_id','')}"
                    for res in residues[:6]
                )
                binding_sites.append({
                    "id": site_id,
                    "description": site_desc,
                    "ligand": ligand_id,
                    "residues": res_list or "N/A",
                })
    except Exception:
        pass

    # Missing residues count
    missing_residues = "N/A"
    try:
        val_r = requests.get(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/{pdb_id.lower()}", timeout=15)
        if val_r.ok:
            molecules = val_r.json().get(pdb_id.lower(), {}).get("molecules", [])
            missing_count = 0
            for mol in molecules:
                for chain in mol.get("chains", []):
                    for res in chain.get("residues", []):
                        if not res.get("observed", True):
                            missing_count += 1
            missing_residues = str(missing_count) if missing_count else "None reported"
    except Exception:
        pass

    # Classification from PDBe summary
    classification = keywords
    try:
        pdbe_r = requests.get(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_id.lower()}", timeout=15)
        if pdbe_r.ok:
            summ = pdbe_r.json().get(pdb_id.lower(), [{}])[0]
            classification = summ.get("compound", [{}])[0].get("molecule", keywords)
    except Exception:
        pass

    return {
        "pdb_id": pdb_id,
        "title": title,
        "classification": classification,
        "keywords": text_keywords,
        "organism": organism,
        "chains": ", ".join(chains) if chains else "N/A",
        "bio_assemblies": str(bio_assemblies),
        "resolution": f"{resolution} Å" if resolution else "N/A",
        "ligands": ligands,
        "binding_sites": binding_sites,
        "missing_residues": missing_residues,
    }


def _build_info_text(info: dict) -> str:
    lines = [
        "=" * 55,
        f"  PROTEIN INFORMATION — {info['pdb_id']}",
        "=" * 55,
        f"Protein Name     : {info['classification']}",
        f"Title            : {info['title']}",
        f"Classification   : {info['keywords']}",
        f"Organism         : {info['organism']}",
        f"Chains           : {info['chains']}",
        f"Resolution       : {info['resolution']}",
        f"Bio. Assemblies  : {info['bio_assemblies']}",
        f"Missing Residues : {info['missing_residues']}",
        "",
        "LIGANDS / CO-FACTORS:",
    ]
    if info["ligands"]:
        for lig in info["ligands"]:
            lines.append(f"  • {lig}")
    else:
        lines.append("  None identified")

    lines += ["", "BINDING SITES:"]
    if info["binding_sites"]:
        for site in info["binding_sites"]:
            lines.append(f"  Site {site['id']}: {site['description']}")
            lines.append(f"    Ligand   : {site['ligand'] or 'N/A'}")
            lines.append(f"    Residues : {site['residues']}")
    else:
        lines.append("  No binding site data available")

    lines.append("=" * 55)
    return "\n".join(lines)


def _ai_analysis(info: dict) -> str:
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        return "[OpenAI API key not set – skipping detailed analysis]"

    client = OpenAI(api_key=api_key)

    ligand_str = ", ".join(info["ligands"]) if info["ligands"] else "none identified"
    site_str = "; ".join(
        f"Site {s['id']} ({s['description']}) residues: {s['residues']}"
        for s in info["binding_sites"]
    ) if info["binding_sites"] else "none available"

    prompt = f"""You are an expert structural biologist and medicinal chemist.
Analyze the following PDB entry and provide a detailed, structured analysis.

PDB ID: {info['pdb_id']}
Title: {info['title']}
Classification: {info['classification']}
Organism: {info['organism']}
Chains: {info['chains']}
Resolution: {info['resolution']}
Ligands/Co-factors: {ligand_str}
Binding Sites: {site_str}
Missing Residues: {info['missing_residues']}

Provide a comprehensive analysis with the following clearly labeled sections:

1. BINDING SITES
   - Describe the binding site(s), their location, and functional significance.

2. KEY RESIDUES
   - Identify catalytic, allosteric, or drug-interacting residues and their roles.

3. INTERACTION TYPES
   - Describe likely protein-ligand interactions (H-bonds, hydrophobic, electrostatic, π-stacking, metal coordination).

4. DESIGN CONSIDERATIONS
   - Summarize structural features important for drug/inhibitor design targeting this protein.

5. LIGAND FEATURES
   - Analyze known ligands or co-factors: size, functional groups, binding mode, pharmacophore features.

6. DRUG DEVELOPMENT INSIGHTS
   - Discuss druggability, selectivity challenges, known inhibitor classes, and therapeutic relevance.

7. ANALYSIS RECOMMENDATIONS
   - Suggest follow-up experiments, computational studies, or structural optimization strategies.

Be specific, concise, and scientifically accurate. Use the structural data provided."""
    response = client.chat.completions.create(
        model="gpt-5-nano-2025-08-07",
        messages=[{"role": "user", "content": prompt}],
    )
    return response.choices[0].message.content


def _analyze_protein(pdb_id: str) -> str:
    pdb_id = pdb_id.strip().upper()
    if not pdb_id:
        return "Please enter a PDB ID."

    try:
        info = _fetch_pdb_info(pdb_id)
    except requests.HTTPError as e:
        if e.response is not None and e.response.status_code == 404:
            return f"PDB ID '{pdb_id}' not found. Please check the ID and try again."
        return f"Error fetching PDB data: {e}\nPlease verify the PDB ID is correct."
    except requests.exceptions.ConnectionError:
        return "Network error: could not reach the PDB server. Please check your internet connection."
    except requests.exceptions.Timeout:
        return "Request timed out. The PDB server may be slow — please try again."
    except Exception as e:
        return f"Unexpected error: {e}"

    info_text = _build_info_text(info)

    try:
        analysis_text = _ai_analysis(info)
    except Exception as e:
        analysis_text = f"[AI analysis failed: {e}]"

    return info_text + "\n\n" + "DETAILED ANALYSIS\n" + "=" * 55 + "\n" + analysis_text


def create_tab():
    with gr.Tab("Protein Analysis"):
        with gr.Row():
            with gr.Column(scale=1):
                pdb_input = gr.Textbox(
                    label="Enter PDB ID",
                    placeholder="e.g., 1AZ5",
                    lines=1,
                )
                analyze_btn = gr.Button(
                    "Analyze Protein",
                    size="lg",
                    elem_classes=["analyze-btn"],
                )

            with gr.Column(scale=2):
                summary_output = gr.Textbox(
                    label="Protein Summary & Detailed Analysis",
                    interactive=False,
                    lines=30,
                    placeholder=(
                        "Results will appear here after clicking 'Analyze Protein'.\n\n"
                        "Includes:\n"
                        "  • Protein Name, Title, Classification\n"
                        "  • Organism, Chains, Biological Assembly\n"
                        "  • Binding Sites, Missing Residues\n\n"
                        "Detailed Analysis:\n"
                        "  • Binding Sites\n"
                        "  • Key Residues\n"
                        "  • Interaction Types\n"
                        "  • Design Considerations\n"
                        "  • Ligand Features\n"
                        "  • Drug Development Insights\n"
                        "  • Analysis Recommendations"
                    ),
                )

        gr.Examples(examples=[["1AZ5"], ["6LU7"], ["1IEP"], ["2HYY"]], inputs=pdb_input)

        analyze_btn.click(_analyze_protein, inputs=[pdb_input], outputs=[summary_output])
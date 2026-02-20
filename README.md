# AI4DrugDesign_Spring2026

A simple Gradio-based application for AI4DrugDesign course.

## Installation

1. **Install uv** (if not already installed):

   **macOS/Linux:**
   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

   **Windows (PowerShell):**
   ```powershell
   powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
   ```

   After installation, restart your terminal or run:
   - macOS/Linux: `source $HOME/.local/bin/env` or `source $HOME/.cargo/env`
   - Windows: Restart PowerShell

   Verify installation:
   ```bash
   uv --version
   ```


2. **Install dependencies**:
   ```bash
   uv sync
   ```
   This will:
   - Create a virtual environment in `.venv/`
   - Install Python 3.12 (specified in `.python-version`)
   - Install all dependencies from `pyproject.toml`

## Running the App

```bash
uv run python app.py
```

To activate the virtual environment manually:
- **macOS/Linux:** `source .venv/bin/activate`
- **Windows (PowerShell):** `.venv\Scripts\Activate.ps1`
- **Windows (CMD):** `.venv\Scripts\activate.bat`



## PREVIOUS RELEASE SPECS
The following is a previous version of the app and the requirements it had:

### User Docs
AI Drug Design Platform 
User Manual 
Auburn University Senior Design Project (COMP 4710) Fall 2025
Table of Contents 
1. Introduction 
2. System Requirements 
3. Getting Started 
4. Platform Features 
a. Protein Analysis 
b. Compound Evaluation 
c. Compound Optimization 
d. Compound Visualization 
e. Docking Preparation 
f. Compound Docking/Screening Export 5. Example Workflows 
6. Interpreting Results 
7. Troubleshooting 
8. Further Resources
Introduction 
The AI For Drug Design Platform is an encompassing tool that integrates computational chemistry with the analytical capabilities of artificial intelligence to aid researchers and users in the drug design, drug analysis, and drug discovery process. This platform, which was developed in Auburn University’s Computer Science and Software Engineering Department in the Samuel Ginn college of Engineering, utilizes Google’s Gemini large language models, RDKit’s molecular modeling capabilities, AutoDock Vina’s Docking component, and an assortment of other tools to provide a useful system for drug creation related tasks. 
Key Features 
➔ Protein target analysis and binding site identification 
➔ Evaluation of the compounds based on target binding likelihood, suitability, and other considerations 
➔ Optimization that suggests structural improvements, goal mapping, and drug-likeness consideration for the compound 
➔ Compound Screening which generates the molecular structure for users to see a 2D representation and also the capability to export a 3D file for an enterprise level screening software 
➔ Docking that allows for the generation of more specific data regarding binding affinity among other things 
Team Members 
➔ Jaelin Robinson 
➔ Drew Black 
➔ Ryan Thomas 
Faculty Oversight 
Dr. Gerry Dozier, Professor, Department of Computer Science and Software Engineering, Auburn University 
Dr. Jakita Thomas, Professor, Department of Computer Science and Software Engineering, Auburn University
Project Sponsor 
Dr. Rajesh Amin, Department of Drug Discovery and Development, Auburn University
System Requirements 
To run the AI Drug Design platform, you will need: ➔ Python 3.10+ 
➔ Internet Connection for API calls 
➔ Required Packages: 
◆ Genai from google 
◆ Io 
◆ Bio.PDB 
◆ Rdkit 
◆ Re 
◆ Os 
◆ Subprocess 
◆ Numpy 
◆ Gradio 
◆ dotenv
Getting Started 
Installation 
1. Install Python 3.10 or higher if not already installed 
2. Clone or download the repository 
3. Install the required dependencies mentioned above. Here is an example: 4. Launch the app: 
5. Copy the url to a browser
Platform Features 
Protein Analysis 
The Protein Analysis tab allows you to explore a variety of protein targets and evaluate based on those targets. 
Analyzing a Protein: 
1. Navigate to the “Protein Analysis” tab 
2. Enter the PDB ID in the “Enter PDB ID” field (e.g. 1AZ5). 
3. Click the “Analyze Protein” Button 
4. Review the analysis in the right panel including: 
a. Information on the protein(Protein Name, Title, Classification, Organism, Chains, Keywords, Biological Assembly, Binding Sites, and Missing Residues) 
b. Detailed Analysis containing 
i. Binding Sites 
ii. Key Residues 
iii. Interaction Types 
iv. Design Considerations 
v. Ligand Features 
vi. Drug Development Insights 
vii. Analysis Recommendations 
Compound Evaluation
This tab allows you to evaluate the compound based on both the target protein and the SMILES compound 
Evaluating a Compound 
1. Navigate to the “Compound Evaluation” tab 
2. Enter a valid PDB ID in the “Target Protein PDB ID” field(e.g., 6LU7) 3. Enter the SMILES Compound in the “Compound SMILES” field (e.g,CC(=O)OC1=CC=CC=C1C)=O)O 
4. Click the “Evaluate Compound” button 
5. Review the Evaluation results, which include: 
a. Likelihood of Target Binding 
b. Potential Activity Against the Target 
c. Pharmacokinetic Considerations 
d. Structural Improvements for Better Target Activity 
e. Potential Off-Target Concerns 
f. Overall Suitability for Further Development 
Compound Optimization
This tab allows you to view possible optimizations that can be applied to a compound 
Optimizing A Compound 
1. Navigate to the “Compound Optimization” Tab 
2. Enter a valid SMILES Compound in the “Compound SMILES” field(e.g., CC(=O)OC1=CC=CC=C1C(=O)O 
3. Enter valid optimization goals in the “Optimization Goals” field (e.g., Increase solubility, reduce lipophilocity) 
4. Click the “Optimize Compound: button 
5. Review the optimization suggestions, which include: 
a. The initial properties of the compound 
b. Suggested Structural Modifications 
c. Optimization Goal Mapping 
d. Drug-Likeness Consideration 
Compound Visualization

This tab allows you to use a 2D visualization of the smiles compound 
Visualizing a Compound 
1. Navigate to the “Compound Visualization: tab 
2. Enter the SMILES compound in the “Compound SMILES” field 3. The platform will then generate a 2D Structure Visualization Docking Preparation 
This tab allows the user to prepare files for docking with AutoDock Vina. Preparing to Dock 
1. Navigate to the “Docking Preparation: tab 
2. Enter the PDB ID in the “Target Protein PDB ID” tab 
3. Enter the SMILES Compound in the “Compound SMILES” tab 4. Click the “Prepare Files for Vina” button 
5. Review the results of the preparation, which will include: 
a. A confirmation that Preparation was successful 
b. Protein PDB saved and Ligand PDB saved 
c. Steps on proceeding with AutoDock Vina 
Molecular Docking
This tab will allow users to perform molecular docking and locate exactly where a compound is in the the 3D grid for the screening software, and then output a file that can be exported to a enterprise level screening software 
Molecular Docking 
1. Navigate to the “Molecular Docking” tab 
2. Enter the PDB ID in “Target PDB ID” field 
3. Enter the SMILES Ligand in the “Ligand SMILES” field 
4. If user is not waare of Binding Box Coordinates they can click the “Auto-Calculate Center from Co-Ligand” to attempt to auto locate the Ligand 
5. Otherwise they can input the coordinates manually in the tabs “Center X”, “Center Y”, and “Center Z” 
6. Click the “Run AutoDock Vina” button 
7. This will output 
a. Vina Binding Affinities 
b. A Downloadable Docked PDBQT file
Example Workflows 
1. Protein Analysis 
1. Navigate to the “Protein Analysis” tab 
2. Enter “1AZ5” in the target name field 
3. Click “Analyze Protein” 
4. Review the analysis of 1AZ5 as a drug target 
2. Compound Evaluation 
1. Navigate to the “Compound Evaluation” tab 
2. Enter “1AZ5” in the target name field 
3. Enter “CC(=O)Oc1ccccc1C(=O)O” in the SMILES field 
4. Click “Evaluate Compound” 
5. Review the details of the evaluation 
3. Compound Optimization 
1. Navigate to the “Compound Optimization” tab 
2. Enter “CC(=O)Oc1ccccc1C(=O)O” in the SMILES field 
3. Enter “increases solubility, reduce lipophilocity” in the optimization goals field 
4. Click “Optimize Compound” 
5. Review the optimation suggestions 
4. Compound Visualization 
1. Navigate to the “Compound Visualization” tab 
2. Enter “CC(=O)Oc1ccccc1C(=O)O” in the SMILES field 
3. Click “Visualize Compound” 
4. View the 2D representation of the compound 
5. Preparation for Docking 
1. Navigate to the “Docking Preparation” tab 
2. Enter “1AZ5” in the target name field 
3. Enter “CC(=O)Oc1ccccc1C(=O)O” in the SMILES field 
4. Click “Prepare Files for Vina” 
5. Verify that files have been prepared
6. Molecular Docking 
1. Navigate to the “Molecular Docking” tab 
2. Enter “1AZ5” in the target name field 
3. Enter “CC(=O)Oc1ccccc1C(=O)O” in the SMILES filed 4. Click “Auto-Calculate Center from Co-Ligand” 5. Verify that coordinates have been found 
6. Click “Run AutoDock Vina” 
7. View the results in “Vina Binding Affinities” 8. Download the Docked PDBQT file
Interpreting Results 
Protein Analysis Results 
- Binding Sites: The primary binding sites for the protein 
- Key Residues: The key residues that appear at the binding sites - Interaction Types: The types of interactions that can occur at the binding sites 
- Design Considerations: considerations or details that could affect the design - Ligand Features: features that the ideal ligand should exhibit - Drug Development Insights: Insights regarding the protein’s use in drug development 
- Analysis Recommendations: Recommendations taken from the analysis Compound Evaluation Results 
- Likelihood of Target Binding: probabilities of successful binding - Potential Activity Against the Target: Activity that could inhibit the target - Pharmacokinetic Considerations: considerations ranging from molecular weight to druglikeness Score 
- Structural Improvements for Better Target Activity: Improvements that can be made that allow for better target activity 
- Potential Off-Target Concerns: Concerns away from the target that should be considered 
- Overall Suitability for Further Development: Likelihood to continue to the next phase of development 
Compound Optimization Results 
- Initial Properties: Initial properties of the smiles compound - Suggested Structural Modifications: suggested modifications that are given that will align with the optimization goals 
- Optimization Goal Mapping: How each modification maps to optimization goals 
- Drug-Likeness Considerations: Consideration that will prioritize the maintenance or improvement of drug-likeness in regards to the optimization goals 
Compound Visualization Results
- Visual Molecular Structure for User Visualization 
Docking Preparation Results 
- Indication of status of preparation 
- Indication of Protein and Ligand PDB save status 
- Instructions on how to proceed with AutoDockVina 
Molecular Docking Results 
- Binding Box coordinates for the compound 
- Vina Binding affinities that include: 
- Table showing mode, affinity, distance from best mode 
- Downloadable file for export to a screening Software 
Protein Analysis 
Based on the PDB ID, this section provides key insights into the proteins characteristics 
● Protein Name: The common or descriptive name of the protein structure associated with the PDB ID 
● Title: The formal title of the PDB entry 
● Classification: A short label describing the general type of molecule or function 
● Organism: The species the protein was derived from 
● Keywords: Descriptive tags that categorize the structure, such as its function, method, or biological role 
● Biological Assembly: The functional form of the protein in nature ● Binding Sites: Locations on the protein where ligands, substrates, inhibitors, or cofactors bind 
● Missing Residues: A list of amino acids not resolved in the structure Compound Evaluation
When assessing a compound against a chosen target, the system provides insights into: 
● Likelihood of Target Binding: Predicts how strongly the compound is expected to interact with the selected protein target 
● Potential Activity: Estimates the compound’s functional effect on the target, such as inhibition or activation 
● Pharmacokinetic Considerations: Provides insights into absorption, distribution, metabolism, and excretion characteristics that may affect drug viability 
● Structural Improvements: Suggests possible modifications to enhance binding affinity, stability, or overall performance 
● Potential Off-target Concerns: Identifies any predicted interactions with unintended proteins that may lead to side effects 
● Overall Suitability: Summarizes the compound’s overall potential as a viable therapeutic candidate based on combined evaluation factors 
Compound Optimization 
The system evaluates the input compound and offers guidance for structural or functional improvements 
● Suggested Structural Modifications: Proposes changes to the compound’s structure to improve target 
● Optimization Goal Mapping:Links each suggested modification to specific optimization objectives 
● Drug-Likeness Consideration:Assesses whether the optimized compound meets key drug-likeness criteria 
Compound Visualization 
From input SMILES shows 2D structure 
Vina Docking 
Outputs Vina binding affinities and downloadable docked PDBQT file for complex molecular docking
Troubleshooting 
Common Issues and Solutions 
Issue: Invalid or incorrect PDB ID 
Solution: The system cannot retrieve the protein structure if the PDB ID is misspelled or does not exist. Double-check spelling 
Issue: Invalid Smiles string 
Solution: The chemical notation (SMILES) entered for the compound is incorrect, malformed, or contains unsupported characters. Make sure input follows formatting rules. 
Issue: API-related errors prevent the system from retrieving or processing data. 
Solution: Check your internet connection, verify your API key is correct and active, and ensure your account has sufficient credits or usage limits. 
Issue: Required dependencies not installed or missing. 
Solution: Install all required software libraries and packages, then restart the program to ensure dependencies are loaded. Python must be 3.10 or higher 
Issue: Gradio interface does not launch or display. 
Solution: Ensure Gradio is correctly copied, check that all dependencies are met, and restart the program or your browser.
Further Resources 
Material and Documentation 
● Protein Data Bank: https://www.rcsb.org/ 
● RDKit Documentation: https://www.rdkit.org/docs/index.html ● AutoDock Vina Documentation: 
https://autodock-vina.readthedocs.io/en/latest/ 
● National Library of Medicine: https://pubchem.ncbi.nlm.nih.gov/ ● Python Dependencies: https://pypi.org/ 
License Information 
This project is licensed under the Apache License 2.0. The Apache License was chosen to provide patent protection in addition to copyright, which is particularly relevant for drug discovery applications. 
Copyright © 2025 Jaelin Robinson, Drew Black, Ryan Thomas

### Developer Docs

AI Drug Design Platform 
Developer Manual 
Auburn University Senior Design Project (COMP 4710) Fall 2025
Table of Contents 
1. Introduction 
2. System Architecture 
3. Development Environment Setup 4. Application Structure 
5. Core Components 
6. API Integration 
7. Data Processing Pipeline 
8. Extending the Platform 
9. Testing Framework 
10.Deployment Guide 
11.Troubleshooting for Developers 12.Contributing Guidelines
1. Introduction 
The AI Drug Design platform is an encompassing tool that integrates computational chemistry with the analytical capabilities of artificial intelligence to aid researchers and users in the drug design, drug analysis, and drug discovery process. 
Project Overview 
The AI Drug Design platform utilizes Google’s Gemini large language models, RDKit’s molecular modeling capabilities, AutoDock Vina’s Docking component, and an assortment of other tools to provide a useful system for drug creation related tasks. The application is built using python and Gradio-based graphical user interface(GUI). 
Technical Stack 
➔ Backend: 3.10+ 
➔ GUI: Gradio 
➔ Molecular Modeling: RDKit 
➔ AI Integration: Google genai 
➔ Data Processing: NumPy 
➔ Visualization: Rdkit.Chem 
➔ Docking: AutoDock VIna 
Team Members 
➔ Jaelin Robinson 
➔ Drew Black 
➔ Ryan Thomas 
Project Oversight 
➔ Faculty Advisors: 
◆ Dr. Gerry Dozier, Department of Computer Science and Software Engineering, Auburn University 
◆ Dr. Jakita Thomas, Department of Computer Science and Software Engineering, Auburn University
➔ Project Sponsor: Dr. Rajesh Amin, Department of Drug Discovery and Development, Auburn University
2. System Architecture 
The AI Drug Design platform follows a modular architecture that integrates one into the other. 
High Level Architecture 
The System contains these layers: 
1. GUI Layer: 
a. Simple GUI based on Gradio library 
b. Separates core functionalities into different tabs: 
i. Protein Analysis 
ii. Compound evaluation 
iii. Compound optimization 
iv. Compound visualization 
v. Docking/Screening exports 
c. Hosted on local and public servers 
2. Core Processing Layer: 
a. Consistent API integration using Google Gemini 
b. PDB database connectivity and data retrieval 
c. Compound property calculation
3. Development Environment Setup 
Prerequisites: 
➔ Python 3.10+ 
➔ Git 
➔ Pip 
Local Development Setup 
1. Clone the repository: 
2. Create and activate a virtual environment: 
3. Install dependencies: 
4. Set up environment variables: 
a. Create a ‘.env’ file in the repository root 
b. Add Your GEMINI API key: ‘GEMINI_API_KEY=your_api_key’ 5. Run the application locally:

4. Application Structure 
The application follows the below file structure: 

Key Files 
➔ Core.py: The core of the AI For Drug Design System that holds most of the computational logic and implementation 
➔ Wrapper.py: The wrapper for the core system that encompassed the UI that integrates these functionalities from Core.py
5. Core Components 
Gradio Interface 
The GUI is built on the Gradio library and provides the following tabs: 
● Protein analysis 
● Compound Evaluation 
● Compound Optimization 
● Vina Docking 
Each tab contains specific components for input, processing, and output display. 
Protein analysis 
Provides insights into potential drug discovery opportunities for a given target protein. 
Key Functions: 
● Identifies main ligand-binding sites. 
● Finds key residues for ligand interaction 
● Describes possible chemical interactions with a given compound 
Compound Evaluation 
Evaluates a compound’s compatibility with a target protein 
Key Functions: 
● Provides a summary of a compounds chemical characteristics ● Estimates likelihood of target binding 
● Highlights potential activity against the target 
● Provides pharmacokinetic considerations 
● Suggests structural improvements for better target activity
Compound Optimization 
Compound optimization proposes targeted chemical modifications to improve drug properties such as binding affinity, solubility, stability, safety, and overall drug-likeness. 
Key functions: 
● Estimates expected property changes following compound modification ● Provides modified SMILES structure 
● Explains rationale behind suggested modifications 
● Optimizes compound to user-input goals 
Vina Docking 
Simulates protein-compound interaction at given binding
6. API Integration 
Gemini API Integration 
The application uses Google’s Gemini API to generate AI-powered analysis and optimization of proteins and chemical compounds. 
Key aspects: 
● Prompt construction for specialized chemistry tasks 
● Real-time model response streaming 
● Support for advanced response formatting using markdown 
PDB Integration 
The application connects to and retrieves data from the Protein Data Bank. Key aspects: 
● PDB ID validation 
● Data retrieval and parsing
7. Data Processing Pipeline 
Molecular Data Processing 
The platform processes molecular data via the following pipeline: 
1. Input Parsing: User-input SMILES string is parsed and validated 2. Structure Generation: Molecular structure data is generated using RDKit 3. Property Calculation: Calculates basic and ADME compound properties 
Protein Data Processing 
Protein data follows this processing pipeline 
1. PDB Data Retrieval: Data is retrieved from PDB 
2. Structure Analysis: Binding pockets and key residues are identified 3. AI-Asisted Analysis: AI models provide insights about the protein 4. Results Compilation: Analysis results are compiled for presentation
8. Extending the Platform 
Adding New Features 
To add new features to the platform; 
1. Identify the appropriate module for the feature 
2. Implement core functionality in said module 
3. Add any necessary UI augmentation/components to the interface 4. Update the workflow so that it incorporates the new feature 5. Test the new functionality 
6. Update documentation that aligns the augmentations 
Integrating New Models or Tools 
To integrate new models or tools; 
1. Add required dependencies 
2. Create adapter functions to integrate the new tool 
3. Modify the affected module to utilize the new tool 4. Update the UI to highlight the new capabilities 
5. Test thoroughly
9. Testing Framework 
To verify the platform’s functionality, you can perform tests by entering our sample inputs directly into the system: 
● Compound Analysis: Input a SMILES string to check compound evaluation and optimization results. 
● Target Analysis: Enter a protein target name to verify correct retrieval and analysis of target information. 
● Protein Analysis: Provide a valid PDB ID to test protein structure visualization and binding site identification. 
By comparing the system’s outputs with expected results, users can confirm that each module is functioning correctly. This approach allows for straightforward verification. 
Test Examples 
The examples contains sample test cases for the platform, including: 
● Sample SMILES strings for testing compound analysis. 
● Sample target names for testing target analysis. 
● Sample PDB IDs for testing protein analysis.
10. Deployment Guide 
The platform can be set up using multiple deployment methods: ● Hugging Face Spaces: Quick and simple deployment with minimal setup. 
● Local Server: Run the platform on a local machine for internal testing or use. 
To deploy the platform on Hugging Face Spaces: 
● Create a new Space on Hugging Face. 
● Upload the core project files 
● Configure the required environment variables in the Space settings. ● Set the Space SDK to Gradio. 
● Launch the Space to make the platform accessible online. 
To deploy on a local server: 
● Set up the development environment 
● Start the application with: `python ‘Wrapper.py` 
● Access the application at the provided local URL 
Environment Variables 
Required environment variables: 
● GEMINI API key: ‘GEMINI_API_KEY=your_api_key
11. Troubleshooting for Developers 
Issue: Missing dependencies or modules. 
Solution: Ensure all required Python packages are installed. 
Issue: API-related errors 
Solution: Check that your internet connection is active, confirm the gemeni_api_key is set correctly in your environment, and ensure your account has sufficient credits 
Issue: RDKit installation problems 
Solution: Ensure all dependencies are met and your system is compatible, or install RDKit using conda: conda install -c conda-forge rdkit 
Issue: Environment variables not loading 
Solution: Store required variables in a .env file in the repository root and ensure your application reads them using a library like python-dotenv. Issue: Low performance or timeouts 
Solution: Check system resources (CPU and memory), make sure all dependencies are up to date, and verify that any API calls are responding correctly.
12. Contributing Guidelines 
We welcome contributions from developers! To ensure smooth collaboration, please follow these guidelines: 
1. Code of Conduct 
● Be respectful, professional, and inclusive in all communications and contributions. 
● Follow project conventions and community standards. 
2. Getting Started 
● Clone the repository: 
● Set up the development environment as described in the setup instructions. ● Install dependencies 
3. Branching and Workflow 
● Use feature branches for new development 
● Commit changes with clear, descriptive messages. 
● Push your branch and create a pull request for review. 
4. Code Style 
● Follow Python best practices (PEP 8). 
● Keep functions and modules modular and well-documented. ● Include comments where necessary for clarity. 
5. Documentation 
● Update README, user manuals, or in-code comments for any new functionality. 
● Include instructions or examples if the change affects usage.

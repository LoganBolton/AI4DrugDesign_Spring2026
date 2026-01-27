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


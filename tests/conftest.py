"""Pytest configuration for running tests from the workspace root."""

from pathlib import Path
import sys


# The project is currently stored below `03_lab_management/oswald` and is not
# installed into the active Conda environment yet. Add its parent directory
# explicitly so `import oswald` works during this pre-packaging stage.
PROJECT_PARENT = Path(__file__).resolve().parents[1].parent
if str(PROJECT_PARENT) not in sys.path:
    sys.path.insert(0, str(PROJECT_PARENT))
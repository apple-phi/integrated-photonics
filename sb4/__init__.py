"""sb4: A Python package for Lumerical FDTD simulation and analysis."""

__version__ = "0.1.0"

# Make key functions available at the package level
from .simulation import run_simulation
from .results import process_and_save_results

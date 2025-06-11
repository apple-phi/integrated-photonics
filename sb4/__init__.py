"""sb4: A Python package for Lumerical FDTD simulation and analysis."""

__version__ = "0.1.0"

import logging
import rich.logging

FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO",
    format=FORMAT,
    datefmt="[%X]",
    handlers=[rich.logging.RichHandler(markup=True)],
)

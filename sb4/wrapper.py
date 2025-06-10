import pathlib
import sys
import logging

logger = logging.getLogger(__name__)
root = pathlib.Path("C:\\Program Files\\Lumerical")

# Recursive search for a python/lumapi.py file in the specified directory
for dirpath, dirnames, filenames in root.walk():
    if "lumapi.py" in filenames:
        sys.path.append(str(dirpath))
        logger.info(f"Found lumapi at: {dirpath / 'lumapi.py'}")
        break
else:
    raise ImportError(f"lumapi.py not found anywhere in {root}")

import lumapi as lumapi_raw  # type: ignore[import-untyped]

lumapi = lumapi_raw  # For compatibility with existing code
u = 1e-6  # Micrometer unit for convenience

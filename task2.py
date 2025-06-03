import sys
import pathlib

root = pathlib.Path("C:\\Program Files\\Lumerical")

# Recursive search for a python/lumapi.py file in the specified directory
for dirpath, dirnames, filenames in root.walk():
    if "lumapi.py" in filenames:
        sys.path.append(str(dirpath))
        print("Found lumapi at:", dirpath / "lumapi.py")
        break
else:
    raise ImportError("lumapi.py not found in the specified directory tree.")

import lumapi  # type: ignore[import-untyped]
import numpy as np
import matplotlib.pyplot as plt

with lumapi.FDTD() as fdtd:
    lambda_range = np.linspace(300e-9, 1100e-9, 500)
    c = 2.99792458e8
    f_range = c / lambda_range
    au_index = fdtd.getfdtdindex(
        "Au (Gold) - CRC", f_range, np.min(f_range), np.max(f_range)
    )  # Use the getfdtdindex command to obtain the correct complex index for gold

    stackRT_result = fdtd.stackrt(
        np.transpose(au_index), np.array([10e-9]), f_range
    )  # Use the stackrt command to calculate the transmission and reflection
# Visualize using matplotlib
fig, ax = plt.subplots()
ax.plot(lambda_range * 1e9, stackRT_result["Ts"], label="Transmission")
ax.set_xlabel("Wavelength [nm]")
ax.set_ylabel("Transmission")
ax.legend()
plt.show()

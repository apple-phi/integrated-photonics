import os
import time
import sys
import pathlib

root = pathlib.Path("C:\\Program Files\\Lumerical")

# Recursive search for a python/lumapi.py file in the specified directory
for dirpath, dirnames, filenames in root.walk():
    if "lumapi.py" in filenames:
        sys.path.append(str(dirpath))
        print("Found lumapi at: '", dirpath / "lumapi.py'")
        break
else:
    raise ImportError("lumapi.py not found in the specified directory tree.")

import lumapi  # type: ignore[import-untyped]
import numpy as np
import matplotlib.pyplot as plt

u = 1e-6

with lumapi.FDTD(hide=False) as fdtd:
    fdtd.addfdtd(
        x=0,
        y=0,
        z=0,
        x_span=35 * u,
        y_span=18 * u,
        z_span=0.5 * u,
        background_material="SiO2 (Glass) - Palik",
    )
    # Parallel waveguides
    wg1 = fdtd.addrect(
        name="wg_1",
        material="Si (Silicon) - Palik",
        x=0,
        y=0.325 * u,
        z=0,
        x_span=10 * u,
        y_span=0.5 * u,
        z_span=0.22 * u,
    )
    wg2 = fdtd.addrect(
        name="wg_2",
        material="Si (Silicon) - Palik",
        x=0,
        y=-0.325 * u,
        z=0,
        x_span=10 * u,
        y_span=0.5 * u,
        z_span=0.22 * u,
    )
    # Input / Output connecting s-bends
    sbend_tr = fdtd.addwaveguide(
        name="sbend_tr",
        material="Si (Silicon) - Palik",
        x=5 * u,
        y=0.325 * u,
        z=0,
        base_width=0.5 * u,
        base_height=0.22 * u,
        base_angle=90,
    )
    sbend_tr.poles = np.array([[0, 0], [5 * u, 0], [5 * u, 5 * u], [10 * u, 5 * u]])
    sbend_br = fdtd.addwaveguide(
        name="sbend_br",
        material="Si (Silicon) - Palik",
        x=5 * u,
        y=-0.325 * u,
        z=0,
        base_width=0.5 * u,
        base_height=0.22 * u,
        base_angle=90,
    )
    sbend_br.poles = (
        np.array(
            [
                [0, 0],
                [5, 0],
                [5, -5],
                [10, -5],
            ]
        )
        * u
    )
    sbend_tl = fdtd.addwaveguide(
        name="sbend_tl",
        material="Si (Silicon) - Palik",
        x=-5 * u,
        y=0.325 * u,
        z=0,
        base_width=0.5 * u,
        base_height=0.22 * u,
        base_angle=90,
    )
    sbend_tl.poles = (
        np.array(
            [
                [0, 0],
                [-5, 0],
                [-5, 5],
                [-10, 5],
            ]
        )
        * u
    )
    sbend_bl = fdtd.addwaveguide(
        name="sbend_bl",
        material="Si (Silicon) - Palik",
        x=-5 * u,
        y=-0.325 * u,
        z=0,
        base_width=0.5 * u,
        base_height=0.22 * u,
        base_angle=90,
    )
    sbend_bl.poles = (
        np.array(
            [
                [0, 0],
                [-5, 0],
                [-5, -5],
                [-10, -5],
            ]
        )
        * u
    )
    # Input / Output straight waveguides
    wg_tr = fdtd.addrect(
        name="wg_tr",
        material="Si (Silicon) - Palik",
        x=17.5 * u,
        y=5.325 * u,
        z=0,
        x_span=5 * u,
        y_span=0.5 * u,
        z_span=0.22 * u,
    )
    wg_br = fdtd.addrect(
        name="wg_br",
        material="Si (Silicon) - Palik",
        x=17.5 * u,
        y=-5.325 * u,
        z=0,
        x_span=5 * u,
        y_span=0.5 * u,
        z_span=0.22 * u,
    )
    wg_tl = fdtd.addrect(
        name="wg_tl",
        material="Si (Silicon) - Palik",
        x=-17.5 * u,
        y=5.325 * u,
        z=0,
        x_span=5 * u,
        y_span=0.5 * u,
        z_span=0.22 * u,
    )
    wg_bl = fdtd.addrect(
        name="wg_bl",
        material="Si (Silicon) - Palik",
        x=-17.5 * u,
        y=-5.325 * u,
        z=0,
        x_span=5 * u,
        y_span=0.5 * u,
        z_span=0.22 * u,
    )

    # Mode source
    source = fdtd.addmode(
        name="mode_source",
        x=-17 * u,  # need to be just inside FDTD box
        y=5.325 * u,
        z=0,
        injection_axis="x-axis",
        mode_selection="fundamental TE mode",
        center_wavelength=1.55 * u,
        wavelength_span=0,
        # wavelength_start=1.55e-6,
        # wavelength_stop=1.55e-6,
    )

    # Mode monitors
    mon_tr = fdtd.addpower(
        name="mon_tr",
        monitor_type="2D X-normal",
        x=17 * u,  # need to be just inside FDTD box
        y=5.325 * u,
        z=0,
        y_span=2 * u,  # No x-span for 2D x-plane
        z_span=0.25 * u,
    )
    mon_br = fdtd.addpower(
        name="mon_br",
        monitor_type="2D X-normal",
        x=17 * u,  # need to be just inside FDTD box
        y=-5.325 * u,
        z=0,
        y_span=2 * u,
        z_span=0.25 * u,
    )
    mon_znorm = fdtd.addpower(
        name="mon_zplane",
        monitor_type="2D Z-normal",
        x=0,
        y=0,
        z=0,
        x_span=34 * u,
        y_span=14 * u,
    )

    # Mode expansion monitors
    me_tr = fdtd.addmodeexpansion(
        name="me_tr",
        mode_selection="fundamental TE mode",
        monitor_type="2D X-normal",
        x=17 * u,  # need to be just inside FDTD box
        y=-5.325 * u,
        z=0,
        y_span=2 * u,
        z_span=0.25 * u,
    )
    fdtd.setexpansion("input", "mon_tr")
    me_br = fdtd.addmodeexpansion(
        name="me_br",
        mode_selection="fundamental TE mode",
        monitor_type="2D X-normal",
        x=17 * u,  # need to be just inside FDTD box
        y=-5.325 * u,
        z=0,
        y_span=2 * u,
        z_span=0.25 * u,
    )
    fdtd.setexpansion("input", "mon_br")

    # E = fdtd.getresult("mon_znorm", "E")
    # fdtd.visualize(E)

    fdtd.save("task2.fsp")
    fdtd.run()

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


def plot_plane(data, title: str):
    """Helper function to plot a 2D plane from the FDTD result."""
    x = data["x"].flatten() * 1e6  # Convert to µm
    y = data["y"].flatten() * 1e6  # Convert to µm
    E = data["E"]  # shape (nx, ny, 1, 1, 3)

    # Collapse the singleton dimensions and compute |E|^2
    Ex = E[:, :, 0, 0, 0]
    Ey = E[:, :, 0, 0, 1]
    Ez = E[:, :, 0, 0, 2]

    Intensity = np.abs(Ex) ** 2 + np.abs(Ey) ** 2 + np.abs(Ez) ** 2

    X, Y = np.meshgrid(x, y, indexing="ij")  # (764×276) grids

    plt.figure(figsize=(6, 4))
    pcm = plt.pcolormesh(X, Y, Intensity, shading="auto")
    plt.xlabel("x (μm)")
    plt.ylabel("y (μm)")
    plt.title(title)
    plt.colorbar(pcm, label="Intensity (a.u.)")
    plt.savefig(
        f"data/task3/{title.replace(' ', '_').replace('(', '_').replace(')', '_').replace('=', '_').lower()}.png"
    )
    plt.show()


with lumapi.FDTD(hide=False) as fdtd:
    fdtd.load("task3.fsp")  # Sim run from ./task3_sim.py

    # print("Available results:")
    # res_mon_tr = fdtd.getresult("mon_tr")
    # print("mon_tr:", res_mon_tr)
    # res_mon_br = fdtd.getresult("mon_br")
    # print("mon_br:", res_mon_br)
    # res_me_tr = fdtd.getresult("me_tr")
    # print("me_tr:", res_me_tr)
    # res_me_br = fdtd.getresult("me_br")
    # print("me_br:", res_me_br)
    # res_mon_zplane = fdtd.getresult("mon_zplane")
    # print("mon_zplane:", res_mon_zplane)

    res_mon_tr_T = fdtd.getresult("mon_tr", "T")
    print("mon_tr total transmitted power:", res_mon_tr_T["T"].item())
    res_mon_br_T = fdtd.getresult("mon_br", "T")
    print("mon_br total transmitted power:", res_mon_br_T["T"].item())

    res_me_tr_exp = fdtd.getresult("me_tr", "expansion for input")
    T_net_tr = res_me_tr_exp["T_net"].item()
    print("me_tr mode expansion T_net:", T_net_tr)
    res_me_br_exp = fdtd.getresult("me_br", "expansion for input")
    T_net_br = res_me_br_exp["T_net"].item()
    print("me_br mode expansion T_net:", T_net_br)

    res = fdtd.getresult("mon_zplane", "E")
    plot_plane(res, f"Z-plane E field intensity (tr={T_net_tr:.4f}, br={T_net_br:.4f})")

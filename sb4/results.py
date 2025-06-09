"""sb4.results: Processes and saves results from Lumerical FDTD simulations."""

import os
import sys
import pathlib
import json

# import csv # Removed, no longer creating summary CSV
import logging

import numpy as np
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

# Ensure lumapi can be imported
try:
    import lumapi  # type: ignore[import-untyped]
except ImportError as e:
    lumerical_install_dir = pathlib.Path(
        os.getenv("LUMERICAL_INSTALL_DIR", "C:\\\\Program Files\\\\Lumerical")
    )
    found_lumapi = False
    for version_dir in lumerical_install_dir.iterdir():
        if version_dir.is_dir():
            py_api_dir = version_dir / "api" / "python"
            if (py_api_dir / "lumapi.py").exists():
                sys.path.append(str(py_api_dir))
                logger.info(f"Found lumapi at: '{py_api_dir / 'lumapi.py'}'")
                import lumapi  # type: ignore[import-untyped]

                found_lumapi = True
                break
            alt_py_api_dir = version_dir / "python"
            if (alt_py_api_dir / "lumapi.py").exists():
                sys.path.append(str(alt_py_api_dir))
                logger.info(f"Found lumapi at: '{alt_py_api_dir / 'lumapi.py'}'")
                import lumapi  # type: ignore[import-untyped]

                found_lumapi = True
                break
    if not found_lumapi:
        logger.error(
            "lumapi.py not found. Please ensure Lumerical is installed and "
            "LUMERICAL_INSTALL_DIR environment variable is set, or lumapi.py is in PYTHONPATH."
        )
        raise ImportError(
            "lumapi.py not found. Please ensure Lumerical is installed and "
            "LUMERICAL_INSTALL_DIR environment variable is set, or lumapi.py is in PYTHONPATH."
        ) from e

u = 1e-6


def plot_plane_parametric(
    fdtd_obj, title_prefix: str, layout_id: str, target_dir: pathlib.Path
):
    """Helper function to plot a 2D plane from the FDTD result and save it."""
    res_zplane = fdtd_obj.getresult("mon_zplane", "E")
    x = res_zplane["x"].flatten() * 1e6  # Convert to µm
    y = res_zplane["y"].flatten() * 1e6  # Convert to µm
    E_field_data = res_zplane["E"]  # shape (nx, ny, 1, 1, 3)

    Ex = E_field_data[:, :, 0, 0, 0]
    Ey = E_field_data[:, :, 0, 0, 1]
    Ez = E_field_data[:, :, 0, 0, 2]

    Intensity = np.abs(Ex) ** 2 + np.abs(Ey) ** 2 + np.abs(Ez) ** 2
    X, Y = np.meshgrid(x, y, indexing="ij")

    plt.figure(figsize=(8, 5))
    pcm = plt.pcolormesh(X, Y, Intensity, shading="auto", cmap="viridis")
    plt.xlabel("x (μm)")
    plt.ylabel("y (μm)")
    plot_title = f"{title_prefix} - {layout_id}"
    plt.title(plot_title)
    plt.colorbar(pcm, label="Intensity (a.u.)")
    plt.axis("equal")

    plot_filename = f"z_plane_intensity_{layout_id}.png"
    save_path = target_dir / plot_filename
    plt.savefig(save_path)
    logger.info(f"Saved z-plane plot to: {save_path}")
    plt.close()


def process_and_save_results(
    sim_filepath: str,  # Full path to the .fsp file
    params_dict: dict,  # Parameters used for this simulation (in meters)
    output_dir: str,  # Directory where results.json and plots will be saved for this layout
    plot_z_plane: bool = False,
):
    """Loads results from a Lumerical simulation, processes, and saves them
    to the specified output directory.
    """
    output_path = pathlib.Path(output_dir)
    # layout_id is the stem of the simulation file (e.g., "wg1w0p5_wg2w0p5_...")
    layout_id = pathlib.Path(sim_filepath).stem

    logger.info(f"Processing results for layout: {layout_id} from file: {sim_filepath}")
    logger.info(f"Results will be saved in: {output_path}")


    results_data = {}
    results_data["parameters"] = params_dict

    try:
        with lumapi.FDTD(hide=True) as fdtd:
            fdtd.load(sim_filepath)

            # Extract power transmission from mode expansion monitors
            res_me_tr_exp = fdtd.getresult("me_tr", "expansion for me_tr")
            t_net_tr = (
                res_me_tr_exp["T_net"].item()
                if "T_net" in res_me_tr_exp
                else float("nan")
            )
            results_data["T_net_tr"] = t_net_tr
            print(f"  me_tr mode expansion T_net: {t_net_tr:.4f}")

            res_me_br_exp = fdtd.getresult("me_br", "expansion for me_br")
            t_net_br = (
                res_me_br_exp["T_net"].item()
                if "T_net" in res_me_br_exp
                else float("nan")
            )
            results_data["T_net_br"] = t_net_br
            print(f"  me_br mode expansion T_net: {t_net_br:.4f}")

            if plot_z_plane:
                plot_title_prefix = (
                    f"Z-plane E-field Intensity (tr={t_net_tr:.4f}, br={t_net_br:.4f})"
                )
                plot_plane_parametric(
                    fdtd,
                    title_prefix=plot_title_prefix,
                    layout_id=layout_id,
                    target_dir=pathlib.Path(output_dir),
                )

    except Exception as e:
        print(f"Error during Lumerical results processing for {layout_id}: {e}")
        results_data["T_net_tr"] = float("nan")  # Indicate error in results
        results_data["T_net_br"] = float("nan")
        results_data["error"] = str(e)

    results_json_path = output_path / "results.json"
    with open(results_json_path, "w", encoding="utf-8") as f_json:
        json.dump(results_data, f_json, indent=4)
    logger.info(f"Saved detailed results to: {results_json_path}")

    # Summary CSV is removed. Each layout folder has its own results.json.

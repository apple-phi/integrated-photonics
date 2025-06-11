"""sb4.results: Processes and saves results from Lumerical FDTD simulations."""

import os
import sys
import pathlib
import json
import logging
import numpy as np
import matplotlib.pyplot as plt

from .wrapper import lumapi, u

logger = logging.getLogger(__name__)


def plot_plane_parametric(fdtd_obj, title_prefix: str, target_dir: pathlib.Path):
    """Helper function to plot a 2D plane from the FDTD result and save it.
    The plot filename is fixed to 'z_plane_intensity.png'.
    """
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
    plot_title = title_prefix  # Removed run_id from title
    plt.title(plot_title)
    plt.colorbar(pcm, label="Intensity (a.u.)")
    plt.gca().set_aspect("auto")  # Ensure plot fills axes
    plt.tight_layout()  # Adjust layout to prevent overlap

    plot_filename = "z_plane_intensity.png"  # Fixed filename
    save_path = target_dir / plot_filename
    plt.savefig(save_path)
    logger.info(f"Saved z-plane plot to: {save_path}")
    plt.close()


def process_and_save_results(
    sim_filepath: str,  # Full path to the .fsp file (e.g., .../sim_HHMM_DDMMYY/simulation.fsp)
    params_dict: dict,  # Parameters used for this simulation (in MICRONS)
    output_dir: str,  # Directory where results.json and plots will be saved (e.g., .../sim_HHMM_DDMMYY)
    plot_z_plane: bool = False,
):
    """Loads results from a Lumerical simulation, processes, and saves them
    to the specified output directory.
    Assumes fixed filenames like 'results.json' and 'z_plane_intensity.png' within output_dir.
    """
    output_path = pathlib.Path(output_dir)
    # The unique identifier for the run (e.g., sim_YYYYMMDD_HHMMSS_ffffff) is the name of the output_dir itself.
    run_id = output_path.name

    logger.info(f"Processing results for run: {run_id} from file: {sim_filepath}")
    logger.info(f"Results will be saved in: {output_path}")

    results_data = {}
    results_data["parameters"] = params_dict  # Already in microns
    results_data["run_id"] = run_id  # Store run_id in results for traceability

    try:
        with lumapi.FDTD(hide=True) as fdtd:
            fdtd.load(sim_filepath)

            # Extract power transmission from mode expansion monitors
            res_me_tr_exp = fdtd.getresult("me_tr", "expansion for me_tr")
            # Use .item() to get scalar from 0-dim array, handle missing key gracefully
            t_net_tr = res_me_tr_exp.get("T_net").item() if isinstance(res_me_tr_exp.get("T_net"), np.ndarray) and res_me_tr_exp.get("T_net").size == 1 else res_me_tr_exp.get("T_net", float("nan"))
            if not isinstance(t_net_tr, (int, float)):
                t_net_tr = float("nan")  # Ensure it's a number
            results_data["T_net_tr"] = t_net_tr
            logger.info(f"  Run {run_id} - me_tr mode expansion T_net: {t_net_tr:.4f}")

            res_me_br_exp = fdtd.getresult("me_br", "expansion for me_br")
            t_net_br = res_me_br_exp.get("T_net").item() if isinstance(res_me_br_exp.get("T_net"), np.ndarray) and res_me_br_exp.get("T_net").size == 1 else res_me_br_exp.get("T_net", float("nan"))
            if not isinstance(t_net_br, (int, float)):
                t_net_br = float("nan")  # Ensure it's a number
            results_data["T_net_br"] = t_net_br
            logger.info(f"  Run {run_id} - me_br mode expansion T_net: {t_net_br:.4f}")

            if plot_z_plane:
                plot_title_prefix = f"Z-plane E-field Intensity (tr={t_net_tr:.4f}, br={t_net_br:.4f})"
                # Pass run_id (folder name) as layout_id for plot title consistency
                plot_plane_parametric(fdtd, title_prefix=plot_title_prefix, target_dir=output_path)

    except Exception as e:
        logger.error(f"Error during Lumerical results processing for run {run_id}: {e}")
        results_data["T_net_tr"] = float("nan")  # Indicate error in results
        results_data["T_net_br"] = float("nan")
        results_data["error"] = str(e)

    results_json_path = output_path / "results.json"  # Fixed filename
    with open(results_json_path, "w", encoding="utf-8") as f_json:
        json.dump(results_data, f_json, indent=4)
    logger.info(f"Saved detailed results for run {run_id} to: {results_json_path}")

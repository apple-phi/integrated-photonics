"""sb4.simulation: Parametric Lumerical FDTD simulation setup."""

import os
import time
import sys
import pathlib
import numpy as np
import logging
import json
from datetime import datetime  # For timestamp generation
from typing import Tuple, Dict, Any, Optional  # For type hinting

from .wrapper import lumapi, u, TE0, SiO2, Si, Xnorm, Ynorm, Znorm  # Importing constants and lumapi wrapper
from .results import process_and_save_results, plot_plane_parametric  # Ensure plot_plane_parametric is imported

logger = logging.getLogger(__name__)  # Use a named logger


def _find_existing_successful_run(output_base_dir: pathlib.Path, params_to_match: Dict[str, Any]) -> Tuple[Optional[pathlib.Path], Optional[Dict[str, Any]]]:
    """Scans for an existing successfully completed run with matching parameters."""
    if not output_base_dir.exists():
        return None, None

    for run_folder in output_base_dir.iterdir():
        if not run_folder.is_dir():
            continue

        results_json_path = run_folder / "results.json"
        if not results_json_path.exists():
            continue

        try:
            with open(results_json_path, "r", encoding="utf-8") as f:
                existing_results = json.load(f)

            # Parameter comparison
            if existing_results.get("parameters") != params_to_match:
                continue

            # Success criteria
            is_successful = (
                "error" not in existing_results
                and "T_net_tr" in existing_results
                and "T_net_br" in existing_results
                and not np.isnan(existing_results.get("T_net_tr", float("nan")))
                and not np.isnan(existing_results.get("T_net_br", float("nan")))
            )

            if is_successful:
                return run_folder, existing_results
            else:
                logger.info(f"Found existing run in {run_folder} with matching parameters, but it was not successful or results are incomplete.")

        except json.JSONDecodeError:
            logger.warning(f"Could not decode results.json in {run_folder}. Skipping this folder.")
        except Exception as e_check:
            logger.warning(f"Error checking results in {run_folder}: {e_check}. Skipping this folder.")

    return None, None


def _handle_existing_run_plotting(run_folder: pathlib.Path, existing_results: Dict[str, Any], hide_fdtd_gui_for_plot: bool = True):
    """Handles plotting for an existing run if the plot is missing."""
    sim_fsp_path = run_folder / "simulation.fsp"
    plot_png_path = run_folder / "z_plane_intensity.png"
    existing_run_id = run_folder.name

    if plot_png_path.exists():
        logger.info(f"Plot already exists for existing run at {plot_png_path}")
        return

    if not sim_fsp_path.exists():
        logger.warning(f"Simulation file {sim_fsp_path} not found for existing run {existing_run_id}. Cannot generate plot.")
        return

    logger.info(f"Plot not found for existing run {existing_run_id}. Generating plot at {plot_png_path}...")
    try:
        with lumapi.FDTD(hide=hide_fdtd_gui_for_plot) as fdtd_plot:
            fdtd_plot.load(str(sim_fsp_path))
            tr_val = existing_results.get("T_net_tr", float("nan"))
            br_val = existing_results.get("T_net_br", float("nan"))
            plot_title_prefix = f"Z-plane E-field Intensity (tr={tr_val:.4f}, br={br_val:.4f})"

            plot_plane_parametric(fdtd_obj=fdtd_plot, title_prefix=plot_title_prefix, target_dir=run_folder)  # Pass folder name as layout_id for plot title consistency
            logger.info(f"Generated plot for existing run at {plot_png_path}")
    except Exception as e_plot:
        logger.warning(f"Could not generate plot for existing run {existing_run_id}: {e_plot}")


def add_rect(fdtd: lumapi.FDTD, name: str, xyz: tuple, span: tuple, material: str = Si):
    """Helper function to add a rectangular waveguide."""
    x, y, z = xyz
    x_span, y_span, z_span = span
    fdtd.addrect(name=name, material=material, x=x, y=y, z=z, x_span=x_span, y_span=y_span, z_span=z_span)


def add_sbend(fdtd: lumapi.FDTD, name: str, wh: tuple, xyz: tuple, poles: list | np.ndarray, material: str = Si):
    """Helper function to add a waveguide S-bend."""
    base_width, base_height = wh
    x, y, z = xyz
    sbend = fdtd.addwaveguide(name=name, material=material, base_width=base_width, base_height=base_height, base_angle=90, x=x, y=y, z=z)
    sbend["poles"] = np.array(poles)
    return sbend


def run_simulation(params: dict, output_base_dir: str = "task3", plot_z_plane_each_run: bool = False, hide_fdtd_gui: bool = False):
    """
    Sets up and runs a Lumerical FDTD simulation for a dual waveguide coupler
    based on the provided parameters. Saves all outputs to a unique, timestamped
    directory under output_base_dir. Checks for existing successful runs with
    identical parameters before starting a new simulation.

    Args:
        params (dict): Dictionary of simulation parameters (in meters).
        output_base_dir (str): Base directory where run-specific folders will be created.
        plot_z_plane_each_run (bool): Whether to generate and save z-plane plot.
        hide_fdtd_gui (bool): Whether to hide the Lumerical FDTD CAD window.
    """

    # 1. Check for existing successful runs with matching parameters
    base_path = pathlib.Path(output_base_dir)
    existing_run_folder, existing_run_results = _find_existing_successful_run(base_path, params)

    if existing_run_folder and existing_run_results:
        existing_run_id = existing_run_folder.name
        logger.info(f"Found existing successful run ({existing_run_id}) with matching parameters in: {existing_run_folder}")
        print(f"Skipping simulation: Found existing successful run ({existing_run_id}) with matching parameters in {existing_run_folder}")

        if plot_z_plane_each_run:
            _handle_existing_run_plotting(existing_run_folder, existing_run_results, hide_fdtd_gui_for_plot=True)
        return  # Skip new simulation

    # 2. If no matching successful run is found, proceed with new simulation
    run_id = f"sim_{datetime.now().strftime('%Y%m%d_%H%M%S_%f')}"
    logger.info(f"Generated new run ID: {run_id}")

    current_run_output_dir = pathlib.Path(output_base_dir) / run_id
    current_run_output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory for this run: {current_run_output_dir}")

    # Use fixed filenames within the run-specific directory
    sim_filepath = current_run_output_dir / "simulation.fsp"
    # results.json will be created by process_and_save_results in current_run_output_dir

    # Removed the block for checking existing results. Each call creates a new run.

    logger.info(f"Starting simulation for run: {run_id} with parameters: {params}")

    wg1_w = params["wg1_width"]
    wg2_w = params["wg2_width"]
    sep = params["separation"]
    L_coup = params["coupling_length"]
    wg_z = params["wg_z_span"]
    sbend_y_total = params["fan_out_y_offset"]
    sbend_x_total = params["sbend_x_extent"]
    center_wl = params["center_wavelength"]

    monitor_sbend_offset = params["monitor_offset_from_sbend"]
    fdtd_padding = params["fdtd_xy_padding"]
    wg_ext_past_fdtd = params["wg_extension_past_fdtd_edge"]

    # --- Derived Geometric Calculations ---
    # Y positions of the center of the parallel waveguides
    wg1_y_cent = sep / 2 + wg1_w / 2
    wg2_y_cent = -(sep / 2 + wg2_w / 2)

    # X coordinates for coupling section
    coupling_x_start = -L_coup / 2
    coupling_x_end = L_coup / 2

    # Absolute Y positions for I/O waveguides (where S-bends connect)
    io_wg_top_y_abs = wg1_y_cent + sbend_y_total
    io_wg_bottom_y_abs = wg2_y_cent - sbend_y_total

    # X-coordinates for S-bend connection points to coupling WGs
    # These are also the points where S-bends start their x-travel
    s_bend_connect_left_x = coupling_x_start
    s_bend_connect_right_x = coupling_x_end

    # X-coordinates for S-bend ends (where they connect to I/O WGs)
    s_bend_end_left_x = s_bend_connect_left_x - sbend_x_total
    s_bend_end_right_x = s_bend_connect_right_x + sbend_x_total

    # Monitor/Source positions (define the 'active' simulation x-region)
    mon_pos_left_x = s_bend_end_left_x - monitor_sbend_offset
    mon_pos_right_x = s_bend_end_right_x + monitor_sbend_offset

    # FDTD region calculation (X-dimension)
    fdtd_cent_x = (mon_pos_right_x + mon_pos_left_x) / 2
    fdtd_active_x_span = mon_pos_right_x - mon_pos_left_x
    fdtd_tot_x_span = fdtd_active_x_span + 2 * fdtd_padding
    fdtd_actual_left_edge = fdtd_cent_x - fdtd_tot_x_span / 2
    fdtd_actual_right_edge = fdtd_cent_x + fdtd_tot_x_span / 2

    # FDTD region calculation (Y-dimension) - to fully contain I/O WGs + padding
    # Structure y_extremes are based on the I/O waveguides outer edges
    structure_min_y = io_wg_bottom_y_abs - wg2_w / 2
    structure_max_y = io_wg_top_y_abs + wg1_w / 2
    fdtd_cent_y = (structure_max_y + structure_min_y) / 2
    fdtd_tot_y_span = (structure_max_y - structure_min_y) + 2 * fdtd_padding
    fdtd_actual_bottom_edge = fdtd_cent_y - fdtd_tot_y_span / 2
    fdtd_actual_top_edge = fdtd_cent_y + fdtd_tot_y_span / 2

    # FDTD Z-dimension
    fdtd_tot_z_span = wg_z + 2 * fdtd_padding  # Assuming wg_z is the core thickness

    # I/O Waveguide definitions
    # Top-Left I/O WG
    wg_tl_io_connect_x = s_bend_end_left_x
    wg_tl_io_outer_edge_x = fdtd_actual_left_edge - wg_ext_past_fdtd
    wg_tl_io_len = wg_tl_io_connect_x - wg_tl_io_outer_edge_x
    wg_tl_io_x_center = (wg_tl_io_connect_x + wg_tl_io_outer_edge_x) / 2

    # Top-Right I/O WG
    wg_tr_io_connect_x = s_bend_end_right_x
    wg_tr_io_outer_edge_x = fdtd_actual_right_edge + wg_ext_past_fdtd
    wg_tr_io_len = wg_tr_io_outer_edge_x - wg_tr_io_connect_x
    wg_tr_io_x_center = (wg_tr_io_connect_x + wg_tr_io_outer_edge_x) / 2

    # Bottom-Left I/O WG
    wg_bl_io_connect_x = s_bend_end_left_x  # Same x-connection as top-left
    wg_bl_io_outer_edge_x = fdtd_actual_left_edge - wg_ext_past_fdtd
    wg_bl_io_len = wg_bl_io_connect_x - wg_bl_io_outer_edge_x
    wg_bl_io_x_center = (wg_bl_io_connect_x + wg_bl_io_outer_edge_x) / 2

    # Bottom-Right I/O WG
    wg_br_io_connect_x = s_bend_end_right_x  # Same x-connection as top-right
    wg_br_io_outer_edge_x = fdtd_actual_right_edge + wg_ext_past_fdtd
    wg_br_io_len = wg_br_io_outer_edge_x - wg_br_io_connect_x
    wg_br_io_x_center = (wg_br_io_connect_x + wg_br_io_outer_edge_x) / 2

    # --- Assertions for geometry and FDTD boundaries ---
    epsilon = 1e-9  # For float comparisons

    # Assert I/O waveguides extend through FDTD boundaries
    assert wg_tl_io_x_center - wg_tl_io_len / 2 < fdtd_actual_left_edge + epsilon, "TL I/O WG does not extend through left FDTD boundary"
    assert wg_tr_io_x_center + wg_tr_io_len / 2 > fdtd_actual_right_edge - epsilon, "TR I/O WG does not extend through right FDTD boundary"
    assert wg_bl_io_x_center - wg_bl_io_len / 2 < fdtd_actual_left_edge + epsilon, "BL I/O WG does not extend through left FDTD boundary"
    assert wg_br_io_x_center + wg_br_io_len / 2 > fdtd_actual_right_edge - epsilon, "BR I/O WG does not extend through right FDTD boundary"

    # Assert S-bends and coupling waveguides are within the 'active' FDTD region (inside PMLs)
    # Effective PML thickness is fdtd_padding (can be refined if Lumerical has specific PML settings)
    active_region_inner_left = fdtd_actual_left_edge + fdtd_padding
    active_region_inner_right = fdtd_actual_right_edge - fdtd_padding

    # S-bend X extents
    # sb_tl_start_x is coupling_x_start, sb_tl_end_x is s_bend_end_left_x
    assert s_bend_end_left_x > active_region_inner_left - epsilon, "S-bend TL ends within left PML or too close"
    assert coupling_x_start < active_region_inner_right + epsilon, "S-bend TL starts too far right (or coupling too long)"
    assert coupling_x_start > active_region_inner_left - epsilon, "Coupling WG (left side) starts in PML"
    assert s_bend_end_right_x < active_region_inner_right + epsilon, "S-bend TR ends within right PML or too close"
    assert coupling_x_end > active_region_inner_left - epsilon, "S-bend TR starts too far left (or coupling too long)"
    assert coupling_x_end < active_region_inner_right + epsilon, "Coupling WG (right side) ends in PML"

    # Coupling waveguide itself
    assert coupling_x_start >= active_region_inner_left - epsilon, "Coupling WG starts inside left PML"
    assert coupling_x_end <= active_region_inner_right + epsilon, "Coupling WG ends inside right PML"

    # Assert source and monitors are within the active FDTD region
    assert mon_pos_left_x >= active_region_inner_left - epsilon, "Source/Left Monitor is in left PML"
    assert mon_pos_left_x <= active_region_inner_right + epsilon, "Source/Left Monitor is too far right"
    assert mon_pos_right_x <= active_region_inner_right + epsilon, "Right Monitor is in right PML"
    assert mon_pos_right_x >= active_region_inner_left - epsilon, "Right Monitor is too far left"

    # S-bend poles (relative to S-bend start points on coupling waveguide)
    # Assumes S-bend x-extent is split half/half for the two horizontal segments
    sb_x_intermed = sbend_x_total / 2
    poles_tr_rel = [[0, 0], [sb_x_intermed, 0], [sb_x_intermed, sbend_y_total], [sbend_x_total, sbend_y_total]]
    poles_tl_rel = [[0, 0], [-sb_x_intermed, 0], [-sb_x_intermed, sbend_y_total], [-sbend_x_total, sbend_y_total]]
    poles_br_rel = [[0, 0], [sb_x_intermed, 0], [sb_x_intermed, -sbend_y_total], [sbend_x_total, -sbend_y_total]]
    poles_bl_rel = [[0, 0], [-sb_x_intermed, 0], [-sb_x_intermed, -sbend_y_total], [-sbend_x_total, -sbend_y_total]]

    try:
        with lumapi.FDTD(hide=hide_fdtd_gui) as fdtd:
            fdtd.newproject()
            fdtd.addfdtd(x=fdtd_cent_x, y=fdtd_cent_y, z=0, x_span=fdtd_tot_x_span, y_span=fdtd_tot_y_span, z_span=fdtd_tot_z_span, background_material=SiO2, mesh_accuracy=2)

            # Coupling waveguides
            add_rect(fdtd, name="wg_1_coupling", xyz=(0, wg1_y_cent, 0), span=(L_coup, wg1_w, wg_z))
            add_rect(fdtd, name="wg_2_coupling", xyz=(0, wg2_y_cent, 0), span=(L_coup, wg2_w, wg_z))

            # I/O Waveguides
            add_rect(fdtd, name="wg_tl_io", xyz=(wg_tl_io_x_center, io_wg_top_y_abs, 0), span=(wg_tl_io_len, wg1_w, wg_z))
            add_rect(fdtd, name="wg_tr_io", xyz=(wg_tr_io_x_center, io_wg_top_y_abs, 0), span=(wg_tr_io_len, wg1_w, wg_z))
            add_rect(fdtd, name="wg_bl_io", xyz=(wg_bl_io_x_center, io_wg_bottom_y_abs, 0), span=(wg_bl_io_len, wg2_w, wg_z))
            add_rect(fdtd, name="wg_br_io", xyz=(wg_br_io_x_center, io_wg_bottom_y_abs, 0), span=(wg_br_io_len, wg2_w, wg_z))

            # S-Bends
            # Start S-bends from the connection points on the coupling waveguides
            add_sbend(fdtd, name="sbend_tl", wh=(wg1_w, wg_z), xyz=(s_bend_connect_left_x, wg1_y_cent, 0), poles=poles_tl_rel)
            add_sbend(fdtd, name="sbend_tr", wh=(wg1_w, wg_z), xyz=(s_bend_connect_right_x, wg1_y_cent, 0), poles=poles_tr_rel)
            add_sbend(fdtd, name="sbend_bl", wh=(wg2_w, wg_z), xyz=(s_bend_connect_left_x, wg2_y_cent, 0), poles=poles_bl_rel)
            add_sbend(fdtd, name="sbend_br", wh=(wg2_w, wg_z), xyz=(s_bend_connect_right_x, wg2_y_cent, 0), poles=poles_br_rel)

            # Mode Source (on wg_tl_io, at monitor_pos_left_x)
            fdtd.addmode(
                name="mode_source",
                x=mon_pos_left_x,
                y=io_wg_top_y_abs,
                z=0,
                injection_axis="x-axis",
                direction="forward",
                mode_selection=TE0,
                center_wavelength=center_wl,
                wavelength_span=0,
                y_span=wg1_w + 2 * u,
                z_span=wg_z + 2 * u,
            )

            # Power Monitors (at monitor_pos_right_x)
            mon_y_span = max(wg1_w, wg2_w) + 3 * u
            mon_z_span = wg_z + 3 * u

            fdtd.addpower(name="mon_tr", monitor_type=Xnorm, x=mon_pos_right_x, y=io_wg_top_y_abs, z=0, y_span=mon_y_span, z_span=mon_z_span)
            fdtd.addpower(name="mon_br", monitor_type=Xnorm, x=mon_pos_right_x, y=io_wg_bottom_y_abs, z=0, y_span=mon_y_span, z_span=mon_z_span)

            # Z-plane field monitor
            # Spans the entire FDTD box x and y dimensions
            fdtd.addpower(name="mon_zplane", monitor_type=Znorm, x=fdtd_cent_x, y=fdtd_cent_y, z=0, x_span=fdtd_active_x_span, y_span=fdtd_tot_y_span)

            # Mode Expansion Monitors (co-located with power monitors)
            fdtd.addmodeexpansion(name="me_tr", mode_selection=TE0, monitor_type=Xnorm, x=mon_pos_right_x, y=io_wg_top_y_abs, z=0, y_span=mon_y_span, z_span=mon_z_span)
            fdtd.setexpansion("me_tr", "mon_tr")
            fdtd.addmodeexpansion(name="me_br", mode_selection=TE0, monitor_type=Xnorm, x=mon_pos_right_x, y=io_wg_bottom_y_abs, z=0, y_span=mon_y_span, z_span=mon_z_span)
            fdtd.setexpansion("me_br", "mon_br")

            logger.info(f"Saving simulation file to: {sim_filepath}")
            fdtd.save(str(sim_filepath))

            logger.info(f"Running simulation for {run_id}...")
            fdtd.run()
            logger.info(f"Simulation {run_id} finished.")

            # --- End of Core Lumerical Simulation Setup and Execution ---

    except Exception as e:
        logger.error(f"Error during Lumerical FDTD setup or run for {run_id}: {e}")
        # Save error information to a results.json in the run directory
        error_results_path = current_run_output_dir / "results.json"
        error_data = {"parameters": params, "error": str(e), "run_id": run_id}
        with open(error_results_path, "w", encoding="utf-8") as f_err:
            json.dump(error_data, f_err, indent=4)
        logger.info(f"Error details saved to {error_results_path}")
        return  # Stop further processing for this run

    # Call process_and_save_results outside the FDTD try/except block
    # It will handle its own Lumerical session for loading results.
    logger.info(f"Proceeding to process and save results for run {run_id}.")
    process_and_save_results(
        sim_filepath=str(sim_filepath),
        params_dict=params,  # Pass original params (in meters)
        output_dir=str(current_run_output_dir),  # Pass the specific dir for this run's results
        plot_z_plane=plot_z_plane_each_run,
    )
    logger.info(f"Finished processing and saving results for run {run_id}.")

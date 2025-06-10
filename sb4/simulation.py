"""sb4.simulation: Parametric Lumerical FDTD simulation setup."""

import os
import time
import sys
import pathlib
import numpy as np
import logging
import json  # Added for checking existing results
from wrapper import lumapi, u
from .results import plot_plane_parametric, process_and_save_results

logger = logging.getLogger(__name__)  # Use a named logger


def create_layout_id(params: dict, u_val: float = 1e-6) -> str:
    """Creates a unique, filesystem-friendly ID from simulation parameters."""
    # Scale back to nominal values for readability in dirname
    # Ensure consistent order of parameters for consistent IDs
    # Using fixed precision for float to string conversion
    param_strs = []
    # Define a specific order for parameters to ensure consistency
    param_keys_ordered = [
        "wg1_width",
        "wg2_width",
        "separation",
        "coupling_length",
        "wg_z_span",
        "fan_out_y_offset",
        "sbend_x_extent",
        "center_wavelength",
    ]
    for key in param_keys_ordered:
        if key in params:
            value = params[key] / u_val
            # Format to a reasonable number of decimal places, avoiding scientific notation for typical values
            if isinstance(value, float):
                # Use a format that handles integers and floats gracefully
                param_strs.append(
                    f"{key.replace('_width', 'w').replace('_length', 'len').replace('separation', 'sep').replace('wavelength', 'wl').replace('_', '').replace('center', 'c').replace('extent', 'ext').replace('offset', 'off')}{value:.3g}".replace(
                        ".", "p"
                    )
                )
            else:
                param_strs.append(
                    f"{key.replace('_width', 'w').replace('_length', 'len').replace('separation', 'sep').replace('wavelength', 'wl').replace('_', '').replace('center', 'c').replace('extent', 'ext').replace('offset', 'off')}{value}".replace(
                        ".", "p"
                    )
                )

    return "_".join(param_strs)


def add_rect(
    fdtd: lumapi.FDTD,
    name: str,
    xyz: tuple,
    span: tuple,
    material: str = "Si (Silicon) - Palik",
):
    """Helper function to add a rectangular waveguide."""
    x, y, z = xyz
    x_span, y_span, z_span = span
    fdtd.addrect(
        name=name,
        material=material,
        x=x,
        y=y,
        z=z,
        x_span=x_span,
        y_span=y_span,
        z_span=z_span,
    )


def add_sbend(
    fdtd: lumapi.FDTD,
    name: str,
    wh: tuple,
    xyz: tuple,
    poles: np.ndarray,
    material: str = "Si (Silicon) - Palik",
):
    """Helper function to add a waveguide S-bend."""
    base_width, base_height = wh
    x, y, z = xyz
    sbend = fdtd.addwaveguide(
        name=name,
        material=material,
        base_width=base_width,
        base_height=base_height,
        base_angle=90,
        x=x,
        y=y,
        z=z,
    )
    sbend["poles"] = poles
    return sbend


def run_simulation(
    params: dict,
    output_base_dir: str = "task3",  # Changed from sim_files_base_dir and results_base_dir
    plot_z_plane_each_run: bool = False,
    hide_fdtd_gui: bool = False,
):
    """
    Sets up and runs a Lumerical FDTD simulation for a dual waveguide coupler
    based on the provided parameters. Saves all outputs to a unique directory
    per layout under output_base_dir. Skips simulation if valid results exist.

    Args:
        params (dict): Dictionary of simulation parameters (in meters).
        output_base_dir (str): Base directory where layout-specific folders will be created.
        plot_z_plane_each_run (bool): Whether to generate and save z-plane plot.
        hide_fdtd_gui (bool): Whether to hide the Lumerical FDTD CAD window.
    """
    layout_id = create_layout_id(params, u_val=u)
    logger.info(f"Generated layout ID: {layout_id}")

    layout_output_dir = pathlib.Path(output_base_dir) / layout_id
    layout_output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory for this layout: {layout_output_dir}")

    sim_filepath = layout_output_dir / f"{layout_id}.fsp"
    results_json_path = layout_output_dir / "results.json"

    # Check if valid results already exist
    if results_json_path.exists():
        try:
            with open(results_json_path, "r", encoding="utf-8") as f:
                existing_results = json.load(f)
            # Define what constitutes "valid" results (e.g., has expected keys, no error field)
            if (
                "T_net_tr" in existing_results
                and "T_net_br" in existing_results
                and "error" not in existing_results
            ):
                logger.info(
                    f"Valid results found at {results_json_path}. Skipping simulation for layout {layout_id}."
                )
                # Optionally, regenerate plot if requested and not present, or if forced
                if plot_z_plane_each_run:
                    plot_filename = f"z_plane_intensity_{layout_id}.png"
                    plot_save_path = layout_output_dir / plot_filename
                    if not plot_save_path.exists():
                        logger.info(
                            f"Plot {plot_filename} not found. Generating from existing simulation file: {sim_filepath}"
                        )
                        if sim_filepath.exists():
                            try:
                                with lumapi.FDTD(
                                    hide=True
                                ) as fdtd:  # Use hide=True for consistency
                                    fdtd.load(str(sim_filepath))
                                    # Ensure params_dict is available for plot_plane_parametric if needed by title
                                    # For now, title prefix is simple
                                    plot_title_prefix = f"Z-plane E-field Intensity (tr={existing_results.get('T_net_tr', float('nan')):.4f}, br={existing_results.get('T_net_br', float('nan')):.4f})"
                                    plot_plane_parametric(
                                        fdtd,
                                        title_prefix=plot_title_prefix,
                                        layout_id=layout_id,  # Use layout_id as run_id for plotting
                                        target_dir=layout_output_dir,
                                    )
                            except Exception as e_plot:
                                logger.warning(
                                    f"Could not generate plot for existing results: {e_plot}"
                                )
                        else:
                            logger.warning(
                                f"Simulation file {sim_filepath} not found. Cannot generate plot."
                            )
                    else:
                        logger.info(f"Plot {plot_filename} already exists.")

                # Ensure results are "processed" in terms of logging, even if loaded from file
                logger.info(
                    f"Loaded existing results for {layout_id}: TR={existing_results.get('T_net_tr', 'N/A')}, BR={existing_results.get('T_net_br', 'N/A')}"
                )
                return  # Exit early as simulation is skipped
            else:
                logger.info(
                    f"Existing results.json at {results_json_path} is incomplete or indicates a previous error. Rerunning simulation."
                )
        except json.JSONDecodeError:
            logger.warning(
                f"Could not decode existing results.json at {results_json_path}. Rerunning simulation."
            )
        except Exception as e_check:
            logger.warning(
                f"Error checking existing results at {results_json_path}: {e_check}. Rerunning simulation."
            )

    logger.info(
        f"Starting simulation for layout: {layout_id} with parameters: {params}"
    )

    # Parameter unpacking
    wg1_w = params.get("wg1_width", 0.5 * u)
    wg2_w = params.get("wg2_width", 0.5 * u)
    sep = params.get("separation", 0.15 * u)
    L_coupling = params.get("coupling_length", 10 * u)
    wg_z = params.get("wg_z_span", 0.22 * u)
    sbend_y_total = params.get("fan_out_y_offset", 5 * u)
    sbend_x_total = params.get("sbend_x_extent", 10 * u)
    center_wl = params.get("center_wavelength", 1.55 * u)

    # New parameters for precise FDTD and I/O WG control
    monitor_sbend_offset = params.get("monitor_offset_from_sbend", 2.5 * u)
    fdtd_padding = params.get("fdtd_xy_padding", 2.5 * u)
    wg_ext_past_fdtd = params.get("wg_extension_past_fdtd_edge", 1.0 * u)

    # --- Derived Geometric Calculations ---
    # Y positions of the center of the parallel waveguides
    wg1_y_center = sep / 2 + wg1_w / 2
    wg2_y_center = -(sep / 2 + wg2_w / 2)

    # X coordinates for coupling section
    coupling_x_start = -L_coupling / 2
    coupling_x_end = L_coupling / 2

    # Absolute Y positions for I/O waveguides (where S-bends connect)
    io_wg_top_y_abs = wg1_y_center + sbend_y_total
    io_wg_bottom_y_abs = wg2_y_center - sbend_y_total

    # X-coordinates for S-bend connection points to coupling WGs
    # These are also the points where S-bends start their x-travel
    s_bend_connect_left_x = coupling_x_start
    s_bend_connect_right_x = coupling_x_end

    # X-coordinates for S-bend ends (where they connect to I/O WGs)
    s_bend_end_left_x = s_bend_connect_left_x - sbend_x_total
    s_bend_end_right_x = s_bend_connect_right_x + sbend_x_total

    # Monitor/Source positions (define the 'active' simulation x-region)
    monitor_pos_left_x = s_bend_end_left_x - monitor_sbend_offset
    monitor_pos_right_x = s_bend_end_right_x + monitor_sbend_offset

    # FDTD region calculation (X-dimension)
    fdtd_center_x = (monitor_pos_right_x + monitor_pos_left_x) / 2
    fdtd_active_x_span = monitor_pos_right_x - monitor_pos_left_x
    fdtd_total_x_span = fdtd_active_x_span + 2 * fdtd_padding
    fdtd_actual_left_edge = fdtd_center_x - fdtd_total_x_span / 2
    fdtd_actual_right_edge = fdtd_center_x + fdtd_total_x_span / 2

    # FDTD region calculation (Y-dimension) - to fully contain I/O WGs + padding
    # Structure y_extremes are based on the I/O waveguides outer edges
    structure_min_y = io_wg_bottom_y_abs - wg2_w / 2
    structure_max_y = io_wg_top_y_abs + wg1_w / 2
    fdtd_center_y = (structure_max_y + structure_min_y) / 2
    fdtd_total_y_span = (structure_max_y - structure_min_y) + 2 * fdtd_padding
    fdtd_actual_bottom_edge = fdtd_center_y - fdtd_total_y_span / 2
    fdtd_actual_top_edge = fdtd_center_y + fdtd_total_y_span / 2

    # FDTD Z-dimension
    fdtd_total_z_span = wg_z + 2 * fdtd_padding  # Assuming wg_z is the core thickness

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
    assert (
        wg_tl_io_x_center - wg_tl_io_len / 2 < fdtd_actual_left_edge + epsilon
    ), "TL I/O WG does not extend through left FDTD boundary"
    assert (
        wg_tr_io_x_center + wg_tr_io_len / 2 > fdtd_actual_right_edge - epsilon
    ), "TR I/O WG does not extend through right FDTD boundary"
    assert (
        wg_bl_io_x_center - wg_bl_io_len / 2 < fdtd_actual_left_edge + epsilon
    ), "BL I/O WG does not extend through left FDTD boundary"
    assert (
        wg_br_io_x_center + wg_br_io_len / 2 > fdtd_actual_right_edge - epsilon
    ), "BR I/O WG does not extend through right FDTD boundary"

    # Assert S-bends and coupling waveguides are within the 'active' FDTD region (inside PMLs)
    # Effective PML thickness is fdtd_padding (can be refined if Lumerical has specific PML settings)
    active_region_inner_left = fdtd_actual_left_edge + fdtd_padding
    active_region_inner_right = fdtd_actual_right_edge - fdtd_padding

    # S-bend X extents
    # sb_tl_start_x is coupling_x_start, sb_tl_end_x is s_bend_end_left_x
    assert (
        s_bend_end_left_x > active_region_inner_left - epsilon
    ), "S-bend TL ends within left PML or too close"
    assert (
        coupling_x_start < active_region_inner_right + epsilon
    ), "S-bend TL starts too far right (or coupling too long)"
    assert (
        coupling_x_start > active_region_inner_left - epsilon
    ), "Coupling WG (left side) starts in PML"

    assert (
        s_bend_end_right_x < active_region_inner_right + epsilon
    ), "S-bend TR ends within right PML or too close"
    assert (
        coupling_x_end > active_region_inner_left - epsilon
    ), "S-bend TR starts too far left (or coupling too long)"
    assert (
        coupling_x_end < active_region_inner_right + epsilon
    ), "Coupling WG (right side) ends in PML"

    # Coupling waveguide itself
    assert (
        coupling_x_start >= active_region_inner_left - epsilon
    ), "Coupling WG starts inside left PML"
    assert (
        coupling_x_end <= active_region_inner_right + epsilon
    ), "Coupling WG ends inside right PML"

    # Assert source and monitors are within the active FDTD region
    assert (
        monitor_pos_left_x >= active_region_inner_left - epsilon
    ), "Source/Left Monitor is in left PML"
    assert (
        monitor_pos_left_x <= active_region_inner_right + epsilon
    ), "Source/Left Monitor is too far right"
    assert (
        monitor_pos_right_x <= active_region_inner_right + epsilon
    ), "Right Monitor is in right PML"
    assert (
        monitor_pos_right_x >= active_region_inner_left - epsilon
    ), "Right Monitor is too far left"

    # S-bend poles (relative to S-bend start points on coupling waveguide)
    # Assumes S-bend x-extent is split half/half for the two horizontal segments
    sb_x_intermediate = sbend_x_total / 2

    poles_tr_rel = np.array(
        [
            [0, 0],
            [sb_x_intermediate, 0],
            [sb_x_intermediate, sbend_y_total],
            [sbend_x_total, sbend_y_total],
        ]
    )
    poles_tl_rel = np.array(
        [
            [0, 0],
            [-sb_x_intermediate, 0],
            [-sb_x_intermediate, sbend_y_total],
            [-sbend_x_total, sbend_y_total],
        ]
    )
    poles_br_rel = np.array(
        [
            [0, 0],
            [sb_x_intermediate, 0],
            [sb_x_intermediate, -sbend_y_total],
            [sbend_x_total, -sbend_y_total],
        ]
    )
    poles_bl_rel = np.array(
        [
            [0, 0],
            [-sb_x_intermediate, 0],
            [-sb_x_intermediate, -sbend_y_total],
            [-sbend_x_total, -sbend_y_total],
        ]
    )

    try:
        with lumapi.FDTD(hide=hide_fdtd_gui) as fdtd:
            fdtd.newproject()
            fdtd.addfdtd(
                x=fdtd_center_x,
                y=fdtd_center_y,
                z=0,
                x_span=fdtd_total_x_span,
                y_span=fdtd_total_y_span,
                z_span=fdtd_total_z_span,
                background_material="SiO2 (Glass) - Palik",
                mesh_accuracy=2,
            )

            # Parallel (Coupling) Waveguides
            fdtd.addrect(
                name="wg_1_coupling",
                material="Si (Silicon) - Palik",
                x=0,
                y=wg1_y_center,
                z=0,
                x_span=L_coupling,
                y_span=wg1_w,
                z_span=wg_z,
            )
            fdtd.addrect(
                name="wg_2_coupling",
                material="Si (Silicon) - Palik",
                x=0,
                y=wg2_y_center,
                z=0,
                x_span=L_coupling,
                y_span=wg2_w,
                z_span=wg_z,
            )

            # Input/Output Straight Waveguides
            fdtd.addrect(
                name="wg_tl_io",
                material="Si (Silicon) - Palik",
                x=wg_tl_io_x_center,
                y=io_wg_top_y_abs,
                z=0,
                x_span=wg_tl_io_len,
                y_span=wg1_w,
                z_span=wg_z,
            )
            fdtd.addrect(
                name="wg_tr_io",
                material="Si (Silicon) - Palik",
                x=wg_tr_io_x_center,
                y=io_wg_top_y_abs,
                z=0,
                x_span=wg_tr_io_len,
                y_span=wg1_w,
                z_span=wg_z,
            )
            fdtd.addrect(
                name="wg_bl_io",
                material="Si (Silicon) - Palik",
                x=wg_bl_io_x_center,
                y=io_wg_bottom_y_abs,
                z=0,
                x_span=wg_bl_io_len,
                y_span=wg2_w,
                z_span=wg_z,
            )
            fdtd.addrect(
                name="wg_br_io",
                material="Si (Silicon) - Palik",
                x=wg_br_io_x_center,
                y=io_wg_bottom_y_abs,
                z=0,
                x_span=wg_br_io_len,
                y_span=wg2_w,
                z_span=wg_z,
            )

            # S-Bends
            # Start S-bends from the connection points on the coupling waveguides
            sbtl = fdtd.addwaveguide(
                name="sbend_tl",
                material="Si (Silicon) - Palik",
                base_width=wg1_w,
                base_height=wg_z,
                base_angle=90,
                x=s_bend_connect_left_x,
                y=wg1_y_center,
                z=0,
            )
            sbtl["poles"] = poles_tl_rel

            sbtr = fdtd.addwaveguide(
                name="sbend_tr",
                material="Si (Silicon) - Palik",
                base_width=wg1_w,
                base_height=wg_z,
                base_angle=90,
                x=s_bend_connect_right_x,
                y=wg1_y_center,
                z=0,
            )
            sbtr["poles"] = poles_tr_rel

            sbbl = fdtd.addwaveguide(
                name="sbend_bl",
                material="Si (Silicon) - Palik",
                base_width=wg2_w,
                base_height=wg_z,
                base_angle=90,
                x=s_bend_connect_left_x,
                y=wg2_y_center,
                z=0,
            )
            sbbl["poles"] = poles_bl_rel

            sbbr = fdtd.addwaveguide(
                name="sbend_br",
                material="Si (Silicon) - Palik",
                base_width=wg2_w,
                base_height=wg_z,
                base_angle=90,
                x=s_bend_connect_right_x,
                y=wg2_y_center,
                z=0,
            )
            sbbr["poles"] = poles_br_rel

            # Mode Source (on wg_tl_io, at monitor_pos_left_x)
            fdtd.addmode(
                name="mode_source",
                x=monitor_pos_left_x,
                y=io_wg_top_y_abs,
                z=0,
                injection_axis="x-axis",
                direction="forward",
                mode_selection="fundamental TE mode",  # Parameterize if needed
                center_wavelength=center_wl,
                wavelength_span=0,
                y_span=wg1_w + 2 * u,
                z_span=wg_z + 2 * u,  # Span to cover mode well
            )

            # Power Monitors (at monitor_pos_right_x)
            monitor_y_span_val = max(wg1_w, wg2_w) + 3 * u
            monitor_z_span_val = wg_z + 3 * u

            fdtd.addpower(
                name="mon_tr",
                monitor_type="2D X-normal",
                x=monitor_pos_right_x,
                y=io_wg_top_y_abs,
                z=0,
                y_span=monitor_y_span_val,
                z_span=monitor_z_span_val,
            )
            fdtd.addpower(
                name="mon_br",
                monitor_type="2D X-normal",
                x=monitor_pos_right_x,
                y=io_wg_bottom_y_abs,
                z=0,
                y_span=monitor_y_span_val,
                z_span=monitor_z_span_val,
            )

            # Z-plane field monitor
            # Spans the entire FDTD box x and y dimensions
            fdtd.addpower(
                name="mon_zplane",
                monitor_type="2D Z-normal",
                x=fdtd_center_x,
                y=fdtd_center_y,
                z=0,
                x_span=fdtd_active_x_span,
                y_span=fdtd_total_y_span,
            )

            # Mode Expansion Monitors (co-located with power monitors)
            me_tr = fdtd.addmodeexpansion(
                name="me_tr",
                mode_selection="fundamental TE mode",
                monitor_type="2D X-normal",
                x=monitor_pos_right_x,
                y=io_wg_top_y_abs,
                z=0,
                y_span=monitor_y_span_val,
                z_span=monitor_z_span_val,
            )
            fdtd.setexpansion("me_tr", "mon_tr")

            me_br = fdtd.addmodeexpansion(
                name="me_br",
                mode_selection="fundamental TE mode",
                monitor_type="2D X-normal",
                x=monitor_pos_right_x,
                y=io_wg_bottom_y_abs,
                z=0,
                y_span=monitor_y_span_val,
                z_span=monitor_z_span_val,
            )
            fdtd.setexpansion("me_br", "mon_br")

            print(f"Saving simulation file to: {sim_filepath}")
            fdtd.save(str(sim_filepath))

            print(f"Running simulation for {layout_id}...")
            fdtd.run()
            logger.info(f"Simulation {layout_id} finished.")

            # --- Process Results ---
            # This call will now happen *after* the FDTD context manager closes,
            # ensuring the .fsp file is fully written and closed by Lumerical.
            # The process_and_save_results function will reopen it.

    except Exception as e:
        logger.error(f"Error during Lumerical FDTD setup or run for {layout_id}: {e}")
        # Save error information to results.json
        error_results = {"parameters": params, "error": str(e)}
        with open(results_json_path, "w", encoding="utf-8") as f_err:
            json.dump(error_results, f_err, indent=4)
        logger.info(f"Error details saved to {results_json_path}")
        return  # Stop further processing for this run

    # Call process_and_save_results outside the FDTD try/except block
    # It will handle its own Lumerical session for loading results.
    logger.info(f"Proceeding to process and save results for {layout_id}.")
    process_and_save_results(
        sim_filepath=str(sim_filepath),
        params_dict=params,  # Pass original params (in meters)
        output_dir=str(
            layout_output_dir
        ),  # Pass the specific dir for this layout's results
        plot_z_plane=plot_z_plane_each_run,
    )
    logger.info(f"Finished processing and saving results for {layout_id}.")

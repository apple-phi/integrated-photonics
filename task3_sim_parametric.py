"""Parametric Lumerical FDTD simulation setup for a dual waveguide coupler."""

import os
import time
import sys
import pathlib
import numpy as np

# Attempt to find lumapi dynamically if direct import fails
try:
    import lumapi  # type: ignore[import-untyped]
except ImportError as e:
    lumerical_install_dir = pathlib.Path(
        os.getenv("LUMERICAL_INSTALL_DIR", "C:\\Program Files\\Lumerical")
    )
    found_lumapi = False
    for version_dir in lumerical_install_dir.iterdir():
        if version_dir.is_dir():
            py_api_dir = version_dir / "api" / "python"
            if (py_api_dir / "lumapi.py").exists():
                sys.path.append(str(py_api_dir))
                print(f"Found lumapi at: '{py_api_dir / 'lumapi.py'}'")
                import lumapi  # type: ignore[import-untyped]

                found_lumapi = True
                break
            alt_py_api_dir = version_dir / "python"
            if (alt_py_api_dir / "lumapi.py").exists():
                sys.path.append(str(alt_py_api_dir))
                print(f"Found lumapi at: '{alt_py_api_dir / 'lumapi.py'}'")
                import lumapi  # type: ignore[import-untyped]

                found_lumapi = True
                break
    if not found_lumapi:
        raise ImportError(
            "lumapi.py not found. Please ensure Lumerical is installed and "
            "LUMERICAL_INSTALL_DIR environment variable is set, or lumapi.py is in PYTHONPATH."
        ) from e

# Import the results processing function
# Assuming task3_res_parametric.py is in the same directory or PYTHONPATH
from task3_res_parametric import process_and_save_results

u = 1e-6  # micrometers


def run_simulation(
    params: dict,
    run_id: str,
    sim_files_base_dir: str = "data/task3_simulations",
    results_base_dir: str = "data/task3_results_parametric",
    plot_z_plane_each_run: bool = False,
    hide_fdtd_gui: bool = True,
):
    """
    Sets up and runs a Lumerical FDTD simulation for a dual waveguide coupler
    based on the provided parameters.

    Args:
        params (dict): Dictionary of simulation parameters including:
            wg1_width (float): Width of the top waveguide (meters).
            wg2_width (float): Width of the bottom waveguide (meters).
            separation (float): Edge-to-edge separation between waveguides (meters).
            coupling_length (float): Length of the parallel coupling section (meters).
            wg_z_span (float): Thickness of the waveguides (meters).
            fan_out_y_offset (float): Y-offset for the S-bend fan-out (meters).
            sbend_x_extent (float): X-extent of the S-bends (meters).
            # wg_io_length (float): Length of the input/output straight sections (meters). REMOVED - NOW DYNAMIC
            center_wavelength (float): Center wavelength for source/monitors (meters).
            monitor_offset_from_sbend (float): X-distance from S-bend end to monitor/source.
            fdtd_xy_padding (float): Padding for FDTD region boundaries (for PMLs + margin).
            wg_extension_past_fdtd_edge (float): How far I/O WGs extend beyond FDTD boundary.
        run_id (str): Unique identifier for this simulation run.
        sim_files_base_dir (str): Base directory to save .fsp simulation files.
        results_base_dir (str): Base directory to save processed results.
        plot_z_plane_each_run (bool): Whether to generate and save z-plane plot for this run.
        hide_fdtd_gui (bool): Whether to hide the Lumerical FDTD CAD window.
    """
    print(f"Starting simulation run: {run_id} with parameters: {params}")

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

    # --- Lumerical FDTD Setup ---
    sim_output_dir = pathlib.Path(sim_files_base_dir) / run_id
    sim_output_dir.mkdir(parents=True, exist_ok=True)
    sim_filepath = str(sim_output_dir / f"simulation_{run_id}.fsp")

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
            fdtd.save(sim_filepath)

            print(f"Running simulation for {run_id}...")
            # fdtd.run() # This would run in the foreground if hide=False
            # For background run or more control, consider fdtd.runsweep() or scripting engine commands
            # For now, a simple run. If it hangs, os.environ["QT_QPA_PLATFORM"] = "offscreen" might be needed before lumapi import.
            run_script = f"""
run;
"""
            fdtd.eval(run_script)  # Using eval to ensure it completes
            # A long sleep might indicate the simulation is running in a separate process
            # and this script might proceed before completion if not handled carefully.
            # Checking for simulation end via Lumerical's job manager or file status might be more robust.
            # For now, assume `eval("run;")` blocks or Lumerical handles it.
            print(f"Simulation {run_id} likely completed or running in background.")

    except Exception as e:
        print(f"Error during Lumerical simulation setup or run for {run_id}: {e}")
        # Optionally, save error information to a file or re-raise
        # For now, we will let process_and_save_results handle missing/corrupt files if run fails
        pass  # Allow processing to be attempted, it will likely fail to load if sim didn't save

    # After simulation (or attempted simulation), process results
    print(f"Proceeding to process results for {run_id}...")
    process_and_save_results(
        sim_filepath=sim_filepath,
        params_dict=params,
        run_id=run_id,
        results_base_dir=results_base_dir,
        plot_z_plane=plot_z_plane_each_run,
    )
    print(f"Finished processing for run: {run_id}")


if __name__ == "__main__":
    print("Running example for task3_sim_parametric.py")

    default_params = {
        "wg1_width": 0.5 * u,
        "wg2_width": 0.5 * u,
        "separation": 0.15 * u,
        "coupling_length": 10 * u,
        "wg_z_span": 0.22 * u,
        "fan_out_y_offset": 5 * u,
        "sbend_x_extent": 10 * u,
        "center_wavelength": 1.55 * u,
        "monitor_offset_from_sbend": 2.5 * u,  # Increased offset reflected here
        "fdtd_xy_padding": 2.5 * u,
        "wg_extension_past_fdtd_edge": 1.0 * u,
    }
    example_run_id = (
        f"w1_{default_params['wg1_width']/u:.2f}_w2_{default_params['wg2_width']/u:.2f}_"
        f"sep_{default_params['separation']/u:.2f}_L_{default_params['coupling_length']/u:.0f}_"
        f"sbX_{default_params['sbend_x_extent']/u:.0f}_sbY_{default_params['fan_out_y_offset']/u:.0f}_"
        f"monOff_{default_params['monitor_offset_from_sbend']/u:.1f}_"
        f"fdtdPad_{default_params['fdtd_xy_padding']/u:.1f}_"
        f"wgExt_{default_params['wg_extension_past_fdtd_edge']/u:.1f}"
    )

    # Ensure the QT_QPA_PLATFORM is set if running headlessly or issues occur
    # os.environ["QT_QPA_PLATFORM"] = "offscreen" # Uncomment if GUI issues

    run_simulation(
        params=default_params,
        run_id=example_run_id,
        sim_files_base_dir="data/task3_simulations_parametric_example",  # Use a distinct dir for example
        results_base_dir="data/task3_results_parametric_example",
        plot_z_plane_each_run=True,
        hide_fdtd_gui=False,  # Show GUI for this example run, set to True for sweeps
    )

    print(
        f"Example simulation run '{example_run_id}' complete. Check output directories."
    )

    # To run a sweep, you would loop through parameter combinations and call run_simulation for each.
    # Example (pseudo-code for a sweep):
    # separations = [0.1*u, 0.15*u, 0.2*u]
    # coupling_lengths = [10*u, 20*u, 30*u]
    # for sep_val in separations:
    #     for len_val in coupling_lengths:
    #         current_params = default_params.copy()
    #         current_params["separation"] = sep_val
    #         current_params["coupling_length"] = len_val
    #         run_id = f"sep{sep_val/u:.2f}_L{len_val/u:.0f}" # ... add other params to id
    #         run_simulation(params=current_params, run_id=run_id, hide_fdtd_gui=True)

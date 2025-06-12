"""sb4: Command-line interface for Lumerical FDTD simulations."""

import rich.pretty
import typer
import pathlib
import json
import logging
from typing_extensions import Annotated
from typing import Optional, List, Dict, Any, Union, Tuple
import questionary
from datetime import datetime
import numpy as np
import itertools
import math  # Added math
import matplotlib.pyplot as plt  # Added matplotlib

from sb4.wrapper import lumapi
from sb4.simulation import run_simulation, u
from sb4.results import process_and_save_results, plot_plane_parametric

logger = logging.getLogger(__name__)

app = typer.Typer(
    name="sb4",
    help="A CLI for running and analyzing Lumerical FDTD simulations for integrated photonics.",
    add_completion=False,
    no_args_is_help=True,
)

DEFAULT_OUTPUT_DIR = pathlib.Path("./data/task3")
SWEEP_MANIFEST_PATH = DEFAULT_OUTPUT_DIR / "sweep_manifest.json"
SWEEP_PLOTS_DIR = DEFAULT_OUTPUT_DIR / "sweep_plots"


# Ensure default directories exist
DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
SWEEP_PLOTS_DIR.mkdir(parents=True, exist_ok=True)


def _are_params_close(params1: Dict[str, float], params2: Dict[str, float], rel_tol: float = 1e-9, abs_tol: float = 1e-12) -> bool:
    """Compares two parameter dictionaries for closeness.
    Considers only keys present in params1.
    """
    if not params1 or not params2:
        return False
    # Check only keys that define a unique simulation instance from run_simulation perspective
    # These are the 11 parameters passed to run_simulation
    relevant_keys = [
        "wg1_width",
        "wg2_width",
        "separation",
        "coupling_length",
        "wg_z_span",
        "fan_out_y_offset",
        "sbend_x_extent",
        "center_wavelength",
        "monitor_offset_from_sbend",
        "fdtd_xy_padding",
        "wg_extension_past_fdtd_edge",
    ]
    for key in relevant_keys:
        val1 = params1.get(key)
        val2 = params2.get(key)
        if val1 is None or val2 is None:  # Should not happen if data is consistent
            return False
        if not math.isclose(val1, val2, rel_tol=rel_tol, abs_tol=abs_tol):
            return False
    return True


def _parse_param_input(input_str: str, param_name: str) -> List[float]:
    """Parses a parameter input string.
    Returns a list of float values.
    If input_str is a float, returns a list with that one float.
    If input_str is 'start:end:steps', returns a list from np.linspace.
    """
    try:
        return [float(input_str)]
    except ValueError:
        parts = input_str.split(":")
        if len(parts) == 3:
            try:
                start = float(parts[0])
                end = float(parts[1])
                steps = int(parts[2])
                if steps < 1:
                    raise ValueError("Number of steps must be at least 1.")
                # np.linspace(start, end, 1) correctly gives [start]
                return list(np.linspace(start, end, steps))
            except ValueError as e:
                raise typer.BadParameter(f"Invalid sweep format for '{param_name}'. Expected float or 'start:end:steps' (e.g., '0.4:0.6:3'). Error: {e}")
        else:
            raise typer.BadParameter(f"Invalid format for '{param_name}'. Expected float or 'start:end:steps'.")


@app.command()
@app.command("r", hidden=True)
def run(
    wg1_width: Annotated[float, typer.Option(help="Width of the top waveguide (microns).")] = 0.5,
    wg2_width: Annotated[float, typer.Option(help="Width of the bottom waveguide (microns).")] = 0.5,
    separation: Annotated[float, typer.Option(help="Edge-to-edge separation (microns).")] = 0.15,
    coupling_length: Annotated[float, typer.Option(help="Length of coupling section (microns).")] = 10.0,
    # Fixed parameters (using their previous default values):
    # wg_z_span: float = 0.22,
    # fan_out_y_offset: float = 5.0,
    # sbend_x_extent: float = 10.0,
    center_wavelength: Annotated[float, typer.Option(help="Center wavelength for the simulation (microns).")] = 1.55,
    # monitor_offset_from_sbend: float = 2.5,
    # fdtd_xy_padding: float = 2.5,
    # wg_extension_past_fdtd_edge: float = 1.0,
    output_dir: Annotated[
        pathlib.Path,
        typer.Option(help="Base directory for simulation outputs (layout-specific folders will be created here)."),
    ] = DEFAULT_OUTPUT_DIR,
    plot_z_plane: Annotated[bool, typer.Option(help="Plot Z-plane intensity after simulation.")] = True,
    show_gui: Annotated[bool, typer.Option(help="Show Lumerical FDTD CAD window during simulation.")] = True,
):
    """
    Run a Lumerical FDTD simulation with specified parameters for waveguide dimensions and coupling.
    Other geometric and simulation parameters are fixed to default values.
    All length parameters are in MICRONS.
    Outputs (.fsp, results.json, plots) are saved in a unique directory under `output_dir`.
    A new timestamped folder (sim_HHMM_DDMMYY) is created for each run.

    User-configurable parameters:
    - wg1_width: Width of the top waveguide (microns).
    - wg2_width: Width of the bottom waveguide (microns).
    - separation: Edge-to-edge separation (microns).
    - coupling_length: Length of coupling section (microns).
    - center_wavelength: Center wavelength for the simulation (microns).

    Fixed parameters (defaults in microns):
    - wg_z_span: 0.22
    - fan_out_y_offset: 5.0
    - sbend_x_extent: 10.0
    - monitor_offset_from_sbend: 2.5
    - fdtd_xy_padding: 2.5
    - wg_extension_past_fdtd_edge: 1.0
    """
    params_microns = {
        "wg1_width": wg1_width,
        "wg2_width": wg2_width,
        "separation": separation,
        "coupling_length": coupling_length,
        "wg_z_span": 0.22,  # Fixed
        "fan_out_y_offset": 5.0,  # Fixed
        "sbend_x_extent": 10.0,  # Fixed
        "center_wavelength": center_wavelength,
        "monitor_offset_from_sbend": 2.5,  # Fixed
        "fdtd_xy_padding": 2.5,  # Fixed
        "wg_extension_past_fdtd_edge": 1.0,  # Fixed
    }
    logger.info("Running simulation with parameters (microns):")
    rich.pretty.pprint(params_microns, max_length=80, expand_all=True)

    # run_simulation now expects parameters in microns for storage and comparison,
    # but will handle conversion to meters internally for Lumerical.
    logger.info(f"Output base directory: {output_dir}")

    run_simulation(params_microns=params_microns, output_base_dir=str(output_dir), plot_z_plane_each_run=plot_z_plane, hide_fdtd_gui=not show_gui)
    logger.info(f"Simulation and processing finished. Check output in: {output_dir}")


@app.command()
@app.command("s", hidden=True)
def sweep(
    wg1_width_str: Annotated[str, typer.Option(help="Top WG width (microns) or 'start:end:steps'.")] = "0.5",
    wg2_width_str: Annotated[str, typer.Option(help="Bottom WG width (microns) or 'start:end:steps'.")] = "0.5",
    separation_str: Annotated[str, typer.Option(help="Edge-to-edge separation (microns) or 'start:end:steps'.")] = "0.15",
    coupling_length_str: Annotated[str, typer.Option(help="Coupling length (microns) or 'start:end:steps'.")] = "10.0",
    center_wavelength_str: Annotated[str, typer.Option(help="Center wavelength (microns) or 'start:end:steps'.")] = "1.55",
    output_dir: Annotated[
        pathlib.Path,
        typer.Option(help="Base directory for simulation outputs."),
    ] = DEFAULT_OUTPUT_DIR,
    plot_z_plane: Annotated[bool, typer.Option(help="Plot Z-plane intensity after each simulation run.")] = True,
    show_gui: Annotated[bool, typer.Option(help="Show Lumerical FDTD CAD window during simulation runs.")] = True,
    match_widths: Annotated[bool, typer.Option(help="If true the waveguide widths will be matched (wg1_width = wg2_width). You must ensure you pass the same parameters for both waveguides.")] = False,
    sweep_description: Annotated[Optional[str], typer.Option(help="Optional description for this sweep.")] = None,
):
    """
    Run a sweep of Lumerical FDTD simulations by varying one or more user-configurable parameters.
    Each parameter can be a fixed value or a sweep range defined as 'start:end:steps'.
    This sweep definition will be recorded in data/task3/sweep_manifest.json.

    User-configurable parameters (provide as float or 'start:end:steps'):
    - wg1_width_str: Top waveguide width (microns).
    - wg2_width_str: Bottom waveguide width (microns).
    - separation_str: Edge-to-edge separation (microns).
    - coupling_length_str: Length of coupling section (microns).
    - center_wavelength_str: Center wavelength for the simulation (microns).

    Fixed parameters (defaults in microns):
    - wg_z_span: 0.22
    - fan_out_y_offset: 5.0
    - sbend_x_extent: 10.0
    - monitor_offset_from_sbend: 2.5
    - fdtd_xy_padding: 2.5
    - wg_extension_past_fdtd_edge: 1.0

    Each combination of swept parameters results in an independent simulation run,
    saved in a unique timestamped folder (sim_HHMM_DDMMYY) under `output_dir`.
    """
    if match_widths and wg1_width_str != wg2_width_str:
        logger.warning("match_widths is enabled but wg1_width_str and wg2_width_str are not the same. wg2_width_str will be set to wg1_width_str.")
        wg2_width_str = wg1_width_str

    current_sweep_config = {
        "params_input": {
            "wg1_width_str": wg1_width_str,
            "wg2_width_str": wg2_width_str,
            "separation_str": separation_str,
            "coupling_length_str": coupling_length_str,
            "center_wavelength_str": center_wavelength_str,
        },
        "match_widths": match_widths,
    }

    # Generate description if not provided
    if not sweep_description:
        swept_params_desc = []
        if len(_parse_param_input(wg1_width_str, "wg1")) > 1:
            swept_params_desc.append("wg1_w" if not match_widths else "matched_w")
        if not match_widths and len(_parse_param_input(wg2_width_str, "wg2")) > 1:
            swept_params_desc.append("wg2_w")
        if len(_parse_param_input(separation_str, "sep")) > 1:
            swept_params_desc.append("sep")
        if len(_parse_param_input(coupling_length_str, "cl")) > 1:
            swept_params_desc.append("cl")
        if len(_parse_param_input(center_wavelength_str, "cwl")) > 1:
            swept_params_desc.append("cwl")
        current_sweep_config["description"] = f"Sweep {', '.join(swept_params_desc) if swept_params_desc else 'single run'}"
    else:
        current_sweep_config["description"] = sweep_description

    # Load existing manifest or initialize
    manifest_data = []
    if SWEEP_MANIFEST_PATH.exists():
        try:
            with open(SWEEP_MANIFEST_PATH, "r", encoding="utf-8") as f:
                manifest_data = json.load(f)
            if not isinstance(manifest_data, list):  # Ensure it's a list
                logger.warning(f"{SWEEP_MANIFEST_PATH} does not contain a list. Initializing new manifest.")
                manifest_data = []
        except json.JSONDecodeError:
            logger.warning(f"Could not decode {SWEEP_MANIFEST_PATH}. Initializing new manifest.")
            manifest_data = []

    # Check for duplicates before adding
    is_duplicate = False
    for existing_sweep in manifest_data:
        # Compare based on the defining string parameters and match_widths
        if all(existing_sweep.get(k) == current_sweep_config[k] for k in ["wg1_width_str", "wg2_width_str", "separation_str", "coupling_length_str", "center_wavelength_str", "match_widths"]):
            is_duplicate = True
            logger.info(
                f"This sweep configuration already exists in {SWEEP_MANIFEST_PATH} (Index {manifest_data.index(existing_sweep)}: '{existing_sweep.get('description')}'). It will not be added again."
            )
            # Optionally, update description if a new one is provided? For now, no.
            break

    if not is_duplicate:
        manifest_data.append(current_sweep_config)
        try:
            with open(SWEEP_MANIFEST_PATH, "w", encoding="utf-8") as f:
                json.dump(manifest_data, f, indent=2)
            logger.info(f"Sweep definition '{current_sweep_config['description']}' saved to {SWEEP_MANIFEST_PATH}")
        except IOError as e:
            logger.error(f"Could not write to sweep manifest {SWEEP_MANIFEST_PATH}: {e}")

    param_inputs = {
        "wg1_width": wg1_width_str,
        "wg2_width": wg2_width_str,
        "separation": separation_str,
        "coupling_length": coupling_length_str,
        "center_wavelength": center_wavelength_str,
    }

    parsed_param_values = {}
    sweep_param_names = []

    for name, input_str in param_inputs.items():
        values = _parse_param_input(input_str, name)
        parsed_param_values[name] = values
        sweep_param_names.append(name)

    value_lists_for_product = [parsed_param_values[name] for name in sweep_param_names]
    all_combinations = list(itertools.product(*value_lists_for_product))

    # Filter combinations if match_widths is True *before* calculating total_runs
    if match_widths:
        all_combinations = [combo for combo in all_combinations if math.isclose(combo[sweep_param_names.index("wg1_width")], combo[sweep_param_names.index("wg2_width")])]

    total_runs = len(all_combinations)

    if total_runs == 0:
        logger.info("No simulation runs to perform based on sweep parameters (possibly due to match_widths filter).")
        return

    logger.info(f"Starting sweep '{current_sweep_config['description']}' with {total_runs} total simulation run(s).")

    fixed_params = {
        "wg_z_span": 0.22,
        "fan_out_y_offset": 5.0,
        "sbend_x_extent": 10.0,
        "monitor_offset_from_sbend": 2.5,
        "fdtd_xy_padding": 2.5,
        "wg_extension_past_fdtd_edge": 1.0,
    }

    actual_runs_performed = 0
    for i, combo_values in enumerate(all_combinations):
        current_params_microns = dict(zip(sweep_param_names, combo_values))
        # This check is now implicitly handled by the pre-filtering of all_combinations
        # if match_widths:
        #     if not math.isclose(current_params_microns["wg1_width"], current_params_microns["wg2_width"]):
        #         # This should not be reached if pre-filtering is correct
        #         logger.debug(f"Skipping combo due to match_widths: {current_params_microns}")
        #         continue

        current_params_microns.update(fixed_params)

        logger.info(f"--- Running Sweep Iteration {actual_runs_performed + 1}/{total_runs} ---")
        logger.info("Parameters for this run (microns):")
        rich.pretty.pprint(current_params_microns)

        try:
            run_simulation(params_microns=current_params_microns, output_base_dir=str(output_dir), plot_z_plane_each_run=plot_z_plane, hide_fdtd_gui=not show_gui)
            logger.info(f"--- Completed Sweep Iteration {actual_runs_performed + 1}/{total_runs} ---")
            actual_runs_performed += 1
        except Exception as e:
            logger.error(f"Error during sweep iteration {actual_runs_performed + 1}/{total_runs} with params: {current_params_microns}")
            logger.error(f"Error details: {e}")
            # typer.echo(f"Error in sweep run {i+1}/{total_runs}. Check logs. Continuing with next run if any.", err=True) # Replaced
            logger.error(f"Error in sweep run {actual_runs_performed + 1}/{total_runs}. Check logs. Continuing with next run if any.")

    logger.info(f"Sweep '{current_sweep_config['description']}' finished. {actual_runs_performed} simulations were attempted/completed.")


@app.command("analyse")
@app.command("a", hidden=True)
def analyse_results(
    output_dir: Annotated[
        pathlib.Path,
        typer.Option(help="Base directory containing layout-specific simulation output folders."),
    ] = DEFAULT_OUTPUT_DIR,
    plot_z_plane: Annotated[
        bool,
        typer.Option(help="Generate/regenerate Z-plane intensity plot during analysis."),
    ] = True,
):
    """
    Analyse existing Lumerical FDTD simulation results.
    Scans the `output_dir` for layout folders (identified by having a `results.json` file).
    Presents a list of available layouts for selection, then displays their results.
    """
    logger.info(f"Scanning for simulation run folders with results.json in: {output_dir}")

    layout_folders: List[pathlib.Path] = []
    for item in output_dir.iterdir():
        if item.is_dir() and (item / "results.json").exists():
            layout_folders.append(item)

    if not layout_folders:
        logger.warning(f"No simulation run folders with results.json found in {output_dir}.")
        # typer.echo(f"No simulation run folders with results.json found in {output_dir}.") # Replaced
        raise typer.Exit()

    layout_folders.sort(key=lambda x: (x / "results.json").stat().st_mtime, reverse=True)

    choices = []
    for folder in layout_folders:
        results_mod_time = datetime.fromtimestamp((folder / "results.json").stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")
        choices.append(questionary.Choice(title=f"{folder.name} (Results: {results_mod_time})", value=str(folder)))

    selected_layout_folder_str = questionary.select(
        "Select a layout folder to analyse:",
        choices=choices,
        use_shortcuts=True,
    ).ask()

    if not selected_layout_folder_str:
        logger.info("No layout folder selected for analysis.")
        raise typer.Exit()

    selected_layout_folder = pathlib.Path(selected_layout_folder_str)
    layout_id = selected_layout_folder.name
    sim_fsp_path = selected_layout_folder / "simulation.fsp"
    results_json_path = selected_layout_folder / "results.json"

    logger.info(f"Analysing layout: {layout_id} in folder: {selected_layout_folder}")

    if not results_json_path.exists():
        logger.error(f"Critical: results.json not found in {selected_layout_folder} after selection. This should not happen.")
        # typer.echo(f"Error: results.json not found in {selected_layout_folder}.") # Replaced
        raise typer.Exit(code=1)

    with open(results_json_path, "r", encoding="utf-8") as f:
        results_data = json.load(f)

    logger.info(f"\\n--- Results for Run: {layout_id} ---")
    rich.pretty.pprint(results_data, max_length=80, expand_all=True)

    params_microns_from_results = results_data.get("parameters")

    if plot_z_plane:
        if params_microns_from_results and sim_fsp_path.exists():
            logger.info(f"Attempting to generate Z-plane plot for {layout_id} using {sim_fsp_path}")
            expected_plot_path = selected_layout_folder / "z_plane_intensity.png"
            if expected_plot_path.exists():
                logger.info(f"Z-plane plot should be available at: {expected_plot_path}")
            else:
                logger.info(f"Z-plane plot was not generated during the initial run or not found at: {expected_plot_path}")
                logger.info("To generate plots, ensure 'plot_z_plane' is enabled during the 'run' command or re-process.")
        elif not sim_fsp_path.exists():
            logger.warning(f"Simulation file {sim_fsp_path} not found. Cannot generate Z-plane plot.")
        else:
            logger.warning("Parameters not found in results.json, cannot reliably generate Z-plane plot during analysis.")

    logger.info(f"--- End of Analysis for {layout_id} ---")


@app.command(name="plot-sweep")
@app.command("ps", hidden=True)
def plot_sweep(
    output_dir: Annotated[pathlib.Path, typer.Option(help="Base directory containing simulation outputs.")] = DEFAULT_OUTPUT_DIR,
):
    """
    Plot results from a previously run sweep defined in sweep_manifest.json.
    Allows selection of a sweep and a result attribute (T_net_br or T_net_tr) to plot.
    Generates a 1D line plot or a 2D heatmap based on sweep dimensionality.
    Plots are saved in data/task3/sweep_plots/.
    """
    if not SWEEP_MANIFEST_PATH.exists():
        logger.error(f"Sweep manifest file not found: {SWEEP_MANIFEST_PATH}")
        raise typer.Exit(code=1)

    try:
        with open(SWEEP_MANIFEST_PATH, "r", encoding="utf-8") as f:
            sweep_manifest: List[Dict] = json.load(f)
        if not sweep_manifest or not isinstance(sweep_manifest, list):
            logger.error(f"Sweep manifest {SWEEP_MANIFEST_PATH} is empty or not a list.")
            raise typer.Exit(code=1)
    except json.JSONDecodeError:
        logger.error(f"Could not decode sweep manifest {SWEEP_MANIFEST_PATH}.")
        raise typer.Exit(code=1)
    except IOError as e:
        logger.error(f"Could not read sweep manifest {SWEEP_MANIFEST_PATH}: {e}")
        raise typer.Exit(code=1)

    # Select sweep
    sweep_choices = [questionary.Choice(title=f"[{i}] {s.get('description', 'No description')} (match_widths={s.get('match_widths')})", value=i) for i, s in enumerate(sweep_manifest)]
    if not sweep_choices:
        logger.error("No sweeps found in the manifest.")
        raise typer.Exit(code=1)

    selected_sweep_idx = questionary.select("Select a sweep to plot:", choices=sweep_choices).ask()

    if selected_sweep_idx is None:
        logger.info("No sweep selected.")
        raise typer.Exit()

    selected_sweep_config = sweep_manifest[selected_sweep_idx]
    logger.info(f"Selected sweep: {selected_sweep_config.get('description')}")

    # Select result attribute
    available_attributes = ["T_net_tr", "T_net_br"]
    result_attribute = questionary.select("Select result attribute to plot:", choices=available_attributes).ask()

    if not result_attribute:
        logger.info("No result attribute selected.")
        raise typer.Exit()

    # --- Gather results for the selected sweep ---
    params = selected_sweep_config["params_input"]
    param_strings_for_sweep = {
        "wg1_width": params["wg1_width"],
        "wg2_width": params["wg2_width"],
        "separation": params["separation"],
        "coupling_length": params["coupling_length"],
        "center_wavelength": params["center_wavelength"],
    }

    parsed_params_for_sweep = {name: _parse_param_input(val, name) for name, val in param_strings_for_sweep.items()}
    sweep_param_names_ordered = list(param_strings_for_sweep.keys())  # Keep consistent order

    value_lists_for_product = [parsed_params_for_sweep[name] for name in sweep_param_names_ordered]
    all_expected_combos_values = list(itertools.product(*value_lists_for_product))

    fixed_sim_params = {
        "wg_z_span": 0.22,
        "fan_out_y_offset": 5.0,
        "sbend_x_extent": 10.0,
        "monitor_offset_from_sbend": 2.5,
        "fdtd_xy_padding": 2.5,
        "wg_extension_past_fdtd_edge": 1.0,
    }

    all_expected_param_dicts_for_sweep = []
    for combo_vals in all_expected_combos_values:
        p_dict = dict(zip(sweep_param_names_ordered, combo_vals))
        if selected_sweep_config["match_widths"]:
            if not math.isclose(p_dict["wg1_width"], p_dict["wg2_width"]):
                continue  # Skip if widths don't match for a match_widths sweep
        p_dict.update(fixed_sim_params)
        all_expected_param_dicts_for_sweep.append(p_dict)

    gathered_results_data = []
    for sim_run_dir in output_dir.iterdir():
        if sim_run_dir.is_dir() and (sim_run_dir / "results.json").exists():
            try:
                with (sim_run_dir / "results.json").open("r", encoding="utf-8") as f_res:
                    run_result_data = json.load(f_res)
                actual_params = run_result_data.get("parameters")
                if actual_params:
                    for expected_params in all_expected_param_dicts_for_sweep:
                        if _are_params_close(expected_params, actual_params):
                            gathered_results_data.append(run_result_data)
                            break  # Found match for this run_result_data
            except Exception as e:
                logger.warning(f"Could not process results.json in {sim_run_dir}: {e}")

    if not gathered_results_data:
        logger.error(f"No simulation results found for the selected sweep '{selected_sweep_config.get('description')}'. Ensure simulations have been run.")
        raise typer.Exit(code=1)

    logger.info(f"Found {len(gathered_results_data)} matching simulation results for the sweep.")

    # --- Determine sweep dimensionality and parameters ---
    swept_param_values = {}  # Stores actual values found in results
    param_names_map = {  # Maps from config string name to actual param name if different
        "wg1_width_str": "wg1_width",
        "wg2_width_str": "wg2_width",
        "separation_str": "separation",
        "coupling_length_str": "coupling_length",
        "center_wavelength_str": "center_wavelength",
    }

    for p_str_name, p_actual_name in param_names_map.items():
        # Use parsed_params_for_sweep to know the original intent of sweep range
        original_values_for_param = parsed_params_for_sweep[p_actual_name]
        if len(original_values_for_param) > 1:
            # Collect actual values from results to handle float precision issues for axis ticks
            actual_vals_from_results = sorted(list(set(res["parameters"][p_actual_name] for res in gathered_results_data if p_actual_name in res["parameters"])))
            # Filter actual_vals_from_results to be close to original_values_for_param
            # This is complex due to floating point. For now, trust the gathered_results_data are correct.
            # Use the unique values from the gathered results for this parameter.
            if actual_vals_from_results:
                swept_param_values[p_actual_name] = actual_vals_from_results

    if selected_sweep_config["match_widths"] and "wg1_width" in swept_param_values:
        # swept_param_values["matched_wg_width"] = swept_param_values.pop("wg1_width")
        # if "wg2_width" in swept_param_values:  # Should be identical if logic is correct
        #     swept_param_values.pop("wg2_width")
        swept_param_values.pop("wg2_width")
        assert "wg1_width" in swept_param_values, "wg1_width should be present if match_widths is True"

    active_swept_dims_names = list(swept_param_values.keys())
    num_dims = len(active_swept_dims_names)

    logger.info(f"Sweep dimensionality: {num_dims}D. Swept parameters: {active_swept_dims_names}")

    # --- Plotting ---
    plt.style.use("seaborn-v0_8-darkgrid")  # Using a seaborn style

    plot_filename_base = f"sweep_idx_{selected_sweep_idx}_{selected_sweep_config.get('description','').replace(' ','_').lower()}_{result_attribute}"

    if num_dims == 0:  # Should ideally be caught by len(original_values_for_param) > 1
        logger.info("Selected sweep has no varying parameters. Plotting the single data point.")
        if gathered_results_data:
            single_res = gathered_results_data[0]
            val = single_res.get(result_attribute, float("nan"))
            logger.info(f"Single point: {result_attribute} = {val:.4f}")
            # No plot to generate, or could plot a single point if desired.
        else:
            logger.error("No data found for single point 'sweep'.")
        raise typer.Exit()

    elif num_dims == 1:
        param_name_x = active_swept_dims_names[0]
        x_values = swept_param_values[param_name_x]  # Already sorted unique values

        # Create a mapping from x_value to result_attribute
        # Handle cases where multiple runs might exist for one x_value (should not happen with current setup)
        # or if some results are missing.
        plot_points = []
        for x_val in x_values:
            # Find corresponding results
            # Need to be careful with float comparisons again
            y_val_for_x = [res[result_attribute] for res in gathered_results_data if math.isclose(res["parameters"][param_name_x], x_val) and result_attribute in res]
            if y_val_for_x:  # Take the first if multiple (should be one)
                plot_points.append((x_val, y_val_for_x[0]))
            else:  # Should not happen if gathered_results_data is complete for the sweep points
                plot_points.append((x_val, np.nan))

        plot_points.sort(key=lambda item: item[0])  # Ensure sorted by x_value
        final_x_coords, final_y_coords = zip(*plot_points)

        plt.figure(figsize=(10, 6))
        plt.plot(final_x_coords, final_y_coords, marker="o", linestyle="-")
        plt.xlabel(f"{param_name_x} (microns)")
        plt.ylabel(result_attribute)
        plt.title(f"{result_attribute} vs {param_name_x}")
        plt.grid(True)
        plot_filepath = SWEEP_PLOTS_DIR / f"{plot_filename_base}_line.png"
        plt.savefig(str(plot_filepath).replace(":", "_"))  # Replace ':' in filename to avoid issues
        logger.info(f"1D line plot saved to: {plot_filepath}")
        plt.close()

    elif num_dims == 2:
        param_name_x = active_swept_dims_names[0]
        param_name_y = active_swept_dims_names[1]

        x_coords = swept_param_values[param_name_x]  # Sorted unique
        y_coords = swept_param_values[param_name_y]  # Sorted unique

        Z_data = np.full((len(y_coords), len(x_coords)), np.nan)

        for res_data_point in gathered_results_data:
            params_of_point = res_data_point["parameters"]
            attr_val = res_data_point.get(result_attribute)
            if attr_val is None:
                continue

            # Find indices, robustly
            try:
                # Find the closest x_coord and y_coord
                x_val = params_of_point[param_name_x]
                y_val = params_of_point[param_name_y]

                # Find index of x_val in x_coords (list of unique sorted values)
                ix = -1
                for idx_search, xc in enumerate(x_coords):
                    if math.isclose(xc, x_val):
                        ix = idx_search
                        break

                iy = -1
                for idx_search, yc in enumerate(y_coords):
                    if math.isclose(yc, y_val):
                        iy = idx_search
                        break

                if ix != -1 and iy != -1:
                    Z_data[iy, ix] = attr_val
                else:
                    logger.debug(f"Could not map point {params_of_point} to grid for heatmap.")
            except (ValueError, KeyError) as e:  # Should not happen if data is clean
                logger.warning(f"Skipping data point for heatmap due to mapping error: {e} for point {params_of_point}")

        plt.figure(figsize=(10, 8))
        # pcolormesh expects X, Y to define cell boundaries if Z is (ny, nx)
        # If X, Y are cell centers, shading='auto' or 'gouraud' can work if Z is (ny, nx)
        # Let's use imshow which is simpler for grid data
        # We need to ensure x_coords and y_coords are monotonic for extent

        # For pcolormesh with X, Y as centers, Z should be (len(y_coords), len(x_coords))
        # X_mesh, Y_mesh = np.meshgrid(x_coords, y_coords) # X_mesh (ny,nx), Y_mesh (ny,nx)
        # plt.pcolormesh(X_mesh, Y_mesh, Z_data, shading='auto', cmap='viridis')

        # Simpler: imshow if coordinates are regularly spaced (or can be treated as such for visualization)
        # extent: left, right, bottom, top
        # We need to handle non-uniform spacing if using imshow with extent.
        # pcolormesh is better for potentially non-uniform grids.
        # If x_coords and y_coords are the coordinates of the *centers* of the cells.
        # Z_data[iy, ix] means Z at (y_coords[iy], x_coords[ix])

        # To use pcolormesh correctly when x_coords, y_coords are centers,
        # we might need to convert them to edges or ensure Z_data matches.
        # Let's try with x_coords, y_coords as direct inputs for axes if they are monotonic.
        # The Z_data is (num_y_points, num_x_points).

        # Create meshgrid for pcolormesh
        X_mesh, Y_mesh = np.meshgrid(x_coords, y_coords)  # X_mesh (len(y_coords), len(x_coords)), Y_mesh also

        if np.all(np.isnan(Z_data)):
            logger.error("All Z_data for heatmap is NaN. Cannot plot.")
            raise typer.Exit(code=1)

        pcm = plt.pcolormesh(X_mesh, Y_mesh, Z_data, shading="gouraud", cmap="viridis", vmin=np.nanmin(Z_data), vmax=np.nanmax(Z_data))
        plt.colorbar(pcm, label=result_attribute)
        plt.xlabel(f"{param_name_x} (microns)")
        plt.ylabel(f"{param_name_y} (microns)")
        plt.title(f"{result_attribute} vs {param_name_x} & {param_name_y}")

        plot_filepath = SWEEP_PLOTS_DIR / f"{plot_filename_base}_heatmap.png"
        plt.savefig(str(plot_filepath).replace(":", "_"))
        logger.info(f"2D heatmap plot saved to: {plot_filepath}")
        plt.close()

    else:  # num_dims > 2
        logger.error(f"Plotting for {num_dims}D sweeps is not implemented.")
        raise NotImplementedError(f"Plotting for {num_dims}D sweeps is not implemented.")


if __name__ == "__main__":
    app()

# python -m sb4 s --coupling-length-str 0:30:15 --center-w

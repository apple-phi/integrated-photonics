"""sb4: Command-line interface for Lumerical FDTD simulations."""

import rich.pretty
import typer
import pathlib
import json
import logging
from typing_extensions import Annotated
from typing import Optional, List, Dict, Any, Union, Tuple  # Added Union, Tuple
import questionary
from datetime import datetime
import numpy as np  # Added numpy
import itertools  # Added itertools

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

DEFAULT_OUTPUT_DIR = pathlib.Path("./data/task3")  # Unified output directory

# Ensure default directory exists
DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


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
):
    """
    Run a sweep of Lumerical FDTD simulations by varying one or more user-configurable parameters.
    Each parameter can be a fixed value or a sweep range defined as 'start:end:steps'.

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
    param_inputs = {
        "wg1_width": wg1_width_str,
        "wg2_width": wg2_width_str,
        "separation": separation_str,
        "coupling_length": coupling_length_str,
        "center_wavelength": center_wavelength_str,
    }

    parsed_param_values = {}
    sweep_param_names = []  # Keep order for itertools.product

    for name, input_str in param_inputs.items():
        values = _parse_param_input(input_str, name)
        parsed_param_values[name] = values
        sweep_param_names.append(name)

    value_lists_for_product = [parsed_param_values[name] for name in sweep_param_names]

    all_combinations = list(itertools.product(*value_lists_for_product))
    total_runs = len(all_combinations)

    if total_runs == 0:
        logger.info("No simulation runs to perform based on sweep parameters.")
        return

    logger.info(f"Starting sweep with {total_runs} total simulation run(s).")
    typer.echo(f"Starting sweep with {total_runs} total simulation run(s).")

    fixed_params = {
        "wg_z_span": 0.22,
        "fan_out_y_offset": 5.0,
        "sbend_x_extent": 10.0,
        "monitor_offset_from_sbend": 2.5,
        "fdtd_xy_padding": 2.5,
        "wg_extension_past_fdtd_edge": 1.0,
    }

    for i, combo_values in enumerate(all_combinations):
        current_params_microns = dict(zip(sweep_param_names, combo_values))
        current_params_microns.update(fixed_params)  # Add fixed parameters

        logger.info(f"--- Starting Sweep Run {i+1}/{total_runs} ---")
        logger.info(f"--- Running Sweep Iteration {i+1}/{total_runs} ---")
        logger.info("Parameters for this run (microns):")
        rich.pretty.pprint(current_params_microns)

        try:
            run_simulation(params_microns=current_params_microns, output_base_dir=str(output_dir), plot_z_plane_each_run=plot_z_plane, hide_fdtd_gui=not show_gui)
            logger.info(f"--- Completed Sweep Run {i+1}/{total_runs} ---")
        except Exception as e:
            logger.error(f"Error during sweep run {i+1}/{total_runs} with params: {current_params_microns}")
            logger.error(f"Error details: {e}")
            typer.echo(f"Error in sweep run {i+1}/{total_runs}. Check logs. Continuing with next run if any.", err=True)

    logger.info("Sweep finished.")


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
        if item.is_dir() and (item / "results.json").exists():  # Check for results.json to identify a valid sim folder
            layout_folders.append(item)

    if not layout_folders:
        logger.warning(f"No simulation run folders with results.json found in {output_dir}.")
        typer.echo(f"No simulation run folders with results.json found in {output_dir}.")
        raise typer.Exit()

    # Sort folders by modification time of their results.json, most recent first
    layout_folders.sort(key=lambda x: (x / "results.json").stat().st_mtime, reverse=True)

    choices = []
    for folder in layout_folders:
        # Display folder name (timestamp) and modification time for clarity
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
    # The layout_id is now the timestamped folder name.
    # Files inside will have fixed names e.g. simulation.fsp, results.json
    layout_id = selected_layout_folder.name
    sim_fsp_path = selected_layout_folder / "simulation.fsp"  # Adjusted to fixed filename
    results_json_path = selected_layout_folder / "results.json"

    logger.info(f"Analysing layout: {layout_id} in folder: {selected_layout_folder}")

    if not results_json_path.exists():
        logger.error(f"Critical: results.json not found in {selected_layout_folder} after selection. This should not happen.")
        typer.echo(f"Error: results.json not found in {selected_layout_folder}.")
        raise typer.Exit(code=1)

    with open(results_json_path, "r", encoding="utf-8") as f:
        results_data = json.load(f)

    typer.echo(f"\n--- Results for Run: {layout_id} ---")
    rich.pretty.pprint(results_data, max_length=80, expand_all=True)

    # Parameters are stored in results_data["parameters"] in units of microns
    params_microns_from_results = results_data.get("parameters")

    if plot_z_plane:
        if params_microns_from_results and sim_fsp_path.exists():
            # Need to pass the FDTD object or path to plot_plane_parametric
            # For simplicity, let's assume plot_plane_parametric can take the sim_fsp_path
            # and handle loading if necessary, or it's called from a context where fdtd object is available.
            # The original plot_plane_parametric in results.py takes an fdtd_obj.
            # This might require a change in how plotting is invoked here or in results.py
            # For now, let's assume we might need to load the simulation to plot.
            # This part of analyse_results might need further refinement based on how plotting is handled.
            logger.info(f"Attempting to generate Z-plane plot for {layout_id} using {sim_fsp_path}")
            # This is a placeholder for how plotting would be re-invoked.
            # The original `plot_plane_parametric` is in `results.py` and expects an FDTD object.
            # We might need a helper in `results.py` that can load and plot.
            # For now, just logging. The actual plotting during 'analyse' might be complex
            # if it requires re-running parts of the simulation or loading large files.
            # The original `process_and_save_results` calls `plot_plane_parametric`.
            # If `plot_z_plane` is true during `run`, it's already generated.
            # This option in `analyse` would be for *re-generating* or generating if missed.
            # Let's print a message indicating where the plot would be if generated during the run.
            expected_plot_path = selected_layout_folder / f"z_plane_intensity.png"  # Adjusted to fixed filename
            if expected_plot_path.exists():
                typer.echo(f"Z-plane plot should be available at: {expected_plot_path}")
            else:
                typer.echo(f"Z-plane plot was not generated during the initial run or not found at: {expected_plot_path}")
                typer.echo("To generate plots, ensure 'plot_z_plane' is enabled during the 'run' command.")
        elif not sim_fsp_path.exists():
            logger.warning(f"Simulation file {sim_fsp_path} not found. Cannot generate Z-plane plot.")
            typer.echo(f"Simulation file {sim_fsp_path} not found. Cannot generate Z-plane plot.")
        else:  # params_microns_from_results is None
            logger.warning("Parameters not found in results.json, cannot reliably generate Z-plane plot during analysis.")
            typer.echo("Parameters not found in results.json, cannot reliably generate Z-plane plot during analysis.")

    typer.echo(f"--- End of Analysis for {layout_id} ---")


if __name__ == "__main__":
    app()

"""sb4: Command-line interface for Lumerical FDTD simulations."""

import rich.pretty
import typer
import pathlib
import json
import logging
from typing_extensions import Annotated
from typing import Optional, List, Dict, Any
import questionary
from datetime import datetime

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


@app.command()
@app.command("r", hidden=True)
def run(
    wg1_width: Annotated[float, typer.Option(help="Width of the top waveguide (microns).")] = 0.5,
    wg2_width: Annotated[float, typer.Option(help="Width of the bottom waveguide (microns).")] = 0.5,
    separation: Annotated[float, typer.Option(help="Edge-to-edge separation (microns).")] = 0.15,
    coupling_length: Annotated[float, typer.Option(help="Length of coupling section (microns).")] = 10.0,
    wg_z_span: Annotated[float, typer.Option(help="Thickness of waveguides (microns).")] = 0.22,
    fan_out_y_offset: Annotated[float, typer.Option(help="S-bend fan-out Y-offset (microns).")] = 5.0,
    sbend_x_extent: Annotated[float, typer.Option(help="S-bend X-extent (microns).")] = 10.0,
    center_wavelength: Annotated[float, typer.Option(help="Center wavelength (microns).")] = 1.55,
    monitor_offset_from_sbend: Annotated[float, typer.Option(help="Monitor X-offset from S-bend end (microns).")] = 2.5,
    fdtd_xy_padding: Annotated[float, typer.Option(help="FDTD XY padding (microns).")] = 2.5,
    wg_extension_past_fdtd_edge: Annotated[float, typer.Option(help="WG extension beyond FDTD edge (microns).")] = 1.0,
    output_dir: Annotated[
        pathlib.Path,
        typer.Option(help="Base directory for simulation outputs (layout-specific folders will be created here)."),
    ] = DEFAULT_OUTPUT_DIR,
    plot_z_plane: Annotated[bool, typer.Option(help="Plot Z-plane intensity after simulation.")] = True,
    show_gui: Annotated[bool, typer.Option(help="Show Lumerical FDTD CAD window during simulation.")] = True,
    # Removed run_id, sim_dir, results_dir as they are now handled internally by layout_id and output_dir
):
    """
    Run a Lumerical FDTD simulation with specified parameters.
    All length parameters are in MICRONS.
    Outputs (.fsp, results.json, plots) are saved in a unique directory under `output_dir` for each parameter set.
    A new timestamped folder is created for each run.
    """
    params_microns = {
        "wg1_width": wg1_width,
        "wg2_width": wg2_width,
        "separation": separation,
        "coupling_length": coupling_length,
        "wg_z_span": wg_z_span,
        "fan_out_y_offset": fan_out_y_offset,
        "sbend_x_extent": sbend_x_extent,
        "center_wavelength": center_wavelength,
        "monitor_offset_from_sbend": monitor_offset_from_sbend,
        "fdtd_xy_padding": fdtd_xy_padding,
        "wg_extension_past_fdtd_edge": wg_extension_past_fdtd_edge,
    }
    logger.info("Running simulation with parameters (microns):")
    rich.pretty.pprint(params_microns, max_length=80, expand_all=True)

    params_meters = {key: value * u for key, value in params_microns.items()}

    # layout_id creation is now handled within run_simulation which creates a timestamped folder
    logger.info(f"Output base directory: {output_dir}")

    # run_simulation will now create a unique timestamped directory inside output_dir
    # and will handle its own logging regarding the specific directory created.
    run_simulation(params=params_meters, output_base_dir=str(output_dir), plot_z_plane_each_run=plot_z_plane, hide_fdtd_gui=not show_gui)
    logger.info(f"Simulation and processing finished. Check output in: {output_dir}")


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

    # Parameters are stored in results_data["parameters"] in units of meters
    params_meters_from_results = results_data.get("parameters")

    if plot_z_plane:
        if params_meters_from_results and sim_fsp_path.exists():
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
        else:  # params_meters_from_results is None
            logger.warning("Parameters not found in results.json, cannot reliably generate Z-plane plot during analysis.")
            typer.echo("Parameters not found in results.json, cannot reliably generate Z-plane plot during analysis.")

    typer.echo(f"--- End of Analysis for {layout_id} ---")


if __name__ == "__main__":
    app()

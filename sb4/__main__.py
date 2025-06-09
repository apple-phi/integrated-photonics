"""sb4: Command-line interface for Lumerical FDTD simulations."""

import rich.logging
import typer
import pathlib
import os
import json
import logging
from typing_extensions import Annotated
from typing import Optional, List, Dict, Any  # Added List, Dict, Any
import questionary
import datetime

from sb4.simulation import run_simulation, u, create_layout_id  # Added create_layout_id
from sb4.results import process_and_save_results, plot_plane_parametric  # Added plot_plane_parametric

FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO",
    format=FORMAT,
    datefmt="[%X]",
    handlers=[rich.logging.RichHandler(markup=True)],
)
logger = logging.getLogger(__name__)

app = typer.Typer(
    name="sb4",
    help="A CLI for running and analyzing Lumerical FDTD simulations for integrated photonics.",
    add_completion=False,
    no_args_is_help=True,
)

DEFAULT_OUTPUT_DIR = pathlib.Path("./task3")  # Unified output directory

# Ensure default directory exists
DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


@app.command()
@app.command("r", hidden=True)
def run(
    wg1_width: Annotated[
        float, typer.Option(help="Width of the top waveguide (microns).")
    ] = 0.5,
    wg2_width: Annotated[
        float, typer.Option(help="Width of the bottom waveguide (microns).")
    ] = 0.5,
    separation: Annotated[
        float, typer.Option(help="Edge-to-edge separation (microns).")
    ] = 0.15,
    coupling_length: Annotated[
        float, typer.Option(help="Length of coupling section (microns).")
    ] = 10.0,
    wg_z_span: Annotated[
        float, typer.Option(help="Thickness of waveguides (microns).")
    ] = 0.22,
    fan_out_y_offset: Annotated[
        float, typer.Option(help="S-bend fan-out Y-offset (microns).")
    ] = 5.0,
    sbend_x_extent: Annotated[
        float, typer.Option(help="S-bend X-extent (microns).")
    ] = 10.0,
    center_wavelength: Annotated[
        float, typer.Option(help="Center wavelength (microns).")
    ] = 1.55,
    monitor_offset_from_sbend: Annotated[
        float, typer.Option(help="Monitor X-offset from S-bend end (microns).")
    ] = 2.5,
    fdtd_xy_padding: Annotated[
        float, typer.Option(help="FDTD XY padding (microns).")
    ] = 2.5,
    wg_extension_past_fdtd_edge: Annotated[
        float, typer.Option(help="WG extension beyond FDTD edge (microns).")
    ] = 1.0,
    output_dir: Annotated[
        pathlib.Path, typer.Option(help="Base directory for simulation outputs (layout-specific folders will be created here).")
    ] = DEFAULT_OUTPUT_DIR,
    plot_z_plane: Annotated[
        bool, typer.Option(help="Plot Z-plane intensity after simulation.")
    ] = False,
    show_gui: Annotated[
        bool, typer.Option(help="Show Lumerical FDTD CAD window during simulation.")
    ] = True,
    # Removed run_id, sim_dir, results_dir as they are now handled internally by layout_id and output_dir
):
    """
    Run a Lumerical FDTD simulation with specified parameters.
    All length parameters are in MICRONS.
    Outputs (.fsp, results.json, plots) are saved in a unique directory under `output_dir` for each parameter set.
    If a simulation for the given parameters already exists and has valid results, it will be skipped.
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
    params_meters = {key: value * u for key, value in params_microns.items()}

    layout_id = create_layout_id(params_meters, u_val=u)
    logger.info(f"Target layout ID: {layout_id}")
    logger.info(f"Parameters (microns): {params_microns}")
    logger.info(f"Output base directory: {output_dir}")

    run_simulation(
        params=params_meters,
        output_base_dir=str(output_dir),
        plot_z_plane_each_run=plot_z_plane,
        hide_fdtd_gui=not show_gui,
    )
    # The run_simulation function now handles skipping and also calls results processing internally.
    logger.info(f"Simulation and processing for layout {layout_id} finished.")


@app.command("analyse")
@app.command("a", hidden=True)
def analyse_results(
    output_dir: Annotated[
        pathlib.Path,
        typer.Option(help="Base directory containing layout-specific simulation output folders."),
    ] = DEFAULT_OUTPUT_DIR,
    plot_z_plane: Annotated[
        bool, typer.Option(help="Generate/regenerate Z-plane intensity plot during analysis.")
    ] = True,
):
    """
    Analyse existing Lumerical FDTD simulation results.
    Scans the `output_dir` for layout folders (identified by having a `results.json` file).
    Presents a list of available layouts for selection, then displays their results.
    """
    logger.info(f"Scanning for layout folders with results.json in: {output_dir}")

    layout_folders: List[pathlib.Path] = []
    for item in output_dir.iterdir():
        if item.is_dir() and (item / "results.json").exists():
            layout_folders.append(item)

    if not layout_folders:
        logger.warning(f"No layout folders with results.json found in {output_dir}.")
        typer.echo(f"No layout folders with results.json found in {output_dir}.")
        raise typer.Exit(code=1)

    # Sort folders by modification time of their results.json, most recent first
    layout_folders.sort(key=lambda x: (x / "results.json").stat().st_mtime, reverse=True)

    choices = []
    for folder in layout_folders:
        results_json_path = folder / "results.json"
        try:
            with open(results_json_path, "r") as f:
                results_data = json.load(f)
            # Display key results for quick identification
            # Using .get for robustness if keys are missing
            tr_val = results_data.get("T_net_through_TL", float('nan'))
            br_val = results_data.get("T_net_coupled_BL", float('nan'))
            # Format to a reasonable number of decimal places
            title = f"{folder.name} (TR: {tr_val:.3f}, BL: {br_val:.3f}, mod: {results_json_path.stat().st_mtime:.0f})"
            if "error" in results_data:
                title += " [ERROR]"
        except Exception as e:
            title = f"{folder.name} (Error reading results: {e}, mod: {results_json_path.stat().st_mtime:.0f})"

        choices.append(questionary.Choice(title=title, value=str(folder)))

    selected_layout_folder_str = questionary.select(
        "Select a layout folder to analyse:",
        choices=choices,
        use_shortcuts=True,
    ).ask()

    if not selected_layout_folder_str:
        logger.info("No folder selected. Exiting analysis.")
        typer.echo("No folder selected. Exiting.")
        raise typer.Exit()

    selected_layout_folder = pathlib.Path(selected_layout_folder_str)
    layout_id = selected_layout_folder.name
    sim_fsp_path = selected_layout_folder / f"{layout_id}.fsp"
    results_json_path = selected_layout_folder / "results.json"

    logger.info(f"Analysing layout: {layout_id} in folder: {selected_layout_folder}")

    if not results_json_path.exists():
        logger.error(f"Critical: results.json not found in {selected_layout_folder}, though it was expected.")
        typer.echo(f"Error: results.json missing in {selected_layout_folder}.")
        raise typer.Exit(code=1)

    with open(results_json_path, "r", encoding="utf-8") as f:
        results_data = json.load(f)

    typer.echo(f"\n--- Results for Layout: {layout_id} ---")
    typer.echo(json.dumps(results_data, indent=2))

    # Parameters are stored in results_data["parameters_meters"]
    params_meters_from_results = results_data.get("parameters_meters")

    if plot_z_plane:
        plot_filename = f"z_plane_intensity_{layout_id}.png"
        plot_save_path = selected_layout_folder / plot_filename

        if plot_save_path.exists():
            logger.info(f"Plot {plot_filename} already exists at {plot_save_path}.")
            typer.echo(f"Plot available at: {plot_save_path}")
        elif sim_fsp_path.exists() and params_meters_from_results:
            logger.info(f"Plot {plot_filename} not found. Generating from: {sim_fsp_path}")
            try:
                # Need lumapi for plotting if we regenerate it here
                import lumapi  # Ensure lumapi is available in this scope if not already
                with lumapi.FDTD(hide=True) as fdtd:
                    fdtd.load(str(sim_fsp_path))
                    # Check if mon_zplane exists in the loaded simulation
                    monitor_list = fdtd.get("monitorlist")
                    if "mon_zplane" not in monitor_list:
                        logger.warning("'mon_zplane' monitor not found in the .fsp file. Cannot generate plot.")
                        typer.echo("Warning: 'mon_zplane' monitor not found in the .fsp file. Cannot generate plot.")
                    else:
                        # Use results from JSON for title consistency
                        tr_val = results_data.get("T_net_through_TL", float('nan'))
                        bl_val = results_data.get("T_net_coupled_BL", float('nan'))
                        plot_title_prefix = (
                            f"Z-plane E-field Intensity (TL={tr_val:.4f}, BL={bl_val:.4f})"
                        )
                        plot_plane_parametric(
                            fdtd_obj=fdtd,
                            title_prefix=plot_title_prefix,
                            layout_id=layout_id,
                            target_dir=selected_layout_folder,
                        )
                        typer.echo(f"Generated plot: {plot_save_path}")

            except ImportError:
                logger.error("Lumerical API (lumapi) not found. Cannot generate plot.")
                typer.echo("Error: Lumerical API not found. Plot generation failed.")
            except Exception as e_plot:
                logger.error(f"Could not generate plot for {layout_id}: {e_plot}")
                typer.echo(f"Error generating plot: {e_plot}")
        elif not params_meters_from_results:
            logger.warning(f"Parameters not found in results.json. Cannot reliably (re)generate plot for {layout_id}.")
            typer.echo("Warning: Parameters missing in results.json, cannot generate plot.")
        else:
            logger.warning(f"Simulation file {sim_fsp_path} not found. Cannot generate plot for {layout_id}.")
            typer.echo(f"Warning: Simulation file {sim_fsp_path.name} missing, cannot generate plot.")

    typer.echo(f"--- End of Analysis for {layout_id} ---")


if __name__ == "__main__":
    app()

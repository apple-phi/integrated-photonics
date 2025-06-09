"""Processes and saves results from Lumerical FDTD simulations."""

import os
import sys
import pathlib
import json
import csv

import numpy as np
import matplotlib.pyplot as plt

# Ensure lumapi can be imported (assuming it's on the path or handled by the caller)
# If lumapi is not found by Python's default search, the dynamic path addition
# from the original scripts should be included here or in the calling script.
# For now, we assume lumapi will be importable.
try:
    import lumapi  # type: ignore[import-untyped]
except ImportError as e:
    # Attempt to find lumapi dynamically if direct import fails
    lumerical_install_dir = pathlib.Path(
        os.getenv("LUMERICAL_INSTALL_DIR", "C:\\Program Files\\Lumerical")
    )
    found_lumapi = False
    # More specific search for lumapi.py within typical Lumerical structures
    for version_dir in lumerical_install_dir.iterdir():
        if version_dir.is_dir():
            py_api_dir = version_dir / "api" / "python"
            if (py_api_dir / "lumapi.py").exists():
                sys.path.append(str(py_api_dir))
                print(f"Found lumapi at: '{py_api_dir / "lumapi.py"}'")
                import lumapi  # type: ignore[import-untyped]

                found_lumapi = True
                break
            # Older Lumerical versions might have it in a different spot (e.g. /python/lumapi.py directly in version_dir)
            alt_py_api_dir = version_dir / "python"
            if (alt_py_api_dir / "lumapi.py").exists():
                sys.path.append(str(alt_py_api_dir))
                print(f"Found lumapi at: '{alt_py_api_dir / "lumapi.py"}'")
                import lumapi  # type: ignore[import-untyped]

                found_lumapi = True
                break
    if not found_lumapi:
        raise ImportError(
            "lumapi.py not found. Please ensure Lumerical is installed and "
            "LUMERICAL_INSTALL_DIR environment variable is set, or lumapi.py is in PYTHONPATH."
        ) from e

u = 1e-6


def plot_plane_parametric(
    fdtd_obj, title_prefix: str, run_id: str, target_dir: pathlib.Path
):
    """Helper function to plot a 2D plane from the FDTD result and save it."""
    res_zplane = fdtd_obj.getresult("mon_zplane", "E")
    x = res_zplane["x"].flatten() * 1e6  # Convert to µm
    y = res_zplane["y"].flatten() * 1e6  # Convert to µm
    E = res_zplane["E"]  # shape (nx, ny, 1, 1, 3)

    Ex = E[:, :, 0, 0, 0]
    Ey = E[:, :, 0, 0, 1]
    Ez = E[:, :, 0, 0, 2]

    Intensity = np.abs(Ex) ** 2 + np.abs(Ey) ** 2 + np.abs(Ez) ** 2
    X, Y = np.meshgrid(x, y, indexing="ij")

    plt.figure(figsize=(8, 5))
    pcm = plt.pcolormesh(
        X, Y, Intensity, shading="auto", cmap="viridis"
    )  # Use a common colormap
    plt.xlabel("x (μm)")
    plt.ylabel("y (μm)")
    plot_title = f"{title_prefix} - {run_id}"
    plt.title(plot_title)
    plt.colorbar(pcm, label="Intensity (a.u.)")
    plt.axis("equal")  # Ensure aspect ratio is maintained

    plot_filename = f"z_plane_intensity_{run_id}.png"
    save_path = target_dir / plot_filename
    plt.savefig(save_path)
    print(f"Saved z-plane plot to: {save_path}")
    plt.close()  # Close the plot to free memory


def process_and_save_results(
    sim_filepath: str,
    params_dict: dict,
    run_id: str,
    results_base_dir: str,
    plot_z_plane: bool = False,
):
    """Loads results from a Lumerical simulation, processes, and saves them."""
    results_base_path = pathlib.Path(results_base_dir)
    current_run_results_dir = results_base_path / run_id
    current_run_results_dir.mkdir(parents=True, exist_ok=True)

    print(f"Processing results for run: {run_id} from file: {sim_filepath}")

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
                    run_id=run_id,
                    target_dir=current_run_results_dir,
                )

    except Exception as e:
        print(f"Error during Lumerical results processing for {run_id}: {e}")
        results_data["T_net_tr"] = float("nan")  # Indicate error in results
        results_data["T_net_br"] = float("nan")
        results_data["error"] = str(e)

    # Save detailed results for this run to JSON
    results_json_path = current_run_results_dir / "results.json"
    with open(results_json_path, "w", encoding="utf-8") as f_json:
        json.dump(results_data, f_json, indent=4)
    print(f"Saved detailed results to: {results_json_path}")

    # Append to summary CSV
    summary_csv_path = results_base_path / "summary_results.csv"
    summary_data_row = {
        "run_id": run_id,
        **params_dict,  # Unpack parameters into the row
        "T_net_tr": results_data.get("T_net_tr", float("nan")),
        "T_net_br": results_data.get("T_net_br", float("nan")),
        "results_path": str(current_run_results_dir),
        "error": results_data.get("error", ""),
    }

    file_exists = summary_csv_path.exists()
    with open(summary_csv_path, "a", newline="", encoding="utf-8") as f_csv:
        # Define fieldnames based on the keys in summary_data_row to ensure order and completeness
        # It's important that params_dict keys are consistent for CSV header
        # For simplicity, assuming params_dict keys are stable. A more robust solution
        # might involve explicitly defining headers or getting them from the first run.
        fieldnames = list(summary_data_row.keys())
        writer = csv.DictWriter(f_csv, fieldnames=fieldnames)
        if (
            not file_exists or summary_csv_path.stat().st_size == 0
        ):  # Check if file is new or empty
            writer.writeheader()
        writer.writerow(summary_data_row)
    print(f"Appended summary to: {summary_csv_path}")


if __name__ == "__main__":
    # Example usage (for testing this script directly)
    print("Running example for task3_res_parametric.py")
    # This example assumes a simulation file exists from a previous run.
    # In a real sweep, task3_sim_parametric.py would call process_and_save_results.

    # Create dummy parameters and paths for the example
    example_run_id = "example_w1_0.5_w2_0.5_sep_0.1_len_10"
    example_params = {
        "wg1_width": 0.5e-6,
        "wg2_width": 0.5e-6,
        "separation": 0.1e-6,
        "coupling_length": 10e-6,
        "wg_z_span": 0.22e-6,  # Added for completeness if needed by CSV header
        "fan_out_y_offset": 5e-6,
        "sbend_x_extent": 5e-6,
        "wg_io_length": 5e-6,
    }
    # IMPORTANT: For this example to run, you need an actual .fsp file.
    # Let's assume one was created by task3_sim_parametric.py in a known location.
    # If you have run task3_sim_parametric.py with default_params, it might be:
    # c:/Users/ln373/Desktop/integrated-photonics/data/task3_simulations/default_run/simulation_default_run.fsp
    # Update this path to an actual .fsp file from your system for testing.
    example_sim_file = pathlib.Path(
        r"C:\Users\ln373\Desktop\integrated-photonics\data\task3_simulations_parametric_example\w1_0.50_w2_0.50_sep_0.15_L_10_sbX_10_sbY_5_monOff_2.5_fdtdPad_2.5_wgExt_1.0\simulation_w1_0.50_w2_0.50_sep_0.15_L_10_sbX_10_sbY_5_monOff_2.5_fdtdPad_2.5_wgExt_1.0.fsp"
    )
    example_results_dir = "data/task3_results_parametric_example"

    if example_sim_file.exists():
        # Create the directory for example results
        pathlib.Path(example_results_dir).mkdir(parents=True, exist_ok=True)
        process_and_save_results(
            sim_filepath=str(example_sim_file),
            params_dict=example_params,
            run_id=example_run_id,
            results_base_dir=example_results_dir,
            plot_z_plane=True,
        )
        print(f"Example processing complete. Check '{example_results_dir}'.")
    else:
        print(
            f"Example simulation file '{example_sim_file}' not found. Skipping example run."
        )
        print(
            "Please update 'example_sim_file' path in task3_res_parametric.py for direct testing."
        )

    # Example of how to read the summary CSV with pandas (if you want to use it for analysis)
    # import pandas as pd
    # summary_df = pd.read_csv(pathlib.Path(example_results_dir) / "summary_results.csv")
    # print("\nPandas DataFrame from summary_results.csv:")
    # print(summary_df.head())

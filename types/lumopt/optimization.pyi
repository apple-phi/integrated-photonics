from _typeshed import Incomplete
from lumopt.lumerical_methods.lumerical_scripts import get_fields as get_fields, get_fields_on_cad as get_fields_on_cad, get_lambda_from_cad as get_lambda_from_cad
from lumopt.utilities.base_script import BaseScript as BaseScript
from lumopt.utilities.fields import FieldsNoInterp as FieldsNoInterp
from lumopt.utilities.gradients import GradientFields as GradientFields
from lumopt.utilities.plotter import Plotter as Plotter
from lumopt.utilities.simulation import Simulation as Simulation
from lumopt.utilities.wavelengths import Wavelengths as Wavelengths

class SuperOptimization:
    plot_history: Incomplete
    fields_on_cad_only: Incomplete
    plotter: Incomplete
    optimizations: Incomplete
    weights: Incomplete
    old_dir: Incomplete
    full_fom_hist: Incomplete
    fom_hist: Incomplete
    params_hist: Incomplete
    grad_hist: Incomplete
    continuation_max_iter: int
    def __init__(self, optimizations: Incomplete | None = None, plot_history: bool = False, fields_on_cad_only: bool = False, weights: Incomplete | None = None) -> None: ...
    def __add__(self, other): ...
    def __del__(self) -> None: ...
    one_forward: bool
    optimizer: Incomplete
    target_fom: Incomplete
    fom_names: Incomplete
    plot_fom_on_log_scale: Incomplete
    last_grad: Incomplete
    def initialize(self, start_params: Incomplete | None = None, bounds: Incomplete | None = None, working_dir: Incomplete | None = None): ...
    def init_plotter(self) -> None: ...
    def plot_fom(self, fomax, paramsax, gradients_ax) -> None: ...
    def plot_gradient(self, fig, ax_fields, ax_gradients) -> None: ...
    workingDir: Incomplete
    def prepare_working_dir(self, working_dir) -> None: ...
    num_threads: Incomplete
    calling_file_name: Incomplete
    base_file_path: Incomplete
    def run(self, working_dir: Incomplete | None = None, num_threads: Incomplete | None = None): ...

class Optimization(SuperOptimization):
    base_script: Incomplete
    wavelengths: Incomplete
    fom: Incomplete
    geometry: Incomplete
    optimizer: Incomplete
    use_var_fdtd: Incomplete
    hide_fdtd_cad: Incomplete
    source_name: Incomplete
    use_deps: bool
    custom_deps: Incomplete
    store_all_simulations: Incomplete
    save_global_index: Incomplete
    unfold_symmetry: Incomplete
    label: Incomplete
    plot_fom_on_log_scale: Incomplete
    calling_file_name: Incomplete
    base_file_path: Incomplete
    def __init__(self, base_script, wavelengths, fom, geometry, optimizer, use_var_fdtd: bool = False, hide_fdtd_cad: bool = False, use_deps: bool = True, plot_history: bool = True, store_all_simulations: bool = True, save_global_index: bool = False, label: Incomplete | None = None, source_name: str = 'source', fields_on_cad_only: bool = False) -> None: ...
    def check_gradient(self, test_params, dx, working_dir: Incomplete | None = None): ...
    def run(self, working_dir: Incomplete | None = None): ...
    def plotting_function(self, params) -> None: ...
    sim: Incomplete
    fom_hist: Incomplete
    def initialize(self, working_dir) -> None: ...
    def save_fields_to_vtk(self, cur_iteration) -> None: ...
    def save_index_to_vtk(self, cur_iteration) -> None: ...
    def make_forward_sim(self, params, iter, co_optimizations: Incomplete | None = None, one_forward: bool = False): ...
    forward_fields_wl: Incomplete
    forward_fields: Incomplete
    forward_fields_iter: Incomplete
    def process_forward_sim(self, iter, co_optimizations: Incomplete | None = None, one_forward: bool = False): ...
    def callable_fom(self, params): ...
    def make_adjoint_sim(self, params, iter, co_optimizations: Incomplete | None = None, one_forward: bool = False): ...
    adjoint_fields: Incomplete
    scaling_factor: Incomplete
    def process_adjoint_sim(self, iter, co_optimizations: Incomplete | None = None, one_forward: bool = False) -> None: ...
    last_grad: Incomplete
    def callable_jac(self, params): ...
    gradient_fields: Incomplete
    gradients: Incomplete
    def calculate_gradients(self): ...
    def plot_gradient(self, fig, ax1, ax2) -> None: ...
    @staticmethod
    def add_index_monitor(sim, monitor_name, wavelengths) -> None: ...
    @staticmethod
    def cross_section_monitor_props(monitor_type): ...
    @staticmethod
    def set_global_wavelength(sim, wavelengths) -> None: ...
    @staticmethod
    def set_source_wavelength(sim, source_name, multi_freq_src, freq_pts) -> None: ...
    @staticmethod
    def set_use_legacy_conformal_interface_detection(sim, flagVal) -> None: ...
    @staticmethod
    def check_simulation_was_successful(sim): ...
    @staticmethod
    def deactivate_all_sources(sim) -> None: ...

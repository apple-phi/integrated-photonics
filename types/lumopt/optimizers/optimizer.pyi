from _typeshed import Incomplete

class Optimizer:
    max_iter: Incomplete
    scaling_factor: Incomplete
    scale_initial_gradient_to: Incomplete
    penalty_fun: Incomplete
    penalty_jac: Incomplete
    penalty_jac_approx: Incomplete
    logging_path: Incomplete
    logfile: Incomplete
    current_fom: Incomplete
    current_gradients: Incomplete
    current_params: Incomplete
    fom_hist: Incomplete
    gradients_hist: Incomplete
    params_hist: Incomplete
    iteration: int
    fom_scaling_factor: int
    fom_calls: int
    def __init__(self, max_iter, scaling_factor: Incomplete | None = None, scale_initial_gradient_to: int = 0, penalty_fun: Incomplete | None = None, penalty_jac: Incomplete | None = None, logging_path: Incomplete | None = None) -> None: ...
    bounds: Incomplete
    scaling_offset: Incomplete
    def initialize(self, start_params, callable_fom, callable_jac, bounds, plotting_function) -> None: ...
    start_point: Incomplete
    def reset_start_params(self, start_params, scale_initial_gradient_to) -> None: ...
    def auto_detect_scaling(self, min_required_rel_change) -> None: ...
    callback: Incomplete
    def define_callback(self, plotting_function) -> None: ...
    def report_writing(self) -> None: ...
    def concurrent_adjoint_solves(self): ...
    @staticmethod
    def create_jac_approx(fom_func): ...
